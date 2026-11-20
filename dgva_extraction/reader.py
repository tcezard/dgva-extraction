"""Database reader, core fields (sprint 3, ticket 13): streams the variant records of one study on one assembly.

It applies the rules written down in sprints 1 and 2 (reference/dgva-oracle-schema.md):

- submitted or remapped: plan step 1.4. Only submitted placements become records; a "conflict" (submitted method with
  a remap score) or an unknown method fails the study.
- W1: a call with no region is written, with no PARENT. Region-to-region links are kept.
- W2: a type with no SO term, a region with no assertion method or an unknown assembly fails the study.
- W3: a placement with no start, or with a start but no chromosome, is not written and goes on the side list.
- confidence intervals: CIPOS*/CIEND* first, then CI_START/CI_STOP (signed in estd219, unsigned elsewhere), then the
  inner and outer coordinates.
- identical placements of one accession (estd208) are written once; the others go on the side list.
- translocation regions with no placement of their own are positioned from their supporting calls (nstd83).
"""
import heapq
import itertools
import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Iterator, List, Optional

from dgva_extraction.keys import SymbolicAlleles
from dgva_extraction.record import ConfidenceInterval, RecordKind, VariantRecord

logger = logging.getLogger(__name__)

NO_POSITION_METHODS = frozenset({'Remap failed', 'Submitted without placement on assembly',
                                 'Submitted on non-INSDC assembly'})
# Studies whose VARIANT_PLACEMENT_CI_START/_CI_STOP hold VCF-style signed values (CI_MINUS <= 0): "Open checks:
# measured" in the reference document. Every other study holds unsigned distances (estd204).
SIGNED_CI_STUDIES = frozenset({'estd219'})
# Columns that must be the same on every row of one placement (rows repeat per parent and per duplicated CI row).
PLACEMENT_COLUMNS = (
    'accession', 'region_or_call', 'so_term', 'so_id', 'exchange_type', 'in_type_table', 'assertion_method',
    'assembly_known', 'placement_method', 'has_remap_score', 'chromosome', 'start_pos', 'stop_pos',
    'outer_start_pos', 'inner_start_pos', 'inner_stop_pos', 'outer_stop_pos', 'ci_start_minus', 'ci_start_plus',
    'ci_stop_minus', 'ci_stop_plus', 'ciposleft', 'ciposright', 'ciendleft', 'ciendright',
)


class ExtractionError(RuntimeError):
    """The study can't be converted as it is in the database (W2, conflicting rows, unknown method...)."""


@dataclass(frozen=True)
class SideListEntry:
    study: str
    accession: str
    kind: str
    reason: str
    assembly: Optional[str] = None
    placement_id: Optional[int] = None
    detail: str = ''

    COLUMNS = ('study', 'assembly', 'accession', 'kind', 'reason', 'placement_id', 'detail')

    def row(self):
        return [self.study, self.assembly or '', self.accession, self.kind, self.reason,
                '' if self.placement_id is None else str(self.placement_id), self.detail]


def write_side_list(path, entries: Iterable[SideListEntry]):
    with open(path, 'w') as handle:
        handle.write('\t'.join(SideListEntry.COLUMNS) + '\n')
        for entry in sorted(entries, key=lambda item: (item.accession, item.reason, item.placement_id or 0)):
            handle.write('\t'.join(value.replace('\t', ' ') for value in entry.row()) + '\n')


def placement_class(method, has_remap_score):
    """The ticket 2 rule (plan step 1.4)."""
    if method in NO_POSITION_METHODS:
        return 'no_position'
    if method is not None and method.startswith('Remapped'):
        return 'remapped'
    if method is None or method.startswith('Submitted'):
        return 'conflict' if has_remap_score else 'submitted'
    return 'other'


def _as_int(value):
    return None if value is None else int(value)


def _interval_from_unsigned(left, right):
    """CIPOSLEFT/RIGHT and unsigned CI_MINUS/CI_PLUS are distances: left 17 means -17."""
    if left is None and right is None:
        return None
    if (left is not None and left < 0) or (right is not None and right < 0):
        raise ExtractionError(f'Unsigned confidence interval with a negative value: {left}, {right}')
    return ConfidenceInterval(None if left is None else -int(left), _as_int(right))


def _interval_from_signed(minus, plus):
    if minus is None and plus is None:
        return None
    if (minus is not None and minus > 0) or (plus is not None and plus < 0):
        raise ExtractionError(f'Signed confidence interval that does not contain 0: {minus}, {plus}')
    return ConfidenceInterval(_as_int(minus), _as_int(plus))


@dataclass
class StudyReader:
    """Reads one study on one submitted assembly. Iterate records(); then read side_list and counts."""
    database: object
    study: str
    assembly: str
    alleles: SymbolicAlleles = field(default_factory=SymbolicAlleles.load)
    signed_ci_studies: frozenset = SIGNED_CI_STUDIES
    side_list: List[SideListEntry] = field(default_factory=list)
    counts: Counter = field(default_factory=Counter)

    @property
    def binds(self):
        return {'study': self.study, 'assembly': self.assembly}

    def records(self) -> Iterator[VariantRecord]:
        """Records in database order (chromosome name, start, accession), translocation regions merged in."""
        placements = self._placement_records(self.database.stream('study_records', self.binds))
        translocations = self._translocation_records(self.database.stream('translocation_regions', self.binds))
        return heapq.merge(placements, translocations, key=lambda record: (record.chromosome, record.start,
                                                                          record.accession))

    # --- Placements ---

    def _placement_records(self, rows) -> Iterator[VariantRecord]:
        seen_at_position = {}
        position = None
        for placement_id, group in itertools.groupby(rows, key=lambda row: row['placement_id']):
            group = list(group)
            self.counts['rows'] += len(group)
            self.counts['placements'] += 1
            row = self._single_placement_row(placement_id, group)
            record = self._placement_record(placement_id, row, group)
            if record is None:
                continue
            # estd208: the same accession placed several times with identical coordinates is written once.
            if (record.chromosome, record.start) != position:
                position = (record.chromosome, record.start)
                seen_at_position = {}
            identity = (record.accession, record.end, record.outer_start, record.inner_start, record.inner_end,
                        record.outer_end, record.ci_pos, record.ci_end, record.placement_method, record.parents)
            if identity in seen_at_position:
                self.counts['duplicate_placements'] += 1
                self._side(record.accession, record.kind.value, 'duplicate_placement', placement_id,
                           f'identical to placement {seen_at_position[identity]}')
                continue
            seen_at_position[identity] = placement_id
            self.counts['records'] += 1
            yield record

    @staticmethod
    def _single_placement_row(placement_id, group):
        first = group[0]
        for other in group[1:]:
            different = [column for column in PLACEMENT_COLUMNS if other.get(column) != first.get(column)]
            if different:
                raise ExtractionError(f'Placement {placement_id} ({first["accession"]}) has rows that disagree on '
                                      f'{", ".join(different)} (duplicated rows in a table without a primary key)')
        return first

    def _placement_record(self, placement_id, row, group) -> Optional[VariantRecord]:
        accession = row['accession']
        kind = row['region_or_call']
        cls = placement_class(row['placement_method'], row['has_remap_score'])
        if cls != 'submitted':
            raise ExtractionError(f'{accession}: placement {placement_id} has method {row["placement_method"]!r} '
                                  f'and remap score flag {row["has_remap_score"]}: class {cls}, decided case by case')
        self._check_w2(accession, kind, row)
        if row['placement_method'] is None:
            self.counts['null_method_placements'] += 1

        start = row['start_pos']
        if start is None:
            start = row['outer_start_pos'] if row['outer_start_pos'] is not None else row['inner_start_pos']
        end = row['stop_pos']
        if end is None:
            end = row['outer_stop_pos'] if row['outer_stop_pos'] is not None else row['inner_stop_pos']
        if start is None:
            self._side(accession, kind, 'no_start', placement_id, 'W3: submitted placement with no start')
            return None
        if row['chromosome'] is None:
            self._side(accession, kind, 'no_chromosome', placement_id,
                       f'W3: submitted placement starting at {start} with no chromosome or contig name')
            return None
        if end is None:
            allele = self.alleles.for_term(row['so_term'])
            if allele is not None and allele.end_is_pos:
                end = start
            else:
                self._side(accession, kind, 'no_end', placement_id, f'submitted placement starting at {start} '
                                                                    f'with no end')
                return None

        ci_pos, ci_end = self._confidence_intervals(accession, row)
        parents = self._parents(accession, group)
        info = {}
        if kind == RecordKind.REGION.value:
            info['ASSERTM'] = row['assertion_method']
        return VariantRecord(
            accession=accession, kind=kind, so_term=row['so_term'], so_id=row['so_id'],
            chromosome=row['chromosome'], start=int(start), end=int(end),
            outer_start=_as_int(row['outer_start_pos']), inner_start=_as_int(row['inner_start_pos']),
            inner_end=_as_int(row['inner_stop_pos']), outer_end=_as_int(row['outer_stop_pos']),
            ci_pos=ci_pos, ci_end=ci_end, parents=parents, placement_method=row['placement_method'],
            placement_id=placement_id, info=info,
        )

    def _check_w2(self, accession, kind, row):
        if not row['assembly_known']:
            raise ExtractionError(f'{accession}: assembly {self.assembly} is not in GENOME_ASSEMBLY (W2)')
        if not row['in_type_table']:
            raise ExtractionError(f'{accession}: not in VARIANT_{kind.upper()} (W2)')
        if not row['so_term'] or not row['so_id']:
            raise ExtractionError(f'{accession}: type {row["exchange_type"]!r} has no SO term (W2)')
        if kind == RecordKind.REGION.value and not row['assertion_method']:
            raise ExtractionError(f'{accession}: region with no assertion method (W2)')

    def _confidence_intervals(self, accession, row):
        """CIPOS*/CIEND* first; then CI_START/CI_STOP; else None (the writer derives them from inner/outer)."""
        try:
            ci_pos = _interval_from_unsigned(row['ciposleft'], row['ciposright'])
            ci_end = _interval_from_unsigned(row['ciendleft'], row['ciendright'])
            if self.study in self.signed_ci_studies:
                from_start = _interval_from_signed(row['ci_start_minus'], row['ci_start_plus'])
                from_stop = _interval_from_signed(row['ci_stop_minus'], row['ci_stop_plus'])
            else:
                from_start = _interval_from_unsigned(row['ci_start_minus'], row['ci_start_plus'])
                from_stop = _interval_from_unsigned(row['ci_stop_minus'], row['ci_stop_plus'])
        except ExtractionError as error:
            raise ExtractionError(f'{accession}: {error}') from error
        if (ci_pos or ci_end) and (from_start or from_stop):
            # Never seen (ticket 8): the two sets come from different studies.
            raise ExtractionError(f'{accession}: has both CIPOS/CIEND and CI_START/CI_STOP rows')
        if from_start or from_stop:
            self.counts['ci_from_ci_start_stop'] += 1
            return from_start, from_stop
        if ci_pos or ci_end:
            self.counts['ci_from_cipos_ciend'] += 1
        return ci_pos, ci_end

    def _parents(self, accession, group):
        parents = set()
        for row in group:
            parent = row['parent_accession']
            if parent is None:
                continue
            if row['parent_withdrawn']:
                self.counts['withdrawn_parents_left_out'] += 1
                continue
            if not row['parent_in_study']:
                raise ExtractionError(f'{accession}: linked to region {parent}, which is not in study {self.study}')
            parents.add(parent)
        if not parents and group[0]['region_or_call'] == RecordKind.CALL.value:
            self.counts['calls_without_region'] += 1
        return tuple(sorted(parents))

    # --- Translocation regions with no placement of their own ---

    def _translocation_records(self, rows) -> Iterator[VariantRecord]:
        records = []
        for accession, group in itertools.groupby(rows, key=lambda row: row['accession']):
            group = sorted(group, key=lambda row: (row['chromosome'] or '', row['start_pos'] or 0,
                                                   row['call_accession'], row['call_placement_id']))
            first = group[0]
            if not first['so_term'] or not first['so_id']:
                raise ExtractionError(f'{accession}: type {first["exchange_type"]!r} has no SO term (W2)')
            if not first['assertion_method']:
                raise ExtractionError(f'{accession}: region with no assertion method (W2)')
            placed = [row for row in group if row['chromosome'] is not None and row['start_pos'] is not None]
            if not placed:
                continue
            anchor = placed[0]
            calls = sorted({row['call_accession'] for row in group})
            self.counts['translocation_regions'] += 1
            self._side(accession, RecordKind.REGION.value, 'positioned_from_calls', None,
                       f'no placement of its own; positioned from {", ".join(calls)}')
            records.append(VariantRecord(
                accession=accession, kind=RecordKind.REGION, so_term=first['so_term'], so_id=first['so_id'],
                chromosome=anchor['chromosome'], start=int(anchor['start_pos']),
                end=int(anchor['stop_pos'] if anchor['stop_pos'] is not None else anchor['start_pos']),
                placement_method=anchor['placement_method'], info={'ASSERTM': first['assertion_method']},
            ))
        records.sort(key=lambda record: (record.chromosome, record.start, record.accession))
        return iter(records)

    def _side(self, accession, kind, reason, placement_id, detail):
        self.counts[f'side_list_{reason}'] += 1
        self.side_list.append(SideListEntry(self.study, accession, kind, reason, self.assembly, placement_id, detail))


def study_side_list(database, study) -> List[SideListEntry]:
    """The variants of a study that are in no VCF: withdrawn, remapped only, no position, no placement."""
    entries = []
    for row in database.stream('study_side_list', {'study': study}):
        detail = ''
        if row['reason'] == 'withdrawn':
            detail = f'{row["withdrawal_reason"]} (withdrawn by {row["withdrawn_by_study"]})'
        elif row['reason'] in ('remapped_only', 'no_position_only'):
            detail = (f'{row["n_placements"]} placement(s): {row["n_remapped"]} remapped, '
                      f'{row["n_no_position"]} with no position')
        entries.append(SideListEntry(study, row['accession'], row['region_or_call'], row['reason'], detail=detail))
    return entries
