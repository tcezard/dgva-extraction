"""VCF writer: turns variant records into a valid, deterministic VCF (sprint 3, ticket 12).

Rules (proposed in ticket 10, see docs/vcf-representation.md):

- one VCF record per DGVa accession and placement, with the accession in the ID column; records at the same position
  are never merged;
- a record whose type has no symbolic allele, or that carries a value with no key in the key table, is refused, and a
  study with a refused record fails: no VCF is produced, and the statistics list every refusal with its reason;
- lines are sorted by contig in assembly-report order, then by position and accession; header lines and INFO keys
  are in a fixed order. The same records give the same bytes.

The writer streams: records must arrive grouped by chromosome and ordered by start within it, as the database query
returns them. Each contig is written to its own temporary file as records arrive (only records sharing one position
are held in memory), and the files are joined in assembly-report order at the end.
"""
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from dgva_extraction.keys import KeyTable, SymbolicAlleles, UnmappedValueError
from dgva_extraction.record import VariantRecord
from dgva_extraction.reference import Contig, Genome, ReferenceLookupError


VCF_VERSION = 'VCFv4.4'
DEFAULT_PLACEMENT_METHOD = 'Submitted genomic'
VALID_BASES = set('ACGTN')
# IUPAC ambiguity codes in the reference are replaced by the first base, in alphabetical order, that they stand for
# (as convertGVFtoVCF does), because VCF allows only A, C, G, T and N in REF.
IUPAC_CODES = {
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
}
HEADER_COLUMNS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
# VCF 4.4: the tandem repeat keys belong to <CNV:TR> records, which need RUS or RUL.
TANDEM_REPEAT_ALLELE = 'CNV:TR'
TANDEM_REPEAT_KEYS = {'RUC', 'RUS', 'RUL', 'RUB', 'RB', 'CIRUC', 'CIRB', 'RUEFF'}


class RecordError(ValueError):
    """The record can't be written. `reason` is a short code used in the statistics."""

    def __init__(self, reason, message):
        super().__init__(message)
        self.reason = reason


class OrderError(RuntimeError):
    """Records did not arrive grouped by contig and ordered by position."""


class StudyFailedError(RuntimeError):
    def __init__(self, statistics):
        self.statistics = statistics
        reasons = ', '.join(f'{reason}: {count}' for reason, count in sorted(statistics.refused.items()))
        super().__init__(f'{sum(statistics.refused.values())} record(s) refused ({reasons}); no VCF was written')


def normalise_reference(bases):
    """Upper case, IUPAC codes replaced; any other character is refused."""
    converted = []
    for base in bases.upper():
        if base in VALID_BASES:
            converted.append(base)
        elif base in IUPAC_CODES:
            converted.append(IUPAC_CODES[base][0])
        else:
            raise RecordError('invalid_reference_base', f'Reference base {base!r} is not a nucleotide or IUPAC code')
    return ''.join(converted)


@dataclass
class VcfLine:
    contig: Contig
    pos: int
    id: str
    ref: str
    alt: str
    info: List[Tuple[str, Optional[str]]]
    format_keys: List[str] = field(default_factory=list)
    sample_values: List[List[str]] = field(default_factory=list)
    symbolic_allele: Optional[str] = None

    def text(self):
        info = ';'.join(key if value is None else f'{key}={value}' for key, value in self.info) or '.'
        columns = [self.contig.vcf_name, str(self.pos), self.id, self.ref, self.alt, '.', '.', info]
        if self.sample_values:
            columns.append(':'.join(self.format_keys) if self.format_keys else '.')
            columns.extend(':'.join(values) if values else '.' for values in self.sample_values)
        return '\t'.join(columns)

    def sort_key(self):
        return self.pos, self.id, self.text()


class VcfLineBuilder:
    """Builds one VCF line from one VariantRecord."""

    def __init__(self, genome: Genome, key_table: KeyTable = None, alleles: SymbolicAlleles = None,
                 samples: Sequence[str] = ()):
        self.genome = genome
        self.key_table = key_table or KeyTable.load()
        self.alleles = alleles or SymbolicAlleles.load()
        self.samples = list(samples)
        if len(set(self.samples)) != len(self.samples):
            raise ValueError('Sample names must be unique')

    def build(self, record: VariantRecord) -> VcfLine:
        try:
            contig = self.genome.contig(record.chromosome)
        except ReferenceLookupError as error:
            raise RecordError('unknown_chromosome', str(error)) from error
        allele_entry = self.alleles.for_term(record.so_term)
        if allele_entry is None:
            raise RecordError('unknown_type', f'{record.accession}: type {record.so_term!r} is not in the allele table')
        self._check_coordinates(record, contig, allele_entry.end_is_pos)

        info: Dict[str, object] = {}
        ci_pos, ci_end = record.confidence_intervals()
        for name, interval in (('CIPOS', ci_pos), ('CIEND', ci_end)):
            if interval is not None and ((interval.lower is not None and interval.lower > 0)
                                         or (interval.upper is not None and interval.upper < 0)):
                raise RecordError('invalid_confidence_interval',
                                  f'{record.accession}: {name} {interval.vcf_values()} does not contain 0')
        if record.sequence:
            if (ci_pos and not ci_pos.is_zero()) or (ci_end and not ci_end.is_zero()):
                raise RecordError('imprecise_with_sequence',
                                  f'{record.accession}: a sequence-resolved allele cannot be imprecise')
            pos, ref, alt = self._sequence_alleles(record, contig, allele_entry.end_is_pos)
            symbolic = None
            is_tandem_repeat = False
        else:
            if allele_entry.allele is None:
                raise RecordError('no_symbolic_allele',
                                  f'{record.accession}: type {record.so_term!r} has no symbolic allele')
            symbolic = allele_entry.allele
            pos = max(record.start - 1, 1)
            ref = normalise_reference(self.genome.fetch(contig, pos, pos))
            alt = f'<{symbolic}>'
            if allele_entry.end_is_pos:
                info['END'] = pos
                # VCF 4.4 requires SVLEN on symbolic alleles; the length of an insertion is unknown without its
                # sequence.
                info['SVLEN'] = [None]
            else:
                info['END'] = record.end
                info['SVLEN'] = [record.end - record.start + 1]
            imprecise = False
            if ci_pos is not None:
                info['CIPOS'] = ci_pos.vcf_values()
                imprecise |= not ci_pos.is_zero()
            if ci_end is not None:
                info['CIEND'] = ci_end.vcf_values()
                imprecise |= not ci_end.is_zero()
            if imprecise:
                info['IMPRECISE'] = True
            if allele_entry.needs_svclaim and 'SVCLAIM' not in record.info:
                if not allele_entry.svclaim:
                    raise RecordError('missing_svclaim',
                                      f'{record.accession}: <{symbolic}> needs SVCLAIM, and type {record.so_term!r} '
                                      f'does not imply one')
                info['SVCLAIM'] = [allele_entry.svclaim]
            is_tandem_repeat = symbolic == TANDEM_REPEAT_ALLELE
            if is_tandem_repeat and not any(key in record.info for key in ('RUS', 'RUL')):
                raise RecordError('missing_repeat_unit',
                                  f'{record.accession}: <{symbolic}> needs RUS or RUL (VCF 4.4)')
        repeat_keys = sorted(set(record.info) & TANDEM_REPEAT_KEYS)
        if repeat_keys and not is_tandem_repeat:
            raise RecordError('repeat_key_without_tandem_repeat',
                              f'{record.accession}: {", ".join(repeat_keys)} can only be written on '
                              f'<{TANDEM_REPEAT_ALLELE}> records (VCF 4.4), not on {alt}')

        if record.so_id:
            info['TYPESOID'] = record.so_id
        if record.parents:
            info['PARENT'] = sorted(set(record.parents))
        if record.placement_method and record.placement_method != DEFAULT_PLACEMENT_METHOD:
            info['PLACEMETHOD'] = record.placement_method
        for key, value in record.info.items():
            if key in info:
                raise RecordError('duplicate_key', f'{record.accession}: INFO {key} is set by the writer itself')
            info[key] = value
        if not self.samples and 'AF' not in info:
            # A VCF with no sample columns needs AF on every line for eva-sub-cli ('allele_frequency' evidence).
            # DGVa has no frequency for most variants: the value is written as missing ('.'), never invented.
            info['AF'] = [None]

        try:
            info_fields = self._format_info(info)
            format_keys, sample_values = self._format_samples(record)
        except UnmappedValueError as error:
            raise RecordError('unmapped_value', f'{record.accession}: {error}') from error
        return VcfLine(contig, pos, record.accession, ref, alt, info_fields, format_keys, sample_values, symbolic)

    def _check_coordinates(self, record, contig, end_is_pos):
        if record.start < 1:
            raise RecordError('invalid_coordinates', f'{record.accession}: start {record.start} is below 1')
        if not end_is_pos and record.end < record.start:
            raise RecordError('invalid_coordinates',
                              f'{record.accession}: end {record.end} is before start {record.start}')
        if max(record.start, record.end) > contig.length:
            raise RecordError('outside_contig', f'{record.accession}: {record.start}-{record.end} is beyond the end '
                                                f'of {contig.vcf_name} ({contig.length})')

    def _sequence_alleles(self, record, contig, is_insertion):
        sequence = record.sequence.upper()
        if not set(sequence) <= VALID_BASES:
            raise RecordError('sequence_not_allele', f'{record.accession}: the sequence is not made of A, C, G, T, N')
        ref_segment = '' if is_insertion else self.genome.fetch(contig, record.start, record.end)
        if record.start > 1:
            pos = record.start - 1
            padding = self.genome.fetch(contig, pos, pos)
            ref, alt = padding + ref_segment, padding + sequence
        else:
            # Nothing precedes position 1: pad with the base that follows the variant instead.
            following = record.start if is_insertion else record.end + 1
            if following > contig.length:
                raise RecordError('outside_contig', f'{record.accession}: no base to pad the allele with')
            pos = 1
            padding = self.genome.fetch(contig, following, following)
            ref, alt = ref_segment + padding, sequence + padding
        ref = normalise_reference(ref)
        if ref == alt:
            raise RecordError('alt_equals_ref', f'{record.accession}: the sequence equals the reference')
        return pos, ref, alt

    def _format_info(self, info):
        keys = sorted((self.key_table.get('INFO', key) for key in info), key=lambda vcf_key: vcf_key.order)
        return [(vcf_key.key, vcf_key.format_value(info[vcf_key.key])) for vcf_key in keys]

    def _format_samples(self, record):
        unknown = set(record.samples) - set(self.samples)
        if unknown:
            raise UnmappedValueError(f'sample(s) {sorted(unknown)} are not in the VCF header')
        if not self.samples:
            return [], []
        used = {key for values in record.samples.values() for key in values}
        if not used:
            # A record with no sample value still needs a FORMAT key: a missing genotype ('.') is written, which
            # states nothing about the samples (no genotype is inferred).
            used = {'GT'}
        format_keys = sorted((self.key_table.get('FORMAT', key) for key in used),
                             key=lambda vcf_key: (vcf_key.key != 'GT', vcf_key.order))
        sample_values = []
        for sample in self.samples:
            values = record.samples.get(sample, {})
            sample_values.append([vcf_key.format_value(values[vcf_key.key]) if vcf_key.key in values else '.'
                                  for vcf_key in format_keys])
        return [vcf_key.key for vcf_key in format_keys], sample_values


@dataclass
class Statistics:
    written_by_kind: Counter = field(default_factory=Counter)
    written_by_type: Counter = field(default_factory=Counter)
    written_by_chromosome: Counter = field(default_factory=Counter)
    refused: Counter = field(default_factory=Counter)
    refused_examples: Dict[str, List[str]] = field(default_factory=lambda: defaultdict(list))
    max_examples: int = 10

    @property
    def written(self):
        return sum(self.written_by_kind.values())

    def add_written(self, record, line):
        self.written_by_kind[record.kind.value] += 1
        self.written_by_type[record.so_term] += 1
        self.written_by_chromosome[line.contig.vcf_name] += 1

    def add_refused(self, record, error):
        self.refused[error.reason] += 1
        if len(self.refused_examples[error.reason]) < self.max_examples:
            self.refused_examples[error.reason].append(str(error))

    def rows(self):
        """(category, value, count) rows, sorted, for a TSV report."""
        rows = []
        for category, counter in (('kind', self.written_by_kind), ('type', self.written_by_type),
                                  ('chromosome', self.written_by_chromosome), ('refused', self.refused)):
            rows.extend((category, value, count) for value, count in sorted(counter.items()))
        return rows

    def write_tsv(self, path):
        with open(path, 'w') as handle:
            handle.write('category\tvalue\tcount\n')
            for category, value, count in self.rows():
                handle.write(f'{category}\t{value}\t{count}\n')
            for reason in sorted(self.refused_examples):
                for example in self.refused_examples[reason]:
                    handle.write(f'refused_example\t{reason}\t{example}\n')


@dataclass
class _ContigPart:
    path: str
    order: int
    last_key: Optional[tuple] = None


class StudyVcfWriter:
    """Writes one VCF (one study and assembly) from a stream of records. Use as a context manager."""

    def __init__(self, output_path, builder: VcfLineBuilder, reference_name=None,
                 meta_lines: Sequence[Tuple[str, str]] = (), source='dgva-extraction'):
        self.output_path = output_path
        self.builder = builder
        self.reference_name = reference_name
        self.meta_lines = list(meta_lines)
        self.source = source
        self.statistics = Statistics()
        self._directory = tempfile.mkdtemp(prefix='.dgva_extraction_', dir=os.path.dirname(output_path) or '.')
        self._parts: Dict[str, _ContigPart] = {}
        self._contigs: Dict[str, Contig] = {}
        self._current: Optional[str] = None
        self._handle = None
        self._buffer: List[Tuple[tuple, str]] = []
        self._used_alleles = set()
        self._used_info = set()
        self._used_format = set()
        self._finished = False

    def write(self, record: VariantRecord):
        try:
            line = self.builder.build(record)
        except RecordError as error:
            self.statistics.add_refused(record, error)
            return None
        self._add(line)
        self.statistics.add_written(record, line)
        if line.symbolic_allele:
            self._used_alleles.add(line.symbolic_allele)
        self._used_info.update(key for key, _ in line.info)
        self._used_format.update(line.format_keys)
        return line

    def refuse(self, record: VariantRecord, reason, message):
        """Record a refusal decided outside the writer (for example by the reader)."""
        self.statistics.add_refused(record, RecordError(reason, message))

    def _add(self, line: VcfLine):
        name = line.contig.vcf_name
        if name != self._current:
            self._switch_contig(line.contig)
        part = self._parts[name]
        key = line.sort_key()
        if part.last_key is not None and key <= part.last_key:
            raise OrderError(f'{line.id} at {name}:{line.pos} arrived after a record that sorts after it')
        if self._buffer and self._buffer[0][0][0] != line.pos:
            if line.pos < self._buffer[0][0][0]:
                raise OrderError(f'{line.id} at {name}:{line.pos} arrived after position {self._buffer[0][0][0]}')
            self._flush()
        self._buffer.append((key, line.text()))

    def _switch_contig(self, contig: Contig):
        self._flush()
        if self._handle is not None:
            self._handle.close()
        name = contig.vcf_name
        if name not in self._parts:
            self._parts[name] = _ContigPart(os.path.join(self._directory, f'{len(self._parts):08d}.part'),
                                            contig.order)
            self._contigs[name] = contig
        self._handle = open(self._parts[name].path, 'a')
        self._current = name

    def _flush(self):
        if not self._buffer:
            return
        self._buffer.sort()
        part = self._parts[self._current]
        for key, text in self._buffer:
            if part.last_key is not None and key <= part.last_key:
                raise OrderError(f'Duplicate or out-of-order record {key[1]} at {self._current}:{key[0]}')
            self._handle.write(text + '\n')
            part.last_key = key
        self._buffer = []

    def header_lines(self):
        lines = [f'##fileformat={VCF_VERSION}', f'##source={self.source}']
        if self.reference_name:
            lines.append(f'##reference={self.reference_name}')
        for name in sorted(self._contigs, key=lambda contig_name: (self._contigs[contig_name].order, contig_name)):
            lines.append(f'##contig=<ID={name},length={self._contigs[name].length}>')
        lines.extend(self.builder.alleles.header_lines(self._used_alleles))
        for field_name, used in (('INFO', self._used_info), ('FORMAT', self._used_format)):
            lines.extend(vcf_key.header_line() for vcf_key in self.builder.key_table.keys(field_name)
                         if vcf_key.key in used)
        lines.extend(f'##{key}={value}' for key, value in self.meta_lines)
        columns = HEADER_COLUMNS + (['FORMAT'] + self.builder.samples if self.builder.samples else [])
        lines.append('\t'.join(columns))
        return lines

    def finish(self, allow_refused=False):
        """Write the VCF. Raises StudyFailedError, and writes nothing, if a record was refused."""
        if self._finished:
            return self.statistics
        self._flush()
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        try:
            if self.statistics.refused and not allow_refused:
                raise StudyFailedError(self.statistics)
            temporary_output = os.path.join(self._directory, 'output.vcf')
            with open(temporary_output, 'w') as output:
                output.write('\n'.join(self.header_lines()) + '\n')
                for name in sorted(self._parts, key=lambda part_name: (self._parts[part_name].order, part_name)):
                    with open(self._parts[name].path) as part:
                        shutil.copyfileobj(part, output)
            os.replace(temporary_output, self.output_path)
        finally:
            self._finished = True
            shutil.rmtree(self._directory, ignore_errors=True)
        return self.statistics

    def abort(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        self._finished = True
        shutil.rmtree(self._directory, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None or not self._finished:
            self.abort()
        return False
