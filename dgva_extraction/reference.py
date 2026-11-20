"""Reference sequence and contig names: an indexed FASTA (samtools .fai) and the NCBI assembly report."""
import os
from dataclasses import dataclass
from typing import Dict, List, Optional


class ReferenceLookupError(ValueError):
    """A contig or a base cannot be found in the reference."""


@dataclass(frozen=True)
class FaiEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


def build_fai(fasta_path):
    """Write <fasta>.fai (samtools format) for an uncompressed FASTA, and return its path."""
    entries = []
    with open(fasta_path, 'rb') as handle:
        name = None
        length = offset = line_bases = line_width = 0
        position = 0
        for line in handle:
            if line.startswith(b'>'):
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                name = line[1:].split()[0].decode()
                length = line_bases = line_width = 0
                offset = position + len(line)
            elif line.strip():
                bases = len(line.rstrip(b'\r\n'))
                if line_bases == 0:
                    line_bases, line_width = bases, len(line)
                length += bases
            position += len(line)
        if name is not None:
            entries.append((name, length, offset, line_bases, line_width))
    fai_path = fasta_path + '.fai'
    with open(fai_path, 'w') as fai:
        for entry in entries:
            fai.write('\t'.join(str(value) for value in entry) + '\n')
    return fai_path


class FastaReference:
    """Random access to an uncompressed FASTA through its .fai index (built if missing)."""

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        fai_path = fasta_path + '.fai'
        if not os.path.exists(fai_path):
            build_fai(fasta_path)
        self.entries: Dict[str, FaiEntry] = {}
        self.names: List[str] = []
        with open(fai_path) as fai:
            for line in fai:
                name, length, offset, line_bases, line_width = line.rstrip('\n').split('\t')[:5]
                self.entries[name] = FaiEntry(name, int(length), int(offset), int(line_bases), int(line_width))
                self.names.append(name)
        self._handle = open(fasta_path, 'rb')

    def __contains__(self, name):
        return name in self.entries

    def length(self, name):
        return self.entries[name].length

    def fetch(self, name, start, end):
        """Bases from start to end, 1-based and inclusive, upper case."""
        entry = self.entries.get(name)
        if entry is None:
            raise ReferenceLookupError(f'Sequence {name!r} is not in {self.fasta_path}')
        if start < 1 or end > entry.length or end < start:
            raise ReferenceLookupError(f'{name}:{start}-{end} is outside the sequence (length {entry.length})')
        first, last = start - 1, end - 1
        byte_start = entry.offset + (first // entry.line_bases) * entry.line_width + first % entry.line_bases
        byte_end = entry.offset + (last // entry.line_bases) * entry.line_width + last % entry.line_bases + 1
        self._handle.seek(byte_start)
        raw = self._handle.read(byte_end - byte_start)
        return raw.replace(b'\n', b'').replace(b'\r', b'').decode().upper()

    def close(self):
        self._handle.close()


@dataclass(frozen=True)
class Contig:
    vcf_name: str
    length: int
    order: int


class Genome:
    """Resolves the chromosome names found in the database to the FASTA, in assembly-report order.

    The VCF CHROM is the name the sequence has in the FASTA, so that the REF check and the header agree. A database
    name is matched first as it is, then through the aliases on its assembly-report row (sequence name, GenBank and
    RefSeq accessions, UCSC name, assigned molecule).
    """

    REPORT_ALIAS_COLUMNS = ('Sequence-Name', 'GenBank-Accn', 'RefSeq-Accn', 'UCSC-style-name', 'Assigned-Molecule')

    def __init__(self, fasta, assembly_report_path=None, assembly_accession=None):
        self.fasta = fasta
        self.assembly_accession = assembly_accession
        self._report_rows: List[Dict[str, str]] = []
        if assembly_report_path:
            self._report_rows = parse_assembly_report(assembly_report_path)
        self._order: Dict[str, int] = {}
        self._aliases: Dict[str, List[str]] = {}
        if self._report_rows:
            for index, row in enumerate(self._report_rows):
                aliases = [row.get(column) for column in self.REPORT_ALIAS_COLUMNS
                           if column != 'Assigned-Molecule' or row.get('Sequence-Role') == 'assembled-molecule']
                aliases = [alias for alias in aliases if alias and alias != 'na']
                fasta_name = next((alias for alias in aliases if alias in fasta), None)
                for alias in aliases:
                    self._aliases.setdefault(alias, aliases)
                if fasta_name is not None:
                    self._order.setdefault(fasta_name, index)
            next_index = len(self._report_rows)
        else:
            next_index = 0
        # Sequences of the FASTA missing from the report come after the report's, in FASTA order.
        for name in fasta.names:
            if name not in self._order:
                self._order[name] = next_index
                next_index += 1
        self._cache: Dict[str, Optional[Contig]] = {}

    def contig(self, name):
        """The Contig for a database chromosome name; raises ReferenceLookupError if it can't be found."""
        if name not in self._cache:
            self._cache[name] = self._resolve(name)
        contig = self._cache[name]
        if contig is None:
            raise ReferenceLookupError(f'Chromosome {name!r} is not in the FASTA or the assembly report')
        return contig

    def _resolve(self, name):
        candidates = [name] + self._aliases.get(name, [])
        if not name.lower().startswith('chr'):
            candidates += [f'chr{name}'] + self._aliases.get(f'chr{name}', [])
        else:
            candidates += [name[3:]] + self._aliases.get(name[3:], [])
        for candidate in candidates:
            if candidate in self.fasta:
                return Contig(candidate, self.fasta.length(candidate), self._order[candidate])
        return None

    def fetch(self, contig, start, end):
        return self.fasta.fetch(contig.vcf_name, start, end)


def parse_assembly_report(path):
    """Rows of an NCBI assembly report, as dicts keyed by the column names of its last comment line."""
    header = None
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('#'):
                header = [column.strip() for column in line.lstrip('#').strip().split('\t')]
                continue
            if not line.strip() or header is None:
                continue
            rows.append(dict(zip(header, line.split('\t'))))
    return rows
