"""The VCF key table and the symbolic alleles (ticket 10), loaded from dgva_extraction/etc."""
import math
import re
from dataclasses import dataclass
from importlib import resources
from typing import Dict, List, Optional, Tuple

import yaml

from dgva_extraction.record import normalise_so_term

VALID_TYPES = {'Integer', 'Float', 'Flag', 'Character', 'String'}
# Characters with a special meaning in INFO and FORMAT values (VCF 4.4, section 1.2: percent encoding). ':' is a
# delimiter in FORMAT values only.
_PERCENT_ENCODED = {'%': '%25', ':': '%3A', ';': '%3B', '=': '%3D', ',': '%2C', '\r': '%0D', '\n': '%0A',
                    '\t': '%09'}
# The symbolic SV allele types VCF 4.4 defines (section 5.3). Only these may have colon-separated subtypes, and only
# these have a length: other symbolic alleles are defined by the file that uses them, with SVLEN=.
RESERVED_SV_TYPES = ('DEL', 'INS', 'DUP', 'INV', 'CNV')
_ALLELE_ID = re.compile(r'[^\s,<>]+')
_PERCENT_PATTERNS = {
    'INFO': re.compile('[' + re.escape(''.join(c for c in _PERCENT_ENCODED if c != ':')) + ']'),
    'FORMAT': re.compile('[' + re.escape(''.join(_PERCENT_ENCODED)) + ']'),
}


class UnmappedValueError(ValueError):
    """A value has no key in the table, or doesn't match its key's Number or Type."""


def _etc(name):
    return resources.files('dgva_extraction').joinpath('etc', name).read_text()


def percent_encode(text, field='FORMAT'):
    return _PERCENT_PATTERNS[field].sub(lambda match: _PERCENT_ENCODED[match.group(0)], text)


@dataclass(frozen=True)
class VcfKey:
    key: str
    field: str
    number: str
    type: str
    description: str
    source: str
    reserved: bool
    status: str
    sprint: int
    order: int

    def header_line(self):
        description = self.description.replace('\\', '\\\\').replace('"', '\\"')
        return f'##{self.field}=<ID={self.key},Number={self.number},Type={self.type},Description="{description}">'

    def format_value(self, value, alt_count=1):
        """Check a typed value against Number and Type, and return its VCF text. Flags return None."""
        if self.type == 'Flag':
            if value is not True:
                raise UnmappedValueError(f'{self.field} {self.key} is a Flag: its value must be True, not {value!r}')
            return None
        values = list(value) if isinstance(value, (list, tuple)) else [value]
        expected = self._expected_count(alt_count)
        if expected is not None and len(values) != expected:
            raise UnmappedValueError(
                f'{self.field} {self.key} has Number={self.number}: expected {expected} value(s), got {values!r}')
        if not values:
            raise UnmappedValueError(f'{self.field} {self.key} has no value')
        return ','.join(self._format_one(item) for item in values)

    def _expected_count(self, alt_count):
        if self.number == 'A':
            return alt_count
        if self.number == 'R':
            return alt_count + 1
        if self.number in ('G', '.'):
            return None
        return int(self.number)

    def _format_one(self, item):
        if item is None:
            return '.'
        if self.type == 'Integer':
            if isinstance(item, bool) or not isinstance(item, int):
                raise UnmappedValueError(f'{self.field} {self.key} is an Integer, not {item!r}')
            return str(item)
        if self.type == 'Float':
            if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
                raise UnmappedValueError(f'{self.field} {self.key} is a Float, not {item!r}')
            return repr(float(item)) if isinstance(item, float) else str(item)
        if not isinstance(item, str) or item == '':
            raise UnmappedValueError(f'{self.field} {self.key} is a {self.type}, not {item!r}')
        if self.type == 'Character' and len(item) != 1:
            raise UnmappedValueError(f'{self.field} {self.key} is a Character, not {item!r}')
        return percent_encode(item, self.field)


class KeyTable:
    """INFO and FORMAT keys, in the fixed order of etc/vcf_keys.yaml."""

    def __init__(self, entries):
        self._keys: Dict[Tuple[str, str], VcfKey] = {}
        for order, entry in enumerate(entries):
            vcf_key = VcfKey(key=entry['key'], field=entry['field'], number=str(entry['number']),
                             type=entry['type'], description=entry['description'], source=entry.get('source', ''),
                             reserved=bool(entry.get('reserved')), status=entry.get('status', 'proposed'),
                             sprint=int(entry.get('sprint', 0)), order=order)
            if vcf_key.field not in ('INFO', 'FORMAT'):
                raise ValueError(f'{vcf_key.key}: field must be INFO or FORMAT')
            if vcf_key.type not in VALID_TYPES:
                raise ValueError(f'{vcf_key.key}: unknown Type {vcf_key.type}')
            if (vcf_key.field, vcf_key.key) in self._keys:
                raise ValueError(f'{vcf_key.field} {vcf_key.key} is listed twice')
            self._keys[(vcf_key.field, vcf_key.key)] = vcf_key

    @classmethod
    def load(cls, text=None):
        return cls(yaml.safe_load(text if text is not None else _etc('vcf_keys.yaml')))

    def get(self, field, key):
        vcf_key = self._keys.get((field, key))
        if vcf_key is None:
            raise UnmappedValueError(f'{field} key {key!r} is not in the key table: the value has no target')
        return vcf_key

    def keys(self, field=None):
        return [vcf_key for vcf_key in self._keys.values() if field is None or vcf_key.field == field]


@dataclass(frozen=True)
class SymbolicAllele:
    so_term: str
    allele: Optional[str]
    end_is_pos: bool
    description: str
    order: int
    svclaim: Optional[str] = None

    def __post_init__(self):
        if self.allele is None:
            return
        if not _ALLELE_ID.fullmatch(self.allele):
            raise ValueError(f'Symbolic allele {self.allele!r} ({self.so_term}): no spaces, commas or angle brackets')
        if ':' in self.allele and not self.is_reserved:
            raise ValueError(f'Symbolic allele {self.allele!r} ({self.so_term}): only the VCF types '
                             f'{", ".join(RESERVED_SV_TYPES)} can have subtypes')

    @property
    def is_reserved(self):
        """True for the SV alleles VCF 4.4 defines (subtypes included); False for ones defined by the file."""
        return bool(self.allele) and self.allele.split(':')[0] in RESERVED_SV_TYPES

    @property
    def needs_svclaim(self):
        """VCF 4.4 requires SVCLAIM on <DEL> and <DUP> alleles, subtypes included."""
        return bool(self.allele) and self.allele.split(':')[0] in ('DEL', 'DUP')

    def header_line(self):
        return f'##ALT=<ID={self.allele},Description="{self.description}">'


class SymbolicAlleles:
    def __init__(self, mapping):
        self._by_term: Dict[str, SymbolicAllele] = {}
        self._allele_order: Dict[str, int] = {}
        for order, (term, entry) in enumerate(mapping.items()):
            allele = SymbolicAllele(normalise_so_term(term), entry.get('allele'), bool(entry.get('end_is_pos')),
                                    entry.get('description', ''), order, entry.get('svclaim'))
            self._by_term[allele.so_term] = allele
            if allele.allele:
                self._allele_order.setdefault(allele.allele, order)

    @classmethod
    def load(cls, text=None):
        return cls(yaml.safe_load(text if text is not None else _etc('symbolic_alleles.yaml')))

    def for_term(self, so_term):
        """The SymbolicAllele for a type, or None if the type is not in the table at all."""
        return self._by_term.get(normalise_so_term(so_term))

    def header_lines(self, used_alleles):
        """##ALT lines for the alleles used, in table order (the first type that maps to each allele names it)."""
        lines: List[str] = []
        for allele in sorted(used_alleles, key=lambda name: self._allele_order[name]):
            first = next(entry for entry in self._by_term.values() if entry.allele == allele)
            lines.append(first.header_line())
        return lines
