"""The variant record: one DGVa accession on one placement, independent of where it was read from."""
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional, Tuple


class RecordKind(str, Enum):
    REGION = 'region'
    CALL = 'call'


@dataclass(frozen=True)
class ConfidenceInterval:
    """Offsets from a position, VCF style: lower <= 0 <= upper. None is an unknown bound (written '.')."""
    lower: Optional[int]
    upper: Optional[int]

    def is_zero(self):
        return self.lower == 0 and self.upper == 0

    def vcf_values(self):
        return [self.lower, self.upper]


def normalise_so_term(term):
    """'copy number gain', 'copy_number_gain' and 'Copy-Number Gain' all become 'copy_number_gain'."""
    return re.sub(r'[\s\-]+', '_', term.strip().lower()) if term else term


@dataclass
class VariantRecord:
    """One VCF record: a DGVa accession on one placement.

    Coordinates are 1-based and inclusive: start is the first affected base and end the last one (the best estimate
    when the variant is imprecise). For an insertion, the sequence is inserted before `start`.
    The confidence intervals are offsets from start and end; when they are not given, they are derived from the
    inner and outer coordinates (see confidence_intervals()).
    """
    accession: str
    kind: RecordKind
    so_term: str
    chromosome: str
    start: int
    end: int
    so_id: Optional[str] = None
    outer_start: Optional[int] = None
    inner_start: Optional[int] = None
    inner_end: Optional[int] = None
    outer_end: Optional[int] = None
    ci_pos: Optional[ConfidenceInterval] = None
    ci_end: Optional[ConfidenceInterval] = None
    sequence: Optional[str] = None
    info: Dict[str, Any] = field(default_factory=dict)
    samples: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    parents: Tuple[str, ...] = ()
    placement_method: Optional[str] = None
    placement_id: Optional[int] = None

    def __post_init__(self):
        self.kind = RecordKind(self.kind)
        self.so_term = normalise_so_term(self.so_term)
        self.parents = tuple(self.parents)

    def confidence_intervals(self):
        """Return (CIPOS, CIEND), each a ConfidenceInterval or None.

        Explicit intervals win. Otherwise: CIPOS = (outer start - start, inner start - start) and
        CIEND = (inner end - end, outer end - end), with '.' for a coordinate that is missing.
        """
        ci_pos = self.ci_pos
        if ci_pos is None and (self.outer_start is not None or self.inner_start is not None):
            ci_pos = ConfidenceInterval(
                None if self.outer_start is None else self.outer_start - self.start,
                None if self.inner_start is None else self.inner_start - self.start)
        ci_end = self.ci_end
        if ci_end is None and (self.inner_end is not None or self.outer_end is not None):
            ci_end = ConfidenceInterval(
                None if self.inner_end is None else self.inner_end - self.end,
                None if self.outer_end is None else self.outer_end - self.end)
        return ci_pos, ci_end
