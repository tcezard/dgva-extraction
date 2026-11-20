"""A call's genotype (FORMAT GT) from its zygosity and the species' ploidy.

DGVa has no genotype table: a call belongs to one sample, and VARIANT_CALL_ZYGOSITY may say whether the sample is
heterozygous, homozygous or hemizygous for it. The ploidy comes from etc/ploidy.yaml:

- Heterozygous: one copy carries the variant (0/1 when diploid);
- Homozygous: every copy carries it (1/1 when diploid);
- Hemizygous: the only copy carries it (1);
- no zygosity: one copy carries it and the others are unknown (1/. when diploid; 1 when haploid).
"""
from importlib import resources

import yaml

ZYGOSITIES = ('Heterozygous', 'Homozygous', 'Hemizygous')


class GenotypeError(ValueError):
    pass


class PloidyMap:
    def __init__(self, mapping):
        self._species = {int(taxonomy): entry for taxonomy, entry in mapping.items()}

    @classmethod
    def load(cls, text=None):
        if text is None:
            text = resources.files('dgva_extraction').joinpath('etc', 'ploidy.yaml').read_text()
        return cls(yaml.safe_load(text))

    def ploidy(self, taxonomy_id, chromosome):
        entry = self._species.get(int(taxonomy_id))
        if entry is None:
            raise GenotypeError(f'Taxonomy {taxonomy_id} is not in etc/ploidy.yaml')
        return int((entry.get('chromosomes') or {}).get(chromosome, entry['ploidy']))


def genotype(zygosity, ploidy):
    """GT for one call of one sample."""
    if ploidy < 1:
        raise GenotypeError(f'Invalid ploidy {ploidy}')
    if zygosity == 'Hemizygous':
        return '1'
    if zygosity == 'Homozygous':
        return '/'.join(['1'] * ploidy)
    if zygosity == 'Heterozygous':
        if ploidy < 2:
            raise GenotypeError('Heterozygous call on a haploid sequence')
        return '/'.join(['0'] * (ploidy - 1) + ['1'])
    if zygosity is None:
        return '/'.join(['1'] + ['.'] * (ploidy - 1))
    raise GenotypeError(f'Unknown zygosity {zygosity!r} (expected one of {", ".join(ZYGOSITIES)})')
