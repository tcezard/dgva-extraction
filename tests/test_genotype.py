import pytest

from dgva_extraction.genotype import GenotypeError, PloidyMap, genotype


@pytest.mark.parametrize('zygosity, ploidy, expected', [
    ('Heterozygous', 2, '0/1'), ('Homozygous', 2, '1/1'), ('Hemizygous', 2, '1'), (None, 2, '1/.'),
    ('Homozygous', 1, '1'), ('Hemizygous', 1, '1'), (None, 1, '1'),
])
def test_genotype_from_zygosity_and_ploidy(zygosity, ploidy, expected):
    assert genotype(zygosity, ploidy) == expected


@pytest.mark.parametrize('zygosity, ploidy', [('Heterozygous', 1), ('heterozygous', 2), ('Unknown', 2), (None, 0)])
def test_impossible_genotypes_are_refused(zygosity, ploidy):
    with pytest.raises(GenotypeError):
        genotype(zygosity, ploidy)


def test_ploidy_map():
    ploidy = PloidyMap.load()
    assert ploidy.ploidy(9606, '1') == 2
    assert ploidy.ploidy('9606', 'MT') == 1
    assert ploidy.ploidy(39947, 'Pt') == 1
    with pytest.raises(GenotypeError):
        ploidy.ploidy(4932, '1')
