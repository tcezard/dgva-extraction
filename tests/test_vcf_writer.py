"""Ticket 12: the writer, from records only (no database).

The edge cases are the checklist of convertGVFtoVCF's tested behaviour (reuse-or-restart review, finding 8): REF lookup
and padding at position 1 and elsewhere, IUPAC codes in the reference, CI bounds including '.', END for each symbolic
allele, symbolic alleles from SO terms, and FORMAT ordering with GT first.
"""
import filecmp
import random
from importlib import resources

import pytest
import yaml

from conftest import CHR1, vcf_validator_errors
from dgva_extraction.keys import SymbolicAlleles
from dgva_extraction.record import Breakend, BreakendSide, ConfidenceInterval, VariantRecord
from dgva_extraction.vcf_writer import (OrderError, RecordError, StudyFailedError, StudyVcfWriter, VcfLineBuilder,
                                        normalise_reference)


def base(sequence, position):
    return sequence[position - 1].upper()


def call(accession, so_term, chromosome, start, end, **kwargs):
    return VariantRecord(accession, 'call', so_term, chromosome, start, end, **kwargs)


def region(accession, so_term, chromosome, start, end, **kwargs):
    return VariantRecord(accession, 'region', so_term, chromosome, start, end, **kwargs)


def data_lines(path):
    with open(path) as handle:
        return [line.rstrip('\n').split('\t') for line in handle if not line.startswith('#')]


def info_dict(info_column):
    return dict(item.split('=', 1) if '=' in item else (item, True) for item in info_column.split(';'))


@pytest.fixture
def builder(genome):
    return VcfLineBuilder(genome)


# --- REF lookup and padding ---

def test_symbolic_allele_is_padded_with_the_preceding_base(builder):
    line = builder.build(call('essv1', 'copy_number_loss', '1', 21, 40))
    assert (line.pos, line.ref, line.alt) == (20, base(CHR1, 20), '<DEL>')
    info = dict(line.info)
    assert info['END'] == '40'
    assert info['SVLEN'] == '20'


def test_symbolic_allele_at_position_1_keeps_pos_1(builder):
    line = builder.build(call('essv1', 'copy_number_loss', '1', 1, 10))
    assert (line.pos, line.ref) == (1, base(CHR1, 1))
    assert dict(line.info)['END'] == '10'


def test_symbolic_allele_at_position_2_is_padded_with_position_1(builder):
    line = builder.build(call('essv1', 'copy_number_loss', '1', 2, 10))
    assert (line.pos, line.ref) == (1, base(CHR1, 1))


def test_sequence_insertion_is_padded_with_the_preceding_base(builder):
    line = builder.build(call('essv1', 'insertion', '1', 30, 30, sequence='ttgac'))
    assert (line.pos, line.ref, line.alt) == (29, base(CHR1, 29), base(CHR1, 29) + 'TTGAC')
    assert dict(line.info).get('END') is None


def test_sequence_insertion_at_position_1_is_padded_with_the_following_base(builder):
    line = builder.build(call('essv1', 'insertion', '1', 1, 1, sequence='GG'))
    assert (line.pos, line.ref, line.alt) == (1, base(CHR1, 1), 'GG' + base(CHR1, 1))


def test_sequence_replacement_covers_the_reference_segment(builder):
    line = builder.build(call('essv1', 'deletion', '1', 11, 14, sequence='A'))
    assert line.pos == 10
    assert line.ref == CHR1[9:14].upper()
    assert line.alt == base(CHR1, 10) + 'A'


def test_sequence_replacement_at_position_1_pads_after_the_segment(builder):
    line = builder.build(call('essv1', 'deletion', '1', 1, 4, sequence='T'))
    assert line.pos == 1
    assert line.ref == CHR1[0:4] + base(CHR1, 5)
    assert line.alt == 'T' + base(CHR1, 5)


def test_lower_case_reference_is_upper_cased(builder):
    line = builder.build(call('essv1', 'copy_number_loss', '1', 122, 125))
    assert line.ref == base(CHR1, 121) == 'A'


# --- IUPAC codes and REF checks ---

def test_iupac_reference_base_is_replaced(builder):
    assert CHR1[60] == 'R'
    line = builder.build(call('essv1', 'copy_number_gain', '1', 62, 70))
    assert (line.pos, line.ref) == (61, 'A')


@pytest.mark.parametrize('bases, expected', [('ACGTN', 'ACGTN'), ('acgt', 'ACGT'), ('RYSWKM', 'ACCAGA'),
                                             ('BDHV', 'CAAA')])
def test_normalise_reference(bases, expected):
    assert normalise_reference(bases) == expected


def test_reference_with_an_invalid_character_is_refused():
    with pytest.raises(RecordError) as error:
        normalise_reference('AC*T')
    assert error.value.reason == 'invalid_reference_base'


def test_sequence_equal_to_reference_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'deletion', '1', 11, 14, sequence=CHR1[10:14]))
    assert error.value.reason == 'alt_equals_ref'


def test_sequence_that_is_not_an_allele_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'insertion', '1', 30, 30, sequence='ACXT'))
    assert error.value.reason == 'sequence_not_allele'


# --- Confidence intervals ---

def test_precise_symbolic_variant_has_no_ci(builder):
    info = dict(builder.build(call('essv1', 'copy_number_loss', '1', 21, 40)).info)
    assert 'IMPRECISE' not in info and 'CIPOS' not in info and 'CIEND' not in info


def test_explicit_cis_make_the_record_imprecise(builder):
    record = call('essv1', 'copy_number_loss', '1', 21, 40, ci_pos=ConfidenceInterval(-7, 8),
                  ci_end=ConfidenceInterval(-3, 4))
    info = dict(builder.build(record).info)
    assert info['CIPOS'] == '-7,8'
    assert info['CIEND'] == '-3,4'
    assert 'IMPRECISE' in info and info['IMPRECISE'] is None


def test_cis_derived_from_inner_and_outer_coordinates(builder):
    record = call('essv1', 'copy_number_gain', '2', 100, 150, outer_start=95, inner_start=103, inner_end=148,
                  outer_end=160)
    info = dict(builder.build(record).info)
    assert info['CIPOS'] == '-5,3'
    assert info['CIEND'] == '-2,10'


def test_missing_ci_bounds_are_written_as_dots(builder):
    record = call('essv1', 'copy_number_gain', '2', 100, 150, outer_start=95, outer_end=160)
    info = dict(builder.build(record).info)
    assert info['CIPOS'] == '-5,.'
    assert info['CIEND'] == '.,10'
    assert 'IMPRECISE' in info


def test_zero_cis_are_written_but_precise(builder):
    record = call('essv1', 'copy_number_gain', '2', 100, 150, ci_pos=ConfidenceInterval(0, 0))
    info = dict(builder.build(record).info)
    assert info['CIPOS'] == '0,0'
    assert 'IMPRECISE' not in info


def test_ci_that_does_not_contain_zero_is_refused(builder):
    record = call('essv1', 'copy_number_gain', '2', 100, 150, ci_pos=ConfidenceInterval(2, 5))
    with pytest.raises(RecordError) as error:
        builder.build(record)
    assert error.value.reason == 'invalid_confidence_interval'


def test_imprecise_sequence_is_refused(builder):
    record = call('essv1', 'insertion', '1', 30, 30, sequence='AC', ci_pos=ConfidenceInterval(-1, 1))
    with pytest.raises(RecordError) as error:
        builder.build(record)
    assert error.value.reason == 'imprecise_with_sequence'


# --- Symbolic alleles from SO terms, and END for each ---

@pytest.mark.parametrize('so_term, allele, end, svlen', [
    ('deletion', '<DEL>', '40', '20'),
    ('copy number loss', '<DEL>', '40', '20'),
    ('alu_deletion', '<DEL:ME:ALU>', '40', '20'),
    ('copy_number_gain', '<DUP>', '40', '20'),
    ('tandem duplication', '<DUP:TANDEM>', '40', '20'),
    ('copy_number_variation', '<CNV>', '40', '20'),
    ('inversion', '<INV>', '40', '20'),
    ('insertion', '<INS>', '20', '.'),
    ('mobile-element insertion', '<INS:ME>', '20', '.'),
    ('line1_insertion', '<INS:ME:LINE1>', '20', '.'),
    ('novel_sequence_insertion', '<INS>', '20', '.'),
    # No VCF 4.4 allele: defined by the file, with no length.
    ('interchromosomal_translocation', '<TRA>', '40', '.'),
    ('translocation', '<TRA>', '40', '.'),
    ('complex_chromosomal_rearrangement', '<COMPLEX_CHROMOSOMAL_REARRANGEMENT>', '40', '.'),
    ('complex_substitution', '<COMPLEX_SUBSTITUTION>', '40', '.'),
    ('sequence_alteration', '<SEQUENCE_ALTERATION>', '40', '.'),
    ('delins', '<DELINS>', '40', '.'),
    ('indel', '<INDEL>', '40', '.'),
    ('copy_number_neutral_loss_of_heterozygosity', '<COPY_NUMBER_NEUTRAL_LOSS_OF_HETEROZYGOSITY>', '40', '.'),
    ('loss_of_heterozygosity', '<LOSS_OF_HETEROZYGOSITY>', '40', '.'),
])
def test_symbolic_allele_and_end_per_type(builder, so_term, allele, end, svlen):
    line = builder.build(call('essv1', so_term, '1', 21, 40))
    info = dict(line.info)
    assert line.alt == allele
    assert info['END'] == end
    assert info['SVLEN'] == svlen


def test_every_type_in_the_allele_table_has_an_allele():
    table = yaml.safe_load(resources.files('dgva_extraction').joinpath('etc', 'symbolic_alleles.yaml').read_text())
    assert [term for term, entry in table.items() if not entry.get('allele')] == []


def test_type_with_no_symbolic_allele_is_refused(genome):
    alleles = SymbolicAlleles.load('undecided_type: {allele: null, end_is_pos: false, description: Not decided}')
    with pytest.raises(RecordError) as error:
        VcfLineBuilder(genome, alleles=alleles).build(call('essv1', 'undecided_type', '1', 21, 40))
    assert error.value.reason == 'no_symbolic_allele'


@pytest.mark.parametrize('allele', ['LOH:CN', 'TWO WORDS', 'A,B', '<CPX>'])
def test_allele_the_file_defines_must_be_a_valid_id(allele):
    with pytest.raises(ValueError):
        SymbolicAlleles.load(f'some_type: {{allele: "{allele}", end_is_pos: false, description: x}}')


def test_type_with_a_sequence_is_written_from_it_whatever_its_allele(builder):
    line = builder.build(call('essv1', 'delins', '1', 21, 22, sequence='TTT'))
    assert not line.alt.startswith('<')


def test_type_not_in_the_allele_table_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'made_up_type', '1', 21, 40))
    assert error.value.reason == 'unknown_type'


@pytest.mark.parametrize('so_term', ['deletion', 'duplication', 'tandem_duplication', 'alu_deletion'])
def test_deletion_or_duplication_without_svclaim_claims_abundance(builder, so_term):
    assert dict(builder.build(call('essv1', so_term, '1', 21, 40)).info)['SVCLAIM'] == 'D'


def test_submitted_svclaim_is_kept(builder):
    assert dict(builder.build(call('essv1', 'deletion', '1', 21, 40, info={'SVCLAIM': ['J']})).info)['SVCLAIM'] == 'J'


def test_copy_number_change_implies_an_abundance_claim(builder):
    assert dict(builder.build(call('essv1', 'copy_number_loss', '1', 21, 40)).info)['SVCLAIM'] == 'D'


def test_tandem_repeat_needs_a_repeat_unit(builder):
    with pytest.raises(RecordError) as error:
        builder.build(region('esv1', 'short_tandem_repeat_variation', '1', 21, 40))
    assert error.value.reason == 'missing_repeat_unit'
    line = builder.build(region('esv1', 'short_tandem_repeat_variation', '1', 21, 40,
                                info={'RUS': ['AC'], 'RUC': [10.0]}))
    assert line.alt == '<CNV:TR>'


def test_repeat_keys_only_on_tandem_repeats(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_gain', '1', 21, 40, info={'RUC': [2.0]}))
    assert error.value.reason == 'repeat_key_without_tandem_repeat'


# --- Coordinates and contigs ---

def test_unknown_chromosome_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '7', 21, 40))
    assert error.value.reason == 'unknown_chromosome'


@pytest.mark.parametrize('name', ['1', 'chr1', 'CM000001.1', 'NC_000001.1'])
def test_chromosome_aliases_resolve_to_the_fasta_name(builder, name):
    assert builder.build(call('essv1', 'copy_number_loss', name, 21, 40)).contig.vcf_name == '1'


def test_end_beyond_the_contig_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '1', 21, len(CHR1) + 1))
    assert error.value.reason == 'outside_contig'


def test_end_before_start_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '1', 21, 20))
    assert error.value.reason == 'invalid_coordinates'


# --- INFO keys, parents, placement method ---

def test_parents_type_and_placement_method(builder):
    record = call('nssv2', 'copy_number_loss', '1', 21, 40, so_id='SO:0001743', parents=('nsv9', 'nsv1', 'nsv9'),
                  placement_method='Submitted cytogenetic')
    info = dict(builder.build(record).info)
    assert info['PARENT_REGION_ACCESSION'] == 'nsv1,nsv9'
    assert info['VC'] == 'SO:0001743'
    assert info['PLACEMENT_METHOD'] == 'Submitted cytogenetic'


def test_default_placement_method_is_not_written(builder):
    record = call('nssv2', 'copy_number_loss', '1', 21, 40, placement_method='Submitted genomic')
    assert 'PLACEMENT_METHOD' not in dict(builder.build(record).info)


def test_info_keys_follow_the_key_table_order(builder):
    record = call('essv1', 'copy_number_loss', '1', 21, 40, info={'HGVS': ['x'], 'AF': [0.25]},
                  ci_pos=ConfidenceInterval(-1, 1), parents=('esv1',), so_id='SO:0001743')
    keys = [key for key, _ in builder.build(record).info]
    assert keys == ['END', 'SVLEN', 'IMPRECISE', 'CIPOS', 'SVCLAIM', 'VC', 'PARENT_REGION_ACCESSION', 'AF', 'HGVS']


def test_value_with_no_key_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '1', 21, 40, info={'NOT_A_KEY': ['x']}))
    assert error.value.reason == 'unmapped_value'


def test_value_with_the_wrong_type_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '1', 21, 40, info={'AF': ['high']}))
    assert error.value.reason == 'unmapped_value'


def test_writer_keys_cannot_be_overridden(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('essv1', 'copy_number_loss', '1', 21, 40, info={'END': 50}))
    assert error.value.reason == 'duplicate_key'


def test_info_strings_are_percent_encoded(builder):
    record = call('essv1', 'copy_number_loss', '1', 21, 40, info={'HGVS': ['NC_1:g.1_2del;x=y,z']})
    assert dict(builder.build(record).info)['HGVS'] == 'NC_1:g.1_2del%3Bx%3Dy%2Cz'


# --- FORMAT ---

def test_format_has_gt_first_and_samples_in_header_order(genome):
    builder = VcfLineBuilder(genome, samples=['S1', 'S2', 'S3'])
    record = call('essv1', 'copy_number_loss', '1', 21, 40,
                  samples={'S3': {'CN': 1.0, 'GT': '0/1'}, 'S1': {'CN': 0.0}})
    line = builder.build(record)
    assert line.format_keys == ['GT', 'CN']
    assert line.sample_values == [['.', '0.0'], ['.', '.'], ['0/1', '1.0']]


def test_record_without_sample_values_gets_a_missing_genotype(genome):
    builder = VcfLineBuilder(genome, samples=['S1'])
    line = builder.build(call('essv1', 'copy_number_loss', '1', 21, 40))
    assert line.format_keys == ['GT']
    assert line.sample_values == [['.']]


def test_unknown_sample_is_refused(genome):
    builder = VcfLineBuilder(genome, samples=['S1'])
    with pytest.raises(RecordError):
        builder.build(call('essv1', 'copy_number_loss', '1', 21, 40, samples={'S9': {'GT': '1'}}))


# --- Whole files ---

def sample_records():
    """One record of each kind of edge case, in database order (contig name, then start)."""
    return [
        region('esv1', 'copy_number_variation', '1', 1, 60, so_id='SO:0001019'),
        call('essv2', 'copy_number_loss', '1', 1, 10, so_id='SO:0001743', parents=('esv1',)),
        call('essv3', 'copy_number_gain', '1', 21, 40, ci_pos=ConfidenceInterval(-5, None),
             ci_end=ConfidenceInterval(None, 3), parents=('esv1',), samples={'S1': {'GT': '0/1', 'CN': 3.0}}),
        call('essv4', 'insertion', '1', 30, 30, sequence='ACGTT', parents=('esv1',)),
        call('essv5', 'alu_insertion', '1', 52, 52, placement_method='Submitted cytogenetic'),
        call('essv6', 'copy_number_loss', '1', 62, 70),
        region('esv7', 'short_tandem_repeat_variation', '1', 101, 120, info={'RUS': ['GGCCAATTGG'], 'RUC': [2.0]}),
        region('esv8', 'inversion', '2', 100, 150, outer_start=95, inner_start=103, inner_end=148, outer_end=160),
        call('essv9', 'deletion', '2', 200, 220, info={'SVCLAIM': ['J']}, so_id='SO:0000159'),
        call('essv10', 'copy_number_gain', '2', 200, 230),
    ]


def write_vcf(genome, path, records, samples=('S1', 'S2')):
    builder = VcfLineBuilder(genome, samples=list(samples))
    with StudyVcfWriter(str(path), builder, reference_name='GCA_TEST.1',
                        meta_lines=[('DGVa_study', 'estd1'), ('DGVa_assembly', 'GCA_TEST.1')]) as writer:
        for record in records:
            writer.write(record)
        return writer.finish()


def test_vcf_from_records_is_valid(genome, tmp_path):
    statistics = write_vcf(genome, tmp_path / 'test.vcf', sample_records())
    assert statistics.written == 10
    assert statistics.written_by_kind == {'call': 7, 'region': 3}
    errors = vcf_validator_errors(str(tmp_path / 'test.vcf'))
    if errors is None:
        pytest.skip('vcf_validator is not installed')
    assert errors == []


def test_alleles_the_file_defines_are_declared_and_valid(genome, tmp_path):
    records = [region('esv1', 'complex_chromosomal_rearrangement', '1', 21, 60),
               call('essv2', 'interchromosomal_translocation', '1', 30, 30, parents=('esv1',)),
               call('essv3', 'sequence_alteration', '1', 41, 50, parents=('esv1',)),
               call('essv4', 'deletion', '1', 52, 58)]
    write_vcf(genome, tmp_path / 'test.vcf', records, samples=())
    text = (tmp_path / 'test.vcf').read_text()
    assert [line for line in text.splitlines() if line.startswith('##ALT')] == [
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=TRA,Description="Translocation; the partner breakpoint is not given">',
        '##ALT=<ID=COMPLEX_CHROMOSOMAL_REARRANGEMENT,Description="Complex chromosomal rearrangement of the region">',
        '##ALT=<ID=SEQUENCE_ALTERATION,Description="Change of the sequence of the region; the kind of change is not given">',
    ]
    errors = vcf_validator_errors(str(tmp_path / 'test.vcf'))
    if errors is None:
        pytest.skip('vcf_validator is not installed')
    assert errors == []


def test_contigs_are_in_assembly_report_order_and_ids_are_accessions(genome, tmp_path):
    write_vcf(genome, tmp_path / 'test.vcf', sample_records())
    lines = data_lines(tmp_path / 'test.vcf')
    assert [line[0] for line in lines] == ['2'] * 3 + ['1'] * 7
    assert [line[2] for line in lines] == ['esv8', 'essv10', 'essv9', 'essv2', 'esv1', 'essv3', 'essv4', 'essv5',
                                          'essv6', 'esv7']
    header = (tmp_path / 'test.vcf').read_text().split('\n')
    assert header.index('##contig=<ID=2,length=320>') < header.index('##contig=<ID=1,length=190>')


def test_records_at_the_same_position_are_not_merged(genome, tmp_path):
    records = [call('essv2', 'copy_number_loss', '2', 100, 150), call('essv1', 'copy_number_loss', '2', 100, 120)]
    write_vcf(genome, tmp_path / 'test.vcf', records)
    lines = data_lines(tmp_path / 'test.vcf')
    assert [(line[1], line[2], info_dict(line[7])['END']) for line in lines] == [('99', 'essv1', '120'),
                                                                                  ('99', 'essv2', '150')]


def test_only_used_keys_and_alleles_are_in_the_header(genome, tmp_path):
    write_vcf(genome, tmp_path / 'test.vcf', [call('essv1', 'copy_number_loss', '2', 100, 150)], samples=())
    header = [line for line in (tmp_path / 'test.vcf').read_text().split('\n') if line.startswith('##')]
    assert [line.split('=<ID=')[1].split(',')[0] for line in header if '=<ID=' in line and 'contig' not in line] == \
        ['DEL', 'END', 'SVLEN', 'SVCLAIM', 'AF']
    assert header[:3] == ['##fileformat=VCFv4.4', '##source=dgva-extraction', '##reference=GCA_TEST.1']


def test_output_is_byte_identical_across_runs(genome, tmp_path):
    write_vcf(genome, tmp_path / 'first.vcf', sample_records())
    write_vcf(genome, tmp_path / 'second.vcf', sample_records())
    assert filecmp.cmp(tmp_path / 'first.vcf', tmp_path / 'second.vcf', shallow=False)


def test_order_of_records_at_one_position_does_not_change_the_output(genome, tmp_path):
    records = [call(f'essv{index}', 'copy_number_loss', '2', 100, 100 + index) for index in range(20)]
    write_vcf(genome, tmp_path / 'first.vcf', records)
    shuffled = list(records)
    random.Random(7).shuffle(shuffled)
    write_vcf(genome, tmp_path / 'second.vcf', shuffled)
    assert filecmp.cmp(tmp_path / 'first.vcf', tmp_path / 'second.vcf', shallow=False)


def test_records_out_of_order_stop_the_writer(genome, tmp_path):
    records = [call('essv1', 'copy_number_loss', '2', 100, 150), call('essv2', 'copy_number_loss', '2', 50, 60)]
    with pytest.raises(OrderError):
        write_vcf(genome, tmp_path / 'test.vcf', records)
    assert not (tmp_path / 'test.vcf').exists()
    assert not any(path.name.startswith('.dgva_extraction_') for path in tmp_path.iterdir())


def test_contig_that_comes_back_later_is_still_sorted(genome, tmp_path):
    # Two database names for one contig ('1' and 'chr1') sort apart in the database.
    records = [call('essv1', 'copy_number_loss', '1', 21, 30), call('essv2', 'copy_number_loss', '2', 21, 30),
               call('essv3', 'copy_number_loss', 'chr1', 41, 50)]
    write_vcf(genome, tmp_path / 'test.vcf', records)
    assert [line[2] for line in data_lines(tmp_path / 'test.vcf')] == ['essv2', 'essv1', 'essv3']


def test_refused_record_fails_the_study_and_writes_nothing(genome, tmp_path):
    records = sample_records() + [call('essv99', 'made_up_type', '2', 300, 310),
                                  call('essv98', 'short_tandem_repeat_variation', '2', 300, 310)]
    with pytest.raises(StudyFailedError) as error:
        write_vcf(genome, tmp_path / 'test.vcf', records)
    assert not (tmp_path / 'test.vcf').exists()
    assert not any(path.name.startswith('.dgva_extraction_') for path in tmp_path.iterdir())
    statistics = error.value.statistics
    assert statistics.refused == {'unknown_type': 1, 'missing_repeat_unit': 1}
    assert statistics.written == 10


def test_statistics_per_kind_type_and_chromosome(genome, tmp_path):
    statistics = write_vcf(genome, tmp_path / 'test.vcf', sample_records())
    assert statistics.written_by_chromosome == {'1': 7, '2': 3}
    assert statistics.written_by_type['copy_number_loss'] == 2
    statistics.write_tsv(tmp_path / 'stats.tsv')
    assert 'kind\tcall\t7' in (tmp_path / 'stats.tsv').read_text()


def test_vcf_without_samples_has_a_missing_af_on_every_line(genome, tmp_path):
    records = [call('essv1', 'copy_number_loss', '2', 100, 150),
               call('essv2', 'copy_number_loss', '2', 160, 170, info={'AF': [0.25]})]
    write_vcf(genome, tmp_path / 'test.vcf', records, samples=())
    assert [info_dict(line[7])['AF'] for line in data_lines(tmp_path / 'test.vcf')] == ['.', '0.25']


def test_vcf_with_samples_has_no_missing_af(genome, tmp_path):
    write_vcf(genome, tmp_path / 'test.vcf', [call('essv1', 'copy_number_loss', '2', 100, 150)], samples=('S1',))
    assert 'AF' not in info_dict(data_lines(tmp_path / 'test.vcf')[0][7])


# --- Ticket 10, point 4: breakends and variants written as several records ---

LEFT, RIGHT = BreakendSide.LEFT, BreakendSide.RIGHT


@pytest.mark.parametrize('side, mate_side, alt', [
    (LEFT, RIGHT, 'T[2:50['), (LEFT, LEFT, 'T]2:50]'), (RIGHT, RIGHT, '[2:50[T'), (RIGHT, LEFT, ']2:50]T'),
])
def test_breakend_notation(builder, side, mate_side, alt):
    # CHR1 position 64 is T (61 R, 62-64 TTG).
    record = call('nssv1', 'interchromosomal_translocation', '1', 63, 63,
                  breakend=Breakend(side, '2', 50, mate_side),
                  info={'EVENT': ['nssv1'], 'EVENTTYPE': ['TRA'], 'MATEID': ['nssv1_2']})
    line = builder.build(record)
    assert (line.pos, line.ref, line.alt) == (63, 'T', alt)
    keys = [key for key, _ in line.info]
    assert 'END' not in keys and 'SVLEN' not in keys and line.symbolic_allele is None


@pytest.mark.parametrize('side, alt', [(LEFT, 'T.'), (RIGHT, '.T')])
def test_single_breakend_notation(builder, side, alt):
    record = call('nsv9', 'complex_chromosomal_rearrangement', '1', 63, 63, vcf_id='nsv9_2',
                  breakend=Breakend(side), info={'EVENT': ['nsv9'], 'EVENTTYPE': ['COMPLEX']})
    line = builder.build(record)
    assert (line.id, line.alt) == ('nsv9_2', alt)
    assert dict(line.info)['EVENT'] == 'nsv9'


def test_breakend_on_an_unknown_mate_chromosome_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('nssv1', 'translocation', '1', 63, 63, breakend=Breakend(LEFT, '7', 5, RIGHT)))
    assert error.value.reason == 'unknown_chromosome'


def test_breakend_spanning_several_bases_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('nssv1', 'translocation', '1', 63, 70, breakend=Breakend(LEFT)))
    assert error.value.reason == 'invalid_breakend'


def test_undeclared_event_type_is_refused(builder):
    with pytest.raises(RecordError) as error:
        builder.build(call('nssv1', 'deletion', '1', 21, 40, info={'EVENT': ['nssv1'], 'EVENTTYPE': ['MADE_UP']}))
    assert error.value.reason == 'unknown_event_type'


def test_variants_written_as_several_records_are_valid(genome, tmp_path):
    records = [
        # A translocation (From chr1:63 '+', To chr2:50 '-'), a complex rearrangement as single breakends, and
        # alternative locations of a CNV.
        call('nssv1', 'interchromosomal_translocation', '1', 63, 63, breakend=Breakend(LEFT, '2', 50, LEFT),
             info={'EVENT': ['nssv1'], 'EVENTTYPE': ['TRA'], 'MATEID': ['nssv1_2'], 'PLACEMENT_RANK': 'BestAvailable'}),
        call('nsv9', 'complex_chromosomal_rearrangement', '1', 100, 100, breakend=Breakend(LEFT),
             info={'EVENT': ['nsv9'], 'EVENTTYPE': ['COMPLEX']}),
        call('nsv9', 'complex_chromosomal_rearrangement', '1', 120, 120, vcf_id='nsv9_2', breakend=Breakend(RIGHT),
             info={'EVENT': ['nsv9'], 'EVENTTYPE': ['COMPLEX']}),
        call('nssv8', 'copy_number_variation', '1', 130, 150, info={'EVENT': ['nssv8'],
                                                                   'PLACEMENT_RANK': 'BestAvailable'}),
        call('nssv1', 'interchromosomal_translocation', '2', 50, 50, vcf_id='nssv1_2',
             breakend=Breakend(LEFT, '1', 63, LEFT),
             info={'EVENT': ['nssv1'], 'EVENTTYPE': ['TRA'], 'MATEID': ['nssv1'], 'PLACEMENT_RANK': 'NotBest'}),
        call('nssv8', 'copy_number_variation', '2', 200, 220, vcf_id='nssv8_2',
             info={'EVENT': ['nssv8'], 'PLACEMENT_RANK': 'NotBest'}),
    ]
    write_vcf(genome, tmp_path / 'test.vcf', records, samples=())
    text = (tmp_path / 'test.vcf').read_text()
    assert '##EVENTTYPE=<ID=COMPLEX,Description="Complex chromosomal rearrangement' in text
    assert '##EVENTTYPE=<ID=TRA' not in text  # reserved by VCF 4.4: not declared
    assert [line[2] for line in data_lines(tmp_path / 'test.vcf')] == ['nssv1_2', 'nssv8_2', 'nssv1', 'nsv9',
                                                                       'nsv9_2', 'nssv8']
    errors = vcf_validator_errors(str(tmp_path / 'test.vcf'))
    if errors is None:
        pytest.skip('vcf_validator is not installed')
    # vcf_validator warns on a single breakend joined on its right ('.T'), whose ALT doesn't start with REF by design.
    assert [error for error in errors if 'do not share the first nucleotide' not in error] == []
