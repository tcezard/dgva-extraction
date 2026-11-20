"""Ticket 13: mapping database rows to records, and the edge cases of sprints 1 and 2, with a fake database."""
from collections import Counter

import pytest

from dgva_extraction.keys import SymbolicAlleles
from dgva_extraction.reader import (ExtractionError, StudyReader, chromosome_sort_key, placement_class,
                                     plan_multi_placements, study_side_list)
from dgva_extraction.record import Breakend, BreakendSide, ConfidenceInterval, RecordKind


def placement_row(placement_id, accession, chromosome='1', start=100, stop=200, **overrides):
    kind = 'region' if accession.startswith(('esv', 'nsv')) else 'call'
    row = {
        'placement_id': placement_id, 'accession': accession, 'region_or_call': kind,
        'so_term': 'copy_number_loss', 'so_id': 'SO:0001743', 'exchange_type': 'copy number loss',
        'in_type_table': 1, 'assertion_method': 'Merging' if kind == 'region' else None, 'assembly_known': 1,
        'placement_method': 'Submitted genomic', 'has_remap_score': 0,
        'chromosome_name': chromosome, 'contig_id': None, 'chromosome': chromosome,
        'start_pos': start, 'stop_pos': stop, 'outer_start_pos': None, 'inner_start_pos': None,
        'inner_stop_pos': None, 'outer_stop_pos': None,
        'ci_start_minus': None, 'ci_start_plus': None, 'ci_stop_minus': None, 'ci_stop_plus': None,
        'ciposleft': None, 'ciposright': None, 'ciendleft': None, 'ciendright': None,
        'parent_accession': None, 'parent_in_study': 0, 'parent_withdrawn': 0,
    }
    row.update(overrides)
    return row


class FakeDatabase:
    def __init__(self, records=(), translocations=(), side_list=(), multi_placements=()):
        self.results = {'study_records': list(records), 'translocation_regions': list(translocations),
                        'study_side_list': list(side_list), 'multi_placement_variants': list(multi_placements)}
        self.calls = []

    def stream(self, name, binds):
        self.calls.append((name, binds))
        return iter(self.results[name])


def read(rows, study='estd1', translocations=(), multi_placements=()):
    reader = StudyReader(FakeDatabase(rows, translocations, multi_placements=multi_placements), study,
                         'GCA_000001405.1')
    return reader, list(reader.records())


@pytest.mark.parametrize('method, remap, expected', [
    (None, 0, 'submitted'), ('Submitted genomic', 0, 'submitted'), ('Submitted cytogenetic', 0, 'submitted'),
    ('Submitted artifact', 0, 'submitted'), ('Remapped', 1, 'remapped'), ('Remapped ad-hoc', 0, 'remapped'),
    ('Remap failed', 0, 'no_position'), ('Submitted without placement on assembly', 0, 'no_position'),
    ('Submitted on non-INSDC assembly', 0, 'no_position'), (None, 1, 'conflict'),
    ('Submitted genomic', 1, 'conflict'), ('Something new', 0, 'other'),
])
def test_placement_class_follows_the_ticket_2_rule(method, remap, expected):
    assert placement_class(method, remap) == expected


def test_reader_queries_with_binds():
    database = FakeDatabase()
    reader = StudyReader(database, 'estd1', 'GCA_1')
    list(reader.records())
    assert ('study_records', {'study': 'estd1', 'assembly': 'GCA_1'}) in database.calls


def test_row_to_record():
    reader, records = read([placement_row(1, 'esv1', so_term='copy number variation', so_id='SO:0001019')])
    record = records[0]
    assert (record.accession, record.kind, record.so_term, record.so_id) == \
        ('esv1', RecordKind.REGION, 'copy_number_variation', 'SO:0001019')
    assert (record.chromosome, record.start, record.end, record.placement_id) == ('1', 100, 200, 1)
    assert record.info == {'REGION_ASSERTION_METHOD': 'Merging'}
    assert record.placement_method == 'Submitted genomic'


def test_calls_point_at_their_regions_and_region_to_region_links_are_kept():
    rows = [
        placement_row(1, 'esv1', parent_accession='esv0', parent_in_study=1),
        placement_row(2, 'essv2', parent_accession='esv1', parent_in_study=1),
        placement_row(2, 'essv2', parent_accession='esv9', parent_in_study=1),
    ]
    reader, records = read(rows)
    assert [record.parents for record in records] == [('esv0',), ('esv1', 'esv9')]
    assert reader.counts['rows'] == 3 and reader.counts['placements'] == 2


def test_w1_call_without_region_is_written_without_parent():
    reader, records = read([placement_row(1, 'essv1')], study='estd236')
    assert records[0].parents == ()
    assert reader.counts['calls_without_region'] == 1


def test_withdrawn_parent_is_left_out():
    reader, records = read([placement_row(1, 'essv1', parent_accession='esv9', parent_in_study=1,
                                          parent_withdrawn=1)])
    assert records[0].parents == ()
    assert reader.counts['withdrawn_parents_left_out'] == 1


def test_parent_from_another_study_fails_the_study():
    with pytest.raises(ExtractionError):
        read([placement_row(1, 'essv1', parent_accession='esv9', parent_in_study=0)])


@pytest.mark.parametrize('overrides, message', [
    ({'so_term': None, 'so_id': None}, 'no SO term'),
    ({'assembly_known': 0}, 'GENOME_ASSEMBLY'),
    ({'in_type_table': 0}, 'not in VARIANT_REGION'),
    ({'assertion_method': None}, 'assertion method'),
])
def test_w2_fails_the_study(overrides, message):
    with pytest.raises(ExtractionError) as error:
        read([placement_row(1, 'esv1', **overrides)])
    assert message in str(error.value)


@pytest.mark.parametrize('method, remap', [('Submitted genomic', 1), (None, 1), ('Unknown method', 0)])
def test_conflict_and_unknown_methods_fail_the_study(method, remap):
    with pytest.raises(ExtractionError):
        read([placement_row(1, 'essv1', placement_method=method, has_remap_score=remap)])


def test_w3_placement_with_no_start_goes_to_the_side_list():
    rows = [placement_row(1, 'esv4011099', start=None, stop=None), placement_row(2, 'esv4011099', start=500,
                                                                                  stop=900)]
    reader, records = read(rows, study='estd236')
    assert [record.placement_id for record in records] == [2]
    assert [(entry.accession, entry.reason, entry.placement_id) for entry in reader.side_list] == \
        [('esv4011099', 'no_start', 1)]


def test_w3_placement_with_a_start_but_no_chromosome_goes_to_the_side_list():
    reader, records = read([placement_row(1, 'essv1', chromosome=None)], study='estd240')
    assert records == []
    assert reader.side_list[0].reason == 'no_chromosome'


def test_contig_name_is_used_when_there_is_no_chromosome():
    reader, records = read([placement_row(1, 'essv1', chromosome='AADC01000001.1', chromosome_name=None,
                                          contig_id='AADC01000001.1', placement_method='Submitted artifact')])
    assert records[0].chromosome == 'AADC01000001.1'


def test_start_and_end_fall_back_to_outer_then_inner():
    rows = [placement_row(1, 'essv1', start=None, stop=None, outer_start_pos=90, inner_start_pos=110,
                          inner_stop_pos=190, outer_stop_pos=210),
            placement_row(2, 'essv2', start=None, stop=None, inner_start_pos=110, inner_stop_pos=190)]
    _, records = read(rows)
    assert (records[0].start, records[0].end) == (90, 210)
    assert records[0].confidence_intervals() == (ConfidenceInterval(0, 20), ConfidenceInterval(-20, 0))
    assert (records[1].start, records[1].end) == (110, 190)


def test_insertion_with_no_end_ends_at_its_start():
    _, records = read([placement_row(1, 'essv1', so_term='insertion', stop=None)])
    assert records[0].end == 100


def test_deletion_with_no_end_goes_to_the_side_list():
    reader, records = read([placement_row(1, 'essv1', stop=None)])
    assert records == [] and reader.side_list[0].reason == 'no_end'


def test_cipos_ciend_are_unsigned_distances():
    _, records = read([placement_row(1, 'essv1', ciposleft=17, ciposright=5, ciendleft=3, ciendright=None)],
                      study='estd214')
    assert records[0].ci_pos == ConfidenceInterval(-17, 5)
    assert records[0].ci_end == ConfidenceInterval(-3, None)


def test_ci_start_stop_are_signed_in_estd219():
    _, records = read([placement_row(1, 'essv1', ci_start_minus=-7, ci_start_plus=8, ci_stop_minus=-2,
                                     ci_stop_plus=4, outer_start_pos=93, inner_start_pos=108)], study='estd219')
    assert records[0].confidence_intervals() == (ConfidenceInterval(-7, 8), ConfidenceInterval(-2, 4))


def test_ci_start_stop_are_unsigned_in_estd204():
    _, records = read([placement_row(1, 'essv1', ci_start_minus=7, ci_start_plus=8, ci_stop_minus=2,
                                     ci_stop_plus=4)], study='estd204')
    assert records[0].confidence_intervals() == (ConfidenceInterval(-7, 8), ConfidenceInterval(-2, 4))


@pytest.mark.parametrize('study, values', [('estd219', (7, 8)), ('estd204', (-7, 8))])
def test_ci_with_an_unexpected_sign_fails_the_study(study, values):
    with pytest.raises(ExtractionError):
        read([placement_row(1, 'essv1', ci_start_minus=values[0], ci_start_plus=values[1])], study=study)


def test_cis_from_inner_and_outer_when_no_ci_table():
    _, records = read([placement_row(1, 'essv1', outer_start_pos=95, inner_start_pos=103, inner_stop_pos=198,
                                     outer_stop_pos=210)])
    assert records[0].ci_pos is None
    assert records[0].confidence_intervals() == (ConfidenceInterval(-5, 3), ConfidenceInterval(-2, 10))


def test_both_ci_sets_on_one_placement_fail_the_study():
    with pytest.raises(ExtractionError):
        read([placement_row(1, 'essv1', ciposleft=1, ciposright=1, ci_start_minus=1, ci_start_plus=1)])


def test_duplicated_ci_rows_that_disagree_fail_the_study():
    with pytest.raises(ExtractionError) as error:
        read([placement_row(1, 'essv1', ciposleft=5, ciposright=5), placement_row(1, 'essv1', ciposleft=6,
                                                                                  ciposright=5)])
    assert 'ciposleft' in str(error.value)


def test_duplicated_rows_that_agree_give_one_record():
    _, records = read([placement_row(1, 'essv1', ciposleft=5, ciposright=5)] * 3)
    assert len(records) == 1


def test_estd208_identical_placements_are_written_once():
    rows = [placement_row(placement_id, 'esv2', inner_start_pos=100, inner_stop_pos=200, start=None, stop=None,
                          placement_method=None) for placement_id in (11, 12, 13, 14, 15)]
    rows.append(placement_row(16, 'esv3', start=100, stop=300))
    reader, records = read(rows, study='estd208')
    assert [(record.accession, record.placement_id) for record in records] == [('esv2', 11), ('esv3', 16)]
    assert [(entry.reason, entry.placement_id) for entry in reader.side_list] == \
        [('duplicate_placement', placement_id) for placement_id in (12, 13, 14, 15)]
    assert reader.counts['null_method_placements'] == 5


def test_same_accession_with_different_coordinates_gives_two_records():
    _, records = read([placement_row(1, 'essv1', stop=200), placement_row(2, 'essv1', stop=250)])
    assert [record.end for record in records] == [200, 250]


def translocation_row(accession, call, chromosome, start, stop, placement_id):
    return {'accession': accession, 'so_term': 'interchromosomal_translocation', 'so_id': 'SO:0002060',
            'exchange_type': 'interchromosomal translocation', 'assertion_method': 'Merging',
            'call_accession': call, 'call_placement_id': placement_id, 'placement_method': 'Submitted genomic',
            'chromosome': chromosome, 'start_pos': start, 'stop_pos': stop}


def test_translocation_region_takes_its_position_from_its_calls():
    rows = [placement_row(1, 'nssv1', chromosome='1', start=50, stop=50),
            placement_row(2, 'nssv2', chromosome='1', start=500, stop=500)]
    translocations = [translocation_row('nsv9', 'nssv3', '2', 300, 300, 3),
                      translocation_row('nsv9', 'nssv2', '1', 500, 500, 2)]
    reader, records = read(rows, study='nstd83', translocations=translocations)
    assert [(record.accession, record.chromosome, record.start) for record in records] == \
        [('nssv1', '1', 50), ('nssv2', '1', 500), ('nsv9', '1', 500)]
    assert records[2].kind == RecordKind.REGION
    assert reader.side_list[0].reason == 'positioned_from_calls'
    assert 'nssv2, nssv3' in reader.side_list[0].detail


def test_study_side_list():
    rows = [
        {'accession': 'esv1', 'region_or_call': 'region', 'reason': 'withdrawn', 'withdrawal_reason': 'Merged',
         'withdrawn_by_study': 'estd222', 'n_placements': 2, 'n_remapped': 0, 'n_no_position': 0},
        {'accession': 'nssv2', 'region_or_call': 'call', 'reason': 'remapped_only', 'withdrawal_reason': None,
         'withdrawn_by_study': None, 'n_placements': 2, 'n_remapped': 2, 'n_no_position': 0},
    ]
    entries = study_side_list(FakeDatabase(side_list=rows), 'nstd35')
    assert [(entry.accession, entry.reason) for entry in entries] == [('esv1', 'withdrawn'),
                                                                      ('nssv2', 'remapped_only')]
    assert 'estd222' in entries[0].detail


# --- Ticket 10, point 4: variants with several submitted placements on one assembly ---

def multi_row(placement_id, accession, chromosome, start, stop=None, so_term='interchromosomal translocation',
              rank='NotBest', order=None, strand=None):
    return {'accession': accession, 'placement_id': placement_id, 'so_term': so_term, 'chromosome': chromosome,
            'start_pos': start, 'stop_pos': start if stop is None else stop, 'placement_rank': rank,
            'breakpoint_order': order, 'strand': strand}


def plan(rows):
    return plan_multi_placements(rows, SymbolicAlleles.load())


@pytest.mark.parametrize('from_strand, to_strand, from_side, to_side', [
    # The four joins of dbVar's VCFs: From record t[p[, t]p], [p[t, ]p]t.
    ('+', '+', BreakendSide.LEFT, BreakendSide.RIGHT),
    ('+', '-', BreakendSide.LEFT, BreakendSide.LEFT),
    ('-', '+', BreakendSide.RIGHT, BreakendSide.RIGHT),
    ('-', '-', BreakendSide.RIGHT, BreakendSide.LEFT),
])
def test_translocation_breakpoints_become_a_breakend_pair(from_strand, to_strand, from_side, to_side):
    plans = plan([multi_row(1, 'nssv1', '5', 500, rank='NotBest', order='To', strand=to_strand),
                  multi_row(2, 'nssv1', '17', 300, rank='BestAvailable', order='From', strand=from_strand)])
    # The best placement keeps the accession; the other gets its own ID. Both carry the accession in EVENT.
    assert (plans[2].vcf_id, plans[1].vcf_id) == ('nssv1', 'nssv1_2')
    assert plans[2].breakend == Breakend(from_side, '5', 500, to_side)
    assert plans[1].breakend == Breakend(to_side, '17', 300, from_side)
    assert plans[2].info == {'EVENT': ['nssv1'], 'EVENTTYPE': ['TRA'], 'MATEID': ['nssv1_2'],
                             'PLACEMENT_RANK': 'BestAvailable'}
    assert plans[1].info['MATEID'] == ['nssv1']


def test_breakend_pair_with_no_breakpoint_order_takes_the_lower_position_as_from():
    counts = Counter()
    plans = plan_multi_placements([multi_row(1, 'nssv1', '17', 300, strand='+', rank='BestAvailable'),
                                   multi_row(2, 'nssv1', '5', 500, strand='-')], SymbolicAlleles.load(), counts)
    # From is chromosome 5 (lower than 17): '-' keeps its right side; To (17) '+' keeps its right side.
    assert plans[2].breakend == Breakend(BreakendSide.RIGHT, '17', 300, BreakendSide.RIGHT)
    assert counts['breakend_pairs_from_lower_position'] == 1


@pytest.mark.parametrize('strands', [('not reported', 'not reported'), (None, None), ('+', None)])
def test_translocation_with_unknown_strand_stays_symbolic_and_linked(strands):
    plans = plan([multi_row(1, 'nssv1', '1', 100, strand=strands[0], rank='BestAvailable'),
                  multi_row(2, 'nssv1', '1', 900, strand=strands[1])])
    assert all(item.breakend is None and item.mode == 'linked' for item in plans.values())
    assert plans[2].info == {'EVENT': ['nssv1'], 'EVENTTYPE': ['TRA'], 'PLACEMENT_RANK': 'NotBest'}


def test_complex_rearrangement_breakpoints_become_single_breakends():
    rows = [multi_row(1, 'nsv9', '1', 100, so_term='complex chromosomal rearrangement', order='From', strand='+'),
            multi_row(2, 'nsv9', '1', 200, so_term='complex chromosomal rearrangement', order='To', strand='+'),
            multi_row(3, 'nsv9', '3', 300, so_term='complex chromosomal rearrangement', order='From', strand='-',
                      rank='BestAvailable')]
    plans = plan(rows)
    assert [plans[i].vcf_id for i in (3, 1, 2)] == ['nsv9', 'nsv9_2', 'nsv9_3']
    assert [plans[i].breakend.side for i in (1, 2, 3)] == [BreakendSide.LEFT, BreakendSide.RIGHT, BreakendSide.RIGHT]
    assert all(plans[i].breakend.is_single and plans[i].info['EVENTTYPE'] == ['COMPLEX'] for i in (1, 2, 3))


def test_inversion_given_as_two_breakpoints_becomes_one_record():
    plans = plan([multi_row(1, 'nssv5', '14', 39513878, so_term='inversion', rank='NotBest'),
                  multi_row(2, 'nssv5', '14', 31269856, so_term='inversion', rank='BestAvailable')])
    assert plans[2].end == 39513878 and plans[2].vcf_id == 'nssv5' and not plans[2].info
    assert plans[1].merged_into == 2


@pytest.mark.parametrize('rows', [
    # b. Alternative locations (nstd116, segmental duplications).
    [multi_row(1, 'nssv8', '8', 7286416, 7366833, so_term='copy number variation', rank='BestAvailable'),
     multi_row(2, 'nssv8', '8', 7669242, 7740180, so_term='copy number variation')],
    # c. Pieces on consecutive contigs (nstd10).
    [multi_row(1, 'nsv4', 'AAEX02001078.1', 1121, 11033, so_term='copy number variation', rank='BestAvailable'),
     multi_row(2, 'nsv4', 'AAEX02001079.1', 1, 321, so_term='copy number variation')],
])
def test_alternative_locations_and_pieces_are_linked_by_event(rows):
    plans = plan(rows)
    assert [plans[1].vcf_id, plans[2].vcf_id] == [rows[0]['accession'], rows[0]['accession'] + '_2']
    assert plans[2].info == {'EVENT': [rows[0]['accession']], 'PLACEMENT_RANK': 'NotBest'}
    assert plans[1].breakend is None and plans[1].end is None


def test_identical_placements_share_a_plan_and_alone_need_none():
    # estd208: five identical placements are one variant on one position: no plan.
    assert plan([multi_row(i, 'esv1', '1', 100, 200, so_term='deletion') for i in range(1, 6)]) == {}
    plans = plan([multi_row(1, 'esv2', '1', 100, 200, so_term='deletion', rank='BestAvailable'),
                  multi_row(2, 'esv2', '1', 100, 200, so_term='deletion'),
                  multi_row(3, 'esv2', '2', 100, 200, so_term='deletion')])
    assert plans[2] is plans[1] and plans[3].vcf_id == 'esv2_2'


def test_no_best_placement_keeps_the_accession_on_the_lowest_placement_id():
    plans = plan([multi_row(7, 'nssv8', '8', 100, 200, so_term='copy number variation', rank=None),
                  multi_row(3, 'nssv8', '8', 500, 600, so_term='copy number variation', rank=None)])
    assert (plans[3].vcf_id, plans[7].vcf_id) == ('nssv8', 'nssv8_2')
    assert 'PLACEMENT_RANK' not in plans[3].info


def test_chromosome_sort_key_is_natural():
    assert sorted(['X', '10', 'chr2', '1', 'MT', 'Y', 'AAEX02001078.1'], key=chromosome_sort_key) == \
        ['1', 'chr2', '10', 'X', 'Y', 'MT', 'AAEX02001078.1']


def test_reader_applies_the_plan_to_the_records():
    rows = [placement_row(1, 'nssv5', '14', 100, 100, so_term='inversion', so_id='SO:1000036'),
            placement_row(2, 'nssv5', '14', 900, 900, so_term='inversion', so_id='SO:1000036'),
            placement_row(3, 'nssv6', '14', 950, 950, so_term='interchromosomal_translocation', so_id='SO:0002060'),
            placement_row(4, 'nssv6', '15', 50, 50, so_term='interchromosomal_translocation', so_id='SO:0002060')]
    multi = [multi_row(1, 'nssv5', '14', 100, so_term='inversion', rank='BestAvailable'),
             multi_row(2, 'nssv5', '14', 900, so_term='inversion'),
             multi_row(3, 'nssv6', '14', 950, order='From', strand='+', rank='BestAvailable'),
             multi_row(4, 'nssv6', '15', 50, order='To', strand='-')]
    reader, records = read(rows, multi_placements=multi)
    assert [(record.id, record.start, record.end) for record in records] == [
        ('nssv5', 100, 900), ('nssv6', 950, 950), ('nssv6_2', 50, 50)]
    assert records[1].breakend == Breakend(BreakendSide.LEFT, '15', 50, BreakendSide.LEFT)
    assert [(entry.accession, entry.reason) for entry in reader.side_list] == [('nssv5', 'merged_breakpoint')]
