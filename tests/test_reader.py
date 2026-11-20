"""Ticket 13: mapping database rows to records, and the edge cases of sprints 1 and 2, with a fake database."""
import pytest

from dgva_extraction.reader import ExtractionError, StudyReader, placement_class, study_side_list
from dgva_extraction.record import ConfidenceInterval, RecordKind


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
    def __init__(self, records=(), translocations=(), side_list=()):
        self.results = {'study_records': list(records), 'translocation_regions': list(translocations),
                        'study_side_list': list(side_list)}
        self.calls = []

    def stream(self, name, binds):
        self.calls.append((name, binds))
        return iter(self.results[name])


def read(rows, study='estd1', translocations=()):
    reader = StudyReader(FakeDatabase(rows, translocations), study, 'GCA_000001405.1')
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
    assert record.info == {'ASSERTM': 'Merging'}
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
