"""Ticket 13.5: rows from the database to a VCF through the CLI, with a fake database (no Oracle)."""
import filecmp

import pytest

from conftest import vcf_validator_errors
from dgva_extraction import cli
from test_reader import FakeDatabase, placement_row


class FakeSession(FakeDatabase):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def sample_study_rows():
    """A small study on chromosomes 1 and 2, in database order, with the sprint 1-2 edge cases."""
    return [
        placement_row(1, 'esv1', '1', 11, 60, so_term='copy_number_variation', so_id='SO:0001019',
                      placement_method=None),
        placement_row(2, 'essv2', '1', 11, 30, placement_method=None, parent_accession='esv1', parent_in_study=1),
        placement_row(3, 'essv3', '1', 21, 60, so_term='copy_number_gain', so_id='SO:0001742',
                      placement_method=None, ciposleft=5, ciposright=5, ciendleft=3, ciendright=10,
                      parent_accession='esv1', parent_in_study=1),
        placement_row(4, 'essv4', '1', 62, 70, placement_method=None),  # W1: no region
        placement_row(5, 'essv5', '1', None, None, placement_method=None),  # W3: no start
        placement_row(6, 'esv6', '2', 100, 150, so_term='inversion', so_id='SO:1000036', placement_method=None,
                      outer_start_pos=95, inner_start_pos=103, inner_stop_pos=148, outer_stop_pos=160),
        placement_row(7, 'esv6', '2', 100, 150, so_term='inversion', so_id='SO:1000036', placement_method=None,
                      outer_start_pos=95, inner_start_pos=103, inner_stop_pos=148, outer_stop_pos=160),
    ]


@pytest.fixture
def run_extract(reference_files, tmp_path, monkeypatch):
    fasta, report = reference_files
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h, port: 1521, service_name: s, user: u, password: p}\n')

    def run(rows, output_dir, study='estd1'):
        monkeypatch.setattr(cli, 'connect', lambda args: FakeSession(rows))
        status = cli.main(['--config', str(config), 'extract', '--study', study, '--assembly', 'GCA_TEST.1',
                           '--fasta', fasta, '--assembly-report', report, '--output-dir', str(output_dir)])
        return status, output_dir / f'{study}_GCA_TEST.1'
    return run


def vcf_ids(path):
    return [line.split('\t')[2] for line in path.read_text().splitlines() if not line.startswith('#')]


def test_extract_writes_regions_and_calls_apart_and_reconciles_the_counts(run_extract, tmp_path):
    status, prefix = run_extract(sample_study_rows(), tmp_path / 'out')
    assert status == 0
    calls = tmp_path / 'out' / 'estd1_GCA_TEST.1.calls.vcf'
    regions = tmp_path / 'out' / 'estd1_GCA_TEST.1.regions.vcf'
    assert vcf_ids(calls) == ['essv2', 'essv3', 'essv4']
    assert vcf_ids(regions) == ['esv6', 'esv1']
    assert '##DGVa_records=calls' in calls.read_text() and '##DGVa_records=regions' in regions.read_text()
    for vcf in (calls, regions):
        assert '##DGVa_placement_method=unknown' in vcf.read_text()
        assert '##DGVa_study=estd1' in vcf.read_text()
        # No sample columns: AF is on every line, as a missing value.
        assert all(';AF=.' in line or '\tAF=.' in line
                   for line in vcf.read_text().splitlines() if not line.startswith('#'))
    side_list = (tmp_path / 'out' / 'estd1_GCA_TEST.1.side_list.tsv').read_text().splitlines()
    assert [line.split('\t')[2:5] for line in side_list[1:]] == [['essv5', 'call', 'no_start'],
                                                                 ['esv6', 'region', 'duplicate_placement']]
    # Records written = submitted placements - placements on the side list.
    written = len(vcf_ids(calls)) + len(vcf_ids(regions))
    assert written == len({row['placement_id'] for row in sample_study_rows()}) - (len(side_list) - 1)
    stats = (tmp_path / 'out' / 'estd1_GCA_TEST.1.stats.tsv').read_text()
    assert 'all\treader\tcalls_without_region\t1' in stats and 'all\trun\tpeak_memory_mb' in stats
    assert 'calls\tkind\tcall\t3' in stats and 'regions\tkind\tregion\t2' in stats
    for vcf in (calls, regions):
        errors = vcf_validator_errors(str(vcf))
        if errors is None:
            pytest.skip('vcf_validator is not installed')
        assert errors == []


def test_study_with_no_region_writes_no_region_vcf(run_extract, tmp_path):
    rows = [row for row in sample_study_rows() if row['region_or_call'] == 'call']
    status, _ = run_extract(rows, tmp_path / 'out')
    assert status == 0
    assert (tmp_path / 'out' / 'estd1_GCA_TEST.1.calls.vcf').exists()
    assert not (tmp_path / 'out' / 'estd1_GCA_TEST.1.regions.vcf').exists()
    assert not any(path.name.startswith('.dgva_extraction_') for path in (tmp_path / 'out').iterdir())


def test_extract_is_deterministic(run_extract, tmp_path):
    run_extract(sample_study_rows(), tmp_path / 'first')
    run_extract(sample_study_rows(), tmp_path / 'second')
    for name in ('estd1_GCA_TEST.1.calls.vcf', 'estd1_GCA_TEST.1.regions.vcf'):
        assert filecmp.cmp(tmp_path / 'first' / name, tmp_path / 'second' / name, shallow=False)


def test_extract_fails_the_study_on_a_refused_record(run_extract, tmp_path):
    rows = sample_study_rows() + [placement_row(8, 'essv8', '2', 200, 210, so_term='sequence_alteration',
                                                so_id='SO:0001059', placement_method=None)]
    status, _ = run_extract(rows, tmp_path / 'out')
    assert status == 1
    # Neither file is written, although only a call was refused.
    assert sorted(path.name for path in (tmp_path / 'out').iterdir()) == ['estd1_GCA_TEST.1.side_list.tsv',
                                                                          'estd1_GCA_TEST.1.stats.tsv']
    stats = (tmp_path / 'out' / 'estd1_GCA_TEST.1.stats.tsv').read_text()
    assert 'calls\trefused\tno_symbolic_allele\t1' in stats
    assert 'essv8' in stats
