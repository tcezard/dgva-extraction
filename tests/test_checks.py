"""The off-hours checks: merged experiment descriptions, calls analyses, validations, and the command (fake database)."""
from dgva_extraction import checks, cli
from test_cli import FakeSession


def experiment(experiment_id, study, submitter_id, experiment_type='Discovery', method='Oligo aCGH',
               analysis='Probe signal intensity', detection='Log2 ratio', platform='Agilent 244K', has_calls=1):
    return {'experiment_id': experiment_id, 'study_accession': study, 'submitter_experiment_id': submitter_id,
            'experiment_type': experiment_type, 'method_type': method, 'analysis_type': analysis,
            'detection_method': detection, 'platform_name': platform, 'has_calls': has_calls}


EXPERIMENTS = [
    # nstd1: two identical discovery experiments (merged), one on another platform, and two validation methods.
    experiment(10, 'nstd1', 2),
    experiment(11, 'nstd1', 1),
    experiment(12, 'nstd1', 3, platform='Illumina 1M'),
    experiment(13, 'nstd1', 4, 'Validation', method='PCR', analysis='Manual observation', detection=None,
               platform=None, has_calls=0),
    experiment(14, 'nstd1', 5, 'Validation', method='FISH', analysis='Manual observation', detection=None,
               platform=None, has_calls=0),
    # An experiment with two platforms: one row per platform.
    experiment(20, 'estd2', 1, platform='A'),
    experiment(20, 'estd2', 1, platform='B'),
]


def test_identical_descriptions_merge_and_ids_follow_the_submitter_ids():
    descriptions = checks.Descriptions.from_rows(EXPERIMENTS)
    assert [descriptions.description_of[i] for i in (11, 10, 12, 13, 14, 20)] == ['D1', 'D1', 'D2', 'V1', 'V2', 'D1']
    assert descriptions.experiments[20].parts['platform_name'] == {'A', 'B'}
    assert descriptions.studies() == ['estd2', 'nstd1']
    assert descriptions.summary()[0].startswith('- Discovery: 4 experiments (4 with calls) in 2 studies; 3 merged')


def test_calls_analyses_add_up_the_experiments_of_one_description():
    descriptions = checks.Descriptions.from_rows(EXPERIMENTS)
    rows = [{'study_accession': 'nstd1', 'assembly_accession': 'GCA_1', 'experiment_id': 10, 'calls': 5},
            {'study_accession': 'nstd1', 'assembly_accession': 'GCA_1', 'experiment_id': 11, 'calls': 7},
            {'study_accession': 'nstd1', 'assembly_accession': 'GCA_1', 'experiment_id': 12, 'calls': 1},
            {'study_accession': 'nstd1', 'assembly_accession': 'GCA_2', 'experiment_id': 12, 'calls': 2},
            {'study_accession': 'estd2', 'assembly_accession': 'GCA_1', 'experiment_id': 20, 'calls': 3},
            {'study_accession': 'estd2', 'assembly_accession': 'GCA_1', 'experiment_id': 99, 'calls': 1}]
    analyses = checks.calls_analyses(rows, descriptions)
    assert analyses[('nstd1', 'GCA_1', 'D1')] == 12
    assert len([key for key in analyses if key[0] == 'nstd1']) == 3
    lines = checks.calls_analyses_summary(analyses, 2)
    assert '**5 calls analyses**' in lines[2]
    assert 'unknown experiment 99' in lines[-1]


def test_validation_summary_finds_several_methods_and_mixed_results():
    descriptions = checks.Descriptions.from_rows(EXPERIMENTS)
    rows = [{'accession': 'nssv1', 'experiment_id': 13, 'validation_result': 'Pass'},
            {'accession': 'nssv1', 'experiment_id': 14, 'validation_result': 'Fail'},
            {'accession': 'nssv2', 'experiment_id': 13, 'validation_result': 'Pass'},
            {'accession': 'nssv2', 'experiment_id': 13, 'validation_result': 'Inconclusive'},
            {'accession': 'nssv3', 'experiment_id': 13, 'validation_result': 'Pass'},
            {'accession': 'nssv3', 'experiment_id': 13, 'validation_result': 'Pass'}]
    summary = checks.validation_summary(rows, descriptions)
    assert (summary.calls, summary.several_descriptions, summary.mixed_results,
            summary.mixed_within_one_description) == (3, 1, 2, 1)
    assert summary.examples[0] == {'accession': 'nssv1', 'study_accession': 'nstd1',
                                   'validations': 'V1:Pass, V2:Fail', 'several_descriptions': 1, 'mixed_results': 1}


def test_checks_command_runs_every_step_and_reports_a_failed_one(tmp_path, monkeypatch):
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h, port: 1521, service_name: s, user: u, password: p}\n')
    session = FakeSession()
    session.results = {
        'study_assembly_list': [{'study_accession': 'nstd1', 'submitted_placements': 3, 'submitted_variants': 3}],
        'so_terms': [{'region_or_call': 'call', 'exchange_type': 'deletion', 'so_term': 'deletion',
                      'so_id': 'SO:0000159', 'variants': 3},
                     {'region_or_call': 'call', 'exchange_type': 'made up', 'so_term': 'made_up', 'so_id': None,
                      'variants': 1}],
        'check_experiment_descriptions': EXPERIMENTS,
        'check_calls_per_experiment': [{'study_accession': 'nstd1', 'assembly_accession': 'GCA_1',
                                        'experiment_id': 10, 'calls': 3, 'placements': 3}],
    }  # No 'check_multiple_validations': that step fails, and is reported.
    monkeypatch.setattr(cli, 'connect', lambda args: session)
    output = tmp_path / 'out'
    assert cli.main(['--config', str(config), 'checks', '--output-dir', str(output)]) == 1
    summary = (output / 'summary.md').read_text()
    for step in cli.CHECK_STEPS:
        assert f'## {step}' in summary
    assert '**differs from** sprint 2' in summary
    assert "call 'made up' (made_up, 1 variants)" in summary
    assert '## validations' in summary and '**Failed:** KeyError' in summary
    assert {path.name for path in output.iterdir()} >= {
        'study_assembly_list.tsv', 'check_types.tsv', 'experiment_descriptions.tsv', 'calls_per_experiment.tsv',
        'calls_analyses.tsv', 'summary.md'}
    # Per-study step over every study with experiments, each queried once.
    assert [binds for name, binds in session.calls if name == 'check_calls_per_experiment'] == [
        {'study': 'estd2'}, {'study': 'nstd1'}]
    assert sum(1 for name, _ in session.calls if name == 'check_experiment_descriptions') == 1
