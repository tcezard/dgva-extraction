"""docs/edge-cases.md must name every way the code refuses a record or leaves one out of the VCF."""
import re
from importlib import resources
from pathlib import Path

EDGE_CASES = Path(__file__).parent.parent / 'docs' / 'edge-cases.md'


def package_text(*parts):
    return resources.files('dgva_extraction').joinpath(*parts).read_text()


def reasons_in_code():
    reasons = set()
    for name in ('vcf_writer.py', 'reader.py'):
        source = package_text(name)
        reasons |= set(re.findall(r"RecordError\('([a-z_]+)'", source))
        reasons |= set(re.findall(r"_side\([^)]*?'([a-z_]+)',", source))
    side_list_sql = package_text('sql', 'study_side_list.sql')
    reasons |= {match for pair in re.findall(r"THEN '([a-z_]+)'|ELSE '([a-z_]+)'", side_list_sql)
                for match in pair if match}
    return reasons


def test_code_reasons_are_found():
    # Guard against the patterns above silently matching nothing.
    reasons = reasons_in_code()
    assert {'no_symbolic_allele', 'unknown_type', 'no_start', 'duplicate_placement', 'withdrawn'} <= reasons


def test_every_reason_is_in_the_edge_case_list():
    listed = EDGE_CASES.read_text()
    missing = sorted(reason for reason in reasons_in_code() if f'`{reason}`' not in listed)
    assert missing == [], f'Add these reasons to docs/edge-cases.md: {missing}'
