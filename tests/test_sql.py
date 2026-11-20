"""Static checks on the packaged SQL (tickets 11 and 13): read-only statements, bind variables only, and portability.

The reader will later point at the PostgreSQL cold-storage copy. The Oracle-only constructs that remain are listed in
ORACLE_ONLY below and in docs/sql-portability.md; any new one fails this test until it is added to both.
"""
import re
from importlib import resources

import pytest

SQL_FILES = sorted(path.name[:-4] for path in resources.files('dgva_extraction').joinpath('sql').iterdir()
                   if path.name.endswith('.sql'))
# Allowed Oracle-only constructs, per file. session_privileges.sql reads Oracle's data dictionary by nature.
ORACLE_ONLY = {
    'session_privileges': {'VARCHAR2', 'DUAL', 'SESSION_PRIVS', 'SESSION_ROLES', 'ALL_TAB_PRIVS', 'ROLE_TAB_PRIVS',
                           'USER'},
}
ORACLE_ONLY_PATTERNS = {
    'VARCHAR2': r'\bVARCHAR2\b', 'DUAL': r'\bDUAL\b', 'NVL': r'\bNVL2?\s*\(', 'DECODE': r'\bDECODE\s*\(',
    'ROWNUM': r'\bROWNUM\b', 'SYSDATE': r'\bSYSDATE\b', 'CONNECT BY': r'\bCONNECT\s+BY\b', 'OUTER (+)': r'\(\+\)',
    'LISTAGG': r'\bLISTAGG\b', 'REGEXP': r'\bREGEXP_\w+', 'TO_CHAR': r'\bTO_CHAR\b', 'MINUS': r'\bMINUS\b',
    'SESSION_PRIVS': r'\bSESSION_PRIVS\b', 'SESSION_ROLES': r'\bSESSION_ROLES\b', 'ALL_TAB_PRIVS': r'\bALL_TAB_PRIVS\b',
    'ROLE_TAB_PRIVS': r'\bROLE_TAB_PRIVS\b', 'USER': r'(?<![_.])\bUSER\b(?!_)', 'DBMS_': r'\bDBMS_\w+',
    'FETCH FIRST': r'\bFETCH\s+FIRST\b',
}
WRITE_WORDS = r'\b(INSERT|UPDATE|DELETE|MERGE|CREATE|ALTER|DROP|TRUNCATE|GRANT|REVOKE|COMMIT|ROLLBACK|BEGIN|' \
              r'DECLARE|EXECUTE|CALL|LOCK)\b'


def sql_text(name):
    return resources.files('dgva_extraction').joinpath('sql', f'{name}.sql').read_text()


def strip_comments_and_strings(sql):
    sql = re.sub(r'--[^\n]*', ' ', sql)
    sql = re.sub(r'/\*.*?\*/', ' ', sql, flags=re.S)
    return re.sub(r"'(?:[^']|'')*'", "''", sql)


def test_there_are_sql_files():
    assert {'session_privileges', 'study_assembly_list', 'study_records', 'study_side_list',
            'translocation_regions', 'so_terms'} <= set(SQL_FILES)


@pytest.mark.parametrize('name', SQL_FILES)
def test_one_read_only_statement(name):
    code = strip_comments_and_strings(sql_text(name)).strip()
    assert re.match(r'^(SELECT|WITH)\b', code, re.I), 'must start with SELECT or WITH'
    assert ';' not in code, 'one statement, with no terminating semicolon'
    assert not re.search(WRITE_WORDS, code, re.I)
    assert not re.search(r'\bFOR\s+UPDATE\b', code, re.I)


@pytest.mark.parametrize('name', SQL_FILES)
def test_binds_match_the_documented_ones(name):
    text = sql_text(name)
    documented = re.search(r'^-- Binds: (.*)$', text, re.M)
    assert documented, 'each file documents its binds'
    used = set(re.findall(r'(?<!:):([a-z_]+)\b', strip_comments_and_strings(text)))
    if documented.group(1).startswith('none'):
        assert used == set()
    else:
        assert used == set(re.findall(r'\b([a-z_]+) \(', documented.group(1)))


@pytest.mark.parametrize('name', SQL_FILES)
def test_oracle_only_constructs_are_listed(name):
    code = strip_comments_and_strings(sql_text(name)).upper()
    found = {label for label, pattern in ORACLE_ONLY_PATTERNS.items() if re.search(pattern, code, re.I)}
    assert found == ORACLE_ONLY.get(name, set())


@pytest.mark.parametrize('name', [name for name in SQL_FILES if name != 'session_privileges'])
def test_tables_are_schema_qualified(name):
    code = strip_comments_and_strings(sql_text(name))
    # Common table expressions ("name AS (") are the only unqualified names allowed.
    ctes = {cte.lower() for cte in re.findall(r'\b([a-z_]\w*)\s+AS\s*\(', code, re.I)}
    for table in re.findall(r'\b(?:FROM|JOIN)\s+([a-z_][\w.]*)', code, re.I):
        assert table.startswith('dgva.') or table.lower() in ctes, table


def test_no_sql_text_is_built_in_python():
    package = resources.files('dgva_extraction')
    for path in package.iterdir():
        if path.name.endswith('.py'):
            source = path.read_text()
            assert not re.search(r'["\']\s*(SELECT|WITH)\b[^"\']*\bFROM\b', source, re.I), \
                f'{path.name} builds SQL in a string'
            assert 'commit()' not in source, f'{path.name} calls commit()'


@pytest.mark.parametrize('name', SQL_FILES)
def test_no_dgva_function_or_package_is_called(name):
    # PUBLIC can execute DGVa packages that write (STAGE_TO_ARCHIVE_*); a SELECT calling one of their functions
    # could write too. Only tables and views of the schema may be named.
    code = strip_comments_and_strings(sql_text(name))
    assert not re.search(r'\bdgva\.\w+\.\w+', code, re.I), 'package member called'
    assert not re.search(r'\bdgva\.\w+\s*\(', code, re.I), 'schema function called'
