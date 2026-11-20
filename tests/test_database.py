"""Ticket 11: the connection module, with a fake connection: streaming, bind variables, no commit, privilege check."""
import pytest

from dgva_extraction.config import ConfigurationError, database_config, load_config
from dgva_extraction.database import DgvaDatabase, WritePrivilegeError, find_write_privileges, load_sql

def privilege(kind, privilege, owner=None, object_name=None, grantee='USER'):
    return {'privilege_kind': kind, 'grantee_kind': 'SESSION' if kind == 'SYSTEM' else grantee,
            'privilege': privilege, 'owner': owner, 'object_name': object_name}


# What the read-only account showed on 2026-11-20: system privileges, plus Oracle's default INHERIT PRIVILEGES.
READ_ONLY_PRIVILEGES = [
    privilege('SYSTEM', 'CREATE SESSION'),
    privilege('SYSTEM', 'SELECT ANY TABLE'),
    privilege('SYSTEM', 'SELECT ANY SEQUENCE'),
    privilege('OBJECT', 'SELECT', 'DGVA', 'VARIANT_STUDY'),
    privilege('OBJECT', 'READ', 'DGVA', 'VARIANT_PLACEMENT'),
    privilege('OBJECT', 'INHERIT PRIVILEGES', 'DGVA', 'DGVA', grantee='PUBLIC'),
]
PUBLIC_PACKAGE_GRANTS = [privilege('OBJECT', 'EXECUTE', 'DGVA', name, grantee='PUBLIC')
                         for name in ('STAGE_TO_ARCHIVE_PKG', 'STAGE_TO_ARCHIVE_VP_PKG')]
PRIVILEGE_COLUMNS = ('privilege_kind', 'grantee_kind', 'privilege', 'owner', 'object_name')


class FakeCursor:
    def __init__(self, connection):
        self.connection = connection
        self.arraysize = 100
        self.prefetchrows = 2
        self.description = None
        self._rows = []
        self.closed = False

    def execute(self, sql, binds=None):
        self.connection.executed.append((sql, binds, self.arraysize))
        columns, rows = self.connection.next_result(sql)
        self.description = [(column.upper(), None) for column in columns]
        self._rows = list(rows)

    def fetchmany(self, size):
        self.connection.fetch_sizes.append(size)
        batch, self._rows = self._rows[:size], self._rows[size:]
        return batch

    def fetchall(self):
        raise AssertionError('Rows must be streamed with fetchmany, never fetched all at once')

    def close(self):
        self.closed = True


class FakeConnection:
    def __init__(self, privileges=READ_ONLY_PRIVILEGES, results=None, autocommit=False):
        self.autocommit = autocommit
        self.privileges = privileges
        self.results = list(results or [])
        self.executed = []
        self.fetch_sizes = []
        self.cursors = []
        self.commits = 0
        self.rollbacks = 0
        self.closed = False

    def cursor(self):
        cursor = FakeCursor(self)
        self.cursors.append(cursor)
        return cursor

    def next_result(self, sql):
        if 'session_privs' in sql:
            return PRIVILEGE_COLUMNS, [tuple(row[column] for column in PRIVILEGE_COLUMNS) for row in self.privileges]
        return self.results.pop(0)

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1

    def close(self):
        self.closed = True


def open_database(connection, arraysize=3):
    database = DgvaDatabase(connection, schema='DGVA', arraysize=arraysize)
    database.check_no_write_privilege()
    return database


def test_rows_are_streamed_in_batches_as_dicts():
    rows = [(index, f'esv{index}') for index in range(7)]
    connection = FakeConnection(results=[(('placement_id', 'accession'), rows)])
    database = open_database(connection, arraysize=3)
    connection.fetch_sizes.clear()
    stream = database.stream('study_records', {'study': 'estd1', 'assembly': 'GCA_1'})
    assert next(stream) == {'placement_id': 0, 'accession': 'esv0'}
    assert connection.fetch_sizes == [3]
    assert len(list(stream)) == 6
    assert connection.fetch_sizes == [3, 3, 3, 3]
    assert connection.cursors[-1].arraysize == 3
    assert connection.cursors[-1].prefetchrows == 4
    assert connection.cursors[-1].closed


def test_values_are_passed_as_bind_variables():
    connection = FakeConnection(results=[(('a',), [])])
    database = open_database(connection)
    list(database.stream('study_records', {'study': 'estd1', 'assembly': 'GCA_1'}))
    sql, binds, _ = connection.executed[-1]
    assert binds == {'study': 'estd1', 'assembly': 'GCA_1'}
    assert 'estd1' not in sql and 'GCA_1' not in sql
    assert ':study' in sql and ':assembly' in sql


@pytest.mark.parametrize('binds', [('estd1',), {'study': ['estd1']}, {'study': object()}])
def test_binds_must_be_a_dict_of_scalars(binds):
    database = open_database(FakeConnection(results=[(('a',), [])]))
    with pytest.raises(TypeError):
        list(database.stream('study_records', binds))


def test_only_packaged_sql_files_can_run():
    database = open_database(FakeConnection())
    with pytest.raises(ValueError):
        list(database.stream('SELECT * FROM dgva.variant_study', {}))
    with pytest.raises(ValueError):
        load_sql('../secrets')


def test_never_commits_and_rolls_back_on_close():
    connection = FakeConnection(results=[(('a',), [(1,)])])
    with open_database(connection) as database:
        list(database.stream('study_records', {'study': 'estd1', 'assembly': 'GCA_1'}))
    assert connection.commits == 0
    assert connection.rollbacks == 1
    assert connection.closed


def test_autocommit_connection_is_refused():
    with pytest.raises(ValueError):
        DgvaDatabase(FakeConnection(autocommit=True))


def test_no_query_runs_before_the_privilege_check():
    database = DgvaDatabase(FakeConnection(results=[(('a',), [])]))
    with pytest.raises(WritePrivilegeError):
        list(database.stream('study_records', {'study': 'estd1', 'assembly': 'GCA_1'}))


def test_privilege_check_runs_with_the_schema_as_a_bind():
    connection = FakeConnection()
    open_database(connection)
    sql, binds, _ = connection.executed[0]
    assert 'session_privs' in sql
    assert binds == {'schema': 'DGVA'}


@pytest.mark.parametrize('write_privilege', [
    privilege('SYSTEM', 'INSERT ANY TABLE'),
    privilege('SYSTEM', 'CREATE TABLE'),
    privilege('SYSTEM', 'DROP ANY TABLE'),
    privilege('OBJECT', 'UPDATE', 'DGVA', 'VARIANT_STUDY'),
    privilege('OBJECT', 'DELETE', 'DGVA', 'VARIANT_PLACEMENT'),
    privilege('OBJECT', 'EXECUTE', 'DGVA', 'STAGE_TO_ARCHIVE_PKG'),
    privilege('OBJECT', 'INSERT', 'DGVA', 'VARIANT_STUDY', grantee='ROLE'),
    privilege('OWNER', 'SCHEMA OWNER', 'DGVA'),
])
def test_privilege_check_stops_on_a_write_privilege_of_the_account(write_privilege):
    connection = FakeConnection(privileges=READ_ONLY_PRIVILEGES + PUBLIC_PACKAGE_GRANTS + [write_privilege])
    with pytest.raises(WritePrivilegeError) as error:
        open_database(connection)
    assert write_privilege['privilege'] in str(error.value)
    assert len(connection.executed) == 1


def test_write_privileges_granted_to_public_are_only_reported(caplog):
    connection = FakeConnection(privileges=READ_ONLY_PRIVILEGES + PUBLIC_PACKAGE_GRANTS)
    with caplog.at_level('INFO', logger='DgvaDatabase'):
        open_database(connection)
    warnings = [record.getMessage() for record in caplog.records if record.levelname == 'WARNING']
    assert len(warnings) == 1
    assert 'EXECUTE on DGVA.STAGE_TO_ARCHIVE_PKG' in warnings[0] and 'PUBLIC' in warnings[0]


def test_inherit_privileges_is_not_a_write_privilege():
    assert find_write_privileges(READ_ONLY_PRIVILEGES) == ([], [])


def test_privilege_check_stops_when_it_returns_nothing():
    with pytest.raises(WritePrivilegeError):
        open_database(FakeConnection(privileges=[]))


def test_find_write_privileges_lists_each_once():
    rows = READ_ONLY_PRIVILEGES + PUBLIC_PACKAGE_GRANTS * 2 + [privilege('SYSTEM', 'UNLIMITED TABLESPACE')] * 2
    assert find_write_privileges(rows) == (
        ['SYSTEM UNLIMITED TABLESPACE'],
        ['OBJECT EXECUTE on DGVA.STAGE_TO_ARCHIVE_PKG', 'OBJECT EXECUTE on DGVA.STAGE_TO_ARCHIVE_VP_PKG'])


def test_connect_closes_the_connection_when_the_check_fails(monkeypatch, tmp_path):
    connection = FakeConnection(privileges=READ_ONLY_PRIVILEGES + [privilege('SYSTEM', 'DELETE ANY TABLE')])
    import oracledb
    monkeypatch.setattr(oracledb, 'connect', lambda **kwargs: connection)
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h, port: 1521, service_name: s, user: u, password: p}\n')
    load_config(str(config))
    with pytest.raises(WritePrivilegeError):
        DgvaDatabase.connect(database_config())
    assert connection.closed and connection.rollbacks == 1 and connection.commits == 0
    assert connection.autocommit is False


def test_config_keeps_the_password_out_of_repr(tmp_path, monkeypatch):
    monkeypatch.delenv('DGVA_PASSWORD', raising=False)
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h, port: 1521, service_name: s, user: u, password: secret}\n')
    load_config(str(config))
    loaded = database_config()
    assert loaded.password == 'secret' and 'secret' not in repr(loaded)
    assert loaded.schema == 'DGVA' and loaded.arraysize == 10000


def test_config_password_from_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv('DGVA_PASSWORD', 'from-env')
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h, port: 1521, service_name: s, user: u}\n')
    load_config(str(config))
    assert database_config().password == 'from-env'


def test_config_reports_missing_values(tmp_path, monkeypatch):
    monkeypatch.delenv('DGVA_PASSWORD', raising=False)
    config = tmp_path / 'config.yaml'
    config.write_text('DGVA: {host: h}\n')
    load_config(str(config))
    with pytest.raises(ConfigurationError) as error:
        database_config()
    assert 'port' in str(error.value) and 'password' in str(error.value)


def test_config_without_a_dgva_section(tmp_path):
    config = tmp_path / 'config.yaml'
    config.write_text('other: {}\n')
    load_config(str(config))
    with pytest.raises(ConfigurationError):
        database_config()


def test_config_found_through_the_environment(tmp_path, monkeypatch):
    config = tmp_path / 'from_env.yaml'
    config.write_text('DGVA: {host: env-host, port: 1521, service_name: s, user: u, password: p}\n')
    monkeypatch.setenv('DGVA_EXTRACTION_CONFIG', str(config))
    load_config(None)
    assert database_config().host == 'env-host'


def test_missing_config_file_is_reported(tmp_path, monkeypatch):
    monkeypatch.delenv('DGVA_EXTRACTION_CONFIG', raising=False)
    monkeypatch.setenv('HOME', str(tmp_path))
    with pytest.raises(FileNotFoundError):
        load_config(str(tmp_path / 'missing.yaml'))


def test_privilege_check_is_logged_under_the_class_name(caplog):
    with caplog.at_level('INFO', logger='DgvaDatabase'):
        open_database(FakeConnection())
    assert [(record.name, record.getMessage()) for record in caplog.records] == \
        [('DgvaDatabase', 'Privilege check passed: the session holds no write privilege of its own')]
