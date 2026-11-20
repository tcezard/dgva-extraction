"""Connection to the DGVa database.

dgva-extraction runs only with a SELECT-only account (sprint 3, "Database safety"). The account is what blocks writes
and DDL; this module adds the guards that catch the wrong credentials being configured:

- at start, it checks that the session holds no write privilege, and stops if it finds one;
- autocommit is off, and nothing here ever commits; the session is rolled back before it is closed;
- SQL comes only from the .sql files shipped in the package (dgva_extraction/sql), and values are passed as bind
  variables, never formatted into the SQL text;
- rows are streamed with a large cursor arraysize, never fetched all at once.

A database error is never swallowed: it propagates and fails the study.
"""
import re
from importlib import resources

from ebi_eva_common_pyutils.logger import AppLogger

from dgva_extraction.config import DEFAULT_ARRAYSIZE, DEFAULT_SCHEMA


SQL_PACKAGE = 'dgva_extraction'
SQL_DIRECTORY = 'sql'
_SQL_NAME = re.compile(r'^[a-z0-9_]+$')

# System privileges that cannot change data or schema objects. Any other system privilege stops the run.
READ_ONLY_SYSTEM_PRIVILEGES = frozenset({
    'CREATE SESSION', 'ALTER SESSION', 'SET CONTAINER',
    'SELECT ANY TABLE', 'READ ANY TABLE', 'SELECT ANY DICTIONARY', 'SELECT ANY SEQUENCE',
})
# Object privileges on the DGVa schema's objects that cannot change data. INHERIT PRIVILEGES on a user is granted to
# PUBLIC by default for every Oracle user: it concerns the privileges that user lends to invoker's-rights code it
# calls, and gives this session no way to write.
READ_ONLY_OBJECT_PRIVILEGES = frozenset({'SELECT', 'READ', 'INHERIT PRIVILEGES', 'INHERIT REMOTE PRIVILEGES'})
PRIVILEGE_QUERY = 'session_privileges'


class WritePrivilegeError(RuntimeError):
    """The session can write. Nothing has been read; the run must stop."""


def load_sql(name):
    """Return the text of dgva_extraction/sql/<name>.sql."""
    if not _SQL_NAME.match(name):
        raise ValueError(f'Invalid SQL file name: {name!r}')
    return resources.files(SQL_PACKAGE).joinpath(SQL_DIRECTORY, f'{name}.sql').read_text()


def find_write_privileges(privilege_rows):
    """Split the write privileges among the rows of session_privileges.sql into those that stop the run and those
    granted to PUBLIC, which are only reported. Returns (blocking, public), each a sorted list of descriptions.

    Each row has privilege_kind (SYSTEM, OBJECT or OWNER), grantee_kind (USER, ROLE, PUBLIC or SESSION), privilege,
    owner and object_name.
    """
    blocking, public = set(), set()
    for row in privilege_rows:
        kind = row['privilege_kind']
        privilege = (row['privilege'] or '').upper()
        if kind == 'SYSTEM' and privilege in READ_ONLY_SYSTEM_PRIVILEGES:
            continue
        if kind == 'OBJECT' and privilege in READ_ONLY_OBJECT_PRIVILEGES:
            continue
        # OWNER rows, unknown kinds and every other privilege count as a write privilege.
        target = f' on {row.get("owner")}.{row.get("object_name")}' if row.get('object_name') else ''
        description = f'{kind} {privilege}{target}'
        if kind == 'OBJECT' and row.get('grantee_kind') == 'PUBLIC':
            public.add(description)
        else:
            blocking.add(description)
    return sorted(blocking), sorted(public)


class DgvaDatabase(AppLogger):
    """A read-only session on the DGVa schema. Use it as a context manager."""

    def __init__(self, connection, schema=DEFAULT_SCHEMA, arraysize=DEFAULT_ARRAYSIZE):
        if getattr(connection, 'autocommit', False):
            raise ValueError('The connection has autocommit on; dgva-extraction never runs with autocommit')
        self.connection = connection
        self.schema = schema
        self.arraysize = arraysize
        self._checked = False

    @classmethod
    def connect(cls, config):
        """Open a connection from a DatabaseConfig and check that it cannot write."""
        import oracledb
        # CLOB columns (sequences) are returned as str, so rows can be streamed without LOB locators.
        oracledb.defaults.fetch_lobs = False
        connection = oracledb.connect(
            user=config.user, password=config.password,
            dsn=oracledb.makedsn(config.host, config.port, service_name=config.service_name),
        )
        connection.autocommit = False
        database = cls(connection, schema=config.schema, arraysize=config.arraysize)
        try:
            database.check_no_write_privilege()
        except BaseException:
            database.close()
            raise
        return database

    def check_no_write_privilege(self):
        """Stop if the session holds any privilege that can change data or objects."""
        rows = list(self._execute(load_sql(PRIVILEGE_QUERY), {'schema': self.schema.upper()}))
        if not rows:
            # A session always holds at least CREATE SESSION: an empty answer means the check did not run.
            raise WritePrivilegeError('The privilege check returned no rows; refusing to continue')
        write_privileges, public_write_privileges = find_write_privileges(rows)
        if write_privileges:
            raise WritePrivilegeError(
                'The DGVa session can write, so dgva-extraction will not run with these credentials. '
                'Use the SELECT-only account. Write privileges found: ' + '; '.join(write_privileges))
        if public_write_privileges:
            # Granted to every account on the database, not to this one: dgva-extraction calls no PL/SQL, and its
            # SQL is checked for calls to DGVa functions and packages (tests/test_sql.py). Ask the DBA to revoke them.
            self.warning('PUBLIC holds privileges that can change DGVa data (any account can use them; '
                         'dgva-extraction never does): ' + '; '.join(public_write_privileges))
        self._checked = True
        self.info('Privilege check passed: the session holds no write privilege of its own')

    def stream(self, sql_name, binds):
        """Run dgva_extraction/sql/<sql_name>.sql with bind variables and yield rows as dicts, lower-case keys."""
        if not self._checked:
            raise WritePrivilegeError('check_no_write_privilege() must pass before any query runs')
        return self._execute(load_sql(sql_name), binds)

    def _execute(self, sql, binds):
        if not isinstance(binds, dict):
            raise TypeError('Bind variables must be passed as a dict')
        for name, value in binds.items():
            if not isinstance(value, (str, int, float, type(None))):
                raise TypeError(f'Bind variable {name!r} must be a scalar, not {type(value).__name__}')
        cursor = self.connection.cursor()
        try:
            cursor.arraysize = self.arraysize
            cursor.prefetchrows = self.arraysize + 1
            cursor.execute(sql, binds)
            columns = [description[0].lower() for description in cursor.description]
            while True:
                rows = cursor.fetchmany(self.arraysize)
                if not rows:
                    break
                for row in rows:
                    yield dict(zip(columns, row))
        finally:
            cursor.close()

    def close(self):
        """Roll back (nothing should be pending) and close. Never commits."""
        if self.connection is None:
            return
        try:
            self.connection.rollback()
        finally:
            self.connection.close()
            self.connection = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
