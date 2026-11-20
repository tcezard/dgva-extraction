# SQL portability: Oracle to the PostgreSQL cold-storage copy (sprint 3, ticket 13.4)

The reader's SQL (`dgva_extraction/sql/`) will later point at the PostgreSQL cold-storage copy (roadmap, Q3). The SQL
uses ANSI joins, common table expressions, `CASE`, `COALESCE` and `NULLIF` only. `tests/test_sql.py` fails if an
Oracle-only construct appears outside the list below.

## Oracle-specific SQL that would need changing

| File | Construct | PostgreSQL change |
|---|---|---|
| `session_privileges.sql` | `SESSION_PRIVS`, `SESSION_ROLES`, `ALL_TAB_PRIVS`, `ROLE_TAB_PRIVS`, `USER`, `DUAL`, `CAST(… AS VARCHAR2(128))` | A different check altogether: `has_table_privilege(current_user, …, 'INSERT,UPDATE,DELETE,TRUNCATE')` over the schema's tables, `pg_roles.rolsuper`, and table ownership in `pg_tables`. The check is needed whichever database is read. |

## Behaviour that differs without a change to the SQL

| What | Oracle | PostgreSQL | Effect on dgva-extraction |
|---|---|---|---|
| Empty strings | `''` is `NULL` | `''` is a value | `NULLIF(chromosome_name, 'NOT_GIVEN')` treats `''` as a name. Oracle never stores `''`, so a copy made with `ora2pg` holds `NULL`s. Check the copy. |
| `ORDER BY` on text | `NLS_SORT` (binary by default) | the column's collation | Contigs may come in another order. The writer copes: it writes each contig to its own file and joins them in assembly-report order. Positions are numbers and sort the same way. |
| `NULL`s in `ORDER BY … ASC` | last | last | Same. |
| Schema prefix `dgva.` | the owner | a schema | Keep the schema named `dgva` in the copy, or set `search_path`. |
| Driver | python-oracledb, `:name` binds | psycopg, `%(name)s` binds | `database.py` would translate `:name` to the driver's style. The `.sql` files keep `:name`. |
| `NUMBER` columns | `int` or `float` in Python | `numeric` comes back as `Decimal` | The reader converts coordinates with `int()`. |
