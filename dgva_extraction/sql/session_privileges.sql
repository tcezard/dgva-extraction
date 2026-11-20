-- Start-up check (sprint 3, "Database safety"): every privilege the current session holds that could matter for
-- writing to the DGVa schema. dgva_extraction.database.find_write_privileges() decides which ones are writes.
--   SYSTEM: every system privilege of the session, direct or through a role (SESSION_PRIVS).
--   OBJECT: privileges on objects owned by the DGVa schema, granted to the user, to PUBLIC, or to an enabled role.
--   OWNER:  the session user is the schema owner, who can always write its own tables.
-- grantee_kind says who holds the privilege: USER (this account), ROLE (one of its enabled roles), PUBLIC (every
-- account on the database) or SESSION (system privileges, which SESSION_PRIVS doesn't attribute). A write privilege
-- held by the account or its roles stops the run; one granted to PUBLIC is logged as a warning, since it is the same
-- for every account and doesn't mean the wrong credentials are configured.
-- Oracle-specific by nature (data dictionary views); a PostgreSQL copy needs its own check.
-- Binds: schema (the owner of the DGVa tables, upper case).
SELECT 'SYSTEM' AS privilege_kind, 'SESSION' AS grantee_kind, sp.privilege,
    CAST(NULL AS VARCHAR2(128)) AS owner, CAST(NULL AS VARCHAR2(128)) AS object_name
FROM session_privs sp
UNION ALL
SELECT 'OBJECT', CASE WHEN tp.grantee = 'PUBLIC' THEN 'PUBLIC' ELSE 'USER' END, tp.privilege, tp.table_schema,
    tp.table_name
FROM all_tab_privs tp
WHERE tp.table_schema = :schema
    AND (tp.grantee = USER OR tp.grantee = 'PUBLIC')
UNION ALL
SELECT 'OBJECT', 'ROLE', rtp.privilege, rtp.owner, rtp.table_name
FROM role_tab_privs rtp
    JOIN session_roles sr ON sr.role = rtp.role
WHERE rtp.owner = :schema
UNION ALL
SELECT 'OWNER', 'USER', 'SCHEMA OWNER', USER, CAST(NULL AS VARCHAR2(128))
FROM dual
WHERE USER = :schema
