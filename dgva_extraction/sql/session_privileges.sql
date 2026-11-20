-- Start-up check (sprint 3, "Database safety"): every privilege the current session holds that could matter for
-- writing to the DGVa schema. dgva_extraction.database.find_write_privileges() decides which ones are writes; any
-- write privilege stops the run.
--   SYSTEM: every system privilege of the session, direct or through a role (SESSION_PRIVS).
--   OBJECT: privileges on objects owned by the DGVa schema, granted to the user, to PUBLIC, or to an enabled role.
--   OWNER:  the session user is the schema owner, who can always write its own tables.
-- Oracle-specific by nature (data dictionary views); a PostgreSQL copy needs its own check.
-- Binds: schema (the owner of the DGVa tables, upper case).
SELECT 'SYSTEM' AS privilege_kind, sp.privilege, CAST(NULL AS VARCHAR2(128)) AS owner,
    CAST(NULL AS VARCHAR2(128)) AS object_name
FROM session_privs sp
UNION ALL
SELECT 'OBJECT', tp.privilege, tp.table_schema, tp.table_name
FROM all_tab_privs tp
WHERE tp.table_schema = :schema
    AND (tp.grantee = USER OR tp.grantee = 'PUBLIC')
UNION ALL
SELECT 'OBJECT', rtp.privilege, rtp.owner, rtp.table_name
FROM role_tab_privs rtp
    JOIN session_roles sr ON sr.role = rtp.role
WHERE rtp.owner = :schema
UNION ALL
SELECT 'OWNER', 'SCHEMA OWNER', USER, CAST(NULL AS VARCHAR2(128))
FROM dual
WHERE USER = :schema
