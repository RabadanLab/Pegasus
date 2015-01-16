/*
    $Id: load_binding_lu.sql 610 2008-12-22 15:54:18Z unsaved $
    Load BINDING Lookup Text Table
*/

\p Creating table BINDING_TMPTXT
CREATE TEMP TEXT TABLE binding_tmptxt (
    id integer,
    name varchar(12)
);

\p Setting text file source
SET TABLE binding_tmptxt SOURCE "binding_lu.ttbl;ignore_first=true;fs=\t"

\p rows in binding_tmptxt:
select count(*) from binding_tmptxt;
\p PRE rows in binding_lu:
select count(*) from binding_lu;

INSERT INTO binding_lu (
    id,
    name
) SELECT
    id,
    name
FROM BINDING_TMPTXT;

commit;

\p POST rows in binding_lu:
select count(*) from binding_lu;
