/*
 * $Id: dsv-altspace.sql 3340 2009-12-14 00:00:49Z unsaved $
 *
 * Tests importing to table in non-default tablespace
 */

CREATE SCHEMA altspace authorization dba;
CREATE TABLE altspace.targtbl (i INT, vc VARCHAR(80));

\m altspace.targtbl.dsv
SELECT COUNT(*) FROM altspace.targtbl;

*if (*? != 2)
    \q Import to alternate tablespace failed
*end if
