/*
 * $Id: atuser.sql 3714 2010-07-22 02:28:38Z unsaved $
 *
 * Test loading other files with @
 */

\i @/tblx.sql

\m @/tblx.dsv

SELECT COUNT(*) FROM tblx;
*if (*? != 2)
    \q Failed to load table deta from @ directory
*end if
