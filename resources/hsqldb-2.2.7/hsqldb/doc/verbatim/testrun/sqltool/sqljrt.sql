/*
 * $Id: sqljrt.sql 3353 2009-12-15 19:52:13Z unsaved $
 *
 * Tests SQL/JRT
 */

create function dehex(VARCHAR(80), INTEGER)
    returns INTEGER
    no sql
    language java
    external name 'CLASSPATH:java.lang.Integer.valueOf'
.;

CALL dehex('12', 16);
*if (*? != 18)
    \q SQL/JRT function failed
*end if
