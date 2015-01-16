/*
 * $Id: sqljrt.sql 3340 2009-12-14 00:00:49Z unsaved $
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
