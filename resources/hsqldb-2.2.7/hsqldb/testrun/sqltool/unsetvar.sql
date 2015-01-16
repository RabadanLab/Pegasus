/*
 * $Id: unsetvar.sql 4809 2011-11-20 21:12:35Z unsaved $
 */

\p *{:unsetvar}

*if (*x != *y)
    \q Two unset variables are not equal
*end if

*x =
*if (*x == *y)
    \q A variable set to '' is equal to an unset variable
*end if

*z =
*if (*x != *z)
    \q Two variables set to '' are not equal
*end if
