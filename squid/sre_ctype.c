/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/* sre_ctype.c
 * 
 * For portability. Some systems have functions tolower, toupper
 * as macros (for instance, MIPS M-2000 RISC/os!)
 * 
 * CVS $Id: sre_ctype.c,v 1.2 2003/10/07 17:00:17 jason Exp $
 */

#include "squidconf.h"

#include <ctype.h>
#include "squid.h"

int
sre_tolower(int c)
{
  if (isupper(c)) return tolower(c);
  else return c;
}

int
sre_toupper(int c)
{
  if (islower(c)) return toupper(c);
  else return c;
}

