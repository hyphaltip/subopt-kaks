/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/* sqdconfig_main.c
 * SRE, Tue Mar  5 15:58:27 2002 [St. Louis]
 * 
 * Small C program designed to print out information on squid's 
 * compile-time configuration options - testsuite scripts can
 * call this to determine what optional stuff is compiled in. 
 * 
 * CVS $Id: sqdconfig_main.c,v 1.1 2003/07/03 13:20:42 jason Exp $
 */


#include <stdio.h>
#include <stdlib.h>
#include "squid.h"

int main(void)
{
#ifdef HAS_64BIT_FILE_OFFSETS
  printf("%-30s true\n", "HAS_64BIT_FILE_OFFSETS");
#else
  printf("%-30s false\n", "HAS_64BIT_FILE_OFFSETS");
#endif
}
