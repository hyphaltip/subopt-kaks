/*****************************************************************
 * @LICENSE@
 *****************************************************************/

#ifndef SQUID_GKI_INCLUDED
#define SQUID_GKI_INCLUDED

/* gki.h
 * SRE, Sat May  1 15:07:22 1999
 * 
 * Declarations of structures, functions for generic key index
 * module: emulation of Perl hashes. See gki.c.
 * 
 * RCS $Id: gki.h,v 1.2 2003/10/07 17:00:17 jason Exp $
 */

/* gki_elem:
 *    key, array index pairs are kept in linked list structures.
 */
struct gki_elem {
  char            *key;
  int              idx;
  struct gki_elem *nxt;
};

/* gki:
 *    a dynamically resized hash structure; 
 *    contains a hash table and associated data
 */
typedef struct {
  struct gki_elem **table;
  
  int primelevel;
  int nhash;
  int nkeys;
} GKI;

GKI *GKIInit(void);
void GKIFree(GKI *hash);
int  GKIHashValue(GKI *hash, char *key);
int  GKIStoreKey(GKI *hash, char *key);
int  GKIKeyIndex(GKI *hash, char *key);
void GKIStatus(GKI *hash);

#endif /* SQUID_GKI_INCLUDED */
