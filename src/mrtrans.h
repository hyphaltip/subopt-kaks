#ifndef _MRTRANS_H_
#define _MRTRANS_H_

/* $Id: mrtrans.h,v 1.4 2004/05/27 02:35:24 jason Exp $ */

/* This header is for the mrtrans.c code which is used to 
   translate from protein to cds alignments when the cds sequences are known.
*/
#include "util.h"

int mrtrans(alignment * protein_align, sequence * cds_seqs, 
	    alignment ** cds_aln,
	    int remove_last_stop);

#define CODON_LENGTH 3
#define STOPCODON1 "TGA"
#define STOPCODON2 "TAA"
#define STOPCODON3 "TAG"

#endif /* _MRTRANS_H_ defined */
