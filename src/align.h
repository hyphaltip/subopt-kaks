#ifndef _ALIGN_H_
#define _ALIGN_H

/* $Id: align.h,v 1.8 2004/07/01 21:29:30 jason Exp $ */

#include "util.h"

/* function prototypes */

int        zuker_align ( const sequence * seq1, 
			 const sequence * seq2, 
			 alignment *** ret_alignments);
int        fill_matrix ( const sequence * seq1, const sequence * seq2,
			 int *** score_matrix,			 
			 whence *** direction_matrix);
int        create_pairwise_zuker_alignment( int i, int j,
				      const sequence * seq1, 
				      const sequence * seq2,
				      const sequence * rev_seq1, 
				      const sequence * rev_seq2,
				      int    ** zuker_matrix,
				      whence ** direction_fmatrix,
				      whence ** direction_rmatrix,
				      alignment * pw_aln);

int        traceback( int i, int j,
		      const sequence * seq1, 
		      const sequence * seq2,
		      int    ** matrix,
		      whence ** direction_matrix,
		      alignment * pw_aln);

int        optimal_align ( const sequence * seq1, 
			   const sequence * seq2,
			   alignment * ret_alignment);

#endif /*_ALIGN_H_ defined */
