#include "align.h"
/* 
 * Author: Jason Stajich
 * $Id: zuker.c,v 1.22 2004/07/01 21:29:31 jason Exp $
 * 
 * Zuker suboptimal alignments
 * Based on algorithm described in
 * Zuker M (1991) "Suboptimal sequence alignment in molecular biology
 *   alignment with error analysis" J Mol Biol 221(2):403-420.
 *
 * This implementation is heavily derived from ideas in 
 * C++ libraries written by Michael Smoot at UVA.
 * 
 */

/* externed in util.h */

int Gapopen = -10;
int Gapext  = -2;
int ** ScoringMatrix;
float MatrixScale;
AlignType Alntype;
int Verbose = 0, rc;
int Upper_bound = INT_MAX, Lower_bound = 0;
int Max_alignments_to_report;
char GapChar = '-', *t;
int Alignment_cutoff_percent = 100;

int 
zuker_align ( const sequence * seq1, 
	      const sequence * seq2,
	      alignment *** ret_alignments)
{
  
  int count =0,i,j,k;
  int ** score_fmatrix, **score_rmatrix, **zuker_matrix;
  whence ** direction_fmatrix, **direction_rmatrix;  
  sequence * rev1, * rev2;
  int  seenscore = 0, possible_score_count =0, num_alignments =0;
  register int idx1,idx2,m,n;
  int optimalscore = -9999; /* arbitrary low number */
  alignment * pw_aln;
	
  m = seq1->length-1;
  n = seq2->length-1;
  
  score_fmatrix = make_int_matrix(seq1->length+1,
				  seq2->length+1); /* N+1 x M+1 for fwd */
  score_rmatrix = make_int_matrix(seq1->length+1,
				  seq2->length+1); /* N+1 x M+1 for rev */
  zuker_matrix = make_int_matrix(seq1->length,
				 seq2->length); /* N x M for sum */
  
  direction_fmatrix = make_whence_matrix(seq1->length+1,
					 seq2->length+1);
  direction_rmatrix = make_whence_matrix(seq1->length+1,
					 seq2->length+1);
  
  rev1 = reverse_seq(seq1);
  rev2 = reverse_seq(seq2);
  
  if( fill_matrix(seq1, seq2, &score_fmatrix, &direction_fmatrix) < 0 ) {
    return -1;
  }
  if(  fill_matrix(rev1, rev2, &score_rmatrix, &direction_rmatrix) < 0 ){
    return -1;
  }
  
  
  if( Verbose > 1 ) {
    print_matrix(score_fmatrix,seq1,seq2,1,"forward");
    print_matrix(score_rmatrix,rev1,rev2,1,"reverse");
  }

  /* Build the Zuker matrix which is fwd+rev */

  for(i=0;i < m;i++) {
    idx1 = toupper(seq1->seqstr[i]) - FIRSTCHAR;
    if( idx1 < 0 || idx1 > 25 ) idx1 = STOPAANUM;
    assert(idx1 >= 0);    
    for(j=0; j < n; j++) {
      /* The numbers are flipped because we're doing fwd+back
       * so are adding numbers from the reverse sequences
       */
      idx2 = toupper(seq2->seqstr[j]) - FIRSTCHAR;      
      if( idx2 < 0 || idx2 > 25 ) idx2 = STOPAANUM;
      assert(idx2 >= 0);    
      seenscore = zuker_matrix[i][j] = 
	ScoringMatrix[idx1][idx2] +
	score_fmatrix[i][j] + 
	score_rmatrix[m - i][n - j];      
      optimalscore = MAX(seenscore,optimalscore);      
    }
  }
  /* optimal score is the maximum score */

  Lower_bound = (int)(optimalscore * (1 - (double)Alignment_cutoff_percent / 100.00));
  if( Verbose > 0 ) {
    fprintf(stderr,"optimal score is %d Lower_bound %d cutoff_percent is %d\n",
	    optimalscore, Lower_bound, Alignment_cutoff_percent);
  }  
  if( Verbose > 1 ) {
    print_matrix(zuker_matrix,seq1,seq2,-1,"zuker");
  }


  /* Now we are going to create the alignments.  To do this we iterate
   * over all the cells in the zuker matrix, if we want to limit the
   * number of alignments we skip those which don't meet the
   * upper/lower boundary requirements
   */
  for( i =0; i < m; i++)
    for( j= 0; j < n; j++ )
      if( zuker_matrix[i][j] >= Lower_bound && 
	  zuker_matrix[i][j] <= Upper_bound )
	possible_score_count++;
  
  for( i =0; i < m; i++ ) 
    for( j = 0; j < n; j++ ) 
      if( zuker_matrix[i][j] >= Lower_bound && 
	  zuker_matrix[i][j] <= Upper_bound &&
	  num_alignments < Max_alignments_to_report) {
	
	if( (pw_aln = calloc( 1,sizeof(alignment))) == NULL ) {
	  fprintf(stderr,"oom pw_aln alloc\n");
	  return -1;
	}
	
	rc = create_pairwise_zuker_alignment(i,j,
				       seq1,seq2,
				       rev1,rev2,
				       zuker_matrix,
				       direction_fmatrix,
				       direction_rmatrix,
				       pw_aln);
	if( rc == 0 ) {
	  (*ret_alignments)[num_alignments++] = pw_aln;	 
	} else {
	  // cleanup_alignment(pw_aln);
	  return -1;
	}
	pw_aln = NULL;
      }

  free(rev1);
  free(rev2);
  cleanup_matrix( (void **)score_fmatrix,    seq1->length +1);
  cleanup_matrix( (void **)score_rmatrix,    seq1->length +1);
  //  fix me soon
  cleanup_matrix( (void **)direction_fmatrix,seq1->length +1);
  cleanup_matrix( (void **)direction_rmatrix,seq1->length +1);
  cleanup_matrix( (void **)zuker_matrix,     seq1->length);
  free(score_fmatrix);
  free(score_rmatrix);
  free(zuker_matrix);
  free(direction_fmatrix);
  free(direction_rmatrix);
  
  return num_alignments;
}

int 
optimal_align ( const sequence * seq1, 
		const sequence * seq2,
		alignment * ret_alignment)
{
  
  int count =0,i,j,k;
  int ** score_matrix;
  whence ** direction_matrix;  
  int  rc, best_loc;
  register int idx1,idx2,m,n;
   
  m = seq1->length-1;
  n = seq2->length-1;
  
  score_matrix = make_int_matrix(seq1->length+2,
				 seq2->length+2);
  
  direction_matrix = make_whence_matrix(seq1->length+2,
					seq2->length+2);
  
  if( (best_loc = fill_matrix(seq1, seq2, &score_matrix, &direction_matrix))
      < 0 ) {
    return -1;
  }
  
  if( Verbose > 1 ) {
    print_matrix(score_matrix,seq1,seq2,1,"forward");
  }

  /* Now we are going to create the alignment.
   */
  
  i = (int)(best_loc / seq1->length);
  j = seq2->length % best_loc;
  if( traceback(i,j,seq1,seq2,
		score_matrix,
		direction_matrix,
		ret_alignment) == 0 ) {
    rc = 1;
  } else {
    rc = -1;
  }
  
  
  cleanup_matrix( (void **)score_matrix,        seq1->length +1);
  cleanup_matrix( (void **)direction_matrix,    seq1->length +1);
  free(score_matrix);
  free(direction_matrix);
  return rc;
}


/* remember that Alntype, Gapext, Gapopen, ScoringMatrix, MatrixScale are all
 * global variables in this implementation
 * non-zero return status indicates failure
 *
 * Based in part on Michael Smoot's AlignmentGenerator.cpp _fill code
 */

int 
fill_matrix ( const sequence * seq1, 
	      const sequence * seq2, 
	      int *** score_matrix,
	      whence *** direction_matrix) 
{
  int i,j; /* i is always down, j is always across */
  register int idx1,idx2;
  int prev, best_loc = 0, best_score = -9999;
  int score, t, e, s, c, beginAcrossGap, tmpAcross, tmpDown,diag;
  int *D, *C, *beginDownGap;
  D = (int *) calloc(seq2->length+1, sizeof(int));
  
  if( D == NULL ) {
    fprintf(stderr,"OOM for D allocate %d\n",seq2->length);    
    goto error;
  }

  C = (int *) calloc(seq2->length+1, sizeof(int));
  if( C == NULL) {
    fprintf(stderr,"OOM for C allocate %d\n",seq2->length);    
    goto error;
  }
  
  if( (beginDownGap = (int *) calloc(seq2->length+1, sizeof(int))) == 0) {
    fprintf(stderr,"OOM for beginDownGap allocate %d\n",seq2->length);    
    goto error;
  }
  
  
  /* e is cost to extend gap
   * c is immedietely left value
   */
  t = e = s = c = beginAcrossGap = 0;

  /* Do the 'gutters' first */   
  (*score_matrix)[0][0] = C[0] = 0;
  (*direction_matrix)[0][0].direction = END_POSITION;
  (*direction_matrix)[0][0].across_index = 0;
  (*direction_matrix)[0][0].down_index = 0;

  t = Gapopen;
  /* do row[0] first */
  for(j = 1; j <= seq2->length; j++ ) {
    if( Alntype == local ) {
      (*direction_matrix)[0][j].direction = END_POSITION;
      (*direction_matrix)[0][j].across_index = 0;
      (*direction_matrix)[0][j].down_index = 0;
    } else if( Alntype == global ) {      
      (*score_matrix)[0][j] = C[j] = t = t + Gapext;
      D[j] = t + Gapopen;
      
      (*direction_matrix)[0][j].direction = EXTENDED_ACROSS;
      (*direction_matrix)[0][j].across_index = 0;
      (*direction_matrix)[0][j].down_index = 0;      
    }  else { 
      fprintf(stderr,"Don't know about the alignment type %d\n",Alntype);
      goto error;
      return -1;
    }
  }

  t = Gapopen;
  for(i = 1; i <= seq1->length; i++ ) {
    if( Alntype == local ) { 
      s = e = c = 0;
      (*score_matrix)[i][0] = 0;
      (*direction_matrix)[i][0].direction = END_POSITION;
      (*direction_matrix)[i][0].across_index = 0;
      (*direction_matrix)[i][0].down_index = 0;      
    } else if( Alntype == global) { 
      s = C[0];
      (*score_matrix)[i][0] = C[0] = c = t = t + Gapext;
      e = t + Gapopen;
      (*direction_matrix)[i][0].direction = EXTENDED_DOWN;      
      (*direction_matrix)[i][0].across_index = 0;
      (*direction_matrix)[i][0].down_index = 0;      
    } else { 
      fprintf(stderr,"Don't know about the alignment type %d\n",Alntype);
      return -1;
    }
    idx1 = toupper(seq1->seqstr[i-1]) - FIRSTCHAR;
    if( idx1 < 0 || idx1 > 25 ) idx1 = STOPAANUM;
    assert(idx1 >= 0);
    beginAcrossGap = 0;    
    for(j = 1; j <= seq2->length; j++ ) {
      idx2 = toupper(seq2->seqstr[j-1]) - FIRSTCHAR;
      if( idx2 < 0 || idx2 > 25 ) idx2 = STOPAANUM;
      
      assert(idx2 >= 0 );
      /* Find the best gap strategy for going ACROSS.
       * If the a gap is created, note where it begins. 
       * c is the value immediately to the left
       * e is the value of extending a gap across
       */
      tmpAcross = MAX(e, c + Gapopen) + Gapext;
      
      if( tmpAcross != (e + Gapext) &&
	  tmpAcross == (c + Gapopen + Gapext) ) {
	beginAcrossGap = j-1;
      }
      e = tmpAcross;
      
      /* Find the best gap strategy for going DOWN.
       * If the a gap is created, note where it begins. 
       * C[j] is the value immediately above
       * D[j] is the value of extending a gap down 
       */
      tmpDown = MAX(D[j], C[j] + Gapopen) + Gapext;
      
      if( tmpDown == (C[j] + Gapopen + Gapext) &&
	  tmpDown != (D[j] + Gapext) )
	beginDownGap[j] = i-1;
      
      D[j] = tmpDown;
      
      diag = s + ScoringMatrix[idx1][idx2];
      
      /* D[j] is the best down
       * e is the best across
       * diag is diag
       */
      c = MAX(D[j], MAX(e, diag) );
      
      /* local (Smith-Waterman)
       *              { 0
       *              { F(i-1,j-1) + s(xi,yi) 
       * F(i,j) = max { F(i-1,j) - gapopen - (#gaps-1 * gapext)
       *              { F(i,j-1) - gapopen - (#gaps-1 * gapext)
       *
       */
      if( Alntype == local ) 
	c = MAX (0, c);
      
      if( Verbose > 1 ) {	
	fprintf(stderr,"diag score is %d, c is %d, D[j] is %d e is %d s is %d C[j] is %d\n",
		diag,c,D[j],e,s,C[j]);
      }
      prev = 0;
      
      if( diag == c ) 
	prev += DIAG;
      
      if( D[j] == c ) 
      {
	if( beginDownGap[j] == i-1 ) {
	  prev += DOWN;
	} else {
	  prev += EXTENDED_DOWN;
	  (*direction_matrix)[i][j].down_index = beginDownGap[j];
	}	
      }
      if( e == c ) {
	if( beginAcrossGap == j-1 ) {
	  prev += ACROSS;
	} else {
	  prev += EXTENDED_ACROSS;
	  (*direction_matrix)[i][j].across_index = beginAcrossGap;
	}
      }
      
      if( Alntype == local && c == 0 )
	prev += END_POSITION;      
      
      (*direction_matrix)[i][j].direction = prev;
      
      s = C[j];
      C[j] = c;
      (*score_matrix)[i][j] = C[j];
      if( c >= best_score ) {
	best_loc = ((i-1) * seq1->length) + j;
	best_score = c;
      }
      
      if(Verbose > 1 ) {
	print_matrix((*score_matrix),seq1,seq2,1,"score");
      }
    }
  }

  free(D);
  free(C);
  free(beginDownGap);

  return best_loc;

  error:
  if( Verbose > 0 ) 
    fprintf(stderr,"Error in pairwise alignment construction\n");
  if( D) 
    free(D);
  if( C) 
    free(C);
  if( beginDownGap )
    free(beginDownGap);
  return -1;
}

/* This is heavily lifted (but ported to C) 
 * from Michael Smoot's C++ ZukerGenerator.cpp code
 */

int 
create_pairwise_zuker_alignment( int i, int j, 
				 const sequence * seq1, 
				 const sequence * seq2,
				 const sequence * rev_seq1, 
				 const sequence * rev_seq2,
				 int    ** zuker_matrix,
				 whence ** dir_fmatrix, 
				 whence ** dir_rmatrix, 
				 alignment * pw_aln)
{
  /* The maximum length the alignment could be */
  int maxlen =   seq1->length + seq2->length +1;
  char * seqstr[NUM_PW_SEQS], * rev_seqstr[NUM_PW_SEQS];
  
  int score = zuker_matrix[i][j];
  int seq_counter = 0;
  
  int bi = i - 1;
  int bj = j - 1;
  int diff,x,dir,k, begin1Index, begin2Index, end1Index, end2Index;

  if( ! pw_aln ) 
    return -1;
  for(k=0;k < NUM_PW_SEQS;k++) {
    if( (seqstr[k]     = (char *)calloc(maxlen, sizeof(char))) == NULL) {
      fprintf(stderr,"OOM for char in pw_aln allocate\n");
      return -1;
    }
    if( (rev_seqstr[k] = (char *)calloc(maxlen, sizeof(char))) == NULL) {
      fprintf(stderr,"OOM for char allocate in pw_aln create\n");      
      return -1;
    }
    
  }
  maxlen--;  
  /* initial starting point */

  /* now we go backwards so to do a prepend we have to
   * do the additions in rev order
   */
  rev_seqstr[0][maxlen - seq_counter]   = seq1->seqstr[i];
  rev_seqstr[1][maxlen - seq_counter++] = seq2->seqstr[j];  
  
  while( bi >= 0 || bj >= 0 ) {
    dir = dir_fmatrix[bi+1][bj+1].direction;
    assert(seq_counter < maxlen);
    if( dir & DIAG ) {
      rev_seqstr[0][maxlen - seq_counter]   = seq1->seqstr[bi--];
      rev_seqstr[1][maxlen - seq_counter++] = seq2->seqstr[bj--];
    } else if( dir & DOWN  ) {
      rev_seqstr[0][maxlen - seq_counter]   = seq1->seqstr[bi--];
      rev_seqstr[1][maxlen - seq_counter++] = '-';
    } else if( dir & ACROSS  ) {
      rev_seqstr[0][maxlen - seq_counter]   = '-';
      rev_seqstr[1][maxlen - seq_counter++] = seq2->seqstr[bj--];
    } else if( dir & EXTENDED_DOWN ) {
      /* we compressed an extended down into a single dir entry 
       * need to untangle.
       */
      diff = bi+1 - dir_fmatrix[bi+1][bj+1].down_index;
      for( x = 0; x < diff; x++ ) {
	rev_seqstr[0][maxlen - seq_counter]   = seq1->seqstr[bi--];
	rev_seqstr[1][maxlen - seq_counter++] = '-';
      }      
    } else if( dir & EXTENDED_ACROSS ) {
      /* we compressed an extended down into a single dir entry 
       * need to untangle.
       */
      diff = bj+1 - dir_fmatrix[bi+1][bj+1].across_index;
      for( x = 0; x < diff; x++ ) {
	rev_seqstr[0][maxlen - seq_counter]   = '-';
	rev_seqstr[1][maxlen - seq_counter++] = seq2->seqstr[bj--];
      }
    } else if ( dir & END_POSITION ) {
      if( Alntype == local )
	break;
      
      if( bj > 0 || bi > 0) {
	fprintf(stderr,
		"Something wrong with fwd direction at bi:%d bj:%d dir: %d", 
		bi,bj,dir);
	return -1;
      } else 
	break;
    } else {
      fprintf(stderr,"Bad fwd direction at %d,%d %d\n",i,j,dir);
      return -1;
    }
  }
  
  /* copy backwards to forwards */
  
  for( k = 0; k < seq_counter; k++ ) {
    seqstr[0][k] = rev_seqstr[0][maxlen - seq_counter + 1 + k];
    seqstr[1][k] = rev_seqstr[1][maxlen - seq_counter + 1 + k];
  }
  begin1Index = bi +1;
  begin2Index = bj +1;
  
  /* now complete the alignment forward */
  i = seq1->length - 2 - i;
  j = seq2->length - 2 - j;
  
  while( i >= 0 || j >= 0 ) {
    dir = dir_rmatrix[i+1][j+1].direction;
    assert(seq_counter < maxlen);
    if( dir & DIAG ) {
      seqstr[0][seq_counter]   = rev_seq1->seqstr[i--];
      seqstr[1][seq_counter++] = rev_seq2->seqstr[j--];
    } else if ( dir & DOWN ) {
      seqstr[0][seq_counter]   = rev_seq1->seqstr[i--];
      seqstr[1][seq_counter++] = '-';
    } else if (dir & ACROSS ) {
      seqstr[0][seq_counter]   = '-';
      seqstr[1][seq_counter++] = rev_seq2->seqstr[j--];    
    } else if( dir & EXTENDED_DOWN ) {
      diff = i+1 - dir_rmatrix[i+1][j+1].down_index;
      for ( x = 0; x < diff; x++ ) {
	seqstr[0][seq_counter]   = rev_seq1->seqstr[i--];
	seqstr[1][seq_counter++] = '-';
      }
    } else if( dir & EXTENDED_ACROSS ) {
      diff = j+1 - dir_rmatrix[i+1][j+1].across_index;
      for ( x = 0; x < diff; x++ ) {
	seqstr[0][seq_counter]   = '-';
	seqstr[1][seq_counter++] = rev_seq2->seqstr[j--];
      }
    } else if( dir & END_POSITION ) {
      if( Alntype == local )
	break;
      
      if( j > 0 || i > 0 ) {
	fprintf(stderr,"Something wrong with rev direction at bi:%d bj:%d dir: %d", bi,bj,dir);
	return -1;	
      } else 
	break;
    } else {
      fprintf(stderr,"Bad rev direction %d,%d\n",i,j);
      return -1;
    }  
  }
  end1Index = seq1->length - i;
  end2Index = seq2->length - j;
  assert(seq_counter < maxlen);
  seqstr[0][seq_counter] = seqstr[1][seq_counter] = '\0';  

  /* memory allocate */

  if( (pw_aln->msa = calloc(1,sizeof(MSA))) == NULL ) {
    fprintf(stderr,"oom pw_aln->msa alloc\n");
    return -1;    
  }  
  
  if( (pw_aln->msa->aseq    = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->aseq alloc\n");
    return -1;
  }
  if( (pw_aln->msa->sqname  = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->sqname alloc\n");
    return -1;
  }
  
  if( (pw_aln->msa->sqdesc  = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->sqdesc alloc\n");
    return -1;
  }
  
  for(k =0; k < NUM_PW_SEQS; k++ ) {
    if( (pw_aln->msa->sqdesc[k] = calloc( SQINFO_DESCLEN, sizeof(char))) 
	== NULL) {
      fprintf(stderr,"oom pw_aln->msa->sqdesc alloc\n");
      return -1;
    } 
    
    if( (pw_aln->msa->sqname[k] = calloc(SQINFO_NAMELEN,sizeof(char)))
	== NULL) {
      fprintf(stderr,"oom sqname[1] allocate\n");
      return -1;
    }  
  }
  
  pw_aln->sqinfo      = (SQINFO *)calloc( NUM_PW_SEQS,
					  sizeof(SQINFO));
  
  /* we'll make copy here to save memory - could take
   * more time, but may be worth it?
   */
  for(k = 0; k < NUM_PW_SEQS; k++ ) {
    pw_aln->msa->aseq[k] = (char *)calloc(seq_counter+1,sizeof(char));
    strncpy(pw_aln->msa->aseq[k], seqstr[k],seq_counter);
    /* free the sequences */
    if( seqstr[k] ) 
      free(seqstr[k]);  
    if( rev_seqstr[k]) 
      free(rev_seqstr[k]);
    seqstr[k] = 0;
  }

  strncpy(pw_aln->msa->sqname[0], seq1->seqname, SQINFO_NAMELEN);
  strncpy(pw_aln->msa->sqname[1], seq2->seqname, SQINFO_NAMELEN);

  pw_aln->msa->nseq      = NUM_PW_SEQS;
  pw_aln->msa->alen      = seq_counter;
  pw_aln->score           = score; /* alignment score */
  pw_aln->sqinfo[0].start = begin1Index;
  pw_aln->sqinfo[0].stop = end1Index;
  strncpy(pw_aln->sqinfo[0].name, seq1->seqname,SQINFO_NAMELEN);
  strncpy(pw_aln->sqinfo[1].name, seq2->seqname,SQINFO_NAMELEN);
  pw_aln->sqinfo[1].flags = pw_aln->sqinfo[0].flags = SQINFO_DESC | SQINFO_ID;
  sprintf(pw_aln->sqinfo[0].desc, "score=%d alignlen=%d start=%d stop=%d",
	  score,seq_counter,begin1Index+1,end1Index+1);
  strncpy(pw_aln->msa->sqdesc[0],pw_aln->sqinfo[0].desc,SQINFO_DESCLEN);
  sprintf(pw_aln->sqinfo[1].desc, "score=%d alignlen=%d start=%d stop=%d",
	  score,seq_counter,begin2Index+1,end2Index+1);
  strncpy(pw_aln->msa->sqdesc[1],pw_aln->sqinfo[1].desc,SQINFO_DESCLEN);

  pw_aln->sqinfo[0].len = seq1->length;
  pw_aln->sqinfo[1].len = seq2->length;

  pw_aln->sqinfo[1].start = begin2Index;
  pw_aln->sqinfo[1].stop = end2Index;

  return 0;  
}


/* This is heavily lifted (but ported to C) 
 * from Michael Smoot's C++ ZukerGenerator.cpp code
 
 * This function only builds an alignment from the forward matrix
 * and forward score
 * i,j must be the end of the alignment path,
 * this will be a traceback
 */

int 
traceback(  int i, int j,
	    const sequence * seq1, 
	    const sequence * seq2,
	    int    ** matrix,
	    whence ** dir_matrix, 
	    alignment * pw_aln)
{
  /* The maximum length the alignment could be */
  int maxlen =   seq1->length + seq2->length +1;
  char * seqstr[NUM_PW_SEQS];
  
  int seq_counter = 0;
  int diff,x,dir,k, begin1Index, begin2Index;
  
  int end1Index = i;
  int end2Index = j;
  
  int score = matrix[i][j];
  assert(i >= 0 && j >= 0);
  
  if( ! pw_aln ) 
    return -1;
  for(k=0;k < NUM_PW_SEQS;k++) {
    if( (seqstr[k]     = (char *)calloc(maxlen, sizeof(char))) == NULL) {
      fprintf(stderr,"OOM for char in pw_aln allocate\n");
      return -1;
    }
    
  }
  maxlen--;
  /* initial starting point */
  i--;
  j--;
  
  while( i >= 0 || j >= 0 ) {
    dir = dir_matrix[i+1][j+1].direction;
    assert(seq_counter < maxlen);
    if( dir & DIAG ) {
      seqstr[0][seq_counter]   = seq1->seqstr[i--];
      seqstr[1][seq_counter++] = seq2->seqstr[j--];
    } else if ( dir & DOWN ) {
      seqstr[0][seq_counter]   = seq1->seqstr[i--];
      seqstr[1][seq_counter++] = '-';
    } else if (dir & ACROSS ) {
      seqstr[0][seq_counter]   = '-';
      seqstr[1][seq_counter++] = seq2->seqstr[j--];    
    } else if( dir & EXTENDED_DOWN ) {
      diff = i+1 - dir_matrix[i+1][j+1].down_index;
      for ( x = 0; x < diff; x++ ) {
	seqstr[0][seq_counter]   = seq1->seqstr[i--];
	seqstr[1][seq_counter++] = '-';
      }
    } else if( dir & EXTENDED_ACROSS ) {
      diff = j+1 - dir_matrix[i+1][j+1].across_index;
      for ( x = 0; x < diff; x++ ) {
	seqstr[0][seq_counter]   = '-';
	seqstr[1][seq_counter++] = seq2->seqstr[j--];
      }
    } else if( dir & END_POSITION ) {
      if( Alntype == local )
	break;
      else {
	//if( j > 0 || i > 0 ) {
	//} else {
	  break;
	  //}
	
      }
      
    } else {
      fprintf(stderr,"Bad fwd direction at %d,%d %d\n",i,j,dir);
      return -1;
    }    
  }
  begin1Index = i+1;
  begin2Index = j+1;
  assert(seq_counter < maxlen);
  seqstr[0][seq_counter] = seqstr[1][seq_counter] = '\0';  

  
  /* memory allocate */

  if( (pw_aln->msa = calloc(1,sizeof(MSA))) == NULL ) {
    fprintf(stderr,"oom pw_aln->msa alloc\n");
    return -1;    
  }  
  
  if( (pw_aln->msa->aseq    = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->aseq alloc\n");
    return -1;
  }
  if( (pw_aln->msa->sqname  = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->sqname alloc\n");
    return -1;
  }
  
  if( (pw_aln->msa->sqdesc  = (char **)calloc( NUM_PW_SEQS, 
					       sizeof(char *))) == NULL){
    fprintf(stderr,"oom pw_aln->msa->sqdesc alloc\n");
    return -1;
  }
  
  for(k =0; k < NUM_PW_SEQS; k++ ) {
    if( (pw_aln->msa->sqdesc[k] = calloc( SQINFO_DESCLEN, sizeof(char))) 
	== NULL) {
      fprintf(stderr,"oom pw_aln->msa->sqdesc alloc\n");
      return -1;
    } 
    
    if( (pw_aln->msa->sqname[k] = calloc(SQINFO_NAMELEN,sizeof(char)))
	== NULL) {
      fprintf(stderr,"oom sqname[1] allocate\n");
      return -1;
    }
    pw_aln->msa->aseq[k] = (char *)calloc(seq_counter+1,sizeof(char));
  
    /* now flipflop */  
    for( i = 1; i <= seq_counter; i++ ) {
      pw_aln->msa->aseq[k][i-1] = seqstr[k][seq_counter-i];
    }    
    
    /* free the sequences */
    if( seqstr[k] ) 
      free(seqstr[k]);  
    seqstr[k] = 0;
  }
  
  pw_aln->sqinfo      = (SQINFO *)calloc( NUM_PW_SEQS,
					  sizeof(SQINFO));  

  strncpy(pw_aln->msa->sqname[0], seq1->seqname, SQINFO_NAMELEN);
  strncpy(pw_aln->msa->sqname[1], seq2->seqname, SQINFO_NAMELEN);

  pw_aln->msa->nseq      = NUM_PW_SEQS;
  pw_aln->msa->alen      = seq_counter;
  pw_aln->score           = score; /* alignment score */
  pw_aln->sqinfo[0].start = begin1Index;
  pw_aln->sqinfo[0].stop = end1Index;
  strncpy(pw_aln->sqinfo[0].name, seq1->seqname,SQINFO_NAMELEN);
  strncpy(pw_aln->sqinfo[1].name, seq2->seqname,SQINFO_NAMELEN);
  pw_aln->sqinfo[1].flags = pw_aln->sqinfo[0].flags = SQINFO_DESC | SQINFO_ID;
  sprintf(pw_aln->sqinfo[0].desc, "score=%d alignlen=%d start=%d stop=%d",
	  score,seq_counter,begin1Index+1,end1Index+1);
  strncpy(pw_aln->msa->sqdesc[0],pw_aln->sqinfo[0].desc,SQINFO_DESCLEN);
  sprintf(pw_aln->sqinfo[1].desc, "score=%d alignlen=%d start=%d stop=%d",
	  score,seq_counter,begin2Index+1,end2Index+1);
  strncpy(pw_aln->msa->sqdesc[1],pw_aln->sqinfo[1].desc,SQINFO_DESCLEN);

  pw_aln->sqinfo[0].len = seq1->length;
  pw_aln->sqinfo[1].len = seq2->length;

  pw_aln->sqinfo[1].start = begin2Index;
  pw_aln->sqinfo[1].stop = end2Index;

  return 0;  
}
