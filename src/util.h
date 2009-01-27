#ifndef _UTIL_H_
#define _UTIL_H_

/* $Id: util.h,v 1.8 2004/07/01 21:29:31 jason Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"
/* other includes above ^^^ */

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define FIRSTCHAR 'A'
#define NUM_PW_SEQS 2
#define STOPAANUM 26
/* Enums */
enum AlignType_e { 
  local, 
  global,
};
typedef enum AlignType_e AlignType;

enum Alphabet_e  { 
  protein, 
  dna, 
  rna,
};

typedef enum Alphabet_e Alphabet;

/* structs */
struct sequence_s
{
  char *    seqstr;
  char *    seqname;
  int       length;
  Alphabet  alphabet;
};
typedef struct sequence_s sequence;


/* this is where we'll keep track of an alignment */
typedef struct alignment_s 
{
  MSA    * msa;
  int    score;       /* alignment score */
  SQINFO *sqinfo;     /* name, id, coord info for each sequence   */
} alignment;


/* Based on ideas from Michael Smoot's code */
typedef struct whence_s
{
  int across_index;
  int down_index;
  int direction; /* direction bit mask */
} whence;

/* for storing the Ka and Ks estimates */

typedef struct pairwise_distances_s 
{
  /* these will are all matricies */
  double ** N;    /* # of nonsynonymous changes */
  double ** dN;   /* nonsynonymous rate */
  double ** S;    /* # of synonomous changes */
  double ** dS;   /* synonymous rate */
  double ** dNdS;   /* Nonsynonymous / synonymous ratio */  
  double ** SEdS; /* standard error for dS */
  double ** SEdN; /* standard error for dN */
  double ** kappa; /* kappa (transversion rate) */ 
  double ** t;    /* t estimates */
} pairwise_distances;


extern  int **     ScoringMatrix;
extern  float      MatrixScale;
extern  int        Gapopen;
extern  int        Gapext;
extern  AlignType  Alntype;
extern  int        Lower_bound;
extern  int        Upper_bound;
extern  int        Max_alignments_to_report;
extern  int        Alignment_cutoff_percent; /* only report alignments 
						which are within N percent of 
						thhe maximum */
extern  int        Verbose;
extern  char       GapChar;

//extern  int        Min_align_len;


/* prototypes */

int    **  make_int_matrix(int x,int y);
double **  make_double_matrix(int x,int y);
whence **  make_whence_matrix(int x,int y);

void       print_matrix( int ** matrix,  
			 const sequence *, const sequence *, 
			 int type, char * label );
void       cleanup_matrix(void ** m,int x);
sequence * reverse_seq( const sequence *);
void       cleanup_alignment(alignment *);

/* direction */
extern const int DIAG, DOWN, ACROSS,EXTENDED_DOWN,EXTENDED_ACROSS,END_POSITION;

#endif /* _UTIL_H_ defined */
