#include "util.h"

/* $Id: util.c,v 1.9 2004/02/27 15:30:41 jason Exp $ */

/**
 * Diagonal (match) constant.
 * A value used for the directional bitmask in the Whence objects.
 */
const int DIAG = 1;

/**
 * Create gap down constant.
 * A value used for the directional bitmask in the Whence objects.
 */
const int DOWN = 2;

/**
 * Create gap across constant.
 * A value used for the directional bitmask in the Whence objects.
 */
const int ACROSS = 4;

/**
 * Gap extended down constant.
 * A value used for the directional bitmask in the Whence objects.
 */
const int EXTENDED_DOWN = 8;

/**
 * Gap extended across constant.
 * A value used for the directional bitmask in the Whence objects.
 */
const int EXTENDED_ACROSS = 16;

/**
 * End position constant.
 * A value used for the directional bitmask in the Whence objects.
 * End position is only to signal the end of local alignments.
 */
const int END_POSITION = 32;

void 
cleanup_alignment( alignment * p) 
{
  int i =0;
  if( ! p ) return;
  if( p->msa ) {    
    for( i = 0; i < p->msa->nseq; i++ ) {
      if( p->msa->aseq && 
	  p->msa->aseq[i] ) free(p->msa->aseq[i]);
      if( p->msa->sqname[i]) free(p->msa->sqname[i]) ;
      if( p->msa->sqdesc[i]) free(p->msa->sqdesc[i]) ;
      p->msa->aseq[i] = p->msa->sqdesc[i] = p->msa->sqname[i] = NULL;
    }
      if( p->msa->aseq)   free(p->msa->aseq);
      if( p->msa->sqname) free(p->msa->sqname);
      if( p->msa->sqdesc) free(p->msa->sqdesc);
      p->msa->sqname = p->msa->sqdesc = p->msa->aseq = NULL;
    
    free(p->msa);
  }  
  p->msa = NULL;
  if( p->sqinfo)  free(p->sqinfo)  ;
  p->sqinfo = NULL;
  free(p);
}

sequence * 
reverse_seq( const sequence * seq ) 
{
  int len,i;
  sequence * newseq;
  if( ! seq ) return 0;  
  len = seq->length;
  newseq = (sequence *) calloc(1,sizeof(sequence));
  if( ! len ) return 0;
  
  newseq->seqstr = (char *)(calloc(len+1,sizeof (char)));
  newseq->length = len;
  for( i =0; i < len; i++ ) {
    newseq->seqstr[i] = seq->seqstr[len-i-1];
  }
  newseq->seqstr[len] = '\0';
  return newseq;  
}

int ** 
make_int_matrix  (int x, int y) 
{
  int i,j;
  int ** matrix;
  if( (matrix = (int **)calloc(x,sizeof(int *))) == NULL ) {
    Die("oom make int matrix\n");
  }
  for(i=0;i<x;i++) {
    matrix[i] = calloc(y,sizeof(int));
    for(j=0;j<y;j++) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

double ** 
make_double_matrix  (int x, int y) 
{
  int i,j;
  double ** matrix;
  if( (matrix = (double **)calloc(x,sizeof(double *))) == NULL ){
    Die("out of memory, make double matrix\n");
  }
  
  for(i=0;i<x;i++) {
    matrix[i] = calloc(y,sizeof(double));
/*    for(j=0;j<y;j++) {
      matrix[i][j] = 0;
    }
*/
  }
  return matrix;
}

whence ** 
make_whence_matrix  (int x, int y) 
{
  int i,j;
  whence ** matrix;
  if( (matrix = (whence **) calloc(x,sizeof(whence *))) == NULL ) {
    Die("oom make whence matrix\n");    
  }
  
  for(i=0;i<x;i++) {
    matrix[i] = calloc(y,sizeof(whence));
    for(j=0;j<y;j++) {
      matrix[i][j].across_index = -1;
      matrix[i][j].down_index = -1;
      matrix[i][j].direction = END_POSITION;      
    }
  }
  return matrix;
}

void print_matrix ( int ** matrix, 
		    const sequence * seq1, const sequence * seq2, 
		    int type, char* label ) {
  int i,j;
  int x = seq1->length;
  int y = seq2->length;
  if( !label )
    label = "matrix";

  if( type < 0 ) {
    x--;
    y--;
  } 
  
  
  printf( "%-9s",label);
    
  for(i=0; i <= y; i++ ) { /* matrix header */
    if( type < 0  ) { 
      if( i < seq2->length ) 
	printf("%5c ",seq2->seqstr[i]);
      else 
	printf("%5s ", "");
    } else if ( i > 0 ) {
      printf("%5c ",seq2->seqstr[i-1]);
    } else { 
      printf("%5s ", "");
    }
  }
  printf("\n%10s","");
  for(i = 0; i <= y; i++ ) { /* number header */
    printf("[%4d]",i);
  }
  printf("\n");
  for(i = 0; i <= x; i++ ){
    if( type < 0 ) {
      if( i < seq1->length )
	printf("[%4d] %c: ",i,seq1->seqstr[i]);
      else 
	printf("[%4d]    ",i);
    } else if( i > 0 ) { 
      printf("[%4d] %c: ",i,seq1->seqstr[i-1]);    
    } else { 
      printf ("[%4d]    ",i);
    }
    for( j = 0; j <= y; j++ ) {
      printf("%5d ", matrix[i][j]); 
    }
    printf("\n");
  }
  printf("\n");
}

void 
cleanup_matrix (void ** m, int x) 
{  
  int i;
  if( ! m ) return;  
  for(i = 0; i < x; i++ ) {
    if( m[i] ) {
      free(m[i]);
    }
    
  }
}
