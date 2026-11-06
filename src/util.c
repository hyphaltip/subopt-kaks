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
  if( ! seq || seq->length <= 0 ) return 0;
  len = seq->length;
  newseq = (sequence *) calloc(1,sizeof(sequence));
  if( ! newseq ) return 0;
  newseq->seqstr = (char *)(calloc(len+1,sizeof (char)));
  if( ! newseq->seqstr ) { free(newseq); return 0; }
  newseq->length   = len;
  newseq->alphabet = seq->alphabet;
  newseq->seqname  = NULL;
  for( i =0; i < len; i++ ) {
    newseq->seqstr[i] = seq->seqstr[len-i-1];
  }
  newseq->seqstr[len] = '\0';
  return newseq;
}

void
free_sequence( sequence * seq )
{
  if( ! seq ) return;
  free(seq->seqstr);
  free(seq->seqname);
  free(seq);
}

int ** 
make_int_matrix  (int x, int y) 
{
  int i;
  int ** matrix;
  if( (matrix = (int **)calloc(x,sizeof(int *))) == NULL ) {
    Die("oom make int matrix\n");
  }
  for(i=0;i<x;i++) {
    matrix[i] = calloc(y,sizeof(int)); /* calloc zero-fills */
  }
  return matrix;
}

double ** 
make_double_matrix  (int x, int y) 
{
  int i;
  double ** matrix;
  if( (matrix = (double **)calloc(x,sizeof(double *))) == NULL ){
    Die("out of memory, make double matrix\n");
  }

  for(i=0;i<x;i++) {
    matrix[i] = calloc(y,sizeof(double));
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

#include <unistd.h> /* access */

/* malloc a copy of a NUL-terminated string. */
static char *
dupstr (const char *s)
{
  size_t n = strlen(s) + 1;
  char *p = (char *) malloc(n);
  if( p ) memcpy(p, s, n);
  return p;
}

/* Join a directory and a filename into a freshly malloc'd "dir/name". */
static char *
join_path (const char *dir, const char *name)
{
  size_t n = strlen(dir) + 1 + strlen(name) + 1;
  char *p = (char *) malloc(n);
  if( p ) snprintf(p, n, "%s/%s", dir, name);
  return p;
}

char *
resolve_matrix_file (const char *spec, const char *matdir)
{
  const char *name = (spec && *spec) ? spec : DEFAULT_MATRIX_NAME;
  const char *dirs[16];
  char root1[PATH_MAX], root2[PATH_MAX];
  char *env, *path;
  int nd = 0, i;

  /* An explicit path is taken verbatim so the caller can report a clear
   * "cannot open <path>" error if it is missing. */
  if( strchr(name, '/') != NULL )
    return dupstr(name);

  /* A bare name that already resolves in the current directory. */
  if( access(name, R_OK) == 0 )
    return dupstr(name);

  /* Otherwise search the matrix directories in priority order. */
  if( matdir && *matdir )
    dirs[nd++] = matdir;
  if( (env = getenv("SUBOPT_MATRIX_DIR")) != NULL && *env )
    dirs[nd++] = env;
  if( (env = getenv("SUBOPTDIR")) != NULL && *env ) {
    snprintf(root1, sizeof root1, "%s/matricies/aa", env);
    dirs[nd++] = root1;
  }
  if( (env = getenv("SUBOPT_DIR")) != NULL && *env ) {
    snprintf(root2, sizeof root2, "%s/matricies/aa", env);
    dirs[nd++] = root2;
  }
  dirs[nd++] = "matricies/aa";
  dirs[nd++] = "../matricies/aa";
  dirs[nd++] = "../../matricies/aa";
#ifdef SUBOPT_DATADIR
  dirs[nd++] = SUBOPT_DATADIR "/matricies/aa";
#endif

  for(i = 0; i < nd; i++) {
    path = join_path(dirs[i], name);
    if( path == NULL )
      return NULL;
    if( access(path, R_OK) == 0 )
      return path;
    free(path);
  }
  return NULL;
}
