/* Author: Jason Stajich <jason.stajichATduke.edu>
 * $Id: yn00_cds_prealigned.c,v 1.3 2004/05/27 21:39:23 jason Exp $
 *
 */

/* The SUBOPTDIR env variable can be used to set the location of
 * of the matricies
 */

#include "align.h"
#include "mrtrans.h"

struct opt_s OPTIONS[] = {
  { "-h",              TRUE, sqdARG_NONE },  
  { "--informat",      FALSE, sqdARG_STRING },
  { "--noheader",      FALSE, sqdARG_NONE},
  { "--output",        FALSE, sqdARG_STRING },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\nUsage: yn00_cds_prealigned [-options] <sequence file>\n\n   Generate the optimal alignments for a pair of CDS sequences using the protein sequence as a template.\n   Available options:\n   -h      : help; print version and usage info\n   -v      : verbose output -- show alignment path states.\n ";
static char experts[] = "\n   --informat    : Format of the input sequence file alignment.\n   --output      : specify an output file instead or writing to STDOUT\n   --noheader    : Do not show the header for the table summary\n\n";

main (int argc, char ** argv ) 
{
  char     *seqfile;            /* name of MSA sequence file     */
  MSAFILE  *msafile;		/* open MSA sequence file        */
  MSA      *msa;                /* store the MSA */
  int       fmt,ofmt=106;	/* format of seqfile         */
  char  *optname;
  char  *optarg, *t;
  int    optind;
  int    seqct = 0;
  int    do_oneline       = 0;
  char   * output_filename = 0;
  int    showheader=1;
  FILE  *ofd, *fd;

  int    i,j,k;
  pairwise_distances pwMLdist, pwNGdist;

  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */

  /* for our purposes this is only pairwise alignments, but
   * would rather do it correctly in case we move to MSA case 
   */
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
  {
    if (strcmp(optname, "--informat") == 0) {
      fmt = String2SeqfileFormat(optarg);
      if (fmt == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    } else if (strcmp(optname, "--outformat") == 0) {
      ofmt = String2SeqfileFormat(optarg);
      if (ofmt == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    } else if (strcmp(optname, "-h") == 0) {
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }  else if(  strcmp(optname, "--output") == 0 ) {
      output_filename = optarg;	  
    } else if( strcmp(optname, "--noheader" ) == 0 ) {
      showheader = 0;
    }      
  }
  if (argc - optind < 1) Die("%s\n", usage);
  
  if( output_filename && strlen(output_filename) != 1 &&
      output_filename[0] != '-') {      
    ofd = fopen(output_filename,"w");
    if( ! ofd ) {
      Die( "could not open file %s",output_filename);
    }
  } else 
    ofd = stdout;
  
  while( optind < argc ) {
    seqfile = argv[optind++];
    /* for starters, lets do this dumb way and read the whole file
       in twice
    */
    if ((msafile = MSAFileOpen(seqfile, fmt, NULL)) == NULL)
      Die("Failed to open sequence file %s for reading", seqfile);

    if( showheader ) {
      fprintf(ofd,"SEQ1\tSEQ2\tdN\tdS\tOMEGA\tN\tS\tkappa\tt\tLENGTH\n");
    }    

    while (msa = MSAFileRead(msafile) )
    {
      if( ! msa->nseq ) continue;      
	
      seqct = msa->nseq;
      pwMLdist.N    = make_double_matrix(seqct,seqct);
      pwMLdist.dN   = make_double_matrix(seqct,seqct);
      pwMLdist.S    = make_double_matrix(seqct,seqct);
      pwMLdist.dS   = make_double_matrix(seqct,seqct);
      pwMLdist.dNdS = make_double_matrix(seqct,seqct);
      pwMLdist.SEdS = make_double_matrix(seqct,seqct);
      pwMLdist.SEdN = make_double_matrix(seqct,seqct);
      pwMLdist.t    = make_double_matrix(seqct,seqct);
      pwMLdist.kappa= make_double_matrix(seqct,seqct);
	
      pwNGdist.dN   = make_double_matrix(seqct,seqct);
      pwNGdist.dS   = make_double_matrix(seqct,seqct);
      pwNGdist.dNdS = make_double_matrix(seqct,seqct);
	
      if( do_kaks_yn00(msa, &pwMLdist,&pwNGdist) < 0 ) {
	fprintf(stderr, "warning: problem with align for %s %s\n",
		msa->sqname[0], msa->sqname[1]);
	continue;
      }
       
      for( i=0; i < msa->nseq; i++ ) {
	/* stip trailing whitespace */
	k = strlen(msa->sqname[i]);
	while( k-- > 0 && msa->sqname[i][k] == ' ');
	msa->sqname[i][k+1] = '\0';

	for( j=i+1;j<msa->nseq;j++ ) {
	  k = strlen(msa->sqname[j]);
	  while( k-- > 0 && msa->sqname[j][k] == ' ');
	  msa->sqname[j][k+1] = '\0';

	  fprintf(ofd,"%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
		  msa->sqname[i],
		  msa->sqname[j],
		  pwMLdist.dN[i][j],
		  pwMLdist.dS[i][j], 
		  pwMLdist.dNdS[i][j],
		  pwMLdist.N[i][j],
		  pwMLdist.S[i][j],
		  pwMLdist.kappa[i][j],
		  pwMLdist.t[i][j],
		  msa->alen);
	}
      }
    
      cleanup_matrix((void **)pwMLdist.N,seqct);
      cleanup_matrix((void **)pwMLdist.dN,seqct);
      cleanup_matrix((void **)pwMLdist.S,seqct);
      
      cleanup_matrix((void **)pwMLdist.dS,seqct);
      
      cleanup_matrix((void **)pwMLdist.SEdS,seqct);
      cleanup_matrix((void **)pwMLdist.SEdN,seqct);
      cleanup_matrix((void **)pwMLdist.t,seqct);
      cleanup_matrix((void **)pwMLdist.dNdS,seqct);
      cleanup_matrix((void **)pwMLdist.kappa,seqct);
      
      cleanup_matrix((void **)pwNGdist.dN,seqct);
      cleanup_matrix((void **)pwNGdist.dS,seqct);
      cleanup_matrix((void **)pwNGdist.dNdS,seqct);
      
      free(pwNGdist.dNdS);
      free(pwNGdist.dN);
      free(pwNGdist.dS);
      
      free(pwMLdist.dNdS);
      free(pwMLdist.dN);
      free(pwMLdist.dS);
      free(pwMLdist.N);
      free(pwMLdist.S);
      free(pwMLdist.SEdS);
      free(pwMLdist.SEdN);
      free(pwMLdist.t);
      free(pwMLdist.kappa);
      
      seqct = 0;
    }
  }

  if( ofd && ofd != stdout )
    fclose(ofd);
  
  return 0;
}
