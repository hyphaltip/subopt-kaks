/* Author: Jason Stajich <jason.stajichATduke.edu>
 * $Id: yn00_zuker_align.c,v 1.16 2004/07/01 21:29:31 jason Exp $
 *
 */

/* The SUBOPTDIR env variable can be used to set the location of
 * of the matricies
 */

#include "align.h"
#include "mrtrans.h"
#include "paml.h"

struct opt_s OPTIONS[] = {
  { "-h",              TRUE, sqdARG_NONE },  
  { "-v",              TRUE, sqdARG_NONE },    
  { "--matrix",        FALSE, sqdARG_STRING },
  { "--global",        FALSE, sqdARG_NONE},  
  { "--optimal",       FALSE, sqdARG_NONE},  
  { "--informat",      FALSE, sqdARG_STRING },
  { "--showalignment", FALSE, sqdARG_NONE},
  { "--outformat",     FALSE, sqdARG_STRING },
  { "--gapopen",       FALSE, sqdARG_STRING },
  { "--gapext",        FALSE, sqdARG_STRING },
  { "--gapchar",       FALSE, sqdARG_STRING },
  { "--upper_bound",   FALSE, sqdARG_STRING },
  { "--max_to_show",   FALSE, sqdARG_STRING },
  { "--lower_bound",   FALSE, sqdARG_STRING },
  { "--sorted",        FALSE, sqdARG_NONE },
  { "--cutoff_percent",FALSE, sqdARG_STRING },
  { "--minlen",        FALSE, sqdARG_STRING },
  { "--separate",      FALSE, sqdARG_NONE },
  { "--noheader",      FALSE, sqdARG_NONE},
  { "--output",        FALSE, sqdARG_STRING },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\nUsage: yn00_zuker_align [-options] <sequence file>\n\n   Generate suboptimal alignments for a pair of sequences.\n   Available options:\n   -h      : help; print version and usage info\n   -v      : verbose output -- show alignment path states.\n ";

static char experts[] = "\n   --optimal     : show the optimal alignment.\n   --global      : generate global rather than the default local alignments\n   --lower_bound : Lower bound for scores to use when generating suboptimal\n                   alignments.\n   --upper_bound : Upper bound for scores to use when generating suboptimal\n                   alignments.\n   --minlen      : minimum length of suboptimal alignments to report.\n                   Default is 50\n   --max_to_show : Maximum number of alignments to show.  Default is 5.\n   --gapopen num : Gap open penalty.  Default is -8.\n   --gapext      : Gap extension penalty.  Default is -2.\n   --informat    : Format of the input sequence file database.\n  --output      : specify an output file instead or writing to STDOUT\n   --matrix      : Path to substitution matrix.\n                   Default is ../matricies/aa/BLOSUM50.\n   --gapchar     : Provide a different gap char from '-'\n   --separate    : Write each alignment in a separate file\n   --sorted      : Return alignments sorted by score.\n   --cutoff_percent: Only show top N percent of alignments.\n   --noheader    : Do not show header in summary table output\n";

static char      Default_output_name[] = "summary.dat";
const int        default_max_alignments = 1000;
const int        default_min_aln_len = 50;
const int        default_cutoff_percent = 10; /* only show alignments 
						 which are  10% of the 
						 maximum 
					      */
const AlignType  default_aln_type   = local;
const char Default_submat[] = "matricies/aa/BLOSUM50";
const char suff[] = "aln";


/* protype */
int comparealn(const void *, const void *);

main (int argc, char ** argv ) 
{
  char     *seqfile;            /* name of sequence file     */
  SQINFO    sqinfo;             /* extra info about sequence */
  SQFILE   *dbfp;		/* open sequence file        */
  int       fmt;	        /* format of seqfile         */
  char     *seq;		/* sequence                  */
  int       type;		/* kAmino, kDNA, kRNA, or kOtherSeq */
  sequence  seqs[2], cds_seqs[2];
  char  *optname;
  char  *optarg;
  int    optind,rc;
  int    be_quiet, do_sorted;
  int    seqct = 0,cdsct = 0;
  int    min_aln_len      = 0;
  int    do_oneline       = 0, aln_count;
  int    separate_files   = 0, provide_ofname = 0;
  char   * output_filename = 0, *submat_file = 0, *t;
  int    showheader=1;
  
  FILE  *ofd, *fd;
  alignment   *cds_aln;
  alignment ** alignments;  /* an array of pairwise alignments */
  
  int    i,j,ik,jk;
  pairwise_distances pwMLdist, pwNGdist;
  /* done initializing and declaring variables */

  Alntype = default_aln_type;

  Max_alignments_to_report = default_max_alignments;
  Alignment_cutoff_percent = default_cutoff_percent;  
  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */
  be_quiet  = FALSE;
  type      =  kOtherSeq;


  /* for our purposes this is only pairwise alignments, but
   * would rather do it correctly in case we move to MSA case 
   */

  pwMLdist.N    = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.dN   = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.S    = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.dS   = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.dNdS = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.SEdS = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.SEdN = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.t    = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwMLdist.kappa= make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);

  pwNGdist.dN   = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwNGdist.dS   = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);
  pwNGdist.dNdS = make_double_matrix(NUM_PW_SEQS,NUM_PW_SEQS);


  // Process command-line arguments  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--matrix") == 0)  submat_file = optarg; 
      else if (strcmp(optname, "--upper_bound")  == 0) 
	Upper_bound  = atoi(optarg); 
      else if (strcmp(optname, "--max_to_show")  == 0)  
	Max_alignments_to_report = atoi(optarg);
      else if (strcmp(optname, "--minlen")  == 0) 
	min_aln_len  = atoi(optarg); 
      else if (strcmp(optname, "--quiet")   == 0)  be_quiet  = TRUE; 
      else if (strcmp(optname, "--gapopen") == 0)  {
	Gapopen = atoi(optarg); 
	if( Gapopen < 0 ) Gapopen *= -1;
	
      } else if (strcmp(optname, "--gapext")  == 0)  {
	Gapext = atoi(optarg); 
	if( Gapext < 0 ) Gapext *= -1;

      } else if (strcmp(optname, "--informat") == 0) {
	fmt = String2SeqfileFormat(optarg);
	if (fmt == SQFILE_UNKNOWN) 
	  Die("unrecognized sequence file format \"%s\"", optarg);
      }  else if( strcmp(optname, "--global") == 0 ) {
	Alntype = global;
      } else if (strcmp(optname, "-h") == 0) {
	puts(usage);
	puts(experts);
        exit(EXIT_SUCCESS);
      } else if ( strcmp(optname, "-v") == 0 ) {
	Verbose = 1;
      } else if ( strcmp(optname, "--gapchar") == 0 ) {
	GapChar = optarg[0];
      } else if(  strcmp(optname, "--separate") == 0 ) {
	separate_files = 1;
      }  else if(  strcmp(optname, "--output") == 0 ) {
	output_filename = optarg;
	provide_ofname    = 1;
      } else if( strcmp(optname, "--cutoff_percent") == 0) {
	Alignment_cutoff_percent = atoi(optarg);
      } else if( strcmp(optname, "--sorted") == 0 ) {
	do_sorted = 1;
      } else if( strcmp(optname, "--noheader" ) == 0 ) {
	showheader = 0;
      }
    }

  if (argc - optind != 1) Die("%s\n", usage);
  seqfile = argv[argc-1];  
  
  if( ! submat_file ) { 
    if( (t = getenv("SUBOPTDIR")) != 0 || 
	(t = getenv("SUBOPT_DIR")) != 0 ) {
      rc = strlen(t)+strlen((void *)Default_submat) + 2;
      submat_file = calloc(rc, sizeof(char));
      sprintf(submat_file, "%s/%s",t,Default_submat);
    } else { 
      submat_file = calloc(strlen((void *)Default_submat) + 24, sizeof(char));
      sprintf(submat_file, "../%s",Default_submat);
    }
  }

  /* open matrix */
  fd = fopen(submat_file, "r");
  if( !fd || 
      ! ParsePAMFile(fd,&ScoringMatrix, &MatrixScale) ) {
    fprintf(stderr, "Cannot parse or open matrix file %s\n",submat_file);
    free(submat_file);
#if DMALLOC
    dmalloc_shutdown();
#endif
    exit(EXIT_SUCCESS);
  }  
  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (fmt == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    fmt = SQFILE_FASTA;
  
  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);
  while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo))
  {
    if( seqct >= 2 ) break; /* only process 2 sequences */
    sqinfo.type = Seqtype(seq);
    if( sqinfo.type == kDNA || sqinfo.type == kRNA ) {
      cds_seqs[cdsct].seqstr = seq;
      seqs[seqct].seqstr = Translate(seq,stdcode1);
      seqs[seqct].seqname = calloc(SQINFO_NAMELEN,sizeof(char));
      cds_seqs[cdsct].seqname = calloc(SQINFO_NAMELEN,sizeof(char));
      strncpy(seqs[seqct].seqname,sqinfo.name,SQINFO_NAMELEN);
      strncpy(cds_seqs[cdsct].seqname,sqinfo.name,SQINFO_NAMELEN);
      
      cds_seqs[cdsct].length = sqinfo.len;
      cds_seqs[cdsct].alphabet = ( sqinfo.type == kDNA ) ? dna : rna;
      seqs[seqct].length = strlen(seqs[seqct].seqstr);
      seqs[seqct].alphabet = protein;
      cdsct++;
      
    } else {
      fprintf(stderr,"Expect CDS sequences (DNA or RNA) not Protein\n");
      goto end;
    }
    seqct++;
    FreeSequence(NULL, &sqinfo);    
  }
  if( seqct != 2 ) {
    fprintf(stderr,"Must have provided a valid file with at least 2 sequences in it");
    goto end;
  }  
  alignments = (alignment **)calloc(Max_alignments_to_report,
				    sizeof(alignment *));
  aln_count = zuker_align(&seqs[0],&seqs[1],&alignments);  
  
  if( aln_count < 0 ) {
    fprintf(stderr,"error in zuker align");
    goto end;
  }
  if( do_sorted ) { 
    qsort(alignments,(size_t)aln_count, sizeof(alignment *),comparealn);    
  }
  if( provide_ofname ) { 
    ofd = fopen(output_filename,"w");
  } else { 
    ofd = fopen(Default_output_name,"w");
  }
  if( showheader ) {
    fprintf(ofd,"SCORE\tdN\tdS\tOMEGA\tN\tS\tkappa\tt\tLENGTH\n");
  }
  
  for( i = 0; i < aln_count; i++ ) {
     rc = mrtrans(alignments[i], cds_seqs, &cds_aln,1);
     
    if( rc != 0  ) { 
      fprintf(stderr, "Could not map the coding sequence to the protein alignemnt for aln %d: %d\n",i,rc);
      goto end;
    }
    /* do KaKs here - from yn00 for now*/
    do_kaks_yn00(cds_aln->msa, &pwMLdist,&pwNGdist);
    for(ik = 0; ik < NUM_PW_SEQS; ik++ ) {	  
      for( jk = ik+1; jk < NUM_PW_SEQS; jk++ ) {
	fprintf(ofd,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
		alignments[i]->score,
		pwMLdist.dN[ik][jk],pwMLdist.dS[ik][jk], 
		pwMLdist.dNdS[ik][jk],
		pwMLdist.N[ik][jk],
		  pwMLdist.S[ik][jk],
		pwMLdist.kappa[ik][jk],
		  pwMLdist.t[ik][jk],
		cds_aln->msa->alen);
      }
    }
    cleanup_alignment(cds_aln);
    cleanup_alignment(alignments[i]);    
  }
  

  if( ofd && ofd != stdout )
    fclose(ofd);      
	  
  end:
  free(submat_file);
  if( alignments) 
    free(alignments);

/*  if( output_filename ) 
    free(output_filename);
*/
  Free2DArray((void **)ScoringMatrix,27);
  for(i =0; i< seqct; i++ ) {
    free(seqs[i].seqstr);
    free(seqs[i].seqname);    
    seqs[i].seqstr = seqs[i].seqname = 0;
  }
  for(i = 0; i < cdsct; i++) {
    free(cds_seqs[i].seqstr);
    free(cds_seqs[i].seqname);    
    cds_seqs[i].seqstr = cds_seqs[i].seqname = 0;
  }
  
  cleanup_matrix((void **)pwMLdist.N,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.dN,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.S,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.dS,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.SEdS,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.SEdN,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.t,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.dNdS,NUM_PW_SEQS);
  cleanup_matrix((void **)pwMLdist.kappa,NUM_PW_SEQS);

  cleanup_matrix((void **)pwNGdist.dN,NUM_PW_SEQS);
  cleanup_matrix((void **)pwNGdist.dS,NUM_PW_SEQS);
  cleanup_matrix((void **)pwNGdist.dNdS,NUM_PW_SEQS);

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
  
#if DMALLOC
  dmalloc_shutdown();
#endif
  return 0;
}

int 
comparealn(const void * b, const void * a) 
{
  alignment *aa, *bb; 
  if( ! a || ! b ) return 0;
  aa = *(alignment **)a, bb = *(alignment **)b;
  return ((aa->score) == (bb->score) ? 0 : (aa->score) >(bb->score) ? 1 : -1 );
}
