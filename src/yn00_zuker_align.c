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
  { "--matrix-dir",    FALSE, sqdARG_STRING },
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

static char experts[] = "\n   --optimal     : show the optimal alignment.\n   --global      : generate global rather than the default local alignments\n   --lower_bound : Lower bound for scores to use when generating suboptimal\n                   alignments.\n   --upper_bound : Upper bound for scores to use when generating suboptimal\n                   alignments.\n   --minlen      : minimum length of suboptimal alignments to report.\n                   Default is 50\n   --max_to_show : Maximum number of alignments to show.  Default is 1000.\n   --gapopen num : Gap open penalty.  Default is -10.\n   --gapext      : Gap extension penalty.  Default is -2.\n   --informat    : Format of the input sequence file database.\n   --output      : specify an output file instead of summary.dat\n   --matrix      : Substitution matrix: a file path, or a bare name (e.g.\n                   BLOSUM62) resolved within a matrix directory. Default BLOSUM62.\n   --matrix-dir  : Directory to search for named matrices (overrides\n                   $SUBOPT_MATRIX_DIR / $SUBOPTDIR).\n   --gapchar     : Provide a different gap char from '-'\n   --separate    : Write each alignment in a separate file\n   --sorted      : Return alignments sorted by score.\n   --cutoff_percent: Only show top N percent of alignments.\n   --noheader    : Do not show header in summary table output\n";

static char      Default_output_name[] = "summary.dat";
const int        default_max_alignments = 1000;
const int        default_min_aln_len = 50;
const int        default_cutoff_percent = 10; /* only show alignments 
						 which are  10% of the 
						 maximum 
					      */
const AlignType  default_aln_type   = local;


/* protype */
int comparealn(const void *, const void *);

int
main (int argc, char ** argv )
{
  char     *seqfile;            /* name of sequence file     */
  SQINFO    sqinfo;             /* extra info about sequence */
  SQFILE   *dbfp;		/* open sequence file        */
  int       fmt;	        /* format of seqfile         */
  char     *seq;		/* sequence                  */
  sequence  seqs[2], cds_seqs[2];
  char  *optname;
  char  *optarg;
  int    optind,rc;
  int    do_sorted = 0;
  int    seqct = 0,cdsct = 0;
  int    min_aln_len      = 0;
  int    aln_count = 0;
  int    separate_files   = 0, provide_ofname = 0;
  char   * output_filename = 0, *submat_file = 0, *matrix_dir = 0;
  int    showheader=1;

  FILE  *ofd = NULL, *fd;
  alignment   *cds_aln = NULL;
  alignment ** alignments = NULL;  /* an array of pairwise alignments */

  int    i,ik,jk;
  pairwise_distances pwMLdist, pwNGdist;
  /* done initializing and declaring variables */

  Alntype = default_aln_type;

  Max_alignments_to_report = default_max_alignments;
  Alignment_cutoff_percent = default_cutoff_percent;
  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */


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
      else if (strcmp(optname, "--matrix-dir") == 0) matrix_dir = optarg;
      else if (strcmp(optname, "--upper_bound")  == 0) 
	Upper_bound  = atoi(optarg); 
      else if (strcmp(optname, "--max_to_show")  == 0)  
	Max_alignments_to_report = atoi(optarg);
      else if (strcmp(optname, "--minlen")  == 0)
	min_aln_len  = atoi(optarg);
      else if (strcmp(optname, "--gapopen") == 0)  {
	/* penalties are stored as negative scores; accept 10 or -10 */
	Gapopen = -abs(atoi(optarg));
      } else if (strcmp(optname, "--gapext")  == 0)  {
	Gapext = -abs(atoi(optarg));

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
  
  {
    char *resolved = resolve_matrix_file(submat_file, matrix_dir);
    if( ! resolved ) {
      fprintf(stderr,
	      "Cannot locate substitution matrix '%s' in any matrix directory.\n"
	      "Use --matrix <file>, --matrix-dir <dir>, or set $SUBOPT_MATRIX_DIR.\n",
	      submat_file ? submat_file : DEFAULT_MATRIX_NAME);
      exit(EXIT_FAILURE);
    }
    submat_file = resolved; /* now heap-owned; freed near end of main */
  }

  /* open matrix */
  fd = fopen(submat_file, "r");
  if( !fd ||
      ! ParsePAMFile(fd,&ScoringMatrix, &MatrixScale) ) {
    fprintf(stderr, "Cannot parse or open matrix file %s\n",submat_file);
    free(submat_file);
    exit(EXIT_FAILURE);
  }
  fclose(fd);
  (void)min_aln_len;    /* parsed for compatibility; not applied */
  (void)separate_files; /* parsed for compatibility; not applied */
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
  SeqfileClose(dbfp);
  if( seqct != 2 ) {
    fprintf(stderr,"Must have provided a valid file with at least 2 sequences in it\n");
    goto end;
  }
  alignments = (alignment **)calloc(Max_alignments_to_report,
				    sizeof(alignment *));
  aln_count = zuker_align(&seqs[0],&seqs[1],&alignments);

  if( aln_count < 0 ) {
    fprintf(stderr,"error in zuker align\n");
    aln_count = 0;
    goto end;
  }
  if( do_sorted ) {
    qsort(alignments,(size_t)aln_count, sizeof(alignment *),comparealn);
  }
  if( ! provide_ofname )
    output_filename = Default_output_name;
  if( strcmp(output_filename, "-") == 0 ) {
    ofd = stdout;
  } else if( (ofd = fopen(output_filename,"w")) == NULL ) {
    fprintf(stderr,"could not open output file %s\n", output_filename);
    goto end;
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
    if( do_kaks_yn00(cds_aln->msa, &pwMLdist,&pwNGdist) < 0 ) {
      fprintf(stderr, "warning: problem with align %d for %s %s\n", i,
	      cds_aln->msa->sqname[0], cds_aln->msa->sqname[1]);
    } else
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
    cds_aln = NULL;
    cleanup_alignment(alignments[i]);
    alignments[i] = NULL;
  }

  end:
  if( ofd && ofd != stdout )
    fclose(ofd);
  cleanup_alignment(cds_aln);
  for( i = 0; i < aln_count; i++ )
    cleanup_alignment(alignments[i]);
  free(submat_file);
  if( alignments)
    free(alignments);

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
