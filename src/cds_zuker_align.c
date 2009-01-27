/* Author: Jason Stajich <jason.stajichATduke.edu>
 * $Id: cds_zuker_align.c,v 1.10 2004/07/01 21:29:30 jason Exp $
 *
 */

/* The SUBOPTDIR env variable can be used to set the location of
 * of the matricies
 */

#include "align.h"
#include "mrtrans.h"
#include <sys/time.h>

struct opt_s OPTIONS[] = {
  { "-h",              TRUE, sqdARG_NONE },  
  { "-v",              TRUE, sqdARG_NONE },    
  { "--matrix",        FALSE, sqdARG_STRING },
  { "--global",        FALSE, sqdARG_NONE},  
  { "--optimal",       FALSE, sqdARG_NONE},  
  { "--informat",      FALSE, sqdARG_STRING },
  { "--outformat",     FALSE, sqdARG_STRING },
  { "--gapopen",       FALSE, sqdARG_STRING },
  { "--gapext",        FALSE, sqdARG_STRING },
  { "--gapchar",       FALSE, sqdARG_STRING },
  { "--upper_bound",   FALSE, sqdARG_STRING },
  { "--max_to_show",   FALSE, sqdARG_STRING },
  { "--lower_bound",   FALSE, sqdARG_STRING },
  { "--minlen",        FALSE, sqdARG_STRING },
  { "--separate",      FALSE, sqdARG_NONE },
  { "--output",        FALSE, sqdARG_STRING },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\nUsage: cds_zuker_align [-options] <sequence file>\n\n   Generate suboptimal alignments for a pair of CDS sequences using the protein sequences as the template.\n   Available options:\n   -h      : help; print version and usage info\n   -v      : verbose output -- show alignment path states.\n ";

static char banner[] = "zuker_align - Generate suboptimal alignments with the Zuker algorithm";

static char experts[] = "\n   --optimal     : show the optimal alignment.\n   --global      : generate global rather than the default local alignments\n   --lower_bound : Lower bound for scores to use when generating suboptimal\n                   alignments.\n   --upper_bound : Upper bound for scores to use when generating suboptimal\n                   alignments.\n   --minlen      : minimum length of suboptimal alignments to report.\n                   Default is 50\n   --max_to_show : Maximum number of alignments to show.  Default is 5.\n   --gapopen num : Gap open penalty.  Default is -10.\n   --gapext      : Gap extension penalty.  Default is -2.\n   --informat    : Format of the input sequence file database.\n   --outformat   : Format of the output sequence alignments.\n   --output      : specify an output file instead or writing to STDOUT\n   --matrix      : Path to substitution matrix.\n                   Default is ../matricies/aa/BLOSUM50.\n   --gapchar     : Provide a different gap char from '-'\n   --separate    : Write each alignment in a separate file\n";

const int        default_cutoff_num = 5;
const int        default_min_aln_len = 50;
const AlignType  default_aln_type   = local;
const char Default_submat[] = "matricies/aa/BLOSUM50";
const char suff[] = "aln";

main (int argc, char ** argv ) 
{
  char     *seqfile;            /* name of sequence file     */
  SQINFO    sqinfo;             /* extra info about sequence */
  SQFILE   *dbfp;		/* open sequence file        */
  int       fmt,ofmt=106;	/* format of seqfile         */
  char     *seq;		/* sequence                  */
  int       type;		/* kAmino, kDNA, kRNA, or kOtherSeq */
  sequence  seqs[2], cds_seqs[2];
  char  *optname;
  char  *optarg, *t, *sfname;
  int    optind;
  int    be_quiet;
  int    seqct = 0,cdsct = 0;
  int    min_aln_len      = 0;
  int    do_oneline       = 0;
  int    separate_files   = 0, provide_ofname = 0;
  char   * output_filename = 0, *submat_file = 0;
  
  FILE  *ofd, *fd;
  alignment   *cds_aln;
  alignment ** alignments;  /* an array of pairwise alignments */

  int    i,j, aln_count, rc;
  struct timeval tp;

  Alntype = default_aln_type;
    
  Max_alignments_to_report = default_cutoff_num;
  
  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */
  be_quiet  = FALSE;
  type      =  kOtherSeq;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--matrix") == 0)  submat_file = optarg; 
      else if (strcmp(optname, "--upper_bound")  == 0) 
	Upper_bound  = atoi(optarg); 
      else if (strcmp(optname, "--lower_bound")  == 0) 
	Lower_bound  = atoi(optarg);
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
      } else if (strcmp(optname, "--outformat") == 0) {
	ofmt = String2SeqfileFormat(optarg);
	if (ofmt == SQFILE_UNKNOWN) 
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
      }
    }

  if (argc - optind != 1) Die("%s\n", usage);
  seqfile = argv[argc-1];
  if( ! output_filename ) { 
    output_filename = calloc(strlen(seqfile)+36,sizeof(char));
    gettimeofday(&tp, NULL);
    sprintf(output_filename,"%s_%ld",seqfile, tp.tv_sec);
  }  

  if( ! submat_file ) { 
    if( (t = getenv("SUBOPTDIR")) != 0 || 
	(t = getenv("SUBOPT_DIR")) != 0 ) {
      fprintf(stderr,"here with %s\n",t);
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
      seqs[seqct].seqname = calloc(strlen(sqinfo.name)+1,sizeof(char));
      cds_seqs[cdsct].seqname = calloc(strlen(sqinfo.name)+1,sizeof(char));
      strcpy(seqs[seqct].seqname,sqinfo.name);
      strcpy(cds_seqs[cdsct].seqname,sqinfo.name);

      cds_seqs[cdsct].length = sqinfo.len;
      cds_seqs[cdsct].alphabet = ( sqinfo.type == kDNA ) ? dna : rna;
      seqs[seqct].length = strlen(seqs[seqct].seqstr);
      seqs[seqct].alphabet = protein;
      cdsct++;
      seqct++;      
    } else {
      fprintf(stderr,"Expect CDS sequences (DNA or RNA) not Protein\n");
      goto end;
    }
    FreeSequence(NULL, &sqinfo);    
  }
  if( seqct != 2 ) {
    fprintf(stderr,"Must have provided a valid file with at least 2 sequences in it");
    goto end;
  }

  alignments = (alignment **)calloc(Max_alignments_to_report,
				    sizeof(alignment *));
  
  aln_count = zuker_align(&seqs[0],&seqs[1],&alignments);
  

  sfname = calloc(strlen(output_filename) + 10,sizeof(char));
  
  for( i = 0; i < aln_count; i++ ) {
    rc = mrtrans(alignments[i], cds_seqs, &cds_aln,0);
    if( rc != 0  ) { 
      fprintf(stderr, "Could not map the coding sequence to the protein alignemnt for aln %d: %d\n",i,rc);
      goto end;
    }
    if( separate_files ) {
      sprintf(sfname,"%s.%d.%s",output_filename,i,suff);
      if( ofd )    // start a new file 
	fclose(ofd);      
      ofd = fopen(sfname,"w");
      if( ! ofd ) {
	fprintf(stderr, "could not open file %s",sfname);
	goto end;
      }
      
    } else if(! ofd && provide_ofname ) { 
      ofd = fopen(output_filename,"w");
      if( ! ofd ) {
	fprintf(stderr, "could not open file %s",output_filename);
	goto end;
      }
    } else if( ! ofd ) { 
      ofd = stdout;
    }

    if( ofmt >= 100 ) { 
      MSAFileWrite(ofd,cds_aln->msa, ofmt,do_oneline);
    } else { 
      for(j=0; j < cds_aln->msa->nseq; j++ ) {
	WriteSeq(ofd, ofmt, 
		 cds_aln->msa->aseq[j],
		 &(cds_aln->sqinfo[j]) );
      }
    }
  }
  
  if( ofd && ofd != stdout )
    fclose(ofd);
  
  end:
  if( sfname) free(sfname) ;
  free(submat_file);
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
  for(i = 0; i < aln_count; i++ ) {
    if( alignments[i] ) {
      cleanup_alignment(alignments[i]);
      alignments[i] = 0;
    }
  }
  
  cleanup_alignment(cds_aln);
#if DMALLOC
  dmalloc_shutdown();
#endif

  return 0;
}
