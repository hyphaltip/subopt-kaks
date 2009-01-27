/* Author: Jason Stajich <jason.stajichATduke.edu>
 * $Id: zuker_align.c,v 1.10 2004/05/27 02:35:09 jason Exp $
 *
 */

#include "align.h"

struct opt_s OPTIONS[] = {
  { "-h",          TRUE, sqdARG_NONE },  
  { "-v",          TRUE, sqdARG_NONE },    
  { "--matrix",   FALSE, sqdARG_STRING },
  { "--global",    FALSE, sqdARG_NONE},  
  { "--optimal",   FALSE, sqdARG_NONE},  
  { "--informat", FALSE, sqdARG_STRING },
  { "--outformat", FALSE, sqdARG_STRING },
  { "--gapopen",  FALSE, sqdARG_STRING },
  { "--gapext",   FALSE, sqdARG_STRING },
  { "--gapchar",   FALSE, sqdARG_STRING },
  { "--upper_bound",   FALSE, sqdARG_STRING },
  { "--max_to_show",   FALSE, sqdARG_STRING },
  { "--lower_bound",   FALSE, sqdARG_STRING },
  { "--minlen",   FALSE, sqdARG_STRING },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\nUsage: suboptaligner [-options] <sequence file>\n\n   Generate suboptimal alignments for a pair of sequences.\n   Available options:\n   -h      : help; print version and usage info\n   -v      : verbose output -- show alignment path states.\n";

static char banner[] = "zuker_align - Generate suboptimal alignments with the Zuker algorithm";

static char experts[] = "\n   --optimal     : show the optimal alignment.\n   --global      : generate global rather than the default local alignments\n   --lower_bound : Lower bound for scores to use when generating suboptimal\n                   alignments.\n   --upper_bound : Upper bound for scores to use when generating suboptimal\n                   alignments.\n   --minlen      : minimum length of suboptimal alignments to report.\n                   Default is 50\n   --max_to_show : Maximum number of alignments to show.  Default is 5.\n   --gapopen     : Gap open penalty.  Default is -10.\n   --gapext      : Gap extension penalty.  Default is -2.\n   --informat    : Format of the input sequence file database.\n   --outformat   : Format of the output sequence alignments.\n   --matrix      : Path to substitution matrix.\n                   Default is ../matricies/aa/BLOSUM50.\n   --gapchar    : Provide a different gap char from '-'\n";

const int        default_cutoff_num = 5;
const int        default_min_aln_len = 50;
const AlignType  default_aln_type   = local;
char * SUBMAT = "../matricies/aa/BLOSUM50";

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
  char  *optarg;
  int    optind;
  int    be_quiet;
  int    seqct = 0, cdsct = 0;
  int    min_aln_len      = 0;
  int    do_oneline       = 0, optimal = 0;
  FILE  *fd;

  alignment ** alignments;  /* an array of pairwise alignments */

  int    i,j, aln_count,len;
  Alntype = default_aln_type;

  Max_alignments_to_report = default_cutoff_num;
  
  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */
  be_quiet  = FALSE;
  type      =  kOtherSeq;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--matrix") == 0)  SUBMAT = optarg; 
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
      } else if ( strcmp(optname,"--optimal") == 0 ) {
	optimal = 1;
      }
    }

  if (argc - optind != 1) Die("%s\n", usage);
  seqfile = argv[argc-1];
  
  /* open matrix */
  fd = fopen(SUBMAT, "r");
  
  if( ! ParsePAMFile(fd,&ScoringMatrix, &MatrixScale) ) {
    fprintf(stderr, "Cannot parse or open matrix file %s\n",SUBMAT);
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
    len = strlen(sqinfo.name)+1;
    if( sqinfo.type == kDNA || sqinfo.type == kRNA ) {
      cds_seqs[cdsct].seqstr = seq;
      seqs[seqct].seqstr = Translate(seq,stdcode1);
      seqs[seqct].seqname = calloc(len,sizeof(char));
      cds_seqs[cdsct].seqname = calloc(len,sizeof(char));
      strncpy(seqs[seqct].seqname,sqinfo.name,len);
      strncpy(cds_seqs[cdsct++].seqname,sqinfo.name,len);

      cds_seqs[cdsct].length = sqinfo.len;
      cds_seqs[cdsct].alphabet = ( sqinfo.type == kDNA ) ? dna : rna;
      seqs[seqct].length = strlen(seqs[seqct].seqstr);
      seqs[seqct].alphabet = protein;

    } else {
      seqs[seqct].seqstr = seq;
      seqs[seqct].seqname = calloc(len,sizeof(char));
      strncpy(seqs[seqct].seqname,sqinfo.name,len);
      seqs[seqct].length = sqinfo.len;
      seqs[seqct].alphabet = protein;
    }
    seqct++;
    FreeSequence(NULL, &sqinfo);    
  }
  if( seqct != 2 ) {
    fprintf(stderr,"Must have provided a valid file with at least 2 sequences in it");
    goto end;
  }

  alignments =(alignment **)calloc(Max_alignments_to_report,
				   sizeof(alignment *));
  aln_count = zuker_align(&seqs[0],&seqs[1],&alignments);  
  
  for( i = 0; i < aln_count; i++ ) {
    if( ofmt >= 100 ) { 
      MSAFileWrite(stdout, alignments[i]->msa,ofmt,do_oneline);
    } else { 
      for(j=0; j < alignments[i]->msa->nseq; j++ ) {
	WriteSeq(stdout, ofmt, 
		 alignments[i]->msa->aseq[j],
		 &(alignments[i]->sqinfo[j]) );
      }
    }
    if( optimal )
      i = aln_count;
    
    /* print_alignment_fasta(alignments[i]); */
  }
  end:
  
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
  return 0;
}
