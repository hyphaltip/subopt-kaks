/* Author: Jason Stajich <jason.stajichATduke.edu>
 * $Id: yn00_cds_optimal.c,v 1.14 2004/07/06 19:46:53 jason Exp $
 *
 */

/* The SUBOPTDIR env variable can be used to set the location of
 * of the matricies
 */

#include "align.h"
#include "mrtrans.h"
#include "paml.h"
#include <sys/time.h> /* used only for getting a reasonably unique filename */

struct opt_s OPTIONS[] = {
  { "-h",              TRUE, sqdARG_NONE },  
  { "-v",              TRUE, sqdARG_NONE },    
  { "--matrix",        FALSE, sqdARG_STRING },
  { "--global",        FALSE, sqdARG_NONE},  
  { "--informat",      FALSE, sqdARG_STRING },
  { "--outformat",     FALSE, sqdARG_STRING },
  { "--gapopen",       FALSE, sqdARG_STRING },
  { "--gapext",        FALSE, sqdARG_STRING },
  { "--gapchar",       FALSE, sqdARG_STRING },
  { "--showtable",     FALSE, sqdARG_NONE },
  { "--noheader",      FALSE, sqdARG_NONE},
  { "--output",        FALSE, sqdARG_STRING },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[] = "\nUsage: yn00_cds_optimal [-options] <sequence file>\n\n   Generate the optimal alignments for a pair of CDS sequences using the protein sequence as a template.\n   Available options:\n   -h      : help; print version and usage info\n   -v      : verbose output -- show alignment path states.\n ";
static char experts[] = "\n   --global      : generate global rather than the default local alignments\n   --gapopen num : Gap open penalty.  Default is -10.\n   --gapext      : Gap extension penalty.  Default is -2.\n   --informat    : Format of the input sequence file database.\n   --outformat   : Format of the output sequence alignments.\n   --output      : specify an output file instead or writing to STDOUT\n   --matrix      : Path to substitution matrix.\n                   Default is ../matricies/aa/BLOSUM50.\n   --gapchar     : Provide a different gap char from '-'\n   --showtable   : Show the summary Ka/Ks table instead of the alignment.\n   --noheader     : Do not show the header for the table summary";

const AlignType  default_aln_type   = local;
const char Default_submat[] = "matricies/aa/BLOSUM50";
const char suff[] = "aln";

main (int argc, char ** argv ) 
{
  char     *seqfile;            /* name of sequence file     */
  SQINFO    sqinfo;             /* extra info about sequence */
  SQFILE   *dbfp;		/* open sequence file        */
  int       fmt,ofmt=106;	/* format of seqfile         */
                                /* 106 is PHYLIP format in SQUID */
  char     *seq;		/* sequence                  */
  int       type;		/* kAmino, kDNA, kRNA, or kOtherSeq */
  sequence  * seqs, * cds_seqs;
  sequence  tmp_seqs[2], tmp_cds_seqs[2];
  char  *optname;
  char  *optarg, *t;
  int    optind;
  int    be_quiet;
  int    seqct = 0,cdsct = 0;
  int    min_aln_len      = 0;
  int    do_oneline       = 0;
  char   * output_filename = 0, *submat_file = 0;
  int    showaln = 1;
  int    showheader=1;
  FILE  *ofd, *fd;
  alignment   *cds_aln;
  alignment * opt_alignment = NULL;  /* place for pairwise alignment */

  int    len,i,j, k, jk,ik,aln_count, rc;
  pairwise_distances pwMLdist, pwNGdist;
  int firsttime = 1;
  
  struct timeval tp;

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
/*
  pwMLdist.N    = pwMLdist.dN   = pwMLdist.S    = 0;
  pwMLdist.dS   = pwMLdist.dNdS = pwMLdist.SEdS = 0;
  pwMLdist.SEdN = pwMLdist.t    = pwMLdist.kappa= 0;
      
  pwNGdist.dN   = pwNGdist.dS   = pwNGdist.dNdS = 0;
*/

  Alntype = default_aln_type;
  
  /* Command line Parse */
  fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */
  be_quiet  = FALSE;
  type      =  kOtherSeq;

  /* for our purposes this is only pairwise alignments, but
   * would rather do it correctly in case we move to MSA case 
   */
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--matrix") == 0)  submat_file = optarg; 
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
      }  else if(  strcmp(optname, "--output") == 0 ) {
	output_filename = optarg;	  
      } else if( strcmp(optname, "--showtable" ) == 0  ) {
	showaln = 0;
      } else if( strcmp(optname, "--noheader" ) == 0 ) {
	showheader = 0;
      }      
    }

  if (argc - optind < 1) Die("%s\n", usage);

  if( ! submat_file ) { 
    if( (t = getenv("SUBOPTDIR")) != 0 || 
	(t = getenv("SUBOPT_DIR")) != 0 ) {
      submat_file = calloc(strlen(t) + 24, sizeof(char));
      sprintf(submat_file, "%s/%s",t,Default_submat);
    } else { 
      submat_file = calloc(strlen((void *)Default_submat) + 24, sizeof(char));
      sprintf(submat_file, "../%s",Default_submat);
    }
  }
  /* open matrix */
  fd = fopen(submat_file, "r");
  
  if( ! ParsePAMFile(fd,&ScoringMatrix, &MatrixScale) ) {
    fprintf(stderr, "Cannot parse or open matrix file %s\n",submat_file);
    free(submat_file);
    exit(EXIT_SUCCESS);
  }
  

  if( output_filename && strlen(output_filename) != 1 &&
      output_filename[0] != '-') {      
    ofd = fopen(output_filename,"w");
    if( ! ofd ) {
      fprintf(stderr, "could not open file %s",output_filename);
      goto end;
    }
  } else 
    ofd = stdout;

  while( optind < argc ) {
    seqfile = argv[optind++];
    
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
      FreeSequence(NULL, &sqinfo);
      seqct++;
    }
    

    cds_seqs = (sequence *)calloc(seqct, sizeof(sequence));
    seqs     = (sequence *)calloc(seqct, sizeof(sequence));
    SeqfileRewind(dbfp);
    seqct=0;

    while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo))
    {
      sqinfo.type = Seqtype(seq);
      if( sqinfo.type == kDNA || sqinfo.type == kRNA ) {

	seqs[seqct].seqstr = Translate(seq,stdcode1);
	/* Let's remove the last codon if it is a stop codon */	
	len = strlen(seqs[seqct].seqstr);
	if( Verbose ) 
	  fprintf(stderr,"seqct is %d length is %d\n",seqct,
		  len);

	if( seqs[seqct].seqstr[len-1] == '*' ) {
	  seqs[seqct].seqstr[len-1] = '\0';
	  seq[strlen(seq) - 3] = '\0';
	}
	cds_seqs[cdsct].seqstr = seq;
	seqs[seqct].seqname = calloc(strlen(sqinfo.name)+1,sizeof(char));
	cds_seqs[cdsct].seqname = calloc(strlen(sqinfo.name)+1,sizeof(char));
	strcpy(seqs[seqct].seqname,sqinfo.name );
	strcpy(cds_seqs[cdsct].seqname,sqinfo.name);	
	cds_seqs[cdsct].length = sqinfo.len;
	cds_seqs[cdsct].alphabet = ( sqinfo.type == kDNA ) ? dna : rna;
	seqs[seqct].length = strlen(seqs[seqct].seqstr);
	
	seqs[seqct].alphabet = protein;
	cdsct++; seqct++;
      } else {
	fprintf(stderr,"Expect CDS sequences (DNA or RNA) not Protein\n");
	goto end;
      }    
      FreeSequence(NULL, &sqinfo);
      if( Verbose && seqct > 3 ) 
	break;
    }
    
    if( seqct < 2 ) {
      fprintf(stderr,"Must have provided a valid file with at least 2 sequences in it");
      goto end;
    }
    
    for( i=0; i  < seqct; i++ ) {
      for(k=i+1; k < seqct; k++ ) {	
	if( (opt_alignment = (alignment *)calloc(1,sizeof(alignment *))) == NULL) {
	  fprintf(stderr,"Could not allocate memory\n");
	  goto end;
	}

	opt_alignment->msa = NULL;
	rc = optimal_align(&seqs[i],&seqs[k],opt_alignment);
  
	if( rc != 1 ) {
	  fprintf(stderr,"Could not make an optimal alignment\n");
	  goto end;
	} else {
	  tmp_cds_seqs[0] = cds_seqs[i];
	  tmp_cds_seqs[1] = cds_seqs[k];
	  rc = mrtrans(opt_alignment, tmp_cds_seqs, &cds_aln,0);
	  if( rc != 0  ) { 
	    fprintf(stderr, "Could not map the coding sequence to the protein alignemnt for aln %d: %d\n",i,rc);
	    goto end;
	  }
	  if( showaln ) {
	    if( ofmt >= 100 ) {
	      MSAFileWrite(ofd,cds_aln->msa, ofmt,do_oneline);
	    } else { 
	      for(j=0; j < cds_aln->msa->nseq; j++ ) {	
		WriteSeq(ofd, ofmt, 
			 cds_aln->msa->aseq[j],
			 &(cds_aln->sqinfo[j]) );
	      }
	    }	    
	  } else {
	    if( showheader && firsttime ) {
	      fprintf(ofd,"SEQ1\tSEQ2\tSCORE\tdN\tdS\tOMEGA\tN\tS\tkappa\tt\tLENGTH\n");
	      firsttime = 0;
	    }
	    if( do_kaks_yn00(cds_aln->msa, &pwMLdist,&pwNGdist) < 0 ) {
	      fprintf(stderr, "warning: problem with align for %s %s\n",
		      cds_aln->msa->sqname[0], cds_aln->msa->sqname[1]);
	      continue;
	    }

	    for(ik = 0; ik < NUM_PW_SEQS; ik++ ) {	  
	      for( jk = ik+1; jk < NUM_PW_SEQS; jk++ ) {
		fprintf(ofd,"%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
			cds_aln->sqinfo[ik].name,
			cds_aln->sqinfo[jk].name,
			opt_alignment->score,
			pwMLdist.dN[ik][jk],pwMLdist.dS[ik][jk], 
			pwMLdist.dNdS[ik][jk],
			pwMLdist.N[ik][jk],
			pwMLdist.S[ik][jk],
			pwMLdist.kappa[ik][jk],
			pwMLdist.t[ik][jk],
			opt_alignment->msa->alen);
	      }
	    }  
	  }
	}
	cleanup_alignment(cds_aln);
	cleanup_alignment(opt_alignment); 
      }
    }
  }
  if( ofd && ofd != stdout )
    fclose(ofd);

  end:
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
  
  return 0;
}
