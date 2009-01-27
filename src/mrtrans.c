/* This is based on ideas in Bill Pearson's mrtrans code 
 * (hence the function name) 
 * I have however written this using Sean Eddy's SQUID
 * library
 * Jason Stajich 2003, Duke University 
 *
 * $Id: mrtrans.c,v 1.11 2004/07/01 21:29:31 jason Exp $
 */


#include "mrtrans.h"

/* for our purposes, we assume that the order of the 
 * CDS sequences is the same as the order of the aligned protein
 * sequences
 * remove_last_stop is a boolean to indicate whether or not
 * we should drop the last codon if it is a stop codon
 */

int mrtrans( alignment * protein_align, 
	     sequence * cds_seqs,
	     alignment ** ret_cds_aln,
	     int remove_last_stop) 
{
  int i,j,k, cds_start_offset, cds_end_offset;
  char * aaptr, * cdsptr, *transptr, *tmpendptr;
  /* this is a 'flushed' alignment */
  int protein_aln_len = protein_align->msa->alen; 
  int cds_aln_len = protein_aln_len * CODON_LENGTH; 
  
  int num_seqs    = protein_align->msa->nseq;
  alignment * cds_aln; /* this is where the new alignment will be
			  stored */
  int offset_end = 0,c=0;
  char * lastcodonptr;
  
  if( (cds_aln = (*ret_cds_aln) = calloc(1, sizeof(alignment))) == NULL ) {
    fprintf(stderr,"could not allocate memory for cds_aln\n");
    return -1;    
  }
  
  cds_aln->msa = (MSA *)calloc(1, sizeof(MSA));
  cds_aln->sqinfo = (SQINFO *)calloc(num_seqs, sizeof(SQINFO));

  if( (cds_aln->msa->aseq    = (char **)calloc( num_seqs, sizeof(char *)))
      == NULL) {
    fprintf(stderr,"could not allocate memory for cds_aln->msa->aseq\n");
    return -1;    
  }
  
  if( (cds_aln->msa->sqname  = (char **)calloc( num_seqs, sizeof(char *)))
      == NULL) {
    fprintf(stderr,"could not allocate memory for cds_aln->msa->sqname\n");
    return -1;
  }
  
  if( (cds_aln->msa->sqdesc  = (char **)calloc( num_seqs, sizeof(char *))) 
      == NULL ) {
    fprintf(stderr,"could not allocate memory for cds_aln->msa->sqdesc\n");
    return -1;
  }  

  cds_aln->msa->nseq = num_seqs;  
  cds_aln->msa->alen = cds_aln_len;
  
  for( k =0; k< num_seqs; k++ ) { 
    /* sanity check that sequence names match up */
    if( strcmp(cds_seqs[k].seqname,protein_align->msa->sqname[k]) != 0 ) { 
	fprintf(stderr,"order of CDS sequences and Protein sequences in the alignment do not match, got '%s' expected '%s'\n",cds_seqs[k].seqname,
		protein_align->msa->sqname[k]);
	return -1;
    }
    cds_start_offset = protein_align->sqinfo[k].start * CODON_LENGTH;
    cds_end_offset   = (protein_align->sqinfo[k].stop-1)  * CODON_LENGTH;
    
    lastcodonptr = cds_seqs[k].seqstr+(cds_end_offset-CODON_LENGTH);
    
    /* Okay this needs to be fixed at some point - to use 
     * codon tables since some sp use these stops as real codons
     */
    if( remove_last_stop && 
	( strncmp(lastcodonptr,STOPCODON1,CODON_LENGTH) == 0 || 
	  strncmp(lastcodonptr,STOPCODON2,CODON_LENGTH) == 0 ||
	  strncmp(lastcodonptr,STOPCODON3,CODON_LENGTH) == 0 ) ) {
      offset_end = CODON_LENGTH;
    }
    
    /* Now let's walk through the sequence and convert protein
     * to the appropriate codon triplet from the CDS sequences
     */

    /* allocate space for the cds alignment */
    if( (cds_aln->msa->aseq[k]    = (char *)calloc(cds_aln_len+1, 
						   sizeof(char))) == NULL) {
      fprintf(stderr,"cannot allocate memory for cds_aln->msa->aseq[%d]",k);
      return -1;    
    }
    
    /* allocate the MSA sequence names */
    if( (cds_aln->msa->sqname[k]  = 
	 (char *)calloc(SQINFO_NAMELEN, sizeof(char))) == NULL){
      fprintf(stderr,"cannot allocate memory for cds_aln->msa->sqname[%d]",k);
      return -1;
    }    
    strncpy(cds_aln->msa->sqname[k], cds_seqs[k].seqname,SQINFO_NAMELEN);
    strncpy(cds_aln->sqinfo[k].name, cds_seqs[k].seqname,SQINFO_NAMELEN);
    
    transptr  = cds_aln->msa->aseq[k]+0; /* point to the beginning of
					  * the mem */
    
    cdsptr    = cds_seqs[k].seqstr+cds_start_offset;
    tmpendptr = cds_seqs[k].seqstr+cds_end_offset;
    c = 0;
    for( aaptr = protein_align->msa->aseq[k]+0; 
	 aaptr < protein_align->msa->aseq[k]+protein_aln_len; 
	 aaptr++ ) {

      if( c > cds_end_offset ) {	
	fprintf(stderr,"seqs are %s,%s\n",
		cds_seqs[0].seqname,
		cds_seqs[1].seqname);
	
	Die("cds and pep aren't matching up!");
      }
      
      if( *aaptr == GapChar ) {
	/* add a triplet */
	for( i =0; i < CODON_LENGTH; i++ )
	  *transptr++ = '-';
	
      } else {
	for( i =0; i < CODON_LENGTH; i++ )
	  *transptr++ = *cdsptr++;
	c++;
      }
    }
    
    /* readjust the length in case we dropped the last stop codon */  
    cds_aln->sqinfo[k].start = cds_start_offset;
    cds_aln->sqinfo[k].stop  = cds_end_offset;
  }
  cds_aln->msa->alen -= offset_end; /* we removed the last codon */

  for(k=0;k<num_seqs; k++ ) { 
    cds_aln->sqinfo[k].stop -= offset_end;
    /* allocate the sequence info for writing this in single seq format */
    cds_aln->sqinfo[k].flags = SQINFO_DESC | SQINFO_ID;
    
    sprintf(cds_aln->sqinfo[k].desc, "score=%d alignlen=%d start=%d stop=%d ",
	    protein_align->score,
	    cds_aln->msa->alen,
	    cds_aln->sqinfo[k].start+1,
	    cds_aln->sqinfo[k].stop+1);
    
    if( (cds_aln->msa->sqdesc[k]  = 
	 (char *)calloc(SQINFO_DESCLEN,sizeof(char))) == NULL){
      fprintf(stderr,"cannot allocate memory for cds_aln->msa->sqname[%d]",k);
      return -1;
    }
    strncpy(cds_aln->msa->sqdesc[k], cds_aln->sqinfo[k].desc, 
	    SQINFO_DESCLEN);    
  }
  return 0;
}

