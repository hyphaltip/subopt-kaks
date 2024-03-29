**VERSION** (unreleased, 2004-02-18)
Suboptimal alignments and Ka/Ks calculations

This src is freely available, but should be considered a personal
communication as I have not properly packaged it or given adequete
credit to the code I am reusing here.  Please contact me if you have
questions or intend to use this for large scale analyses so I can be
sure it is performing properly for you. 

To build this:
```
% cd squid
% ./configure 
% make
% cd ../src
% make
% make install
% cd ..

# set the env variable SUBOPTDIR 
% export SUBOPTDIR=`pwd`

# do a test alignment
% ./bin/yn00_cds_optimal --showtable tst/F35G8.1_CBG16160.cds.fas > table

#Run on a pre-aligned file alignment format 
% yn00_cds_prealigned tst/F35G8.1_CBG16160.cds.aln
```
Should build on OSX and Linux fine.  

SQUID is Sean Eddy's library.  Much of the code for alignment is based
on design ideas from Michael Smoot.  The dN,dS estimation is reusing
Z.Yang's yn00 implementation from his PAML package with some
modification to excise it from the PAML data structures.

To use it, I suggest using yn00_cds_optimal which will calculate an
optimal alignment (local alignment by default) and generate the table
of alignments using the --showtable option.  

Helpful hints:

  You need to pass in a file which contains two CDS sequences, introns
  removed, and with the stop codon removed.  I realize the removing of
  the stop codon is annoying, but this is the requirement of PAML and
  I have not added code in these applications yet to detect and strip
  the last stop codon from a CDS.

  The --showtable will print out a tab delimited summary of Ka and Ks
  and Omega as calculated by Yang's YN00 implementation. 

  You can specify the substitution matrix with --matrix option - pass
  in the whole path to the matrix file.

  -h will print help options for an application

  This is not a particularly robust implementation for some things so
  if you get segfaults let me know, I am finding slight differences in
  memory management between OSX and linux and have not been testing
  this on solaris at all at this point.

```
Usage: yn00_cds_optimal [-options] <sequence file>

   Generate the optimal alignments for a pair of CDS sequences using the protein sequence as a template.
   Available options:
   -h      : help; print version and usage info
   -v      : verbose output -- show alignment path states.
 

   --global      : generate global rather than the default local alignments
   --gapopen num : Gap open penalty.  Default is -8.
   --gapext      : Gap extension penalty.  Default is -2.
   --informat    : Format of the input sequence file database.
   --outformat   : Format of the output sequence alignments.
   --output      : specify an output file instead or writing to STDOUT
   --matrix      : Path to substitution matrix.
                   Default is $SUBOPTDIR/matricies/aa/BLOSUM50 or
                    ../matricies/aa/BLOSUM50 in the absence of
                    SUBOPTDIR env variable
   --gapchar     : Provide a different gap char from '-'
   --showtable   : Show the summary Ka/Ks table instead of the alignment.
   --noheader    : Don't display a header for the table
 ```

Jason Stajich
Duke University

