#!/usr/bin/env perl
# Write one FASTA file per ortholog pair, ready for yn00_cds_optimal.
#
#   make_ortho_pairfiles.pl --pairs orthologs.tsv \
#       --db species1.cds.fa --db species2.cds.fa [--outdir pairfiles] \
#       [--index /path/to/index]
#
# --pairs  : two whitespace-separated sequence IDs per line; lines starting
#            with '#' or '/' are skipped
# --db     : FASTA file(s) to index (repeatable)
# --outdir : output directory (default: pairfiles)
# --index  : Bio::Index::Fasta index file (default: <outdir>/seqs.idx)
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Bio::Index::Fasta;
use Bio::SeqIO;

my @dbs;
my $pairfile;
my $outdir = 'pairfiles';
my $index;
GetOptions(
    'db=s'     => \@dbs,
    'pairs=s'  => \$pairfile,
    'outdir=s' => \$outdir,
    'index=s'  => \$index,
) or die "bad options\n";
die "--pairs and at least one --db are required\n" unless $pairfile && @dbs;

mkdir $outdir unless -d $outdir;
$index ||= File::Spec->catfile( $outdir, 'seqs.idx' );

my $db = Bio::Index::Fasta->new( -filename => $index, -write_flag => 1 );
$db->make_index(@dbs);

open( my $orth, '<', $pairfile ) || die "cannot open $pairfile: $!";
my $n = 0;
while (<$orth>) {
    next if /^\s*$/ || /^\s*[#\/]/;
    my ( $a, $b ) = split;
    my $out = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => '>' . File::Spec->catfile( $outdir, "${a}_$b.fas" )
    );
    for my $id ( $a, $b ) {
        my $seq = $db->get_Seq_by_acc($id) || $db->get_Seq_by_id($id)
            || die "sequence $id not found in @dbs\n";
        $out->write_seq($seq);
    }
    $n++;
}
print STDERR "wrote $n pair files to $outdir\n";
