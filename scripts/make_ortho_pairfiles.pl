#!/usr/bin/perl -w
use strict;

use Bio::Index::Fasta;
use Bio::SeqIO;
my $basedir = '/Users/jason/Projects/briggsae/briggsae_synteny/jason';
my $db = new Bio::Index::Fasta(-filename => '/tmp/wormpep.idx',
			       -write_flag=>1);
$db->make_index("$basedir/2_00/cb25.hybrid.pep2.00",
		 "$basedir/2_00/ws77.pep2.00_longest_splices");
my $outdir    = 'pairfiles';
my $orthofile = "$basedir/ortholog_analysis/input_files/orthologs-2.00";

open(ORTH, $orthofile) || die($!);

while(<ORTH>) {
    next if( /^\// || /^\s+/ || /^\s*\#/ );
    my ($ce,$cb) = split;
    my $out = new Bio::SeqIO(-file => ">$outdir/$ce\_$cb.fas");
    for my $g ( $ce,$cb) {
	$out->write_seq($db->get_Seq_by_acc($g) || die($g));
    }
}
