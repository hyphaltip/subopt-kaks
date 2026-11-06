#!/usr/bin/env perl
# Run yn00_zuker_align on every CDS pair file in a directory.
#
#   run-dir.pl [--dir DIR] [--ext fas] [--iters 10000] [--outdir DIR]
#              [--bin /path/to/yn00_zuker_align] [-- extra yn00_zuker_align args]
#
# Each DIR/<name>.<ext> produces <outdir>/<name>.dat.
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $dir    = '.';
my $ext    = 'fas';
my $iters  = 10_000;
my $outdir;
my $bin    = 'yn00_zuker_align';
GetOptions(
    'dir=s'    => \$dir,
    'ext=s'    => \$ext,
    'iters=i'  => \$iters,
    'outdir=s' => \$outdir,
    'bin=s'    => \$bin,
) or die "bad options\n";
$outdir ||= $dir;
mkdir $outdir unless -d $outdir;

opendir( my $dh, $dir ) || die "cannot open $dir: $!";
for my $file ( sort readdir($dh) ) {
    next unless $file =~ /^(.+)\.\Q$ext\E$/;
    my $base = $1;
    my @cmd  = (
        $bin, '--max_to_show', $iters,
        '--output', File::Spec->catfile( $outdir, "$base.dat" ),
        @ARGV, File::Spec->catfile( $dir, $file )
    );
    system(@cmd) == 0 or warn "failed: @cmd\n";
}
closedir($dh);
