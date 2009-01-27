#!/usr/bin/perl -w

my $dir = shift || ".";
my $ext = shift || "fas";
my $iters = 10_000;

opendir(DIR, $dir) || die($!);
foreach my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.$ext$/;
    `yn00_zuker_align --max_to_show $iters --output $1.dat $dir/$file`;
}
