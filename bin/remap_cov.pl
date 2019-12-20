#!/usr/bin/env perl

# remap_cov.pl
# Martin Kilbinger 2008

use warnings;
use Fatal qw/ open /;

if ($#ARGV!=1) {
    print STDERR "Usage: remap_cov.pl remap.dat cov\n";
    exit 2;
}

$remapdat = $ARGV[0];
open(REMAP, "$remapdat");

@remap = split(" ", $_ = <REMAP>);
chomp @remap;
close REMAP;

$cov = $ARGV[1];
open(COV, "$cov");

@header = split(" ", <COV>);
print $#remap+1, " $header[1] ", $#remap+1, " $header[3]\n";

@mean   = split(" ", <COV>);
foreach $i (0 .. $#remap) {
    print "$mean[$remap[$i]] ";
}
print "\n";


$l = 0;
while (<COV>) {
    @line = split(" ", $_);
    foreach $c (0 .. $#line) {
	    $cov[$l][$c] = $line[$c];
    }
    $l++;
}

foreach $i (0 .. $#remap) {
    foreach $j (0 .. $#remap) {
	    #print STDERR "*** $i $j *** $remap[$i] $remap[$j]\n";
	    print "$cov[$remap[$i]][$remap[$j]] ";
    }
    print "\n";
}

