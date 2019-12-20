#!/usr/bin/env perl -w

# MK 11/2007

use Getopt::Std;
%options=();
getopts("m:", \%options);

if ($#ARGV!=0) {
    print STDERR "Usage: corr_temp [-m imax] chain\n";
    print STDERR "  -m imax: maximum index for correlation function,\n";
    print STDERR "           default: 0.5*chainlenght\n";
    exit 2;
}


### read chain ###
open(CHAIN, "$ARGV[0]") or die "could not open file $ARGV[0]: $!";
$i = 0;
while (<CHAIN>) {
 
    next if /#/;
    @F = split(" ", $_);

    if ($i==0) {
	$npar = $#F-1;
    }

    foreach $c (0 .. $npar-1) {
	$a[$i][$c] = $F[$c+2];
    }

    $i++;
	
}
$nchain = $i;

close CHAIN;


### mean ###
foreach $c (0 .. $npar-1) {
    foreach $i (0 .. $nchain-1) {
	$mean[$c] += $a[$i][$c];
    }
    $mean[$c] /= $nchain;
}


### maximal correlation length ###
if (defined $options{m}) {
    $imax = $options{m};
} else {
    $imax = $nchain/2;
}


### calculate correlation ###
foreach $t (0 .. $imax) {

    print STDERR "$t ";

    foreach $c (0 .. $npar-1) {
	$corr[$t][$c] = 0.0;
	foreach $i (0 .. $nchain-1-$t) {
	    $corr[$t][$c] += ($a[$i][$c] - $mean[$c])*($a[$i+$t][$c] - $mean[$c]);
	}
	if ($nchain-$t>0) { $corr[$t][$c] /= ($nchain-$t); }
	else { $corr[$t][$c] = 0.0; }
    }

}

print STDERR "\n";


# print correlation

foreach $t (0 .. $imax) {

    printf "%6d ", $t;
    foreach $c (0 ..$npar-1) {
	printf " % .4f", $corr[$t][$c]/$corr[0][$c];
    }
    printf "\n";

}


