#!/usr/bin/perl -w

# Rescales the evidence if the prior changes. This is only correct if the posterior
# is close to zero near the (old and new) prior boundaries.
# The prior is p = 1/V with V = Product[th_min-th_min].
# To rescale the posterior P from p1 to p2:
# P2 = P1 * p2/p1 = P1 * V1/V2.

use Getopt::Std;
%options=();
getopts("c:C:e:h", \%options);

if ($#ARGV!=1 || defined $options{h}) {
    print STDERR "Usage: rescale_evidence.pl [OPTIONS] dir1 dir2\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -c       Original config file (default dir1/config_pmc)\n";
    print STDERR "  -C       Rescaled config file (default dir2/config_pmc)\n";
    print STDERR "  -e       Original evidence file (default dir1/evidence)\n";
    print STDERR "  -h       This message\n";
    exit 1;
}

$dir1 = $ARGV[0];
$dir2 = $ARGV[1];

$config1 = defined $options{c} ? $options{c} : "$dir1/config_pmc";
$config2 = defined $options{C} ? $options{C} : "$dir2/config_pmc";
$evi1    = defined $options{e} ? $options{e} : "$dir1/evidence";

($npar1, $vol1, $ndata1) = prior_volume($config1);
($npar2, $vol2, $ndata2) = prior_volume($config2);

die "Inconsistent number of parameters ($npar1, $npar2)" if $npar1!=$npar2;
die "Inconsistent number of data sets ($ndata1, $ndata2)" if $ndata1!=$ndata2;

$f = ($vol1/$vol2)**$ndata1;

print STDERR "Correction: V1=$vol1, V2=$vol2, f=vol1/vol2=$f\n";

open(EVI1, "$evi1") or die "Error for $evi1: $!";

while (<EVI1>) {
    (print && next) if /#/;

    @F = split(" ", $_);

    printf "%d % g % g % g\n", $F[0], $F[1]+log($f)/log(10.0), $F[2]+log($f), $F[3]*$f;
}
close EVI1;


sub prior_volume {

    my ($name) = @_;

    open(CON, "$name") or die "Error for $name: $!";
    while (<CON>) {

	next if /#/;

	@F = split(" ", $_);

	$ndata = $F[1] if /ndata/;
	$npar  = $F[1] if /npar/;

	if (/^min\b/) {
	    foreach $i (0 ..$npar-1) {
		$min[$i] = $F[$i+1];
	    }
	}
	if (/^max\b/) {
	    foreach $i (0 ..$npar-1) {
		$max[$i] = $F[$i+1];
	    }
	}
    }
    close CON;

    $vol = 1;
    foreach $i (0 .. $npar-1) {
	$vol *= ($max[$i] - $min[$i]);
	#print "$i $vol $max[$i] $min[$i]\n";
    }

    ($npar, $vol, $ndata);
}



