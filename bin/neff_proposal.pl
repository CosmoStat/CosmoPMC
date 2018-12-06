#!/usr/bin/perl -w

use Fatal qw/ open /;

# neff_proposal.pl
# Martin Kilbinger 2009
# Calculates the effective number of proposal components from the proposal weights,
# in analogy to the effective sample size, WK09 (20)

usage() if $#ARGV!=0;
usage() if $ARGV[0] eq "-h";

open(PROP, "$ARGV[0]");

($Ncomp, $Ndim) = split(" ", <PROP>);


$neff = 0;
for $n (1 .. $Ncomp) {

  $_ = <PROP>;
  chomp;
  $w[$n] = $_;
  <PROP>;   # Component header
  <PROP>;   # Component mean
  for (1 .. $Ndim) {  # Component covariance
    <PROP>;
  }
  $neff += $w[$n]*$w[$n];
}

print "# Ncomp Neff\n";
printf " %d      %.3f\n", $Ncomp, 1.0/$neff;

close PROP;


sub usage {
  print STDERR "Usage: neff_proposal.pl PROP\n";
  print STDERR "   Calculates the effective number of components for the mix_mvdens file 'PROP'\n";
  exit 1;
}
