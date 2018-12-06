#!/usr/bin/perl -w

use Fatal qw/ open /;
use Getopt::Std;

%options = ();
getopts("h");

usage(0) if defined $options{h};
usage(1) if $#ARGV!=0;


open(PMCSIM, "$ARGV[0]");

$evi = $n = 0;
while (<PMCSIM>) {
  next if /#/;

  @F = split(" ", $_);
  $logw = $F[0];
  $w    = exp($logw);
  $evi += $w;
  $n ++;
}
close PMCSIM;

$evi = $evi/$n;

printf "# log10(E) ln(E) E\n";
printf "%g %g %g\n", log10($evi), log($evi), $evi;

sub log10 {
  my ($x) = @_;

  return log($x)/log(10.0);
}

sub usage {
  my ($ex) = @_;

  print STDERR "Usage: evidence.pl [OPTIONS] SAMPLE\n";
  print STDERR "OPTIONS:\n";
  print STDERR "   -h         This message\n";
  print STDERR "SAMPLE        PMC sample file\n";

  exit $ex if defined $ex;
}
