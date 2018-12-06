#!/usr/bin/perl -w

use Getopt::Std;

### Command line options
%options=();

getopts("h", \%options);

usage() if $#ARGV==-1;
usage() if defined $options{h};


$header = <>;
@F = split(" ", $header);
die "Input file is not in mvdens format" if $#F!=3;
$ndim = $F[0];

$mean = <>;

print $header;
print $mean;

foreach $i (0 .. $ndim-1) {
  @F = split(" ", <>);
  foreach $j (0 .. $ndim-1) {
    print "0.0 " if $i!=$j;
    print "$F[$i] " if $i==$j;
  }
  print "\n";
}

sub usage {
  print STDERR "Usage: diag_mvdens.pl IN\n";
  print STDERR "    Prints the mvdens file 'IN' with the covariance replaced by its diagonal.\n";
  exit 1;
}
