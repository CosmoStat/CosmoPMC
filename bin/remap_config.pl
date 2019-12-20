#!/usr/bin/env perl -w

# remap_config.pl
# Martin Kilbinger 2009
# Creates a new config file according to a parameter remapping.
# See remap.sh
# Usage: remap_config.pl remap.dat config [npar [n_ded]]


use Fatal qw/open/;


# Read remap.dat file
my $remapdat = $ARGV[0];
open(my $remap_fh, "$remapdat");
my @remap = split(" ", $_ = <$remap_fh>);
chomp @remap;
close $remap_fh;


# Set number of parameters, and number of decuded parameters
if ($#ARGV==1) {
  $npar  = $#remap+1;
  $n_ded = 0;
} elsif ($#ARGV==2) {
  $npar = $ARGV[2];
  $n_ded = 0;
} elsif ($#ARGV==3) {
  $npar  = $ARGV[2];
  $n_ded = $ARGV[3];
}


# Remap config file, write to stdout
open(my $con_fh, "$ARGV[1]");
while (<$con_fh>) {

  s/(npar\s*)(\d+)/$1$npar/;
  s/(n_ded\s*)(\d+)/$1$n_ded/;

  if (/^(spar\s*)/ || /^(min\s*)/ || /^(max\s*)/) {

    print $1;
    my @F = split(" ", $_);
    foreach my $p (0 .. $#remap) {
      print "$F[$remap[$p]+1] ";
    }
    print "\n";

  } elsif (/^(indprior\s*)/) {

    print $1;
    my @F = split(" ", $_);
    foreach my $p (0 .. $#remap) {
      my $val = defined $F[$remap[$p]+1] ? $F[$remap[$p]+1] : 0;
      print "$val ";
    }
    print "\n";

  } else {

    print;

  }

}
close $con_fh;

