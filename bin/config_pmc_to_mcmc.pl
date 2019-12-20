#!/usr/bin/env perl -w

# Martin Kilbinger 2008-2010
#
# Creates a cosmo_mcmc configuration file using the parameter and data
# sections of a pmc config file. To be used with cosmo_mcmc to calculate
# the Fisher matrix for the initial pmc proposal.
#
# Usage: config_pmc_to_mcmc.pl config_pmc > config_mcmc


use Fatal qw/ open /;
use Getopt::Std;

### Command line options
%options=();

getopts("rf:h", \%options);

usage() if defined $options{h};

die "Options 'r' and 'f' cannot both be defined" if defined $options{r} && defined $options{f};

# Minimum parameters
while (<>) {
    print;
    last if /^min/;
}
@Fmin = split(" ", $_);

# Maximum parameters
do {
    $_ = <>;
} while (length($_)<=1);
@Fmax = split(" ", $_);
print;

while (<>) {
    last if /nsamples/;
    s/PMC/MCMC/;
    print;
}
print "nchain          0\n";
print "ncov            0\n";
print "fburnin         0.0\n";
print "ndecorr         0\n";
print "fudge           2.4\n";
print "sinitial        Fisher\n";
if (defined $options{r}) {
  print "sstart          ran\n";
} else {
  print "sstart          fid\n";
  # Fiducial parameters. Deduced parameters are ignored by cosmo_mcmc.
  print "fid              ";
  if (defined $options{f}) {
    print "$options{f}\n";
  } else {
    foreach $i (1 .. $#Fmax) {
      printf " %g", ($Fmax[$i]+$Fmin[$i])/2.0;
    }
    print "\n";
  }
}
print "nbinhist        64\n";

sub usage {
  print STDERR "Usage: config_pmc_to_mcmc.pl [OPTIONS] <config_pmc>\n";
  print STDERR "\nOPTIONS:\n";
  print STDERR "  -r            Start cosmo_mcmc with random point\n";
  print STDERR "  -h            This message\n";
  exit 1;
}
