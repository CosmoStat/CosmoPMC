#!/usr/bin/env perl -w

# config_pmc_to_max_and_fish.pl
# Martin Kilbinger 2010
#
# Creates max_post or go_fishing configuration files using the parameter
# and data sections of a pmc config file. Useful for the calculation of
# the Fisher matrix for the initial PMC proposal.



use Fatal qw/ open /;
use Getopt::Std;

### Command line options
%options=();

getopts("c:rf:t:p:MFdh", \%options);

usage() if defined $options{h};
usage() if !(defined $options{M} || defined $options{F});
usage() if defined $options{F} && !(defined $options{f} || defined $options{p});

$def_r = defined $options{r};
$def_f = defined $options{f};
$def_p = defined $options{p};
$sum_rfp = $def_r + $def_f + $def_p;

die "Only one option out of '-r', '-f' and '-p' can be defined" if $sum_rfp>1;
die "Options '-M' and '-F' cannot both be defined" if defined $options{M} && defined $options{F};

$config_pmc = defined $options{c} ? $options{c} : "config_pmc";
$tolerance  = defined $options{t} ? $options{t} : 0.01;

# Minimum parameters
open(CONFIG, $config_pmc);
while (<CONFIG>) {
    print;
    last if /^min/;
}
@Fmin = split(" ", $_);

# Maximum parameters
do {
    $_ = <CONFIG>;
} while (length($_)<=1);
@Fmax = split(" ", $_);
print;

while (<CONFIG>) {
#    last if /nsamples/;
    last if /nchain/ or /nsamples/;
    s/(\#.*)PMC/$1Maximum search (max_post)/ if defined $options{M};
    s/(\#.*)PMC/$1Fisher matrix (go_fishing)/ if defined $options{F};
    print;
}
close CONFIG;

if (defined $options{M}) {

  if (defined $options{r}) {

    # Random starting parameter
    print "sstart          ran\n";

  } else {

    # Fiducial parameters. Deduced parameters are ignored.
    print "sstart          fid\n";
    print "fid              ";

    my $fid = undef;
    # User-specified value
    $fid = $options{f} if defined $options{f};
    # Value from file
    $fid = fid_from_file($options{p}) if defined $options{p};

    if (defined $fid) {
        @par = split(" ", $fid);
	die "Number of parameters (option -f or -p) has to be $#Fmax, not ", $#par+1 if $#par+1 != $#Fmax;
	print "$fid\n";
    } else {
      # Central value
      foreach $i (1 .. $#Fmax) {
	printf " %g", ($Fmax[$i]+$Fmin[$i])/2.0;
      }
      print "\n";
    }

  }

  print "tolerance       $tolerance\n";

} elsif (defined $options{F}) {

  print "fid              ";
  if (defined $options{f}) {
    # User-specified value
    print "$options{f}\n";
  } elsif (defined $options{p}) {
    # Value from file
    $fid = fid_from_file($options{p});
    print "$fid\n";
  }

  print "sinitial        Fisher";
  print "_diag" if defined $options{d};
  print "\n";

}


# Input is either plain text file (with one row of mean values), a 'maxlogP' file or
# a 'mean' file.
sub fid_from_file {
    my ($name) = @_;

    open(IN, $name);
    @F = split(" ", $_ = <IN>);
    if ($_ =~ /max log P/) {
	   $fid = "";
	   foreach $i (10 .. $#F) {
	       $fid = join(" ", $fid, $F[$i]);
	   }
    } elsif ($_ =~ /mean/) {
	   $fid = "";
	   while (<IN>) {
	       next if /#/;
	       @F = split(" ", $_);
	       $fid = join(" ", $fid, $F[2]);
	   }
    } else {
	   $fid = $_;
    }
    close IN;

    return $fid;
}


sub usage {
  print STDERR "Usage: config_pmc_to_max_and_fish.pl [OPTIONS]\n";
  print STDERR "OPTIONS:\n";
  print STDERR "  -M            Create config file for maximum search (max_post)\n";
  print STDERR "  -F            Create config file for Fisher matrix (go_fishing)\n";
  print STDERR "  -c CONFIG     Input PMC config file CONFIG (default: 'config_pmc')\n";
  print STDERR "  -r            Random starting point (for maximum search)\n";
  print STDERR "  -f FID        Fiducial starting point FID. FID is a white-space\n";
  print STDERR "                 separated list in quotes, e.g. '0.25 0.75'\n";
  print STDERR "  -p FILE       Fidcucial parameter from FILE (e.g. 'maxlogP')\n";
  print STDERR "  -t TOLERANCE  Tolerance for maximum-search (default: 0.01)\n";
  print STDERR "  -d            Calculate only diagonal of Fisher matrix (go_fishing)\n";
  print STDERR "  -h            This message\n";
  print STDERR "One of '-M' or '-F' is obligatory\n";
  print STDERR "The default starting point for maximum search is (max-min)/2\n";
  print STDERR "For Fisher matrix ('-F'), a fiducial parameter has to be indicated with '-f FID'\n";
  print STDERR " or '-p FILE'\n";
  exit 1;
}
