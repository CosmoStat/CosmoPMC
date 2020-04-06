#!/usr/bin/env perl

# mean2eps.pl
# Martin Kilbinger 2012
# cosmo_pmc preprocessing script
# Creates an eps file with a table containing the parameter means
# and errors from a 'mean' file


use warnings;
use Cwd;
use Fatal qw/ open close unlink /;
use Getopt::Std;
use File::Basename;


my %options = ();
getopts("vc:o:P:h", \%options);


usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;


my $cwd  = cwd;
set_ENV_COSMOPMC($cwd);
die "Environment variable \$COSMOPMC = '$ENV{COSMOPMC}' not valid"
  unless -s "$ENV{COSMOPMC}/bin/cosmo_pmc.pl";


my $config  = defined $options{c} ? $options{c} : "config_pmc";
my $verbose = defined $options{v} ? 1 : 0;
my $outbase = defined $options{o} ? $options{o} : "$ARGV[0]";


### Main program ###

run_own_command("meanvar2tab.pl -s 1 -p 2 -e  -c $config $ARGV[0] > $outbase.tab", $verbose);
run_own_command("tab2tex.pl -s 1.25 $outbase.tab > $outbase" . "_in.tex", $verbose);
run_own_command("txt2tex.pl $outbase" . "_in.tex > $outbase.tex", $verbose);
my $base = basename("$outbase.tex");
my $dir  = dirname("$outbase.tex");
chdir "$dir";
run_own_command("lde.sh $base", $verbose);
chdir "$cwd";


### Subroutines ###

sub set_ENV_COSMOPMC {
  my ($cwd) = @_;

  # Copy from command argument
  if (defined $options{P}) {
    $ENV{COSMOPMC} = $options{P};
    return;
  }

  # Environment variable defined (in shell)
  return if defined $ENV{COSMOPMC};

  # Use cwd
  if (-e "$cwd/bin/cosmo_pmc.pl") {
    $ENV{COSMOPMC} = $cwd;
    return;
  }

  die "Set environment variable '\$COSMOPMC' or use option '-P PATH'";
}

# Runs CosmoPMC program
sub run_own_command {
  my ($command, $verbose) = @_;
  run_command("$ENV{COSMOPMC}/bin/$command", $verbose);
}

# New: Shift exit error code by 8 bits 
sub run_command {
  my ($command, $verbose) = @_;

  my $cwd = cwd;
  print STDERR "Running '$command' in '$cwd'\n" if $verbose;
  `$command`;
  my $ex = $? >> 8;
  die "Last command '$command' unsuccessful (exit code $ex)\n" if $ex;
}


sub usage {
  my ($ex) = @_;

  print STDERR "Usage: mean2eps.pl [OPTIONS] MEAN\n";
  print STDERR "OPTIONS:\n";
  print STDERR "    MEAN           File containing mean and confidence levels (output of\n";
  print STDERR "                    'cosmo_pmc' or 'histograms_sample'\n";
  print STDERR "    -c CONFIG      Uses config file CONFIG (default: 'config_pmc')\n";
  print STDERR "    -P PATH        Use PATH as CosmoPMC directory (default: environment\n";
  print STDERR "                    variable \$COSMOPMC)\n";
  print STDERR "    -o BASE        Outname BASE (default: <MEAN>)\n";
  print STDERR "    -v             Verbose\n";
  print STDERR "    -h             This message\n";

  exit $ex if defined $ex;
}
