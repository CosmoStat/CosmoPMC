#!/usr/bin/env perl

# get_spar.pl
# Martin Kilbinger 2008
# 2012: Option -p added. New default now: print input string, not p<i>
# Prints the names of parameters readable by various (script) languages.
# The input parameters are strings defined in par.h

use warnings;
use Fatal qw/open/;
use Getopt::Std;
use Cwd;
use File::Basename;

getopt("c:i:P:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV==-1;
usage(2) if $#ARGV==0 && !defined $options{c};
usage(3) if $#ARGV>0 && defined $options{c};

$mylang = $ARGV[0];

my $cwd = cwd;
my $path_bin = dirname(__FILE__);

open(IN, "$path_bin/spar.txt");

$ilang = -1;

while (<IN>) {

  @F = split("&", $_);

  if (/#/ && $ilang==-1) {
    foreach $i (1 .. $#F) {
      $lang = trim($F[$i]);
      if ($lang eq $mylang) {
	$ilang = $i;
	last;
      }
    }
  }

  die "Language $mylang not found in spar.txt" if $ilang==-1;

  #print "($F[0]) ($F[$ilang])\n";
  $thispar = trim($F[0]);
  $par{$thispar} = trim($F[$ilang]);
}
close IN;

if (defined $options{c}) {
  open(CONFIG, $options{c});
  while (<CONFIG>) {
    @F = split(" ", $_);
    next if $#F==-1;
    if ($F[0] eq "spar") {
      @allmypar = @F[1..$#F];
      last;
    }
  }
  close CONFIG;
} else {
  @allmypar = @ARGV[1..$#ARGV];
}

$i = 0;
while ($mypar = shift @allmypar) {
  if (defined $options{i} && $options{i}!=$i) {
    # Don't print this parameter
  } else {
     if (defined $par{$mypar}) {
       print "$par{$mypar}";
     } else {
       #print "p$i";
        print "$mypar";
     }
  }
  print "&" if ! defined $options{i};
  $i++;
}
print "\n";


sub usage {
  my ($ex) = @_;

  print STDERR "Usage:\n";
  print STDERR "  get_spar.pl [OPTIONS] LANG PAR1 [PAR2 [...]]\n";
  print STDERR "  get_spar.pl -c CONFIG LANG\n";
  print STDERR "OPTIONS:\n";
  print STDERR "   -c CONFIG          Configuration file CONFIG (default 'config_pmc')\n";
  print STDERR "   -i INDEX	          Returns only par[INDEX]\n";
  print STDERR "   -p                 Print 'p<i> for unknown parameters instead of input string\n";
  print STDERR "   LANG               One of 'yorick', 'gnuplot', 'TeX', 'R'.\n";
  print STDERR "                       More languages can be defined in spar.txt\n";
  print STDERR "   PAR1 ...           Parameter strings\n";

  exit $ex if defined $ex;
}

# Removes whitespace from the start and end of the string
sub trim {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}
