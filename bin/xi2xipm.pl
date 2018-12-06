#!/usr/bin/perl -w

# xi2xipm.pl
# Martin Kilbinger 2009
# Combines xi+ and xi- into one data vector file
# 04/20??: Default name correctly 'XI.pm' (was 'XI.out')
# 06/2014: Renamed from xi2xi.pl to xi2xipm.pl 


use Fatal qw/ open /;
use Getopt::Std;


### Command line options
%options = ();
getopts("o:v:wkKh", \%options);

usage(0) if defined $options{h};
usage(2) if $#ARGV<0;

my $kappa = (defined $options{k} or defined $options{K}) ? 1 : 0;

die "variance for 'kappa' not yet implemented" if defined $options{v} and $kappa;


my $tol_frac = 0.01;


$outname = defined $options{o} ? $options{o} : $ARGV[0] . ".pm";


# Read input file
foreach $fc (0 .. $#ARGV) {
    open(my $in_fh, "$ARGV[$fc]");
    $i = 0;
    while (<$in_fh>) {
	    next if /#/;
	    ($theta[$fc][$i], $xip[$fc][$i], $xim[$fc][$i], $xix, $w, $sqrtD, $sqrtDcor[$fc][$i], $npair) = split(" ", $_);
	    $i ++;
    }
    close $in_fh;
    $n[$fc] = $i;
}

# Check files lengths
foreach $z (0 .. $#ARGV-1) {
    my $z1 = $z + 1;
    die "Files '$ARGV[$z]' and '$ARGV[$z1]' have different lengths ($n[$z], $n[$z1])" if $n[$z] != $n[$z+1];
}

# Write xi+ and xi- to file
open(my $out_fh, ">$outname");

# xi+
foreach $i (0 .. $n[0] - 1) {

    # Check angular scales consistency
    foreach $z (0 .. $#ARGV-1) {
	    my $z1 = $z + 1;
	    die "Files '$ARGV[$z]' and '$ARGV[$z1]' have incompatible angular scales ($theta[$z][$i], $theta[$z1][$i])"
	      if ! defined $options{w} and abs($theta[$z][$i] - $theta[$z1][$i])/$theta[$z][$i] > $tol_frac; 
    }

    printf {$out_fh} "%.5f ", $theta[0][$i];
    foreach $z (0 .. $#ARGV) {
        my $val = $kappa ? $xip[$z][$i] / 2 : $xip[$z][$i]; 
	    printf {$out_fh} " % .8e", $val;
    }
    printf {$out_fh} "\n";
}

# xi-
if (! defined $options{K}) {
    foreach $i (0 .. $n[0] - 1) {
        printf {$out_fh} "%.5f ", $theta[0][$i];
        foreach $z (0 .. $#ARGV) {
            my $val = $kappa ? $xip[$z][$i] / 2 : $xim[$z][$i];
	        printf {$out_fh} " % .8e", $val;
        }
    printf {$out_fh} "\n";
    }
}
close $out_fh;


# Write variance file (variance = sqrt(D)^2)
if (defined $options{v}) {
  open(my $out_fh, ">$options{v}");
  foreach $i (0 .. $n[0] - 1) {
      printf {$out_fh} "%.5f ", $theta[0][$i];
      foreach $z (0 .. $#ARGV) {
	    printf {$out_fh} " % .8en", $sqrtDcor[$z][$i]*$sqrtDcor[$z][$i];
      }
  }
  close $out_fh;
}

$xix = $w = $npair = $sqrtD = 0;

sub usage {
  my ($ex) = @_;

  print STDERR "xi2xipm.pl [OPTIONS] XI [XI2 [...]]\n";
  print STDERR "Transforms a athena 'xi' format file into cosmo_pmc 'xipm' format file\n";
  print STDERR "OPTIONS:\n";
  print STDERR "  -o NAME        Output file NAME (default XI.pm)\n";
  print STDERR "  -v VAR_OUT     Ouput of variance (Poisson-term) to VAR_OUT\n";
  print STDERR "  -k             Interpret xi as kappa correlation:\n";
  print STDERR "                  Only use 'xi+' and divide by 2.\n";
  print STDERR "  -K             Same as '-k', but create 'xiE' format output\n";
  print STDERR "  -w             Do not check for compatible angular scales (e.g. weighted scales\n";
  print STDERR "                  from 'athena -w' need not be the same)\n";
  print STDERR "  -h             This message\n";
  print STDERR "  XI [XI2 [...]] Input file (athena 'xi' format), more than one file for tomography\n";

  exit $ex if defined $ex;
}
