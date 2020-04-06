#!/usr/bin/env perl

# xip_xim2xi.pl
# Martin Kilbinger 2014
#
# Combined separated xi+ and xi- files, and writes N_corr 'xi' output files.
#
# 04/2015:  Added input, output tomo types as separate options

  
use warnings;
use Fatal qw/ open /;
use Getopt::Std;


### Command line options
%options = ();
getopts("o:s:S:vh", \%options);

usage(0) if defined $options{h};
usage(2) if $#ARGV!=1;


my $tol_frac = 0.01;

my $outbase  = defined $options{o} ? $options{o} : "xi";
my $verbose  = defined $options{v} ? 1 : 0;
my $tomo_out = defined $options{s} ? $options{s} : "all";
my $tomo_in  = defined $options{S} ? $options{S} : "all";

die "Invalid output tomography type $tomo_out (option -s)" unless $tomo_out eq "all" or $tomo_out eq "auto_only" or $tomo_out eq "cross_only";
die "Invalid input tomography type $tomo_in (option -S)" unless $tomo_in eq "all" or $tomo_in eq "auto_only" or $tomo_in eq "cross_only";
die "Input tomography type \'auto_only\' implies the same output type" if $tomo_in eq "auto_only" and !($tomo_out eq "auto_only");
die "Input tomography type \'cross_only\' implies the same output type" if $tomo_in eq "cross_only" and !($tomo_out eq "cross_only");



# Read input files (xi+, xi-)
foreach $fc (0 .. $#ARGV) {
    open(my $in_fh, "$ARGV[$fc]");
    $i = 0;
    while (<$in_fh>) {
	    next if /#/;
        my @F = split(" ", $_);

        if ($i == 0 and $fc == 0) {
            $Nzcorr = $#F;
        } else {
            die "Nzcorr inconsistent in file '$ARGV[$fc]', line $i" if $Nzcorr != $#F;
        }

        $theta[$fc][$i] = $F[0];
        foreach $c (1 ..$#F) {
            $xi[$fc][$i][$c-1] = $F[$c];
        }
	    $i ++;
    }
    close $in_fh;
    $n[$fc] = $i;
}

# Get number of redshif bins
#print STDERR "$Nzcorr + 1 columns found\n" if $verbose;
my $nzbin = -1;
if ($tomo_in eq "all") {
    my $dnzbin = 0.5 * (sqrt(8 * $Nzcorr + 1) - 1);
    use integer;
    $nzbin =  int($dnzbin);
    no integer;
    die "Invalid number of redshift correlations $Nzcorr" if abs($nzbin - $dnzbin) > 0.1;
    print "$nzbin redshift bins found ($nzbin * " . ($nzbin + 1) . " / 2 = $Nzcorr redshift correlations)\n" if $verbose;
} elsif ($tomo_in eq "auto_only") {
    $nzbin = $Nzcorr;
    print "$nzbin redshift bins found, autocorrelation only\n" if $verbose;
} elsif ($tomo_in eq "cross_only") {
    die "cross correlation only on input not yet implemented";
}

# Check files lengths
foreach $z (0 .. $#ARGV-1) {
    my $z1 = $z + 1;
    die "Files '$ARGV[$z]' and '$ARGV[$z1]' have different lengths ($n[$z], $n[$z1])" if $n[$z] != $n[$z+1];
}

# Write xi+ and xi- to Nzcorr files
my $c = -1;
foreach $i_bin (0 .. $nzbin - 1) {
    foreach $j_bin ($i_bin .. $nzbin - 1) {

        next if $i_bin != $j_bin and $tomo_in eq "auto_only";
        next if $i_bin == $j_bin and $tomo_in eq "cross_only";


        $c ++;

        next if $i_bin != $j_bin and $tomo_out eq "auto_only" and $tomo_in eq "all";
        next if $i_bin == $j_bin and $tomo_out eq "cross_only" and $tomo_in eq "all";

        my $outname = "$outbase" . "_" . $i_bin . "_" . $j_bin;
        print "Writing file $outname\n" if $verbose;
        open(my $out_fh, ">$outname");

        # Header
        print {$out_fh} "# theta[arcmin]     xi_p            xi_m\n";

        foreach $i (0 .. $n[0] - 1) {

            # Check angular scales consistency
            foreach $z (0 .. $#ARGV-1) {
                my $z1 = $z + 1;
                die "Files '$ARGV[$z]' and '$ARGV[$z1]' have incompatible angular scales ($theta[$z][$i], $theta[$z1][$i])"
                    if ! defined $options{w} and abs($theta[$z][$i] - $theta[$z1][$i])/$theta[$z][$i] > $tol_frac; 
            }

            # Write theta, xip, xim
            printf {$out_fh} "%.5f ", $theta[0][$i];
            printf {$out_fh} " % .8e % .8e", $xi[0][$i][$c], $xi[1][$i][$c];
            printf {$out_fh} "\n";
        }

        close $out_fh;

    }
}


sub usage {
  my ($ex) = @_;

  print STDERR "xip_xim2xi.pl [OPTIONS] XI_P XI_M\n";
  print STDERR "Transforms xi_p and xi_m files (lensingdemo output) to Nzcorr files of athena 'xi' format\n";
  print STDERR "OPTIONS:\n";
  print STDERR "  -o BASE        Output file name BASE (default 'xi')\n";
  print STDERR "  -s TOMO        For output: Tomography type TOMO, one in 'all', 'auto_only', 'cross_only'\n";
  print STDERR "  -S TOMO        On input: Tomography type TOMO, one in 'all', 'auto_only', 'cross_only'\n";
  print STDERR "  -v             Verbose\n";
  print STDERR "  -h             This message\n";
  print STDERR "  XI_P XI_M      Input files for xi+ and xi-\n";

  exit $ex if defined $ex;
}

