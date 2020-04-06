#!/usr/bin/env perl

# Martin Kilbinger 2008
# Calculates the Bayes factor for two models. The file evidence should be of the format
#   iter log(sum_i w_i) ...
# where the second column is log of the sum of all (PMC) weights (= log of the evidence).

use warnings;
use Fatal qw/ open /;
use Getopt::Std;


%options=();
getopts("shli:f:", \%options);


usage() if defined $options{h};
usage() if $#ARGV!=1;


if (defined $options{i}) {
  @my_iter = split(" ", $options{i});
  $my_iter[1] = $my_iter[0] if $#my_iter==0;
}

if (defined $options{f}) {
  @filename = split(" ", $options{f});
} else {
  $filename[0] = "evidence";
}
$filename[1] = $filename[0] if $#filename==0;


if (defined $options{l}) {
  goto laplace;
}

# Read evidences
open(EVI, "$ARGV[0]/$filename[0]");
while (<EVI>) {
  next if /#/;
  @E = split(" ", $_);
  $iter = $E[0];
  $loge1[$iter] = $E[1];
}
close EVI;

open(EVI, "$ARGV[1]/$filename[1]");
while (<EVI>) {
  next if /#/;
  @E = split(" ", $_);
  $iter = $E[0];
  $loge2[$iter] = $E[1];
}
close EVI;


if ($#loge1!=$#loge2 && ! defined $options{i}) {
  print STDERR "Error: files have different numbers of iterations, use '-i' option\n";
  exit 2;
}


### Output

# Header
if (! defined $options{s}) {
    print "===================\n";
    print "Model 1: $ARGV[0]\n";
    print "Model 2: $ARGV[1]\n";
    print "===================\n";
    print "iter         B_12       1/B_12   |log10 B_12|   Jeffrey's scale                  |ln B_12|  Modified Jeffrey's scale\n";
    print "--------------------------------------------------------------------------------------------------------------------\n";
} else {
  print "# iter       B_12       1/B_12   |log10 B_12|   Jeffrey's scale                  |ln B_12|  Modified Jeffrey's scale\n";
}


if (defined $options{i}) {
  die "Iteration given by '-i' too large" if $#loge1<$my_iter[0] || $#loge2<$my_iter[1];
  evidence($my_iter[0], $my_iter[1], $loge1[$my_iter[0]], $loge2[$my_iter[1]]);
} else {
  if (! defined $options{s}) {
    foreach $i (0 .. $iter) {
      evidence($i, $i, $loge1[$i], $loge2[$i]);
    }
  } else {
    # Last iteration
    evidence($iter, $iter, $loge1[$iter], $loge2[$iter]);
  }
}


laplace:
if (defined $options{l}) {
  # Evidence from Laplace approximation
    print "--------------------------------------------------------------------------------------------------------------------\n";
    open(EVI1, "$ARGV[0]/evidence_fisher");
    open(EVI2, "$ARGV[1]/evidence_fisher");

    <EVI1>;
    <EVI2>;
    @E1 = split(" ", <EVI1>);
    @E2 = split(" ", <EVI2>);

    evidence($E1[0], $E2[0], $E1[1], $E2[1]);

    close EVI1;
    close EVI2;
}


sub evidence {

    my ($iter1, $iter2, $logE1, $logE2) = @_;
    $log_diff   = ($logE1 - $logE2);
    $B_12       = 10**($log_diff);
    $log_B_12   = $log_diff;
    $ln_B_12    = log($B_12);

    if ($iter1==$iter2) {
      printf "  % d %12.3g %11.3g %14.3g   ",
	$iter1, $B_12, 1/$B_12, abs($log_B_12);
    } else {
      printf "  % d,% d %10.3g %11.3g %14.3g   ",
	$iter1, $iter2, $B_12, 1/$B_12, abs($log_B_12);
    }

    # Log10-based Jeffrey scale
    if (abs($log_B_12)<=0.5) {
	print "Weak ";
    } elsif (abs($log_B_12)<=1) {
	print "Substantial ";
    } elsif (abs($log_B_12)<=2) {
	print "Strong ";
    } else {
	print "Decisive ";
    }

    printf "evidence against Model %d   ", $log_B_12<0?1:2;


    # Ln-based (modified) Jeffrey scale
    printf "%10.3g  ", abs($ln_B_12);

    if (abs($ln_B_12)<=1) {
	print "Inconclusive ";
    } elsif (abs($ln_B_12)<=2.5) {
	print "Weak ";
    } elsif (abs($ln_B_12)<=5) {
	print "Moderate ";
    } else {
	print "Strong ";
    }

    printf "evidence";
    printf " against Model %d   ", $ln_B_12<0?1:2 unless abs($ln_B_12)<=1;

    #print "** $logE1 $logE2 **";

    print "\n";

}

sub usage {
    print STDERR "Usage: bayes_factor.pl [OPTIONS] DIR1 DIR2\n";
    print STDERR "Calculates the Bayes factor between models. The corresponding\n";
    print STDERR " evidence files (from PMC) have to be in the directories DIR1 and DIR2\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -i 'ITER1 [ITER2]'  Use iteration ITER1 for DIR1 and ITER2 for DIR2\n";
    print STDERR "                       (default: all iterations)\n";
    print STDERR "  -f 'EVI1 [EVI2]'    Use files DIR1/EVI1 and DIR2/EVI2 (default: 'evidence')\n";
    print STDERR "  -s                  Short output, last iteration only\n";
    print STDERR "  -l                  Laplace approx. from Fisher matrix (denoted with iter=-1)\n";
    print STDERR "  -h                  This message\n";
    exit 1;
}
