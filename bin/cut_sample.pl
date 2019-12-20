#!/usr/bin/env perl -w

use Fatal qw/ open /;
use Getopt::Std;
use Cwd;

%options = ();
getopts("c:", \%options);

$sample = $ARGV[0];

### Read config file
$configname = "config_pmc" if -e "config_pmc";
$configname = $options{c} if defined $options{c};
die "No configuration file found" unless defined $configname;

open(CONFIG, "$configname");
while (<CONFIG>) {
    @F = split(" ", $_);
    next if /^#/;
    next if $#F==-1;

    if ($F[0] eq "npar") {
	$npar = $F[1];
    }
    if ($F[0] eq "n_ded") {
	$n_ded = $F[1];
    }
    if ($F[0] eq "spar") {
	die "Entry 'npar' not found in configuration file" if ! defined $npar;
	die "Entry 'n_ded' not found in configuration file" if ! defined $n_ded;
	foreach $i (0 .. $npar+$n_ded-1) {
	  #$spar[$i] = $F[$i+1];
	}
    }
    if ($F[0] eq "min") {
	die "Entry 'npar' not found in configuration file" if ! defined $npar;
	die "Entry 'n_ded' not found in configuration file" if ! defined $n_ded;
	foreach $i (0 .. $npar+$n_ded-1) {
	    $min[$i] = $F[$i+1];
	}
    }
    if ($F[0] eq "max") {
	foreach $i (0 .. $npar+$n_ded-1) {
	    $max[$i] = $F[$i+1];
	}
    }
}
close CONFIG;

open(IN, "$sample");
open(OUT, ">$sample" . ".out");
while (<IN>) {

  if (/#/) {
    print OUT;
    next;
  }

  @F = split(" ", $_);
  $outside = 0;
  foreach $i (0 .. $#F-2) {
    if ($F[$i+2]<$min[$i] || $F[$i+2]>$max[$i]) {
      $outside = 1;
      last;
    }
  }
  if ($outside==0) {
    print OUT;
  }
}

close IN;
close OUT;
