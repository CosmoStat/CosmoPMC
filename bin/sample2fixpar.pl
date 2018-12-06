#!/usr/bin/perl -w


# sample2fixpar.pl
# Martin Kilbinger 2008
# Marginalizes PMC sample/MCM chain over a parameter with flat prior.


if ($#ARGV!=3 || join(@ARGV) =~ /-h/) {
  print STDERR "Usage: sample2fixpar.pl SAMPLE_IN COL MIN MAX\n";
  print STDERR "    SAMPLE_IN          Input sample (PMC simulation or MCM chain)\n";
  print STDERR "    COL                Column number of fixed parameter\n";
  print STDERR "                        (Note that par #i is in column i+2)\n";
  print STDERR "    MIN, MAX           Minimum and maximum values for fixed parameter\n";
  exit 1;
}

$sample_in = $ARGV[0];
$col       = $ARGV[1];
$par_min   = $ARGV[2];
$par_max   = $ARGV[3];

open(IN, "$sample_in") or die "Could not open file $sample_in: $!";
while (<IN>) {

    (print && next) if /#/;
    @F = split(" ", $_);
    if ($F[$col]>=$par_min && $F[$col]<=$par_max) {
        print;
    }
}

close IN;
