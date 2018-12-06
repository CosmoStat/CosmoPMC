#!/usr/bin/perl -w

use Fatal qw/ open /;

if ($#ARGV!=2) {
    usage();
    exit(1);
}

open(IN, "$ARGV[0]");
$zmin = $ARGV[1];
$zmax = $ARGV[2];
while (<IN>) {
    (print && next) if /#/;
    ($z, $n) = split(" ", $_);
    next if ($z<$zmin);
    if ($z==$zmax) {
	print "$z 0\n";
	exit 0;
    }
    print;
}
close IN;
$n = 0;



sub usage {
    print STDERR "Usage: cut_nofz_hist.pl file zmin zmax\n";
}
