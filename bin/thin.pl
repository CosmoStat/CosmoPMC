#!/usr/bin/perl -w

use Getopt::Std;
%options=();
getopts("t:",\%options);

if ($#ARGV!=0) {
    print STDERR "Usage: thin.pl [-t every] chain\n";
    print STDERR "   Default: every=5\n";
    exit 1;
}

if (defined $options{t}) { $thin = $options{t}; }
else { $thin = 5; }

open(CHAIN, "$ARGV[0]") or die "Could not open chain file $ARGV[0]\n";

$i = 0;
while (<CHAIN>) {
    if (/#/) {print; next;}
    print if $i++%$thin==0;
}
