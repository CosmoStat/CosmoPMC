#!/usr/bin/env perl

# add_par_from_prior.pl
# Martin Kilbinger 2012
# Adds a new parameter to a PMC simulation or sample file.
# The parameter is drawn from a prior distribution.

use warnings;
use Fatal qw/open opendir/;
use Getopt::Std;

my %options = ();
getopts("c:o:p:P:C:s:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;

my $out        = defined $options{o} ? $options{o} : "$ARGV[0].out";
my $prior_dist = defined $options{p} ? $options{p} : "Flat";
my $prior_par  = defined $options{P} ? $options{P} : "-1 1";
my $par_str    = defined $options{s} ? $options{s} : "p0";
my $par_col    = defined $options{C} ? $options{C} : -1;


open(my $in_fh, "$ARGV[0]");
open(my $out_fh, ">$out");


# Header

# npar, n_ded
$_ = <$in_fh>;
my ($npar, $n_ded) = $_ =~ /# npar = (\d*), n_ded = (\d*)/;
$par_col = $npar + 1 if $par_col == -1;
print {$out_fh} "# npar = ", $npar+1, ", n_ded = $n_ded\n";

# weight chi2 param[i]
<$in_fh>;
printf {$out_fh} "# %14s %15s", "weight", "chi2";
for my $i (0 .. $npar) {
    printf {$out_fh} "%16s", "param[$i]";
}
for my $i (0 .. $n_ded-1) {
    printf {$out_fh} "%16s", "param_ded[$i]";
}
print {$out_fh} "\n";

# Parameter names
@F = split(" ", <$in_fh>);
printf {$out_fh} "# %30s", " ";
for my $i (1 .. $#F) {
    printf {$out_fh} "%16s", $par_str if $i-1 == $par_col;
    printf {$out_fh} "%16s", $F[$i];
}
printf {$out_fh} "%16s", $par_str if $#F+1 <= $par_col;
print {$out_fh} "\n";


# Sample points
while (<$in_fh>) {

    my $par_val = random($prior_dist, $prior_par);

    @F = split(" ", $_);

    printf {$out_fh} "%16.9g %15d", $F[0], $F[1];
    for my $i (2 .. $#F) {
	printf {$out_fh} "%16.9g", $par_val if $i-2 == $par_col;
	printf {$out_fh} "%16.9g", $F[$i];
    }
    printf {$out_fh} "%16.9g", $par_val if $#F <= $par_col;
    print {$out_fh} "\n";

}
close $out_fh;
close $in_fh;

sub random {
    my ($dist, $par) = @_;

    my @PAR = split(" ", $par);

    my $r = 1.0e10;

    if ($dist eq "Flat") {

	die "Number of parameters for 'Flat' has to be 2, not ", $#PAR+1, if $#PAR!=1;
	$r = rand ($PAR[1] - $PAR[0]);
	$r = $r + $PAR[0];

    } elsif ($dist eq "Gauss") {

	$r = 0.0;
	
    } else {
	die  "Unknown prior distribution, has to be one in 'Flat', 'Gaussian'";
    }

    return $r;
}




sub usage {
    my ($ex) = @_;

    print STDERR "Usage: add_par_from_prior.pl [OPTIONS] sample\n";
    print STDERR "Adds a new random parameter to a PMC sample file, drawn under a distribution\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -o OUT      Output sample file OUT (default: '<sample>.out'\n";
    print STDERR "  -p DIST     Prior distribution, DIST one of 'Flat' (default), 'Gauss'\n";
    print STDERR "  -P ARG      Prior arguments (white-spaced list if more than one). For DIST =\n";
    print STDERR "               Flat:  ARG = 'min max' (defaut '-1 1')\n";
    print STDERR "               Gauss: ARG = 'mean sigma'\n";
    print STDERR "  -C COL      Column COL of new parameter (default: last)\n";
    print STDERR "  -s STR      Name string STR of new parameter\n";
    print STDERR "  -h          This message\n";

    exit $ex if defined $ex;
}

