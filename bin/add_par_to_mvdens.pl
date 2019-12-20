#!/usr/bin/env perl -w

# add_par_to_mvdens.pl
# Martin Kilbinger 2011
#  July 2012:   added mixmvdens support


use Getopt::Std;
use Fatal qw/ open /;

my %options = ();

getopts("xc:m:v:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;

my $par_mean = defined $options{m} ? $options{m} : 0.0;
my $par_var  = defined $options{v} ? $options{v} : 1.0;


# Read (mix)mvdens file
open($in_fh, "$ARGV[0]");

if (defined $options{x}) {

    my ($ncomp, $ndim) = split(" ", <$in_fh>);
    $ndim ++;
    print "$ncomp $ndim\n";

    foreach my $n (1 .. $ncomp) {
        my ($weight) = split(" ", <$in_fh>);
        print "$weight\n";
        do_mvdens_component($in_fh, %options);
    }

} else {

    do_mvdens_component($in_fh, %options);

}

close $in_fh;


sub do_mvdens_component {
    my ($in_fh, %options)  = @_;

    my ($mean_s, $cov_s, $ndim, $df, $B, $c) = read_mvdens($in_fh);
    my @mean = @{$mean_s};
    my @cov  = @{$cov_s};

    my $col = defined $options{c} ? $options{c} : $ndim;
    die "Column $col (option '-c') larger than number of parameters = $ndim\n" if $col > $ndim;


    # Update parameters
    $ndim ++; $B ++;
    print "$ndim $df $B $c\n";

    splice(@mean, $col, 0, $par_mean);
    foreach $i (0 .. $ndim-1) {
        print "$mean[$i] ";
    }
    print "\n";

    foreach $i (0 .. $ndim) {

        if ($i == $col) {
	        foreach $j (0 .. $ndim-1) {
	            print $i == $j ? "$par_var " : "0.0 ";
	        }
	        print "\n";
        }
        last if $i == $ndim - 1;
    
        @row = split(" ", $cov[$i]);
        splice(@row, $col, 0, 0.0);
        foreach $j (0 .. $ndim-1) {
	        print "$row[$j] ";
        }
        print "\n";
    }
}

sub read_mvdens {
    my ($in_sh) = @_;

    my @F = split(" ", <$in_fh>);
    die "Four columns expected, found ", $#F+1, " on line '@F'" if $#F != 3 and $#F != 1;
    die "Four columns expected, found one. If file is in 'mix_mvdens' format, use the '-x' option",
      if $#F == 1;

    my ($ndim, $df, $B, $c) = @F;
    my @mean = split(" ", <$in_fh>);
    my @cov = ();

    foreach my $i (0 .. $ndim - 1) {
        $cov[$i] = <$in_fh>;
    }

    return (\@mean, \@cov, $ndim, $df, $B, $c);

}



sub usage {
    my ($ex) = @_;

    print STDERR "add_par_to_mvdens.pl (MIX)MVDENS [OPTIONS]\n";
    print STDERR "Adds a parameter to a (mix)mvdens file (e.g. Fisher matrix, PMC proposal)\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -c COL     Adds parameter in column and row COL (default: last column)\n";
    print STDERR "  -m VAL     Parameter mean VAL (default 0)\n";
    print STDERR "  -v VAL     Parameter variance VAL (default 1)\n";
    print STDERR "  -x         File is in 'mixmvdens' format\n";
    print STDERR "   FILE      File name\n";
    print STDERR "  -h         This message\n";

    exit $ex if defined $ex;
}
