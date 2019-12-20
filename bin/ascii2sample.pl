#!/usr/bin/env perl -w

# ascii2sample.pl
# Martin Kilbinger 2013
# Transforms an MCM chain or sample file in raw ascii format
# into a PMC-compatible sample.


use Getopt::Std;

my %options = ();
getopts("c:w:t:ps:N:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;

@columns = defined $options{c} ? split(" ", $options{c}) : ();
foreach my $col (@columns) {
    die "Weight column ($options{w}) can not be contained in parameter column list (@columns)"
        if defined $options{w} and $col == $options{w};
}

die "Options '-p' implies '-w 0'" if defined $options{p} and defined $options{w} and $options{w} != 0;
die "Option '-c' is not compatible with option '-p'" if defined $options{p} and defined $options{c};

my $c_w = (defined $options{w} ? $options{w} : (defined $options{p} ? 0 : -1));
my $offset = ! defined $options{p} ? 0 : 2;
my $wtype = defined $options{t} ? $options{t} : "LOG";
my $n_head = defined $options{N} ? $options{N} : 0;

die "Wrong weight type (options 't') $wtype" unless $wtype eq "LOG" or $wtype eq "LIN";


open(my $in_fh, "$ARGV[0]");
open(my $out_fh, ">$ARGV[0].sample");

my $i = 0;
while (<$in_fh>) {
    next if /#/;

    my @F = split(" ", $_);

    if ($i < $n_head) {
        $i ++;
        next;
    }

    # First line: count parameters, write header
    if ($i == $n_head) {
        if (@columns) {
            $npar = $#columns + 1;
        } else {
            $npar = $#F +1 - $offset;
            @columns = ($offset .. $#F);
        }

        print {$out_fh} "# npar = $npar, n_ded = 0\n";
        printf {$out_fh} "# %14s %15s", "weight", "chi2";
        for my $i (0 .. $npar-1) {
            printf {$out_fh} "%16s", "param[$i]";
        }
        print {$out_fh} "\n";
        if (defined $options{s}) {
            my @names = split(" ", $options{s});
            printf {$out_fh} "# %14s %15s", "", "";
            for my $i (0 .. $npar-1) {
                printf {$out_fh} "%16s", "$names[$i]";
            }
            print {$out_fh} "\n";
        }
    }

    my $weight = $c_w == -1 ? 1 : $F[$c_w];
    #print {$out_fh} "$weight ";
    if ($wtype eq "LIN") {
        next if $weight == 0;
        $weight = log($weight);
    }
    printf {$out_fh} "% 16.8f ", $weight; # This is log(weight)
    my $comp = defined $options{p} ? "$F[1] " : "0 ";
    printf {$out_fh} "%16s ", $comp;     # PMC component or MCMC chi2

    foreach my $col (@columns) {
       #print {$out_fh} "$F[$col] " unless defined $options{w} and $col == $options{w};
       printf {$out_fh} "% 14.9f  ", $F[$col] unless defined $options{w} and $col == $options{w};
    }
    print {$out_fh} "\n";

    $i ++;
}
close $in_fh;
close $out_fh;


sub usage {
    my ($ex) = @_;

    print STDERR "ascii2sample.pl [OPTIONS] FILE\n";
    print STDERR "Transforms an ascii file, e.g. a MCM chain, into a PMC-compatible sample file\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  FILE           Input file\n";
    print STDERR "  -w COL         Use column COL as log(weight) (default: zero log(weight), or COL=0 if '-p' is given)\n";
    print STDERR "  -t WTYPE       Weight type WTYPE = LOG (logarithmic; default), or LIN (linear)\n";   
    print STDERR "  -c 'C0 C1 ...' Use columns C0, C1, ... as parameters (default: all input columns\n";
    print STDERR "                  that not specified as weight columns)\n";
    print STDERR "  -s NAMES       String of parameter names to be written in header\n";
    print STDERR "  -p             File is in PMC sample format, columns: log(weight) comp param_0 param_1 ...\n";
    print STDERR "  -N             N header lines to ignore (default: 0)\n)";
    print STDERR "  -h             This message\n";
    print STDERR "The specification of a weight column implies option '-c'\n";

    exit $ex if defined $ex;
}


