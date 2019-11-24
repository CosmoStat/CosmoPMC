#!/usr/bin/env perl

# proposal_means.pl
# Martin Kilbinger 2008
# Reads in proposal files (mix_mvdens) and plots the mean of 
# each component as a function of iteration.


use Getopt::Std;
use File::Basename;
use File::Which;

%options=();
getopts("d:c:niIP:h", \%options);

usage() if defined $options{h};

$dir  = defined $options{d} ? $options{d} : ".";
$inverted = defined $options{i} ? 1 : 0;
$Inverted = defined $options{I} ? 1 : 0;

$base = "proposal";
$out  = "proposal_mean";

my $path_bin = dirname(__FILE__);

# Read proposal files
$iter = 0;
while (1) {
    $name = "$dir/iter" . "_" . "$iter" . "/$base";
    if (! -s $name) {
        die "No proposal file '$name' found" if $iter==0;
        last;
    }

    open(PROP, "$name") or die "Could not open file $name: $!";
    ($ncomp, $ndim) = split(" ", <PROP>);
    foreach $icomp (0 .. $ncomp-1) {
        $weight[$iter][$icomp] = <PROP>;
        $header = <PROP>;
        @tmp = split(" ", <PROP>);
        foreach $idim (0 .. $ndim-1) {
            $mean[$iter][$icomp][$idim] = $tmp[$idim];
        }
        foreach $idim (0 .. $ndim-1) {
            @var = split(" ", <PROP>);
        }
    }

    $iter++;
}
$niter = $iter-1;


foreach $idim (0 .. $ndim-1) {
    foreach $icomp (0 .. $ncomp-1) {
        $valid[$idim][$icomp] = 0;
    }
}

# Write mean values to files
foreach $idim (0 .. $ndim-1) {
    $name = sprintf("%s_%02d", $out, $idim);
    open(OUT, ">$name") or die "Could not create file $name: $!";
    foreach $iter (0 .. $niter) {
        printf OUT "%3d", $iter;
        foreach $icomp (0 .. $ncomp-1) {
            if ($weight[$iter][$icomp]>0) {
                printf OUT " % .6f", $mean[$iter][$icomp][$idim];
                $valid[$idim][$icomp] = 1;
            } else {
                printf OUT " %9s", "-";
            }
        }
        print OUT "\n";
    }
    close OUT;
}



# Read config file for parameter type
$configname = "$dir/config_pmc" if -e "$dir/config_pmc";
$configname = $options{c} if defined $options{c};
die "No configuration file found" unless defined $configname;

$output = qx($path_bin/get_spar.pl -c $configname gnuplot);
@parlist = split("&", $output);

# Create .gnu file and call gnuplot
$gname = "$out" . "." . "gnu";
open (GNU, ">$gname") or die "Could not create file $gname: $!";
print GNU "unset logs\n";
print GNU "unset grid\n";
if (!$Inverted) {
    if (!$inverted) {
        print GNU "set xtics 1\n";
    } elsif ($inverted) {
        print GNU "set ytics 1\n";
        print GNU "set xtics rotate\n";
    }
} elsif ($Inverted) {
    print GNU "unset xtics\n";
    print GNU "unset ytics\n";
    if (!$inverted) {
        print GNU "set xtics 1\n";
        print GNU "set y2tics 1 autofreq\n";
    } elsif ($inverted) {
        print GNU "set ytics 1\n";
        print GNU "set x2tics rotate\n";
    }
}
#print GNU "set term post eps enhanced color 32\n";
print GNU "set term pdfcairo dashed enhanced font 'Arial' fontscale 0.5\n";

foreach $idim (0 .. $ndim-1) {
    if (!$Inverted) {
        if (!$inverted) {
            print GNU "set xlabel 'iteration'\n";
            print GNU "set ylabel '$parlist[$idim]'\n";
        } elsif ($inverted) {
            print GNU "set ylabel 'iteration'\n";
            print GNU "set xlabel '$parlist[$idim]'\n";
        }
    } elsif ($Inverted) {
        if (!$inverted) {
            print GNU "set xlabel 'iteration'\n";
            print GNU "set y2label '$parlist[$idim]'\n";
        } elsif ($inverted) {
            print GNU "set ylabel 'iteration'\n";
            print GNU "set x2label '$parlist[$idim]'\n";
        }

    }
    $name = sprintf("%s_%02d", $out, $idim);
    $name_out = $name if !$inverted;
    $name_out = $name . "_inv" if  $inverted;
    #print GNU "set output '$name_out.eps'\n";
    print GNU "set output '$name_out.pdf'\n";
    print GNU "pl ";
    foreach $icomp (0 .. $ncomp-1) {
        if ($valid[$idim][$icomp]!=0) {
            print GNU "'$name' u 1:", $icomp+2, " w lp lt $icomp lw 3 t ''" if !$inverted;
            print GNU "'$name' u ", $icomp+2, ":(\$1) w lp lt $icomp lw 3 t ''" if $inverted;
            print GNU ", " if $icomp<$ncomp-1;
        }
    }
    print GNU "\n";
    print GNU "set output\n";
}
close GNU;

if (! defined $options{n}) {
    `gnuplot $gname`;
    `$path_bin/allps2tex.pl -f pdf -t "Proposal means" > all_means.tex`;
     exit_if_not_cmd("pdflatex");
    `pdflatex -halt-on-error all_means`;
    `rm -f all_means.log all_means.dvi all_means.aux`;

}

$#var = 0;
$weight = 0;
$header = 0;

# Exits if command does not exist
sub exit_if_not_cmd {

  my ($cmd) = @_;

  $res = which("$cmd");

  die "$cmd not found in search path" unless defined $res;
}


sub usage {
    print STDERR "Usage: proposal_mean.pl [OPTIONS]\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -d DIR         Directory DIR containing the sub-directories 'iter_*'\n";
    print STDERR "                  with the proposal files (default '.')\n";
    print STDERR "  -c CONFIG      Configuration file CONFIG (default 'DIR/config_pmc')\n";
    print STDERR "  -n             No plotting, only creates '.gnu' file\n";
    print STDERR "  -i             x- and y-axes inverted\n";
    print STDERR "  -I             x- and y-labels on top/right\n";
    print STDERR "                  variable \$COSMOPMC)\n";
    print STDERR "  -h             This message\n";
    exit 0;
}
