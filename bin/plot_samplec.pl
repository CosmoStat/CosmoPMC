#!/usr/bin/env perl

use Fatal qw/ open /;
use warnings;
use Getopt::Std;
use Cwd;

### Command line options
%options=();
getopts("c:t:e:s:h", \%options);

usage() if $options{h};
usage() if $#ARGV!=0;

$title = defined $options{t} ? $options{t} : " ";
$every = defined $options{e} ? $options{e} : 5;
$size  = defined $options{s} ? $options{s} : 0.1;

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

# Parameter strings
$output = qx(get_spar.pl -c $configname yorick);
@parlist = split("&", $output);

### Create yorick file
open(YOUT, ">samplec.i");

print YOUT "include, \"likeli.i\"\n";
print YOUT "include, \"stuff.i\"\n\n";
print YOUT "read_table, \"$sample\", sample\n";
print YOUT "comp_min = sample(2,min)\n";
print YOUT "comp_max = sample(2,max)\n";
print YOUT "pltitle_height = 20\n\n";

$Ntot = $npar + $n_ded;

for ($i=0; $i<$Ntot-1; $i++) {
    for ($j=$i+1; $j<$Ntot; $j++) {

	$suff = "$i" . "_" . "$j";
	print YOUT "window, display=\"\", hcp=\"samplec_$suff\"\n";
	print YOUT "for (c=comp_max,cc=1; c>=comp_min; c--,cc++) {\n";
	print YOUT "  w = where(sample(2,)==c)\n";
	print YOUT "  if (!is_array(w)) continue\n";
	print YOUT "  plmk, sample($j+3,w(::$every)), sample($i+3,w(::$every)), color=color_gnu(cc%Ncolor_gnu), msize=$size, width=11\n";
	print YOUT "}\n";
	print YOUT "limits, $min[$i], $max[$i]\n" if defined $min[$i];
	print YOUT "range, $min[$j], $max[$j]\n" if defined $min[$j];
	printf YOUT "xytitles, \"%s\", \"%s\", [-0.01, 0.015]\n", $parlist[$i], $parlist[$j];
	print YOUT "pltitle, \"$title\"\n";
	print YOUT "eps, \"samplec_$suff\"\n";
	print YOUT "\n\n";

    }
}

print YOUT "quit\n";
close YOUT;

### Call yorick, tex et al.
print STDERR "Calling yorick...\n";
`yorick -batch samplec.i`;

$pwd = cwd;
`all_vs_all.pl -b samplec -e ps -t $pwd > all_samplec.tex`;
print STDERR "Calling latex...\n";
`ldp.sh all_samplec.tex -q`;
`rm -f all_samplec.{log,dvi,aux}`;



sub usage {
    print STDERR "Usage: plot_samplec.pl [Options] sample\n";
    print STDERR "\nOptions:\n";
    print STDERR "  -t TITLE       Title string (default empty)\n";
    print STDERR "  -e EVERY       Plot only one in EVERY point (default 5)\n";
    print STDERR "  -s SIZE        Point size (default 0.1)\n";
    print STDERR "  -h             This message\n";
    exit 1;
}
