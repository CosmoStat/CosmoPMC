#!/usr/bin/env perl -w

# all_vs_all.pl
# Martin Kilbinger 2009
# Creates a triangle plot of all parameters vs. all parameters from file
# base_i_j.ending.


use Getopt::Std;
%options=();
getopts("b:e:l:s:t:dB:S:o:h", \%options);

usage(0) if defined $options{h};

if (!defined $options{b}) {
    print STDERR "all_vs_all.pl: Must specify '-b BASE' option\n";
    print STDERR "all_vs_all.pl: Type 'all_vs_all.pl -h' for more information\n"; 
    exit 1;
}


$base   = $options{b};
$ending = defined $options{e} ? $options{e} : "eps";
$lbase  = defined $options{l} ? $options{l} : "";
$start  = defined $options{s} ? $options{s} : 0;
$title  = defined $options{t} ? $options{t} : "";
$bbstr  = defined $options{B} ? ", bb=$options{B}" : "";
$str    = defined $options{S} ? ", $options{S}" : "";

$version = 1.2;
$author  = "Martin Kilbinger 2008-2010";
$email   = "martin.kilbinger\@universe-cluster.de";

$whole   = 17.5;

$DIR = ".";

opendir(DIR, $DIR) || die "cannot open $DIR: $!";
@files = grep {/($base)_\d*_\d*.($ending)/} readdir(DIR);
@files = sort @files;

$nplot = $#files;

use integer;
#$ncol = $nplot > 0 ? sqrt(2*$nplot+0.25) - 0.5 : 1;
$ncol = sqrt(2*$nplot+0.25) - 0.5 if $nplot > 0;
$ncol = 1 if $nplot==0;
no integer;

if (defined $ncol) {
  $width = $whole/$ncol if !defined $options{l};
  $width = $whole/($ncol+1) if defined $options{l};
} else {
  $width = $whole;
  $ncol  = 0;
}

if (defined $options{o}) {
    open($fh_out, ">$options{o}");
} else {
    $fh_out = *STDOUT;
}

# tex header & preambel
print {$fh_out} "% Created by allps2tex.pl v", $version, " ($author, $email)\n\n";
print {$fh_out} "\\documentclass[a4paper]{article}\n\n";
#print {$fh_out} "\\usepackage{epsfig,a4wide}\n";
print {$fh_out} "\\usepackage{epsfig}\n";
print {$fh_out} "\\pagestyle{empty}\n";
print {$fh_out} "\\setlength{\\parindent}{0ex}\n\n";
print {$fh_out} "\\newlength{\\pswidth}\n";
print {$fh_out} "\\setlength{\\pswidth}{$width cm}\n\n";
print {$fh_out} "\\addtolength{\\textwidth}{3cm}\n\n";
print {$fh_out} "\\begin{document}\n\n";

if (defined $options{t}) {
    $title =~ s/_/\\_/g;
    $title =~ s/(.*MonteCarloTest\/)//g;
    print {$fh_out} "\\section*{$title}\n\n";
}

if (! defined $options{d}) {
    foreach $j ($start .. $start+$ncol) {
	print {$fh_out} "\\hspace*{-2cm}%\n";
	foreach $i ($start .. $j-1) {
	    ($name = "$base $i $j.$ending") =~ s/ /_/g;
	    print {$fh_out} "\\includegraphics[width=\\pswidth$bbstr$str]{$DIR/$name}\n";
	}
	if (defined $options{l}) {
	    ($name = "$lbase $j.$ending") =~ s/ /_/g;
	    print {$fh_out} "\\includegraphics[width=\\pswidth$bbstr$str]{$DIR/$name}\n";
	}
	print {$fh_out} "\n";
    }
} else {
    foreach $i ($start .. $start+$ncol-1) {
	print {$fh_out} "\\hspace*{-2cm}%\n";
	foreach $j ($start .. $i-1) {
	    ($name = "$base $i $i.$ending") =~ s/ /_/g;
	    print {$fh_out} "\\phantom{\\includegraphics[width=\\pswidth$bbstr]{$DIR/$name}}\n";
	}
	foreach $j ($i .. $start+$ncol-1) {
	    ($name = "$base $i $j.$ending") =~ s/ /_/g;
	    print {$fh_out} "\\includegraphics[width=\\pswidth$bbstr]{$DIR/$name}\n";
	}
	print {$fh_out} "\n";
    }
}

print {$fh_out} "\n\\end{document}\n";


#########

sub usage {
    my ($ex) = @_;

    print STDERR "Usage all_vs_all.pl [OPTIONS]\n";
    print STDERR "Creates a tex file for a triangle plot with all combinations BASE_i_j.ENDING\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -b BASE          Basename of files (required)\n";
    print STDERR "  -e ENDING        Default is 'eps'\n";
    print STDERR "  -l 1DBASE        If 1d-plots are present\n";
    print STDERR "  -s START         Start index (default: 0)\n";
    print STDERR "  -t TITLE         Print title string (default: none)\n";
    print STDERR "  -d               Include diagonal (i=j)\n";
    print STDERR "  -B BOUNDINGBOX   'b1 b2 b3 b4'\n";
    print STDERR "  -s STR           Add 'STR' to includegraphics options\n";
    print STDERR "  -o FNAME         Output filename (default: stdout)\n";

    exit $ex if defined $ex;
    exit 2;
}
