#!/usr/bin/env perl -w

# allps2tex.pl
# Martin Kilbinger 2004

# Creates a latex file on stdout with all [e]ps[i] files except the file 'all.ps'
# included. 

# v1.2 2008
# v1.3 2010
# v1.4 2012
# v1.5 2016, article -> extarticle, to allow for more fontsizes (new option -s)
# v1.6 2016, does not include <output>.<suff> (output file created by this script), creates LaTeX problem.


use Getopt::Std;
use File::Basename;

%options = ();
getopts("t:w:f:r:B:s:nH:o:h", \%options);

$version = 1.5;

usage() if defined $options{h};

$width  = defined $options{w} ? $options{w} : 8;
$format = defined $options{f} ? $options{f} : "eps";
$angle  = defined $options{r} ? $options{r} : 0;
$bbstr  = defined $options{B} ? ", bb=$options{B}" : "";
$fsize  = defined $options{s} ? $options{s} : 11;

if (defined $options{o}) {
    open($out_fh, ">$options{o}");
} else {
    $out_fh = *STDOUT;
}



# tex header & preambel
print {$out_fh} "% Created by allps2tex.pl v", $version, " (Martin Kilbinger 2004, martin.kilbinger\@cea.fr)\n\n";
#printf {$out_fh} "\\documentclass[a4paper,%dpt]{article}\n\n", $fsize;
printf {$out_fh} "\\documentclass[a4paper,%dpt]{extarticle}\n\n", $fsize;
print {$out_fh} "\\usepackage{epsfig,a4wide}\n";
print {$out_fh} "\\pagestyle{empty}\n";
print {$out_fh} "\\setlength{\\parindent}{0ex}\n\n";
print {$out_fh} "\\newlength{\\pswidth}\n";
print {$out_fh} "\\setlength{\\pswidth}{$width cm}\n\n";
print {$out_fh} "\\begin{document}\n\n";

if (defined $options{t}) {
  my $title = $options{t};

  # I didn't find an easier way for the following
  # (replace '_' with '\_' but leave '\_' alone).

  # Mask underscore with backslash
  $title =~ s/_/\\_/g;
  # Remove double backslashes before underscore
  $title =~ s/\\\\_/\\_/g;
  print {$out_fh} "\\section*{$title}\n\n";
}

if (defined $options{H}) {
    print {$out_fh} "\\input{$options{H}}\n\n";
}

die "Width has to be smaller or equal 17" if $width > 17;
use integer;
$num_col = 17/$width;   # number of columns of images
no integer;

@ALLDIR = $#ARGV > -1 ? @ARGV[0..$#ARGV] : (".");

# process all dirs
while ($DIR = shift @ALLDIR) {

  # look for ??.ps files
  opendir(DIR, $DIR) || die "cannot open $DIR: $!";

  @file = grep {/.*\.e?psi?\z/ && !/all.*\.e?ps/} readdir(DIR) if $format eq "eps" or $format eq "ps";
  @file = grep {/.*\.pdf\z/ && !/all.*\.pdf/} readdir(DIR) if $format eq "pdf";
  @file = grep {/.*\.png\z/} readdir(DIR) if $format eq "png";
  @file = sort(@file);

  foreach $i (0 .. $#file) {

    ($name = $file[$i]) =~ s/(\.$format)\z//g;

    if (($DIR eq ".") and (defined $options{o}) and (basename($options{o}, ".tex") eq $name)) {
        next;
    }

    if (defined $options{n}) {
	    $n = "$DIR/$name";
	    $n =~ s/_/\\_/g;   # '_' -> '\_'
        $n =~ s/\.\///g;   # './' -> ''
	    $n = "$DIR/$n" unless $DIR eq ".";
    }

    print {$out_fh} "\\begin{minipage}{\\pswidth}\n";

    print {$out_fh} "\\centerline{$n}\n\n" if defined $options{n};

    print {$out_fh} "\\includegraphics[width=\\textwidth$bbstr,type=$format,ext=.$format,read=.$format,angle=$angle]{$DIR/$name}\n";
    print {$out_fh} "\\end{minipage}\n";

    print {$out_fh} "\n" if ($i+1)%$num_col==0;

  }

}

print {$out_fh} "\n\\end{document}\n";

sub usage {
    print STDERR "Usage: allps2tex.pl [OPTIONS] [DIR1 [DIR2 [...]]]\n";
    print STDERR "Writes a latex file to stdout with all [e]ps[i] files except 'all.ps' in given directories.'\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -t TITLE     Title string (default none)\n";
    print STDERR "  -w WIDTH     Width of individual postscript files in cm (default 8)\n";
    print STDERR "  -f [eps|pdf] Format of input files (default 'ps', includes eps, ps, epsi)\n";
    print STDERR "  -r ANGLE     Rotate figues by ANGLE degrees (default 0)\n";
    print STDERR "  -B BBOX      'b1 b2 b3 b4'\n";
    print STDERR "  -s FSIZE     Font size [pt], default FSIZE=11\n";
    print STDERR "  -n           Write file name in tex file after including figure\n";
    print STDERR "  -o FNAME     Output filename (default: stdout)\n";
    print STDERR "  -H FNAME     Include header from file FNAME\n";
    print STDERR "  -h           This message\n";
    print STDERR "  DIR1 ...     Directory list (default '.')\n";
    exit 1;
}

