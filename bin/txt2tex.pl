#!/usr/bin/env perl

use warnings;
$version = 0.1;

# tex header & preambel
print "% Created by txt2tex.pl v", $version, " (MK 2009)\n\n";
print "\\documentclass[a4paper]{article}\n\n";
print "\\usepackage{epsfig,a4wide}\n";
print "\\usepackage{amsmath}\n";
print "\\usepackage{multirow}\n"; 
print "\\pagestyle{empty}\n";
print "\\setlength{\\parindent}{0ex}\n\n";
print "\\begin{document}\n\n";

while (<>) {
    print;
}

print "\n\\end{document}\n";


