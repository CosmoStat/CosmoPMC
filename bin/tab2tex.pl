#!/usr/bin/perl -w
#
# tab2tex.pl
# Martin Kilbinger
# Writes an ascii table to tex format
# The first row with leading comment character ('#')
# contains the column description

# Temp version! For ValP


use Getopt::Std;
%options=();
getopts("abs:S:ml:L:t:r:h",\%options);

usage(0) if defined $options{h};

$bare    = defined $options{b} ? 1 : 0;
$array   = defined $options{a} ? 1 : 0;
$stretch = defined $options{s} ? $options{s} : 1.0;
$fsize   = $options{S} if defined $options{S};
$math    = defined $options{m} ? 1 : 0;
$lines_r = defined $options{l} ? "$options{l}" : "a";
$lines_c = defined $options{L} ? "$options{L}" : "a";
$chr     = defined $options{r} ? $options{r} : " ";


print "{\\Large \\textbf{$options{t}}}\\bigskip\n\n" if defined $options{t};

print "{\\$fsize\n" if defined $fsize;
print "\\renewcommand{\\arraystretch}{$stretch}\n";

if (! $array ) {
  if (! $bare ) {
    print "\\begin{tabular}";
  } 
} else {
  if (! $bare ) {
    print "\\begin{equation}\n";
    print "\\left(\n";
  }
  print "\\begin{array}";
}

$zeile = 0;

while (<>) {

    #@a = split(' ', $_, 9999);
    @a = split(/\s+\s+/, $_, 9999);

    # First line (header) = table description
    if ($zeile==0) {
        $Ntab = /#/ ? $#a-1 : $#a;  # Number of columns. One less if header line (starting with '#')

        if (! $bare) {
            print "{";
            for ($j=0;$j<$Ntab;$j++) {
                if (! $array && ! ($lines_c =~ "n")) {
                    print "|" unless $array;
                }
                print "l";
            }
            if (! $array && ! ($lines_c =~ "n")) {
                print"|";
            }
            print "}";
            if (! $array && ! ($lines_r =~ "n")) {
                # Double line before header
	            print "\\hline\\hline";
            }
            print "\n";
            print"\\rule[-3mm]{0em}{8mm}" if /#/;
        }

	    $zeile ++ unless /#/;  # Increase line number if no header line
    }

    my @list_multi = ();

    for ($i=0; $i<=$#a-1-!$zeile;$i++) {
        $spalte = $a[$i+!$zeile];
        $spalte =~ s/$chr/ /g unless $spalte =~ /\$/;
        $spalte =~ s/e(.{1})0?(\d{1,2})/\\cdot10\^{$1$2}/g if $spalte =~ /\$/;
        print "\$" if $math;
        print "$spalte";
        print "\$" if $math;
        print "\t";
        print " &" unless $i==$#a-1-!$zeile;

        if ($spalte =~ /multirow/) {
            push @list_multi, $i+1;
        }
    }

    if ($zeile!=0) {
        for ($i=$#a; $i<$Ntab; $i++) {
            print " &";
        }
    }
    print "\\\\";
    if ( !$array && ( ($lines_r =~ "a" || ($zeile == 0 && $lines_r =~ "h") ) ||
                !$lines_r =~ "n" ) ) {
        #print "\\hline";
        for ($i=1; $i<=$Ntab; $i++) {
            print "\\cline{$i-$i}" unless $i ~~ @list_multi;
        }
# After-header line
        print "\\hline" if !$zeile;
    }
    print "\n";

    # Add space to first body row if after-header line was drawn
    print"\\rule[-1mm]{0em}{5mm}" if ($lines_r =~ "h" && !$zeile && ! $array);

    $zeile++;
}

if (! $array ) {

  # Final line
  print "\\hline\n" if $lines_r =~ "h";

  if (! $bare ) {
    print "\\end{tabular}\n";
    #print "\\end{center}\n";
  }
} else {
  print "\\end{array}\n"; 
  if (! $bare ) {
    print "\\right)\n";
    print "\\end{equation}\n";
  }
}

print "}\n" if defined $fsize;


sub usage {
  my ($ex) = @_;

  print STDERR "Usage: tab2tex.pl [OPTIONS] file\n";
  print STDERR "\nOPTIONS:\n";
  print STDERR "   -a           Produce tex array, not tex table\n";
  print STDERR "   -b           Bare output, no table/array header\n";
  print STDERR "   -s STRETCH   Set arraystretch to STRETCH\n";
  print STDERR "   -S SIZE      LaTeX font size, e.g. 'small'\n";
  print STDERR "   -m           Add '\$' around entries (tex inline math mode)\n";
  print STDERR "   -l MODE      Print vertical lines between rows according to MODE;\n";
  print STDERR "                 a    all lines (default)\n";
  print STDERR "                 n    no lines\n";
  print STDERR "                 h    header lines\n";
  print STDERR "   -L MODE      Print horizontal lines between columns according to MODE:\n";
  print STDERR "                 a    all lines (default)\n";
  print STDERR "                 n    no lines\n";
  print STDERR "   -t TITLE     Title string TITLE\n";
  print STDERR "   -r CHR       Replace character CHR with space in output (unless in math mode)\n";
  print STDERR "   -h           This message\n";

  exit $ex if $ex>=0;
}
