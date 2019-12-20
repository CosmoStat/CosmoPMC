#!/usr/bin/env perl

# Creates a histogram from a data file.
# See histogram.pl -h for options.
#
# Martin Kilbinger 2007
#  v2.2   2012
#  v2.3   04/2014 (added column name, read from header)


use warnings;
use lib("$ENV{'HOME'}/perl");
use stuff;
use Fatal qw/ open /;
use Getopt::Std;


%options = ();
getopts("cb:N:nlLm:M:w:o:Hfqh", \%options);

$def_b_N = (defined $options{b}) + (defined $options{N});


usage() if $#ARGV==0 && $#ARGV==1;
usage() if $def_b_N!=1;
usage() if defined $options{h};

$version = 2.3;

if (defined $options{o}) {
  open($out_fh, ">$options{o}");
} else {
  $out_fh = *STDOUT;
}

if (defined $options{f}) {
  $options{l} = 1;
  $options{H} = 1;

  print {$out_fh} "# hist\n";
}

# Header
if (! defined $options{H}) {
  print {$out_fh} "# Created with histogram.pl";
  print {$out_fh} " -c" if defined $options{c};
  print {$out_fh} " -b $options{b}" if defined $options{b};
  print {$out_fh} " -N $options{N}" if defined $options{N};
  print {$out_fh} " -n" if defined $options{n};
  print {$out_fh} " -l" if defined $options{l};
  print {$out_fh} " -L" if defined $options{L};
  print {$out_fh} " -m $options{m}" if defined $options{m};
  print {$out_fh} " -M $options{M}" if defined $options{M};
  print {$out_fh} " -w $options{w}" if defined $options{w};
  print {$out_fh} " -o $options{o}" if defined $options{o};
  print {$out_fh} " -H" if defined $options{H};
  print {$out_fh} " -f" if defined $options{f};
  print {$out_fh} " -q" if defined $options{q};
  print {$out_fh} " $ARGV[0] ";
  print {$out_fh} " $ARGV[1]" if $#ARGV==1;
  print {$out_fh} "\n";
  print {$out_fh} "# histogram.pl (MK) version v$version\n";
}

open(FILE, $ARGV[0]);
$col_in  = $#ARGV==1 ? $ARGV[1] : 0;

# Reading header if present
$header = <FILE>;
if (! ($header =~ /#/)) {
  # No header found
  $header = "";
  seek FILE, 0, 0;
} else {
   %column = %{get_columns($header)};
}

# Getting column number
if ($col_in =~ m/[^0-9]/) {

  # Input column contains non-number characters
  print STDERR "Warning: Option '-c' not necessary for non-number column name\n" if (defined $options{c} and ! defined $options{q});

  die "No header found" unless %column;
  $col = $column{$col_in};
  die "Column name \"$col_in\" not found in header" unless defined $col;
  print STDERR "Creating histogram for column with name \"$col_in\", number $col\n" unless defined $options{q};

} else {

  # Input column is alpha-numeric
  if (defined $options{c}) {

    # Force input column to be name
    die "No header found" unless %column;
    $col = $column{$col_in};
    die "Column name \"$col_in\" not found in header" unless defined $col;
    print STDERR "Creating histogram for column with name \"$col_in\", number $col\n" unless defined $options{q};

  } else {

    # Input column is number
    $col = $col_in;
    print STDERR "Creating histogram for column with number $col\n" unless defined $options{q};

  }


}


# Looking for min and max
if (!defined $options{m} || !defined $options{M}) {

  $first = 1;
  while (<FILE>) {

    next if /#/;

    @F = split(" ", $_);
    next if $#F<$col;

    $x = $F[$col];
    if (defined $options{L}) {
      if ($x==0) { next; }
      $x = log($x);
    }
    if ($first==1) {
      $min = $x;
      $max = $x;
      $first = 0;
    }
    $max = $x if $x>$max;
    $min = $x if $x<$min;

  }

}

$min = $options{m} if defined $options{m};
$max = $options{M} if defined $options{M};

if (!defined $min || !defined $max) {
  die "Error: min or max could not be calculated. File $ARGV[0] may be empty, or does contain less than ", $col+1, " columns.";
}

print STDERR "min, max = $min $max\n" unless defined $options{q};

# Number of bins and bin width
if (defined $options{N}) {
    $binN  = $options{N};
    $delta = ($max-$min)/($binN-1);
} elsif (defined $options{b}) {
    $delta = $options{b};
    $binN  = ($max-$min)/$delta + 1;
}

if ($delta==0) {
  print STDERR "Error: min=max, exiting" if !defined $options{q};
  exit 3;
}

# Count occurences in bins
seek(FILE, 0, 0) or die "error: $!";
$sum = 0.0;
while (<FILE>) {
    next if /#/;
    @F = split(" ", $_);
    next if $#F<$col;
    $x = $F[$col];
    if (defined $options{L}) {
       if ($x==0) { next; }
       $x = log($x);
    }
    $ii = ($x-$min)/$delta;
    $i = int($ii);


    # This may happen if -m or -M are given
    next if $i<0;
    next if $i>$binN-1;

    # Add weight or one if no weights defined
    $inc = defined $options{w} ? $F[$options{w}] : 1.0;
    $hist[$i] += $inc;
    $sum      += $inc;
}
close FILE;

die "File $ARGV[0] may be empty, contain only zeros in requested column, or contain less than ", $col+1, " columns." if $sum == 0;

# Normalise histogram
if (defined $options{n}) {
    foreach $i (0 ..$#hist) {
	if (defined $hist[$i]) {
	    $hist[$i] /= $sum;
	} else {
	    $hist[$i] = 0.0;
	}
    }
}

# Output
$add = defined $options{l} ? 0.0 : $delta/2;
foreach $i (0 .. $binN-1) {
    $x = $i*$delta + $min + $add;
    $h = defined $hist[$i] ? $hist[$i] : 0.0;
    printf {$out_fh} "%10.7g %12.5g\n", $x, $h;
}
printf {$out_fh} "%10.7g %12.5g\n", $x+$delta, 0 if defined $options{l};

# end

sub usage {
    print STDERR "Usage: histogram.pl OPTIONS FILE [COL]\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  FILE             File with list of numbers\n";
    print STDERR "  COL              Column number or name (default 0)\n";
    print STDERR "  -c               Force COL to be name (useful if column name is alphanumeric)\n";
    print STDERR "  -b binwidth      Width of bins\n";
    print STDERR "  -N Nbin          Number of bins\n";
    print STDERR "  -n               Normalized histogram\n";
    print STDERR "  -l               Left corner rather than bin center in output\n";
    print STDERR "  -L               Logarithmic bins\n";
    print STDERR "  -m min, -M max   Minimum and maximum (by default determined from data)\n";
    print STDERR "  -w COL           Weights in column COL\n";
    print STDERR "  -o OUT           Output file (default stdout)\n";
    print STDERR "  -H               No header is printed\n";
    print STDERR "  -f               Output in 'nicaea' histogram format. Sets options '-l -H'\n";
    print STDERR "  -q               Quiet mode\n";
    print STDERR "  -h               This message\n";
    print STDERR "  Either the binwidth or number of bins have to be given\n";
    exit(1);
}
