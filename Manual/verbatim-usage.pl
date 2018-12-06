#!/usr/bin/perl -w

# verbatim-usage.pl
# Martin Kilbinger 2010

# Usage: verbatim-usage.pl COMMAND


$cmd = $ARGV[0];

print STDERR "\\begin{verbatim}\n";

`$cmd -h >&2`;

print STDERR "\\end{verbatim}\n";

