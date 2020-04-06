#!/usr/bin/env perl

# Returns 1 if a file is empty or contains only comment lines
# (with '#' at beginning of line)
use warnings;

if ($#ARGV!=0) {
  die "Usage: empty.pl file";
}

open(IN, "$ARGV[0]") or exit 1;

while (<IN>) {
  if (!/^#/) {
    exit 0;
  }
}
exit 1;

