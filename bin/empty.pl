#!/usr/bin/env perl -w

# Returns 1 if a file is empty or contains only comment lines
# (with '#' at beginning of line)

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

