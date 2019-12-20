#!/usr/bin/env perl

# Martin Kilbinger 2008
# Shortens dashes and gaps for dotted style for
# yorick ps output files
use warnings;

while (<>) {

	s/82.5 \] \[ 4.5 61.5 \] \[ 82.5 39.0 4.5 39.0/42.5 \] \[ 4.5 31.5 \] [ 42.5 39.0 4.5 39.0/;
	s/82.5 39.0 4.5 39.0 4.5 39.0/42.5 39.0 4.5 39.0 4.5 39.0/;
	print;
}
