#!/usr/bin/perl -w

use Fatal qw/open/;

while (@ARGV) {

   $last = "";
   $count = 0;
   $all = 0;
   $ncomment = 0;

   open(IN, "$ARGV[0]");
   while (<IN>) {
       if (/#/) { $ncomment++; next; }
       $count++ if $_ ne $last;
       $all++;
       $last = $_;
   }
   close IN;

   printf "%d %d %.5f\n", $count, $all, $count/$all;
   printf STDERR "$ncomment comment lines\n";

   shift;
}
