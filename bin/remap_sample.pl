#!/usr/bin/perl -w

# remap_sample.pl
# Martin Kilbinger 2012
# Called by remap.sh
# Usage: remap_sample.pl remap.dat sample


use Fatal qw/open/;


# Read remap.dat file
my $remapdat = $ARGV[0];
open(my $remap_fh, "$remapdat");
my @remap = split(" ", $_ = <$remap_fh>);
chomp @remap;
close $remap_fh;
my $npar = $#remap + 1;


# Remap sample file
open(my $pmc_fh, "$ARGV[1]");

# Read header
my $h = <$pmc_fh>;
($nded) = $h =~ /n_ded\s*=\s*(\d*)/;
$h = <$pmc_fh>;
$h = <$pmc_fh>;
@parname = split(" ", $h);

# Deduced parameters are not in the name list (cosmo_pmc <= v1.11)
foreach my $i (0 .. $nded - 1) {
  push(@parname, "?");
}

# Write new header
print "# npar = $npar, n_ded = 0\n";
print "#         weight            chi2";
foreach my $i (0 .. $npar - 1) {
  print "        param[", $i, "]";
}
print "\n";
print "#                               ";
foreach my $i (0 .. $npar - 1) {
  printf "%16s", $parname[$remap[$i] + 1];
  #print "$i $remap[$i]   ";
}
print "\n";
                
             
# Read and write body
while (<$pmc_fh>) {
  @line = split(" ", $_);

  # weight, -logL
  printf "%16.9g%16.9g", $line[0], $line[1];
  for my $i (0 .. $npar - 1) {
    printf "%16.9g", $line[$remap[$i] + 2];
  }
  printf "\n";
}

close $pmc_fh;
