#!/usr/bin/perl

use Fatal qw/ open opendir /;
use Getopt::Std;
use Term::ANSIColor;


my %options=();
getopts("d:q:i:o:h", \%options);

usage(0) if defined $options{h};


my $quiet    = defined $options{q} ? 1 : 0;
my $depth    = defined $options{d} ? $options{d} : 1;
my $filename = defined $options{i} ? $options{i} : "pmcsim";
my $outname  = defined $options{o} ? $options{o} : "$filename.cat";
my $DIR      = ".";


if ($depth == 0) {
    # Get '.'
    opendir(my $dir_fh, $DIR);
    @subdir = grep {/^\.$/ && -d "$DIR/$_"} readdir($dir_fh);
    closedir $dir_fh;
} else {
    my @parent = ($DIR);
    foreach (1 .. $depth) {
        my @new = read_subdirs(@parent);
        @parent = @new;
    }
    @subdir = @parent;
}

open(my $fh_out, ">$outname");
my $num = 0;
print STDERR "Looking for $filename:\n" unless $quiet;
foreach my $s (0 .. $#subdir) {

   # Check whether files exist in this subdir
   my @FILE = glob("$subdir[$s]/$filename");
   foreach my $file (@FILE) {

      print STDERR "$file " unless $quiet;
      if (! -s "$file") {
         print STDERR (color("red"), "not found or empty", color("reset"), "\n") if ! $quiet;
         next;
      } else {
         print STDERR (color("green"), "ok", color("reset"), "\n") if ! $quiet;
         $num++;
      }

      open(my $fh, "$file");
      while (<$fh>) {
        if (/#/) {
           if ($num == 1) {
              # Print header of first file
              print {$fh_out} $_;
           }
        } else {
           # Print body
           print {$fh_out} $_;
        }
      }
      close $fh;

   }

}
close $fh_out;


# From stuff.pm
sub read_subdirs {
    my @DIR = @_;

    my @subdirs = ();
    # Parent directories
    @DIR = (".") if $#DIR == -1;

    my $dir_fh = undef;
    foreach $dir (@DIR) {
        opendir($dir_fh, $dir);
        # Alle Unterverzeichnisse
        my @to_add = grep {-d "$dir/$_" && ! /^\.{1,2}$/} readdir($dir_fh);
        foreach $t (0 .. $#to_add) {
            $to_add[$t] = "$dir/$to_add[$t]";
        }
        #print "$dir: @to_add\n";
        closedir($dir_fh);
        @subdirs = (@subdirs, @to_add);
    }

    @subdirs = sort(@subdirs);
    return @subdirs;
}


sub usage {
    my ($ex)  = @_;

    print STDERR "Usage: pmcsim_cat.pl [OPTIONS]\n";
    print STDERR "Concatenates PMC simulations (files 'pmcsim').\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -d DEPTH   Descend to DEPTH subdirectories. Default value is 1. No descend is 0.\n";
    print STDERR "  -i IN      Input file IN (default: 'pmcsim')\n";
    print STDERR "  -o OUT     Output file OUT (default: '<pmcsim>.cat')\n";
    #print STDERR "  -c CONF    Config file CONF (default: not considered)\n";
    print STDERR "  -q         Quiet mode\n";
    print STDERR "  -h         This message\n";

    exit $ex if defined $ex;
}

