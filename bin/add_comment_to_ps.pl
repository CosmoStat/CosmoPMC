#!/usr/bin/env perl

# add_comment_to_ps.pl
# Martin Kilbinger 2011


use Fatal qw/ open /;
use warnings;
use Cwd;
use Getopt::Std;


my %options = ();
getopts("o:l:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;


open(my $ps_fh, "$ARGV[0]");

my $host = `uname -n`;
chomp($host);
my $pwd  = cwd;
my $user = $ENV{USER};

my $lang = defined $options{l} ? $options{l} : 'g';
if ($lang eq 'g') {
    $str = "%%Creator";
} elsif ($lang eq 'y') {
    $str = "%%For";
} else {
    die "Wrong language $lang";
}

my $out_fh = *STDOUT;
open($out_fh, ">$options{o}") if defined $options{o};

my $found = 0;

while (<>) {

    if (/$str/ and $found == 0) {
        print "$str: $user\@$host:$pwd\n";
	    $found = 1;
    }

    print;

}

close $ps_fh;
close $out_fh if defined $options{o};

sub usage {
    my ($ex) = @_;

    print STDERR "Usage:add_comment_to_ps.pl [OPTIONS] FILE\n";
    print STDERR "Adds username, hostname and current directory as\n";
    print STDERR " comment to a ps file.\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  FILE        Input ps file\n";
    print STDERR "  -l LANG     Language LANG:\n";
    print STDERR "               'g': gnuplot (default)\n";
    print STDERR "               'y': yorick\n";
    print STDERR "  -o OUT      Output file OUT (default stdout)\n";
    print STDERR "  -h          This message\n";

    exit $ex if defined $ex;
}


