#!/usr/bin/perl -w


# Martin Kilbinger 2010
# Replaces the obsolete script 'fisher_to_meanvar.pl'.


use Fatal qw/ open close unlink /;
use Getopt::Std;
use File::Basename;


%options = ();
getopts("nxmkh", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;
usage(2) if defined $options{n} and defined $options{m};

$fname = $ARGV[0];


open(my $out_fh, ">fishtmp.i");
print {$out_fh} "include, \"stuff.i\"\n";

if (defined $options{x}) {

  print {$out_fh} "mixmvd = mix_mvdens_read(\"$fname\")\n";

  print {$out_fh} "for (i=1; i<=mixmvd.ncomp; i++) {\n";

  print {$out_fh} "  write, format=\"Component %d, weight %.4f\\n\", i, (*(mixmvd.weight))(i)\n";

  print {$out_fh} "  mvd = (*(mixmvd.mvd))(i)\n";
  print_mean($out_fh);

  print {$out_fh} "  cov = *(mvd.var)\n";
  print_cov($out_fh);

  print {$out_fh} "}\n";

} else {

  print {$out_fh} "mvd = mvdens_read(\"$fname\")\n";
  print_mean($out_fh);

  print {$out_fh} "cov = *(mvd.var)\n";
  print_cov($out_fh);

}

print {$out_fh} "quit\n";

close $out_fh;

system("yorick -batch fishtmp.i");

unlink "fishtmp.i" unless defined $options{k};


sub print_mean {
    my ($out_fh) = @_;
    print {$out_fh} "  write, *(mvd.mean), format=\"%8.4f\", linesize=10000\n";
    print {$out_fh} "  writeNL\n";
}

sub print_cov {
  my ($out_fh) = @_;

  if (! defined $options{n} and ! defined $options{m}) {
    print {$out_fh} "c = corr_coeff(cov)\n";
    print {$out_fh} "w = where(abs(c)>1)\n";
    print {$out_fh} "if (is_array(w)) {\n";
    print {$out_fh} "   cov = make_diag(array(0.0, dimsof(cov)(2)))\n";
    print {$out_fh} "} else {\n";
    print {$out_fh} "   cov = LUsolve(cov)\n";
    print {$out_fh} "}";
  }

  if (defined $options{m}) {
    print {$out_fh} "write, 1/sqrt(diag(cov)), format=\"%8.4f\", linesize=10000\n";
  } else {
    print {$out_fh} "write, sqrt(diag(cov)), format=\"%8.4f\", linesize=10000\n";
  }
  print {$out_fh} "writeNL\n";
}

sub usage {
  my ($ex) = @_;

  print STDERR "fisher_to_meanvar.pl [OPTIONS] file\n";
  print STDERR "OPTIONS:\n";
  print STDERR "    -n             No inverse\n";
  print STDERR "    -m             Marginal errors (don't invert matrix)\n";
  print STDERR "    -x             mixmvdens format (default: mvdens format)\n";
  print STDERR "    -k             Keep temporary file 'fishtmp.i'\n";
  print STDERR "    -h             This message\n";
  print STDERR "Options '-m' and '-n' exclude each other\n";

  exit $ex if defined $ex;
}
