#!/usr/bin/perl -w


# Name: fisher_to_meanvar.pl 
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 2010


use Fatal qw/ open close unlink /;
use Getopt::Std;
use Cwd;
use File::Basename;


%options = ();
getopts("nxmkP:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV != 0;
usage(2) if defined $options{n} and defined $options{m};

$fname = $ARGV[0];

$cwd = cwd;
set_ENV_COSMOPMC($cwd);

open(my $out_fh, ">fishtmp.i");
print {$out_fh} "include, \"$ENV{COSMOPMC}/yorick/stuff.i\"\n";

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

sub set_ENV_COSMOPMC {
  my ($cwd) = @_;

  # Copy from command argument
  if (defined $options{P}) {
    $ENV{COSMOPMC} = $options{P};
    return;
  }

  # Environment variable defined (in shell)
  return if defined $ENV{COSMOPMC} and $ENV{COSMOPMC} ne "";

  # Use cwd
  if (-e "$cwd/bin/cosmo_pmc.pl") {
    $ENV{COSMOPMC} = $cwd;
    return;
  }

  die "Set environment variable '\$COSMOPMC' or use option '-P PATH'";
}


sub usage {
  my ($ex) = @_;

  print STDERR "fisher_to_meanvar.pl [OPTIONS] file\n";
  print STDERR "OPTIONS:\n";
  print STDERR "    -n             No inverse\n";
  print STDERR "    -m             Marginal errors (don't invert matrix)\n";
  print STDERR "    -x             mixmvdens format (default: mvdens format)\n";
  print STDERR "    -k             Keep temporary file 'fishtmp.i'\n";
  print STDERR "    -P PATH        Use PATH as CosmoPMC directory (default: environment\n";
  print STDERR "                    variable \$COSMOPMC)\n";
  print STDERR "    -h             This message\n";
  print STDERR "Options '-m' and '-n' exclude each other\n";

  exit $ex if defined $ex;
}
