#!/usr/bin/perl -w

# proposal_var.pl
# Martin Kilbinger 2010
# Plots the variance of proposal files as function
# of iteration.

use Getopt::Std;
use Fatal qw/ open close unlink /;
use Cwd;

getopts("c:d:P:h", \%options);

usage(0) if defined $options{h};

$path_flag  = defined $options{P} ? "-P $options{P}" : "";
$dir        = defined $options{d} ? $options{d} : ".";
my $config  = defined  $options{c} ? $options{c} : "$dir/config_pmc";
my ($niter, $ncomp, $ndim) = get_info($config);
my $out = "proposal_var";


my $cwd = cwd;
set_ENV_COSMOPMC($cwd);
$yok = check_for_yorick();
die "'yorick' not found" if $yok != 0;

# Get variances with yorick calls
foreach $iter (0 .. $niter - 1) {

  `$ENV{COSMOPMC}/bin/fisher_to_meanvar.pl -x -n $dir/iter_$iter/proposal > itmp.tmp`;
  die "Proposal not valid" unless -s "itmp.tmp";
  open(my $in_fh, "itmp.tmp");
  foreach $icomp (0 .. $ncomp - 1) {

    # Header with weight
    $_ = <$in_fh>;
    @F = split(" ", $_);
    $weight[$iter][$icomp] = $F[3];

    $_ = <$in_fh>;		# Means

    # Variances
    $_ = <$in_fh>;
    my @F = split(" ", $_);
    foreach $idim (0 .. $ndim - 1) {
      $var[$iter][$icomp][$idim] = $F[$idim];
    }
  }

  close $in_fh;

}

foreach $idim (0 .. $ndim - 1) {
  foreach $icomp (0 .. $ncomp - 1) {
    $valid[$idim][$icomp] = 0;
  }
}

# Write variances to files
foreach $idim (0 .. $ndim - 1) {
  $name = sprintf("%s_%02d", $out, $idim);
  open(my $out_fh, ">$name");
  foreach $iter (0 .. $niter - 1) {
    printf {$out_fh} "%3d", $iter;
    foreach $icomp (0 .. $ncomp - 1) {
      if ($weight[$iter][$icomp]>0) {
	printf {$out_fh} " % .6f", $var[$iter][$icomp][$idim];
	$valid[$idim][$icomp] = 1;
      } else {
	printf {$out_fh} " %9s", "-";
      }
    }
    print {$out_fh} "\n";
  }
  close $out_fh;
}


# Read config file for parameter type
$cmd = "$ENV{COSMOPMC}/bin/get_spar.pl -c $config $path_flag gnuplot";
print "*** Running $cmd***\n";
$output = qx($cmd);
@parlist = split("&", $output);

# Create .gnu file and call gnuplot
$gname = "$out" . "." . "gnu";
open (GNU, ">$gname") or die "Could not create file $gname: $!";
print GNU "unset logs\n";
print GNU "unset grid\n";
print GNU "set xtics 1\n";
#print GNU "set term post eps enhanced color 32\n";
print GNU "set term pdfcairo dashed enhanced font 'Arial' fontscale 0.5\n";

foreach $idim (0 .. $ndim-1) {
  print GNU "set xlabel 'iteration'\n";
  print GNU "set ylabel '{/Symbol s}($parlist[$idim])'\n";
  $name = sprintf("%s_%02d", $out, $idim);
  $name_out = $name;
  print GNU "set output '$name_out.pdf'\n";
  print GNU "pl ";
  foreach $icomp (0 .. $ncomp-1) {
    if ($valid[$idim][$icomp]!=0) {
      print GNU "'$name' u 1:", $icomp+2, " w lp lt $icomp lw 3 t ''";
      print GNU ", " if $icomp<$ncomp-1;
    }
  }
  print GNU "\n";
  print GNU "set output\n";
}
close GNU;

#if (! defined $options{n}) {
`gnuplot $gname`;
`$ENV{COSMOPMC}/bin/allps2tex.pl -f pdf -t "Proposal variances" > all_vars.tex`;
#`$ENV{COSMOPMC}/bin/ldp.sh all_vars -q`;
`pdflatex all_vars -q`;
`rm -f all_vars.log all_vars.dvi all_vars.aux`;
#}


### Subroutines ###

sub get_info {
  my ($conf) = @_;
  open(my $conf_fh, "$conf");
  my $nit   = -1;
  my $ncomp = -1;
  my $ndim  = -1;

  while (<$conf_fh>) {
    if (/^niter/) {
      my @F = split(" ", $_);
      $nit = $F[1];
    }
    if (/^ncomp/) {
      my @F = split(" ", $_);
      $ncomp = $F[1];
    }
    if (/^npar/) {
      my @F = split(" ", $_);
      $ndim = $F[1];
    }
  }
  close $conf_fh;

  return ($nit, $ncomp, $ndim);
}

sub set_ENV_COSMOPMC {
  my ($cwd) = @_;

  # Environment variable defined (in shell)
  return if defined $ENV{COSMOPMC};

  # Copy from command argument
  if (defined $options{P}) {
    $ENV{COSMOPMC} = $options{P};
    return;
  }

  # Use cwd
  if (-e "$cwd/bin/cosmo_pmc.pl") {
    $ENV{COSMOPMC} = $cwd;
    return;
  }

  die "Set environment variable '\$COSMOPMC' or use option '-P PATH'";
}

# Returns if yorick exists and seems to work well
sub check_for_yorick {
  open(my $tmp_fh, ">cytmp.i");
  print {$tmp_fh} "write, 2 * 3\nquit\n";
  close $tmp_fh;
  $res = qx(yorick -batch cytmp.i);
  $ex = $?;
  chomp($res);
  unlink "cytmp.i";

  if ($res ne "" and $res == 6 and $ex == 0) { return 0; }
  return 1;
}

sub usage {
  my ($ex) = @_;

  print STDERR "Usage: proposal_var.pl [OPTIONS]\n";
  print STDERR "OPTIONS:\n";
  print STDERR "  -d DIR         Directory DIR containing the sub-directories 'iter_*'\n";
  print STDERR "                  with the proposal files (default '.')\n";
  print STDERR "  -c CONFIG      Configuration file CONFIG (default 'DIR/config_pmc')\n";
  print STDERR "  -P PATH        Use PATH as CosmoPMC root directory (default: environment\n";
  print STDERR "                  variable \$COSMOPMC)\n";
  print STDERR "  -h             This message\n";

  exit $ex unless $ex < 0;
}
