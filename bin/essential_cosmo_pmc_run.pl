#!/usr/bin/env perl

# essential_cosmo_pmc_run.pl
# Martin Kilbinger 2010-2011

use warnings;

use Cwd;
use Fatal qw/ open close unlink /;
use Getopt::Std;

my %options = ();
getopts("c:P:kvh", \%options);


usage(0) if defined $options{h};
$verbose = defined $options{v} ? 1 : 0;
$config  = defined $options{c} ? $options{c} : "config_pmc";

my $cwd  = cwd;

set_ENV_COSMOPMC($cwd);
die "Environment variable \$COSMOPMC = '$ENV{COSMOPMC}' not valid"
  unless -s "$ENV{COSMOPMC}/bin/cosmo_pmc.pl";


### Main program ###

my $iter_dir = "";
if (-s "deduced") {
   $iter_dir = "deduced";
} else {
   my $niter     = get_niter($config);
   my $last_iter = $niter - 1;
   $iter_dir = "iter_$last_iter";
}

perplexity_split("perplexity");
my %plot = ();
$plot{"key"} = "left";
plot_file("perp.txt", %plot);

$plot{"key"} = "none";
$plot{"ylabel"} = "effective number of components";
plot_file("enc", %plot);


# Mean and confidence intervals
run_own_command("mean2eps.pl -c $config $iter_dir/mean");


$title = qx(basename $cwd);
chomp($title);
$title =~ s/_/\\_/g;
unlink "diagnostics.eps" if -e "diagnostics.eps";
unlink "diagnostics.ps" if -e "diagnostics.ps";
run_own_command("allps2tex.pl -f eps -w 8 -t \"$title\" . > diagnostics.tex");
run_own_command("ldp.sh diagnostics.tex");

run_sys_command("pstopdf diagnostics.ps");
run_sys_command("pstopdf proposal_means/all_means.ps proposal_means/all_means.pdf");

# Paste pdf files together
$cmd = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=essential.pdf diagnostics.pdf proposal_means/all_means.pdf";
if (-s "proposal_vars/all_vars.ps") {
  run_sys_command("pstopdf proposal_vars/all_vars.ps proposal_vars/all_vars.pdf");
  $cmd .= " proposal_vars/all_vars.pdf";
}
if (-s "$iter_dir/all_contour2d.pdf") {
  $cmd .= " $iter_dir/all_contour2d.pdf";
}
run_sys_command($cmd);


# Clean up
unlink <diagnostics.aux diagnostics.log diagnostics.eps diagnostics.dvi>;

chdir "$iter_dir";
unlink <mean.log mean.dvi mean.aux>;
unlink <mean.tab mean_in.tex mean.tex mean.eps> unless defined $options{k};
chdir "$cwd";

unlink <fig_perp.txt.eps fig_enc.eps> unless defined $options{k};
unlink <perp.txt perp.txt.gnu enc.gnu> unless defined $options{k};
unlink "diagnostics.tex" unless defined $options{k};



### End main program ###


### Subroutines ###

sub get_niter {
  my ($config) = @_;
  open(my $conf_fh, "$config");
  my $niter = -1;

  while (<$conf_fh>) {
    if (/^niter/) {
      my @F = split(" ", $_);
      $niter = $F[1];
      last;
    }
  }
  close $conf_fh;

  return $niter;
}

# Write perplexity and ess/nsample to a file, for plotting with 'plot_file'
sub perplexity_split {
  my ($perp_name) = @_;

  open(my $perp_fh, "$perp_name");
  open(my $pout_fh, ">perp.txt");
  print {$pout_fh} "# nsample[10^4] perplexity ess/nsample\n";
  $last = 0;
  while (<$perp_fh>) {
    next if /#/;
    @F = split(" ", $_);
    $ns = $F[1] / 10000;
    $this_ns = $F[1] - $last;
    print {$pout_fh} "$ns $F[2] ", $F[3]/$this_ns, "\n";
    $last = $F[1];
  }
  close $perp_fh;
  close $pout_fh;
}


# Calls gnuplot to plot file with name $name. Creates output
# file 'fig_$name.eps'
# From test_suite_cosmo_pmc.pl
sub plot_file {
  my ($name, %plot) = @_;

  # Get number of columns
  open(my $in_fh, "$name");
  my @header = ();
  my $ncol = -1;
  while (<$in_fh>) {

    my @F = split(" ", $_);
    if (/#/) {
      @header = @F;
      next;
    }
    die "Column numbers inconsistent ($ncol, $#F)" if $ncol!=-1 && $ncol!=$#F;
    $ncol = $#F;

  }
  $ncol ++;
  close $in_fh;

  map(s/\'/\'\'/g, @header) if $#header > -1;

  # Create gnu file
  open(my $gnu_fh, ">$name.gnu");
  print {$gnu_fh} "set output 'fig_$name.eps'\n";
  print {$gnu_fh} "set term post eps enhanced color 'Times-Roman' 32\n";
  print {$gnu_fh} "set logs x\n" if has_value($plot{"logx"}, 1);
  print {$gnu_fh} "set logs y\n" if has_value($plot{"logy"}, 1);
  if (defined $plot{"key"}) {
    $key = $plot{"key"};
    print {$gnu_fh} $key eq "none" ? "unset key\n" : "set key $key\n";
  }

  if (defined $plot{"xlabel"}) {
    my $xlabel = $plot{"xlabel"};
    print {$gnu_fh} "set xlabel '$xlabel'\n";
  } else {
    print {$gnu_fh} "set xlabel '$header[1]'\n" if $#header > -1;
  }

  if (defined $plot{"ylabel"}) {
    my $ylabel = $plot{"ylabel"};
    print {$gnu_fh} "set ylabel '$ylabel'\n";
  }

  print {$gnu_fh} "pl '$name' ";
  foreach $n (2 .. $ncol) {
    print {$gnu_fh} "'' " unless $n==2;
    my $title = $#header > -1 ? "$header[$n]" : "";
    print {$gnu_fh} "u 1:(\$$n) w lp t '$title'";
    print {$gnu_fh} ", " unless $n==$ncol;
  }
  print {$gnu_fh} "\n";
  print {$gnu_fh} "set output\n";
  close $gnu_fh;

  run_sys_command("gnuplot $name.gnu");
}

# Returns 1 if $var is defined and $var = $val
sub has_value {
  my ($var, $val) = @_;
  return 0 if ! defined $var;
  return $var == $val;
}

# New: Shift exit error code by 8 bits 
sub run_command {
  my ($command) = @_;

  my $cwd = cwd;
  print STDERR "Running '$command' in '$cwd'\n" if $verbose;
  `$command`;
  my $ex = $? >> 8;
  die "Last command '$command' unsuccessful (exit code $ex)\n" if $ex;
}

# Runs system program
sub run_sys_command {
  my ($command) = @_;
  run_command("$command");
}

# Runs CosmoPMC program
sub run_own_command {
  my ($command) = @_;
  run_command("$ENV{COSMOPMC}/bin/$command");
}

sub usage {
  my ($ex) = @_;

  print STDERR "Usage: essential_cosmo_pmc.pl [OPTIONS]\n";
  print STDERR "OPTIONS:\n";
  print STDERR "    -c CONFIG      Uses config file CONFIG (default: 'config_pmc')\n";
  print STDERR "    -P PATH        Use PATH as CosmoPMC directory (default: environment\n";
  print STDERR "                    variable \$COSMOPMC)\n";
  print STDERR "    -k             Keep temporary files\n";
  print STDERR "    -v             Verbose\n";
  print STDERR "    -h             This message\n";

  exit $ex if defined $ex;
}

sub set_ENV_COSMOPMC {
  my ($cwd) = @_;

  # Copy from command argument
  if (defined $options{P}) {
    $ENV{COSMOPMC} = $options{P};
    return;
  }

  # Environment variable defined (in shell)
  return if defined $ENV{COSMOPMC};

  # Use cwd
  if (-e "$cwd/bin/cosmo_pmc.pl") {
    $ENV{COSMOPMC} = $cwd;
    return;
  }

  die "Set environment variable '\$COSMOPMC' or use option '-P PATH'";
}
