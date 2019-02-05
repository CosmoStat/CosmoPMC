#!/usr/bin/perl -w

# cosmo_pmc.pl
# Martin Kilbinger 2010 - 2011

use Fatal qw/ open /;
use Getopt::Std;
use Cwd;

### Command line options
%options=();
getopts("p:c:f:rm:dDaA:n:S:s:P:eO:qh", \%options);

usage() if defined $options{h};

$def_r = defined $options{r};
$def_f = defined $options{f};
$sum_rf = $def_r + $def_f;
die "Only one option out of '-r' and '-f' can be defined" if $sum_rf>1;

$ncpu        = defined $options{n} ? $options{n} : 1;
$random_flag = defined $options{r} ? "-r" : "";
$max_method  = defined $options{m} ? $options{m} : "a";
$diag_flag   = defined $options{d} ? "-d" : "";
$f_pos_flag  = defined $options{D} ? "" : "-f";
$adaptive    = defined $options{a} ? "-a" : "";
$seed        = defined $options{s} ? $options{s} : -1;
$fid_flag    = defined $options{f} ? "-f \"$options{f}\"" : "";
$plot_mode   = defined $options{p} ? "$options{p}" : "R";
$opt_plot    = defined $options{O} ? $options{O} : "";
$quiet       = defined $options{q} ? "-q" : "";

$mpi_cmd     = $ncpu == 1 ? "" : "mpirun -n $ncpu ";

$cwd = cwd;
my $title = qx(basename $cwd);
chomp($title);

set_ENV_COSMOPMC($cwd);

if (!($max_method eq "c") && !($max_method eq "a")) {
  print STDERR "Wrong maximum-search method (option '-m')\n";
  usage();
  exit 2;
}

# Default file names
$fisher      = "fisher";
$config_fish = "config_fish";
$maxlogP     = "maxlogP";
$config_max  = "config_max";

$config_pmc  = defined $options{c} ? $options{c} : "config_pmc";


if (defined $options{A}) {
  $def_answer = $options{A};
  die "Value of option '-A' has to be 'y' or 'n'" if !($def_answer eq 'y' || $def_answer eq 'n');
} else {
  $def_answer = "";
}

if (defined $options{S}) {
  die "Value of flag '-S' has to be 'M' or 'F'" if !($options{S} eq "M") && !($options{S} eq "F");
} else {
  $options{S} = "";
}

die "PMC config file $config_pmc not found (use '-c CONFIG' open)" if ! -e $config_pmc;

$ex_pmcsim = existence("iter_0/pmcsim", "Delete 'iter_*/pmcsim' and rerun PMC", "iter_*/pmcsim");

if (!$ex_pmcsim) {

  open(CONFIG, "$config_pmc");
  while ($_ = <CONFIG>) {
    if (/^sinitial/) {
      @F = split(" ", $_);
      $sinitial = $F[1];
      last;
    }
  }
  close CONFIG;
  $want_fisher = ($sinitial =~ /fisher/ ? 1 : 0);

  if ($want_fisher) {
    $ex_fisher = existence("$fisher", "Delete and recalculate Fisher matrix '$fisher'", "$fisher");

    if (!$ex_fisher) {
      $ex_conf_fish = existence("$config_fish", "Delete and recreate '$config_fish' from '$config_pmc' and '$maxlogP'",
				               "$config_fish");

      if (!$ex_conf_fish) {
	    $ex_maxlogP = existence("$maxlogP", "Delete '$maxlogP' and recalculate maximum posterior", "$maxlogP");

	    if (!$ex_maxlogP) {
	      $ex_conf_max = existence("$config_max", "Delete and recreate '$config_max' from '$config_pmc'", "$config_max");

	      if (!$ex_conf_max) {
	        runpr("$ENV{COSMOPMC}/bin/config_pmc_to_max_and_fish.pl -M $random_flag $fid_flag -c $config_pmc > $config_max");
	      }

	      runpr("$ENV{COSMOPMC}/exec/max_post -t -m $max_method -s $seed");

	      goto end if $options{S} eq "M";
	    }

	    runpr("$ENV{COSMOPMC}/bin/config_pmc_to_max_and_fish.pl -F -p $maxlogP -c $config_pmc $diag_flag > $config_fish");
      }

	  goto end if $options{S} eq "M";

      runpr("$mpi_cmd$ENV{COSMOPMC}/exec/go_fishing $adaptive $f_pos_flag", 1);
      if (! -e "fisher") {
	    print STDERR "Calculation of Fisher matrix failed. Stopping cosmo_pmc.pl\n";
	    exit 1;
      } else {
	    print STDERR "(If an error occured just now, it might have been caused by 'mpirun' and not 'go_fishing')\n";
      }

      goto end if $options{S} eq "F";
    }
  }

}


# Save previous files
rename "perplexity", "perplexity.prev" if -s "perplexity";

# Run PMC
if (!($plot_mode =~ "o")) {
  # do not run cosmo_pmc if only plot
  runpr("$mpi_cmd$ENV{COSMOPMC}/exec/cosmo_pmc -c $config_pmc $quiet -s $seed");
}

### Plotting ###

# Read final iteration from config file
open(CONFIG, "$config_pmc");
while ($_ = <CONFIG>) {
  if (/^niter/) {
    @F = split(" ", $_);
    $last_iter = $F[1] - 1;
    last;
  }
}
close CONFIG;
$dir = "iter_" . $last_iter;

# Contour plots
chdir "$dir";
print "Creating plots in directory '$dir'\n";
if ($plot_mode =~ "y") {
  # yorick
  runpr("$ENV{COSMOPMC}/bin/plot_contour2d.pl -c ../$config_pmc -1 m1t -o pdf -T \"$title\" $quiet $opt_plot");
}
if ($plot_mode =~ "R") {
  # R
  runpr("Rscript $ENV{COSMOPMC}/bin/plot_confidence.R pmcsim -c ../$config_pmc -W $opt_plot -P $ENV{COSMOPMC}");
}
chdir "..";

# Weight plots
#runpr("$ENV{COSMOPMC}/plot_pmc_final.sh -c $config_pmc");

# Proposal means
mkdir "proposal_means";
chdir "proposal_means";
runpr("$ENV{COSMOPMC}/bin/proposal_mean.pl -d .. -c ../$config_pmc");
chdir "..";

# Proposal variances
mkdir "proposal_vars";
chdir "proposal_vars";
runpr("$ENV{COSMOPMC}/bin/proposal_var.pl -d .. -c ../$config_pmc", 1);
chdir "..";

runpr("$ENV{COSMOPMC}/bin/essential_cosmo_pmc_run.pl -c $config_pmc") if defined $options{e};

end:
print "cosmo_pmc.pl finished\n";


### Subroutines ###

# Returns 0 if '$name' does not exist (after possibly removing file).
sub existence {
  my ($name, $msg, $to_remove) = @_;

  my $exists = 1;

  if (-z $name) {
    print "Warning: '$name' is an empty file.\n";
  }

  if (-e $name) {
    print "$msg? [n] ";
    if ($def_answer eq '') {
      $my_answer = <STDIN>;
      chomp($my_answer);
    } else {
      print "$def_answer";
      $my_answer = $def_answer;
    }

    if ($my_answer eq "y") {
      runpr("rm -fr $to_remove");
      $exists = 0;
    }
  } else {
    $exists = 0;
  }

  return $exists;
}

sub runpr {
  my ($command, $no_exit) = @_;
  print STDERR "*** Running $command ***\n";
  system($command);

  $no_exit = 0 if ! defined $no_exit;
  if ($?) {
    if (! $no_exit) {
      print STDERR "Last command returned error code $?. Stopping cosmo_pmc.pl\n";
      exit $?;
    } else {
      print STDERR "Last command returned error code $?, continuing...\n";
    }
  }
}


sub usage {

  print STDERR "Usage: cosmo_pmc.pl [OPTIONS]";
  print STDERR "\nOPTIONS:\n";
  print STDERR "   -n NCPU              Run PMC in parallel on NPCU cpus using 'mpirun' (default: 1)\n";
  print STDERR "   -c CONFIG            Configuration file for PMC (default: config_pmc)\n";
  print STDERR "   -f FID               Fiducial starting point FID. FID is a white-space\n";
  print STDERR "                         separated list in quotes, e.g. '0.25 0.75'\n";
  print STDERR "   -r                   Random starting point for maximum search\n";
  print STDERR "                         (default: (max-min)/2)\n";
  print STDERR "   -m [c|a]             Maximum-search method: 'c' (cg), 'a' (amoeba)\n";
  print STDERR "   -d                   Calculate only diagonal of Fisher matrix\n";
  print STDERR "   -D                   Do not force Fisher matrix F to be positiv. If F is negative,\n";
  print STDERR "                         script exits with an error\n";
  print STDERR "   -a                   Adaptive numerical differentiation for Fisher matrix\n";
  print STDERR "   -s SEED              Use SEED for random number generator. If SEED=-1 (default)\n";
  print STDERR "                         the current time is used as seed.\n";
  print STDERR "   -S [M|F]             Stops after maximum search ('M') or Fisher matrix ('F')\n";
  print STDERR "   -A [y|n]             Default answer to all questions on stdin\n";
  print STDERR "   -P PATH              Use PATH as CosmoPMC directory (default: environment\n";
  print STDERR "                         variable \$COSMOPMC)\n";
  print STDERR "   -e                   Create 'essential' plots\n";
  print STDERR "   -p PRO               Plotting of marginalized posterior (1d and 2d):\n";
  print STDERR "                         PRO = 'R' (R; default), 'y' (yorick+perl), 'n' (none),\n";
  print STDERR "                         'o' (only). Letters can be combined, e.g. \'yRo\'.\n";
  print STDERR "                         Combinations of letters are possible, e.g. 'yR' or 'oy'\n";
  print STDERR "   -M MULT              Output sample MULT times input (default 1).\n";
  print STDERR "                         Valid if plotting script is 'R'\n";
  print STDERR "   -O OPT               Pass options OPT to 'plot_contour2d.pl'\n";
  print STDERR "   -q                   Quiet mode\n";
  print STDERR "   -h                   This message\n";
  exit 1;
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

