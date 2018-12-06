#!/usr/bin/perl -w

# cosmo_mcmc.pl
# Jean coupon. Adapted from cosmo_pmc.pl (Martin Kilbinger 2010 - 2011)

use Fatal qw/ open /;
use Getopt::Std;
use Cwd;

### Command line options
%options=();
getopts("p:c:o:f:rm:M:dDaA:n:S:s:P:eO:qh", \%options);

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
$plot_mode   = defined $options{p} ? "$options{p}" : "y";
$opt_plot    = defined $options{O} ? $options{O} : "";
$mult        = defined $options{M} ? $options{M} : 2;
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

$config_mcmc  = defined $options{c} ? $options{c} : "config_mcmc";
$output_dir   = defined $options{o} ? $options{o} : ".";

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

die "MCMC config file $config_mcmc not found (use '-c CONFIG' open)" if ! -e $config_mcmc;

$ex_chain = existence("chain.fin", "Delete 'chain.fin' and rerun MCMC", "chain.fin");
$ex_chain = existence("chain.acc", "Delete 'chain.acc' and rerun MCMC", "chain.acc");
if (!$ex_chain) {

  open(CONFIG, "$config_mcmc");
  while ($_ = <CONFIG>) {
    if (/^sinitial/) {
      @F = split(" ", $_);
      $sinitial = $F[1];
      last;
    }
  }
  close CONFIG;
  $want_fisher = ($sinitial =~ /Fisher/ ? 1 : 0);

  if ($want_fisher) {
    $ex_fisher = existence("$fisher", "Delete and recalculate Fisher matrix '$fisher'", "$fisher");

    if (!$ex_fisher) {
      $ex_conf_fish = existence("$config_fish", "Delete and recreate '$config_fish' from '$config_mcmc' and '$maxlogP'",
				               "$config_fish");

      if (!$ex_conf_fish) {
	    $ex_maxlogP = existence("$maxlogP", "Delete '$maxlogP' and recalculate maximum posterior", "$maxlogP");

	    if (!$ex_maxlogP) {
	      $ex_conf_max = existence("$config_max", "Delete and recreate '$config_max' from '$config_mcmc'", "$config_max");

	      if (!$ex_conf_max) {
	        runpr("$ENV{COSMOPMC}/bin/config_pmc_to_max_and_fish.pl -M $random_flag $fid_flag -c $config_mcmc > $config_max");
	      }

	      runpr("$ENV{COSMOPMC}/exec/max_post -t -m $max_method -s $seed");

	      goto end if $options{S} eq "M";
	    }

	    runpr("$ENV{COSMOPMC}/bin/config_pmc_to_max_and_fish.pl -F -p $maxlogP -c $config_mcmc $diag_flag > $config_fish");
      }

	  goto end if $options{S} eq "M";

      runpr("$mpi_cmd$ENV{COSMOPMC}/exec/go_fishing $adaptive $f_pos_flag", 1);
      if (! -e "fisher") {
	    print STDERR "Calculation of Fisher matrix failed. Stopping cosmo_mcmc.pl\n";
	    exit 1;
      } else {
	    print STDERR "(If an error occured just now, it might have been caused by 'mpirun' and not 'go_fishing')\n";
      }

      goto end if $options{S} eq "F";
    }
  }
}


##### REPLACE #######

# Save previous files
# rename "perplexity", "perplexity.prev" if -s "perplexity";

######### WITH ########
# Save old log file before overwriting
#if (-e log_mcmc) then
#    set name = `find . -maxdepth 1 -name log_mcmc -printf "%CY-%Cm-%Cd-%CT\n"`
#    mv log_mcmc log_mcmc-$name
#endif


# Run MCMC
if (!($plot_mode =~ "o")) {
  # do not run cosmo_mcmc if only plot
  runpr("$mpi_cmd$ENV{COSMOPMC}/exec/cosmo_mcmc -c $config_mcmc -o $output_dir -s $seed");
}

### Plotting ###

# Contour plots
chdir "$output_dir";
print "Creating plots in directory '$output_dir'\n";
if ($plot_mode =~ "y") {
  # yorick
  runpr("$ENV{COSMOPMC}/bin/plot_contour2d.pl -c ../$config_mcmc -p 0 -1 m1t -o pdf -T \"$title\" $quiet $opt_plot");
}
if ($plot_mode =~ "R") {
  # R THIS WILL PROBABLY NOT WORK
  runpr("Rscript $ENV{COSMOPMC}/R/sample_from_pmcsimu.R pmcsim -i chain.fin -M $mult") unless -e "sample";
  runpr("Rscript /Users/coupon/source/ecosstat/trunk/bin/plot_confidence.R sample -c ../config_mcmc")
}
chdir "..";

end:
print "cosmo_mcmc.pl finished\n";


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
      print STDERR "Last command returned error code $?. Stopping cosmo_mcmc.pl\n";
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
  print STDERR "   -o DIR               Output directory (default: current)\n";
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
  print STDERR "   -p PRO               Plotting scripts: 'y' (yorick; default), 'R' (R), 'n' (none)\n";
  print STDERR "                         or 'o' (only)\n";
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

