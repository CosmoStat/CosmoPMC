#!/usr/bin/perl -w

# meanvar2tabl.pl
# Martin Kilbinger 2008-2009
#
# Reads a 'mean' file and prints latex-readable output.
# To transform into a latex table, use 'tab2tex.pl' as
# second step.

use Getopt::Std;

%options=();
getopts("s:p:eht:c:S:P:", \%options);

usage() if $#ARGV==-1;

my @conf = ("-1", "68.3", "95.5", "99.7");

usage() if defined $options{h};
$sigma     = defined $options{s} ? $options{s} : 1;
$precision = defined $options{p} ? $options{p} : 3;
$t         = defined $options{t} ? $options{t} : "Mean\$\\pm\$$conf[$sigma]" . "\\%" . "cl." x @ARGV;
$sep       = defined $options{S} ? $options{S} : " ";
@title     = split($sep, $t);


@sig    = (1,0,0);
$sig[0] = 0 if !($sigma =~ /1/);
$sig[1] = 1 if $sigma =~ /2/;
$sig[2] = 1 if $sigma =~ /3/;


usage() if !($sigma =~ /[123]/);


# Read config file for parameter type
$configname = "config_pmc" if -e "config_pmc";
$configname = $options{c} if defined $options{c};
die "No configuration file found" unless defined $configname;

$output = qx($ENV{COSMOPMC}/bin/get_spar.pl -c $configname TeX);
@parlist = split("&", $output);

# Header
print "# Parameter$sep";
foreach $j (0 .. $#ARGV) {
    printf " %39s", "$title[$j]$sep";
}
print "\n";

# Read mean file(s)
foreach $j (0 .. $#ARGV) {
    open(MEAN, "$ARGV[$j]") or die "Could not open file '$ARGV[$i]': $!";
    $l = 0;
    while (<MEAN>) {
	next if /#/;
	$line[$j][$l] = $_;
	$l++;
    }
    
    #print "File $j, nl $l\n";

    if ($j==0) {
	$npar = $l-1;
    } else {
	die "Number of parameters (", $l-1, ") inconsistent (with $npar) in file '$ARGV[$j]'\n"
	    if $npar != $l-1;
    }
}

# Dummy parameter names
if (! defined @parlist || $#parlist<=0) {
    foreach $p (0 .. $npar) {
	$parlist[$p] = "p_\{$p\}";
    }
}

# Write table
foreach $p (0 .. $npar) {

    foreach $j (0 .. $#ARGV) {

	#print STDERR "*** $p $j $line[$j][$p] ***\n";
	if ($line[$j][$p] =~ /#/) {
	    print STDERR "Skipping parameter $p in file $j\n";
	    print " -";
	    next;
	}

	@F = split(" ", $line[$j][$p]);

	# For backward compability
	($mean, $sigp, $sigm) = ($F[1], $F[2*$sigma], $F[2*$sigma+1]) if $#F==7;
	($mean, $sigp, $sigm) = ($F[2], $F[2*$sigma+1], $F[2*$sigma+2]) if $#F==8;

	if (!($mean eq "-" && $sigp eq "-" && $sigm eq "-")) {
	    if (defined $options{e}) {
		if (abs($sigp)>1) { 
		    $i = 1; $s = 0;
		    while ($i*abs($sigp)>1) { $i /= 10; $s++; }
		    $myprec = $s - 2 + $precision;
		} else {
		    $i = 1; $s = 0;
		    while ($i*abs($sigp)<1) { $i *= 10; $s++; }
		    $myprec = $s - 1 + $precision;
		}
	    } else {
		    $myprec = $precision;
	    }

	    $fs = "%." . $myprec . "f";

	    $pm_equal = int($sigp*10**$myprec+0.5)==int($sigm*10**$myprec+0.5);

	    if ($pm_equal!=1) {
		    $format1 = " \t\$$fs^{+$fs}_{-$fs}\$ $sep";
	    } else {
		    $format1 = " \t\$$fs\\pm{$fs}\$ $sep"; 
	    }

	    if ($j==0) {
	            my $str = sprintf("\$%s\$ $sep", $parlist[$p]);
		    printf "%-30s", $str;
		    #printf "\$%s\$ $sep", $parlist[$p];
	    }

	    if ($pm_equal!=1) {
		    printf $format1, $mean, $sigp, $sigm;
	    } else {
		    printf $format1, $mean, $sigp;
	    }
	} else {
	    printf " -";
	}
    }

    print "\n";
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

sub usage {
    print STDERR "Usage: meanvar2tab.pl [OPTIONS] file [file2 [...]]\n";
    print STDERR "\nOptions:\n";
    print STDERR "  -s {123}     68% (1), 95% (2) or 99.7% (3) errors (default = 1)\n";
    print STDERR "  -p PREC      Output with PREC digits (\'\%PREC\' format string)\n";
    print STDERR "  -e           Error(s) written to PREC significant digits (use -p PREC)\n";
    print STDERR "  -c CONFIG    Uses config file CONFIG (default: 'config_pmc')\n";
    print STDERR "  -t TITLE     Title (table heading) TITLE is string list with entries according\n";
    print STDERR "               to the number of input files\n";
    print STDERR "  -S SEP       Use SEP as input separator for TITLE list (default: white space)\n";
    print STDERR "  -P PATH      Use PATH as CosmoPMC directory (default: environment\n";
    print STDERR "                variable \$COSMOPMC)\n";
    print STDERR "  -h           This message\n";
    exit(2);
}
