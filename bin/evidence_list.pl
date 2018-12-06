#!/usr/bin/perl -w

# Martin Kilbinger 2009

use Fatal qw/ open /;
use Getopt::Std;

%options=();
getopts("r:k:s:S:nLlh", \%options);

if (defined $options{h} || $#ARGV==-1) {
    usage();
    exit -1;
}

$fname_evi = defined $options{L} ? "evidence_fisher" : "evidence";

@DIRS = @ARGV;
$dmin = -1;
$logEmin = -1.0e30;
$d = 0;
while ($dir = shift @DIRS) {

    open(EVI, "$dir/$fname_evi");
    while (<EVI>) {
	next if /#/;
	@F = split(" ", $_);
    }
    close EVI;

    $logE[$d] = $F[1];
    $lnE[$d]  = $F[2];
    $E[$d]    = $F[3];
    $dmin     = $d if $logEmin<$logE[$d];

    # Get number of parameters from config file
    if (defined $options{n}) {
	open(CONF, "$dir/config_pmc");
	while (<CONF>) {
	    if (/npar/) {
		@F = split(" ", $_);
		$npar[$d] = $F[1];
		last;
	    }
	}
	close CONF;
    }

    $d++;
}

$D = $d;

$Sep = " ";
$Sep = $options{S} if defined $options{S};

# Print header
printf "# %20s$Sep ", "Model";
printf "%7s$Sep ", "\$n_{\\textrm{par}}\$" if defined $options{n};
printf "%11s$Sep %11s$Sep %11s$Sep\n", "\$\\log{E}\$", "\$\\ln{E}\$", "\$E\$";

if (defined $options{k}) {
  $sep = " ";
  $sep = $options{s} if defined $options{s};
  @keys = split($sep, $options{k});
} else {
  @keys = @ARGV;
  foreach $d (0 .. $D-1) {
    $keys[$d] =~ s/_/\\_/;
  }
}

# Maximum key length
$len = 0;
foreach $d (0 .. $D-1) {
    $l = length $keys[$d];
    $len = $l if $len<$l;
}


if (defined $options{r}) {
    if ($options{r}==-1) {
	$indmin = $dmin;
    } else {
	$indmin = $options{r};
    }
    $logmin = $logE[$indmin];
    $lnmin  = $lnE[$indmin];
} else {
    $logmin = $lnmin = 0.0;
}


foreach $d (0 .. $D-1) {
    $format = sprintf("%%%ds$Sep", $len+2);
    printf $format, $keys[$d];
    printf "%7d$Sep", $npar[$d] if defined $options{n};
    $ratio = exp($lnmin)==0 ? -1 : $E[$d]/exp($lnmin);
    printf "% 11.3f$Sep% 11.3f$Sep% 10.3g$Sep\n", $logE[$d]-$logmin, $lnE[$d]-$lnmin, $ratio;
}


sub usage {
    print STDERR "Usage: evidence_list.pl [OPTIONS] DIR1 [DIR2 [...]]\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -r N            Subtract log(E) from DIRN (default: no subtraction)\n";
    print STDERR "                  For N=-1 subtract log(E_min)\n";
    print STDERR "  -k KEY          Use KEY (string list) instead of\n";
    print STDERR "                  directory names (default)\n";
    print STDERR "  -s SEP          Use SEP as input separator for KEY list\n";
    print STDERR "  -S SEP          Use SEP as output separator\n";
    print STDERR "                  (default for both: white-space)\n";
    print STDERR "  -n              Write number of model parameters\n";
    print STDERR "  -L              Use Laplace approximation (reading file 'evidence_fisher')\n";
    print STDERR "  -h              This message\n";
}
