#!/usr/bin/perl -w

# MK 6/2008
# Cuts a MCM chain (chain.fin) in pieces i=1..niter of length i*nsample. Writes the chains in
# directories iter_<iter>/chain.fin.


use Getopt::Std;
%options=();
getopts("n:h", \%options);
usage() if defined $options{h};

$nsamples = defined $options{n} ? $options{n} : 10000;

open(IN, "chain.fin") or die "Could not open file chain.fin: $!";

$laufend = 0;
$iter    = 0;
$header = "";
while (<IN>) {
    if (/#/) {
	$header = $header . $_;
	next;
    }

    if ($laufend==0) {
	print STDERR "Creating dir iter_$iter\n";
	mkdir "iter_$iter";
	open(OUT, ">iter_$iter/chain.tmp") or die "Could not create file iter_$iter/chain.tmp: $!";
    }

    print OUT $_;
    $laufend++;

    if ($laufend==$nsamples) {
	$nsamples[$iter] = $laufend;
	close OUT;
	$laufend = 0;
	$iter++;
    }
}

# Last chain
$nsamples[$iter] = $laufend;
close OUT;

$fsfinal = $nsamples[$iter]/$nsamples[0];
$itermax = $iter;

# Delete old chains if any
`rm -f iter_*/chain.fin`;

open(HEAD, ">header") or die "Could not create file header: $!";
print HEAD $header;
close HEAD;

foreach $iter (0 .. $itermax) {
    $tmp    = "iter_$iter/chain.tmp";
    if ($iter>0) {
	$iterm1 = $iter - 1;
	$pre    = "iter_$iterm1/chain.noh";
    } else {
	$pre = "";
    }
    $final  = "iter_$iter/chain.noh";
    runpr("cat $pre $tmp > $final");
}
foreach $iter (0 .. $itermax) {
    runpr("cat header iter_$iter/chain.noh > iter_$iter/chain.fin");
}

`rm -f iter_*/chain.tmp iter_*/chain.noh header`;

# Return value for mkmc_mean_var.sh
print "$nsamples[0] $fsfinal $itermax\n";
exit(0);


##############

sub usage {
    print STDERR "Usage cut_mkmc.pl [-n nsamples]\n"
}

sub runpr {
    my ($command) = @_;
    print STDERR "*** $command ***\n";
    $res = `$command`;
    return $res;
}
