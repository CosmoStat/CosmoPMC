#!/usr/bin/env perl

# add_ded.pl
# Martin Kilbinger 2010
# Adds deduced parameters to a PMC simulation file
#  2012:    Read parameter indices from config file
#  08/2013: w(a) = w_0 + w_1 (1 - a) as function of a,
#           to determine pivot w_p


use warnings;
use Fatal qw/open/;
use Getopt::Std;


my %options = ();
getopts("c:n:h", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV < 1;

$sample  = $ARGV[0];
$newspar = $ARGV[1];
$nnew    = 1;


if ($newspar eq "Sigma") {
    die "Usage: add_ded.pl sample \"Sigma\" alpha Omega_m_fid\n" if $#ARGV!=3;
} elsif ($newspar eq "wa") {
    die "Usage: add_ded.pl sample \"wa\" [Na]" if $#ARGV!=1 and $#ARGV!=2;
} else {
    die "Usage: add_ded.pl sample \"$newspar\"\n" if $#ARGV!=1;
}


# Get parameter list from config file
my $conf = defined $options{c} ? $options{c} : "config_pmc";
open(my $c_fh, "$conf");
while (<$c_fh>) {

  next if /^#/;
  @F = split(" ", $_);
  next if $#F==-1;

  if ($F[0] eq "npar") {
    $npar = $F[1];
  }
  if ($F[0] eq "n_ded") {
    $n_ded = $F[1];
  }
  if ($F[0] eq "spar") {
    die "Entry 'npar' not found in configuration file" if ! defined $npar;
    die "Entry 'n_ded' not found in configuration file" if ! defined $n_ded;
    foreach $i (0 .. $npar+$n_ded-1) {
      $spar{$F[$i+1]} = $i;
    }
  }

}
close $c_fh;



if ($newspar eq "Sigma") {

  die "Omega_m not found in parameter list" unless defined $spar{Omega_m};
  die "sigma_8 not found in parameter list" unless defined $spar{sigma_8};

  $alpha  = $ARGV[2];
  $Om0    = $ARGV[3];
  $iOm    = $spar{Omega_m};
  $is8    = $spar{sigma_8};

  $newmin = 0.5;
  $newmax = 1.0;

} elsif ($newspar eq "q_acc") {

  die "Omega_m not found in parameter list" unless defined $spar{Omega_m};
  die "Omega_de not found in parameter list" unless defined $spar{Omega_de};

  $iOm = $spar{Omega_m};
  $iOL = $spar{Omega_de};

  $newmin = -1.0;
  $newmax =  1.0;

} elsif ($newspar eq "Omega_K") {

  if (defined $spar{omega_K} && defined $spar{h_100}) {
    $iok = $spar{omega_K};
    $ih  = $spar{h_100};
  } elsif (defined $spar{Omega_m} && defined $spar{Omega_de}) {
    $iOm = $spar{Omega_m};
    $iOL = $spar{Omega_de};
  } else {
    die "Either (omega_K, h_100) or (Omega_m, Omega_de) have to be in parameter list";
  }

  $newmin = -1.0;
  $newmax =  1.0;

} elsif ($newspar eq "Omega_de") {

  die "Omega_m not found in parameter list" unless defined $spar{Omega_m};
  die "Omega_K not found in parameter list" unless defined $spar{Omega_K};

  $iOm  = $spar{Omega_m};
  $iOK  = $spar{Omega_K};

  $newmin = 0.0;
  $newmax = 2.0;

} elsif ($newspar eq "Omega_m") {

  die "omega_m not found in parameter list" unless defined $spar{omega_m};
  die "h_100 not found in parameter list" unless defined $spar{h_100};

  $iom = $spar{omega_m};
  $ih  = $spar{h_100};

  $newmin = 0.0;
  $newmax = 2.0;
    
} elsif ($newspar eq "Omega_b") {
    
  die "omega_b not found in parameter list" unless defined $spar{omega_b};
  die "h_100 not found in parameter list" unless defined $spar{h_100};
    
  $iob = $spar{omega_b};
  $ih  = $spar{h_100};
    
  $newmin = 0.0;
  $newmax = 0.2;

} elsif ($newspar eq "Omega_c") {

    die "Omega_m not found in parameter list" unless defined $spar{Omega_m};
    die "Omega_b not found in parameter list" unless defined $spar{Omega_b};

    $iOm = $spar{Omega_m};
    $iOb = $spar{Omega_b};
    $ih  = $spar{h_100};

    $newmin = 0.0;
    $newmax = 2.0;

} elsif ($newspar eq "omega_b") {

    $iOb = $spar{Omega_b};
    $ih  = $spar{h_100};

    $newmin = 0.0;
    $newmax = 0.1;

} elsif ($newspar eq "omega_c") {

    $iOc = $spar{Omega_c};
    $ih  = $spar{h_100};

    $newmin = 0.0;
    $newmax = 1.0;

} elsif ($newspar eq "wa") {

    die "w_0_de not found in parameter list" unless defined $spar{w_0_de};
    die "w_1_de not found in parameter list" unless defined $spar{w_1_de};

    $iw0 = $spar{w_0_de};
    $iw1 = $spar{w_1_de};

    $nnew = defined $ARGV[2] ? $ARGV[2] : 5;
    my $zmax = 2;

    $newmin = -10.0;
    $newmax =  10.0;

    @Z_LIST = ();
    @A_LIST = ();
    my $dz = $zmax/($nnew-1);
    for (my $i=0; $i<$nnew; $i++) { 
        $Z_LIST[$i] = $zmax - $dz * $i;
        $A_LIST[$i] = 1.0 / ($Z_LIST[$i] + 1.0);
    }
    open(my $out_fh, ">wa_az_list.txt");
    print {$out_fh} "z = ";
    foreach $i (0 .. $#Z_LIST) { printf {$out_fh} " %.3f", $Z_LIST[$i]; }
    print {$out_fh} "\na = ";
    foreach $i (0 .. $#A_LIST) { printf {$out_fh} " %.3f", $A_LIST[$i]; }
    print {$out_fh} "\n";
    close $out_fh;

} else {

   die "Deduced parameter \'$newspar\' not supported";

}

open(CONFIG, "$conf");
open(OUT, ">$conf+ded");
while (<CONFIG>) {

    # Add $nnew deduced parameters to config file
    s/(spar.*)/$1. ("\t$newspar" x $nnew)/e;
    s/(min.*)/$1. ("\t$newmin" x $nnew)/e;
    s/(max.*)/$1. ("\t$newmax" x $nnew)/e;
    s/n_ded.*(\d+)/my $x = $1+$nnew; "n_ded\t\t" . $x/e;
    print OUT;

}
close CONFIG;
close OUT;

open(SAMPLE, "$sample");
while (<SAMPLE>) {

    chomp $_;

    if (/#/) {

        s/n_ded.*=.*(\d+)/$x = $1; $y=$1+$nnew; "n_ded = " . $y;/e;
        print;

        if (/weight\s*chi2/) {
            foreach my $i (0 .. $nnew-1) {
                print "    param_ded[", $x+$i, "]";
            }
            print "\n";
            $_ = <SAMPLE>;
            chomp $_;

            print $_;
            foreach my $i (0 .. $nnew-1) {
                printf "%16s", $newspar;
            }
        }
        print "\n";
        next;

    }

    my @F = split(" ", $_);
    my @X = ();

    if ($newspar eq "Sigma") {
        $Om = $F[$iOm+2]; 
        $s8 = $F[$is8+2];
        $X[0]  = $s8*($Om/$Om0)**$alpha;
    } elsif ($newspar eq "q_acc") {
        $Om = $F[$iOm+2];
        $OL = $F[$iOL+2];
        $X[0]  = $Om/2.0 - $OL;
    } elsif ($newspar eq "Omega_K") {
        if (defined $spar{omega_K} && defined $spar{h_100}) {
            $ok = $F[$iok+2];
            $h  = $F[$ih+2];
            $X[0]  = $ok/$h/$h;
        } elsif (defined $spar{Omega_m} && defined $spar{Omega_de}) {
            $Om = $F[$iOm+2];
            $OL = $F[$iOL+2];
            $X[0]  = 1.0 - $Om - $OL;
        }
    } elsif ($newspar eq "Omega_de") {
        $Om = $F[$iOm+2];
        $OK = $F[$iOK+2];
        $X[0]  = 1.0 - $Om - $OK;
    } elsif ($newspar eq "Omega_m") {
        $om = $F[$iom+2];
        $h  = $F[$ih+2];
        $X[0]  = $om/$h/$h;
    } elsif ($newspar eq "Omega_b") {
        $ob = $F[$iob+2];
        $h  = $F[$ih+2];
        $X[0]  = $ob/$h/$h;
    } elsif ($newspar eq "Omega_c") {
        $Om = $F[$iOm+2];
        $Ob = $F[$iOb+2];
        $X[0] = $Om - $Ob;
    } elsif ($newspar eq "omega_b") {
        $Ob = $F[$iOb+2];
        $h  = $F[$ih+2];
        $X[0] = $Ob*$h*$h;
    } elsif ($newspar eq "omega_c") {
        $Oc = $F[$iOc+2];
        $h  = $F[$ih+2];
        $X[0] = $Oc*$h*$h;
    } elsif ($newspar eq "wa") {
        foreach my $i (0 .. $#A_LIST) {
            # MKDEBUG: Valid only for sde_param = linder
            my $wa = $F[$iw0+2] + $F[$iw1+2] * (1 - $A_LIST[$i]);
            $X[$i] = $wa;
        }
    }

  print;
  foreach my $i (0 .. $nnew-1) {
      printf " % 15.9g", $X[$i];
  }
  print "\n";

}
close SAMPLE;

sub usage {
  my ($ex) = @_;

  print STDERR "Usage: add_ded.pl [OPTIONS] SAMPLE SPAR [PARAMS]\n";
  print STDERR "OPTIONS:\n";
  print STDERR "   SAMPLE      PMC simulation file, e.g. 'pmcsim'\n";
  print STDERR "   SPAR        Deduced parameter name, one in\n";
  print STDERR "                [Sigma|q_acc|Omega_K|Omega_de|Omega_m|Omega_b|Omega_c|omega_b|omega_c|wa]\n";
  print STDERR "   PARAMS      Further options:\n";
  print STDERR "                SPAR             PARAMS\n";
  print STDERR "                Sigma            alpha Omega_m_fid\n";
  print STDERR "                wa               Na (default 5)\n";
  print STDERR "   -c CONF     Config file CONF (default 'config_pmc')\n";
  print STDERR "   -h          This message\n";

  exit $ex if defined $ex;
}

