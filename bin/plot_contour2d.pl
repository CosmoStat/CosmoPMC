#!/usr/bin/env perl -w

# plot_contour2d.pl
# Martin Kilbinger 2008-2010
# Version 1.2

use Fatal qw/ open /;
use Getopt::Std;
use Cwd;
use File::Basename;

@command =  " $0 @ARGV\n";

### Command line options
%options=();

getopts("c:p:f:t:T:n1:i:s:Srg:G:CN:F:kK:y:w:o:bm:P:qh", \%options);

# Default parameters
$width_def     = 4;

usage() if defined $options{h};

$str = join " ", @command;
open(LOG, ">log_plot"); print LOG $str; close LOG;

# Command line arguments
$do_covar      = defined $options{f} ? $options{f} : 0;
$do_proposal   = defined $options{p} ? $options{p} : 0;
$title         = defined $options{t} ? $options{t} : " ";
$TITLE         = defined $options{T} ? $options{T} : "";
$noshade       = defined $options{n} ? 1 : 0;
$sigma	       = defined $options{s} ? $options{s} : 3;
$square        = defined $options{r} ? 1 : 0;
$smooth_2d_tmp = defined $options{g} ? $options{g} : 0;
$smooth_1d_tmp = defined $options{G} ? $options{G} : $smooth_2d_tmp;
$farb_schema   = defined $options{F} ? $options{F} : 0;
$width         = defined $options{w} ? $options{w} : $width_def;  # Line width
$key           = defined $options{k} ? 1 : 0;
$title_height  = defined $options{y} ? $options{y} : 24;
$output_format = defined $options{o} ? $options{o} : "pdf";
$quiet         = defined $options{q} ? 1 : 0;
$path_flag     = defined $options{P} ? "-P $options{P}" : "";

my $cwd = cwd;
my $path_bin = dirname(__FILE__);
my $yor_inc = "$path_bin/../yorick";

$title = " " if $title eq "";

# Default directory
$ARGV[0] = "." if $#ARGV==-1;

# Create list of smooth factors, 1d and 2d
@smooth_fac_2d = split(" ", $smooth_2d_tmp);
@smooth_fac_1d = split(" ", $smooth_1d_tmp);

if ($#smooth_fac_2d!=0) {
    die "Directory list and 2d-smooth factor list (option -g) have different lengths" 
	if $#ARGV!=$#smooth_fac_2d;
} else {
    foreach $d (0 ..$#ARGV-1) { push @smooth_fac_2d, $smooth_2d_tmp; }
}
if ($#smooth_fac_1d!=0) {
    die "Directory list and 1d-smooth factor list (option -g) have different lengths" 
	if $#ARGV!=$#smooth_fac_1d;
} else {
    foreach $d (0 ..$#ARGV-1) { push @smooth_fac_1d, $smooth_1d_tmp; }
}

# Create list of strings for 'key' in panels
if ($key!=0) {

   # Text for key
   if (! defined $options{K}) {
      if ($#ARGV==-1) {
         print STDERR "Directory is default (\".\"), disabling key flag (-k)\n" unless $quiet;
         $key = 0;
      } else {
        foreach $d (0 .. $#ARGV) {
          $key_str[$d] = $ARGV[$d];
        }
      }
   } else {
      @key_str = split(" ", $options{K});
   }
   die "Number of key-string entries $#key_str differs from number of directories $#ARGV" if $#key_str!=$#ARGV;

   # Replace "_" -> "!_" for yorick
   foreach $d (0 .. $#key_str) {
     $key_str[$d] =~ s/_/\!_/g;
   }

}

die "Output format (Option -o) has to be 'eps' or 'pdf'" if ! ($output_format eq "eps" || $output_format eq "pdf");


# Give warning if shading will hide contours
if (!$noshade) {
  if ($#ARGV>0) {
    print STDERR "Warning: Only one set of contours will be displayed without noshade [-n] option.\n";
  } else {
    $neg = 0;
    foreach $d (0 .. $#smooth_fac_2d) {
      if ($smooth_fac_2d[$d]<0) {
	$neg = 1;
	last;
      }
      if ($neg==1) {
	print STDERR "Warning: Unsmoothed contours will not be displayed without noshade [-n] option.\n";
      }
    }
  }
}


# Normalisation for 1d-posteriors
$norm1d = 0; # Maximum to 1, default
$norm1d = 0 if defined $options{N} && $options{N} eq "m"; # Maximum to 1
$norm1d = 1 if defined $options{N} && $options{N} eq "i"; # Integral to 1


usage() if $sigma<1 || $sigma>3;
usage() if $do_covar!=0 && $do_covar!=1 && $do_covar!=2;

if ($farb_schema == 0) {
   @color = ("", "Blue", "Darkgreen2", "Red", "Orange", "Violett", "Magenta");
} elsif ($farb_schema == 1) {
   @color = ("", "Orange", "Darkgreen2", "Magenta", "Blue", "Red", "Cyan");
} elsif ($farb_schema == 2) {
   @color = ("", "Red", "Red", "Blue", "Blue", "Violett", "Red");
} else {
   die "Color scheme (flag F) unknown\n";
}


# For smoothing
$FFT_scale = 2;


$dir0 = $ARGV[0];

### Read config file
$configname = "config_pmc" if -e "config_pmc";
$configname = $options{c} if defined $options{c};
die "No configuration file found. Use option '-c CONFIG_FILE'." unless defined $configname;

open(CONFIG, "$configname");
while (<CONFIG>) {
    @F = split(" ", $_);
    next if /^#/;
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
	        $spar[$i] = $F[$i+1];
	    }
    }
    if ($F[0] eq "min") {
	    die "Entry 'npar' not found in configuration file" if ! defined $npar;
	    die "Entry 'n_ded' not found in configuration file" if ! defined $n_ded;
	    foreach $i (0 .. $npar+$n_ded-1) {
	        $min[$i] = $F[$i+1];
	    }
    }
    if ($F[0] eq "max") {
	    foreach $i (0 .. $npar+$n_ded-1) {
	        $max[$i] = $F[$i+1];
	    }
    }
}
close CONFIG;


# Get yorick parameter names
$cmd = "$path_bin/get_spar.pl -c $configname yorick $path_flag";
$output = qx($cmd);
die "'get_spar.pl' didn't finish successfully" unless $?==0;
@parlist = split("&", $output);
if ($#parlist<$npar+$n_ded) {
   print STDERR "Warning: parameter list smaller than number of parameters\n" unless $quiet;
   foreach $i (0 .. $npar+$n_ded-1) {
     $parlist[$i] = "p_$i";
  }
}

### Create yorick file
open(YOUT, ">plot_contour2d.i");

print YOUT "include, \"$yor_inc/likeli.i\"\n";
print YOUT "include, \"$yor_inc/stuff.i\"\n\n";

# Read covariance or Fisher matrix
if ($do_covar==1) {
    print YOUT "covar = mvdens_read(\"" . $dir0 . "/covar+ded.fin\")\n";
}
if ($do_covar==2) {
    print YOUT "covar = mvdens_read(\"" . $dir0 . "/fisher\")\n";
    print YOUT "inv = LUsolve(*(covar.var))\n";
    print YOUT "covar.var = &inv\n";
}

# Read PMC proposal
if ($do_proposal > 0) {
    if (-e "proposal_0") {
	$pdir = ".";
    } elsif (-e "../proposal_0") {
	$pdir = "..";
    } else {
	print STDERR "File 'proposal_0' not found, continuing without plotting proposal\n" unless $quiet;
	$do_proposal = 0;
    }
}
if ($do_proposal > 0) {
    if (!defined $options{i}) {
	$niter = 1;
	print STDERR "Setting niter to $niter (use -i flag to change default)\n" unless $quiet;
    } else {
	$niter = $options{i};
    }
    print YOUT "mix_mvd = read_all_proposals($niter, dir=\"$pdir\")\n";
}


### 2d marginal contours
print YOUT "write, format=\"%s\\n\", \"2d plots\"\n" unless $quiet;
$i = 0;
while ($i < $npar+$n_ded-1) {
    $j = $i + 1;
    while ($j < $npar+$n_ded) {

	$suff = "$i" . "_" . "$j";
	print YOUT "window, display=\"\", hcp=\"contour2d_$suff\"\n";

	print YOUT "pltitle_height = $title_height\n";
	print YOUT "limits, $min[$i], $max[$i]\n" if defined $min[$i];
	print YOUT "range, $min[$j], $max[$j]\n" if defined $min[$j];

	### Loop over directories
	$d = 0;
	while ($dir = $ARGV[$d]) {

	    #print STDERR "dir = $dir\n";
	    print YOUT "imin = read_chin(\"$dir/chi2_$suff\", 2, prior=0, quiet=1)\n";
	    print YOUT "chi2 = matrix(chinall(3,), n=nel)\n";

	    if ($smooth_fac_2d[$d]<=0) {
		# Plot unsmoothed contours
		contour($smooth_fac_2d[$d]<0);
	    }
	    if ($smooth_fac_2d[$d]!=0) {
		if (defined $options{C}) {
		    # Smoothing using covariance
		    print YOUT "ind = [$i+1, $j+1]\n";
		    print YOUT "cov = mvdens_read(\"" . $dir . "/covar+ded.fin\")\n";
		    print YOUT "var  = (*cov.var)(ind,ind)\n";
		    print YOUT "var(1,1) = var(1,1) * dimsof(chi2)(2)^2/($max[$i]-($min[$i]))/($max[$i]-($min[$i]))\n";
		    print YOUT "var(2,2) = var(2,2) * dimsof(chi2)(2)^2/($max[$j]-($min[$j]))/($max[$j]-($min[$j]))\n";
		    print YOUT "var(1,2) = var(1,2) * dimsof(chi2)(2)^2/($max[$i]-($min[$i]))/($max[$j]-($min[$j]))\n";
		    print YOUT "var(2,1) = var(2,1) * dimsof(chi2)(2)^2/($max[$i]-($min[$i]))/($max[$j]-($min[$j]))\n";
		    print YOUT "Sigma = var/abs($smooth_fac_2d[$d])\n";
		    print YOUT "chi2 = smooth2d(chi2, dimsof(chi2)(2:), Sigma=Sigma, scale=$FFT_scale)\n";
		} else {
		    # Isotropic Gaussian smoothing
		    print YOUT
			"chi2 = smooth2d(chi2, dimsof(chi2)(2:), factor=abs($smooth_fac_2d[$d]), scale=$FFT_scale)\n";
		}
		# Plot smoothed contours
		contour(0);
	    }

	    $d++;
	}

	printf YOUT "xytitles, \"%s\", \"%s\", [-0.01, 0.015]\n", $parlist[$i], $parlist[$j];
	print YOUT "pltitle, \"$title\"\n";

	printf YOUT "limits, %f, %f, %f, %f\n", $min[$i], $max[$i], $min[$j], $max[$j] if defined $min[$j];
	print YOUT "l = limits()\n";

	if ($square==1) {
	    print YOUT "dx = l(2)-l(1)\n";
	    print YOUT "dy = l(4)-l(3)\n";
	    print YOUT "if (dy>dx) {\n";
	    print YOUT "   r = (dy-dx)/2\n";
	    print YOUT "   l(2) += r\n";
	    print YOUT "   l(1) -= r\n";
	    print YOUT "} else {\n";
	    print YOUT "   r = (dx-dy)/2\n";
	    print YOUT "   l(4) += r\n";
	    print YOUT "   l(3) -= r\n";
	    print YOUT "}\n";
	    print YOUT "limits, l(1), l(2), l(3), l(4)\n";
	}

	if ($do_covar != 0) {
	    # Chain covariance ellipses
	    print YOUT "ind = [$i+1, $j+1]\n";
	    print YOUT "mean = (*covar.mean)(ind)\n";
	    print YOUT "var  = (*covar.var)(ind,ind)\n";
	    #print YOUT "[$i, $j]\n";

	    print YOUT "for (s=dimsof(sig)(2); s>=1; s--) {\n";
	    print YOUT "   sigma = sqrt(-2*log(1.0-erf(s/sqrt(2))))\n";
	    print YOUT "   Z = ellipse(mean, var, sigma=sigma, N=100)\n";
	    print YOUT "   draw_ellipse, Z, type=2\n";
	    #print YOUT "   draw_ellipse, Z, type=2, with_ab=(s==1?1:0), mean=mean, sigma=sigma\n";
	    print YOUT "}\n";

	    print YOUT "x = l(2)-0.05*(l(2)-l(1))\n";
	    print YOUT "y = l(4)-0.05*(l(4)-l(3))\n";
	    print YOUT "str = swrite\(format=\"[ci] %.4f \(%.4f,%.4f\)\", eigenvector\(2,1\)/eigenvector\(1,1\), eigenvector\(1,1\), eigenvector\(2,1\)\)\n";
	    print YOUT "plt, str, x, y, justify=\"RT\", tosys=1\n";
	    print YOUT "x = l(2)-0.05*(l(2)-l(1))\n";
	    print YOUT "y = l(4)-0.1*(l(4)-l(3))\n";
	    print YOUT "str = swrite\(format=\"[ci] %.4f \(%.4f,%.4f\)\", eigenvector\(2,2\)/eigenvector\(1,2\), eigenvector\(1,2\), eigenvector\(2,2\)\)\n";
	    print YOUT "plt, str, x, y, justify=\"RT\", tosys=1, color=Red\n";
	}

        # Proposal not for deduced parameters
	if (($do_proposal & 1) && $j<$npar) {
	    print YOUT "plot_mix_mvdens_mean, mix_mvd, [$i+1, $j+1], msize=5, connect=1\n";
	}
	if (($do_proposal & 2) && $j<$npar) {
	    print YOUT "plot_mix_mvdens_covar, mix_mvd\($niter-1\), [$i+1, $j+1]\n";
	}

	if (defined $options{m}) {
	  @mpar = split(" ", $options{m});
      die "Marked data point (option -m) has different dimension (", $#mpar + 1, ") than given in config file (", $npar + $n_ded, ")"
        if $#mpar + 1 != $npar + $n_ded;
	  print YOUT "plmk, marker=4, ", $mpar[$j], ", ", $mpar[$i], ", width=10, msize=0.5, color=Red\n";
	}
	print YOUT "$output_format, \"contour2d_$suff\"\n";

	# density plot
	print YOUT "winkill, 0\n";
	print YOUT "window, display=\"\", hcp=\"density_$suff\"\n";
	print YOUT "palette, \"heat.gp\"\n";
	printf YOUT "fma; pli, -chi2, %f, %f, %f, %f\n", $min[$i], $min[$j], $max[$i], $max[$j] if defined $min[$j];
	print YOUT "xytitles, \"$parlist[$i]\", \"$parlist[$j]\", [-0.01, 0.015]\n";
	print YOUT "pltitle, \"$title\"\n";
	print YOUT "$output_format, \"density_$suff\"\n";

	print YOUT "\n";
	$j++;
	#if ($j==4) { print YOUT "quit\n"; }

    }
    $i++;
}


### 1d marginal posteriors
if (defined $options{1}) {

    $do_mean = $options{1} =~ "m" ? 1 : 0; 
    $sig1    = $options{1} =~ "1" ? 1 : 0;
    $sig2    = $options{1} =~ "2" ? 2 : 0;
    $sig3    = $options{1} =~ "3" ? 3 : 0;
    $text    = $options{1} =~ "t" ? 1 : 0;
    if ($options{1} =~ "n") {
       $do_mean = $sig1 = $sig2 = $sig3 = $text = 0;
    }
    $do_sigma = "[" . $sig1 . ", " . $sig2 . ", " . $sig3 . "]";
    $text = 0 if $sig1==0;

    print YOUT "write, format=\"%s\\n\", \"1d plots\"\n" unless $quiet;
 
    print YOUT "pmaxall = array(double, $npar+$n_ded, $#ARGV+1)\n";

    $i = 0;
    while ($i < $npar+$n_ded) {

	print YOUT "window, display=\"\", hcp=\"likeli1d_$i\"\n";

	### Loop over directories
	$d = 0;
	while ($dir = $ARGV[$d]) {

	    print YOUT "imin = read_chin(\"$dir/chi2_$i\", 1, prior=0, quiet=1)\n";
	    print YOUT "meansig = array(double, 3)\n";
	    print YOUT "chi2 = chinall(2,)\n";

	    if ($smooth_fac_1d[$d]<=0) {
	        plot_1d(1);
	    }
	    if ($smooth_fac_1d[$d]!=0) {
	        print YOUT "chi2 = smooth1d(chi2, dimsof(chi2)(2), abs($smooth_fac_1d[$d]), scale=$FFT_scale)\n";
	        plot_1d(0);
	    }

	    printf YOUT "limits, %f, %f\n", $min[$i], $max[$i] if defined $min[$i];
	    print YOUT "range, 0, 1\n" if $norm1d==0;
	    print YOUT "range, 0\n" if $norm1d==1;
	    print YOUT "l = limits()\n";
            if ($text == 1) {
	       print YOUT "x = l(2)-0.05*(l(2)-l(1))\n";
	       print YOUT "str = swrite(format=\"%.3f+%.3f-%.3f\", meansig(1), meansig(1)-meansig(2), meansig(3)-meansig(1))\n";
	       print YOUT " plt, str, x, l(4)-0.05*($d+1), justify=\"RT\", tosys=1, color=$color[$d+1]\n";
            }

	    if ($do_proposal>0 && -e "prop_likeli1d_$i") {
		print YOUT "read_table, \"prop_likeli1d_$i\", like\n";
		print YOUT "plg, like(2,)/like(2,max), like(1,), marks=0\n";
	    }

	    printf YOUT "limits, %f, %f\n", $min[$i], $max[$i] if defined $min[$i];
	    $d++;

	}

	print YOUT "pltitle_height = $title_height\n";
	print YOUT "xytitles, \"$parlist[$i]\", \"posterior\", [-0.01, 0.025]\n";
	print YOUT "pltitle, \"$title\"\n";
	print YOUT "$output_format, \"likeli1d_$i\"\n";
	print YOUT "\n\n";

	$i++;

    }

    # Write maximum posterior values to file "max_post"
    print YOUT "f = open(\"max_post\", \"w\")\n";
    print YOUT "spar = array(string, $npar+$n_ded)\n";
    foreach $i (1 .. $npar+$n_ded) {
      print YOUT "spar($i) = \"$spar[$i-1]\"\n";
    }
    print YOUT "for (i=1; i<=$npar+$n_ded; i++) {\n";
    print YOUT "   write, f, format=\"%2d %15s\", i-1, spar(i)\n";
    print YOUT "   for (d=1; d<=$#ARGV+1; d++) {\n";
    print YOUT "     write, f, format=\" % 10.5g\", pmaxall(i,d)\n";
    print YOUT "   }\n";
    print YOUT "   writeNL, f\n";
    print YOUT "}\n";
    print YOUT "close, f\n";

}


print YOUT "quit\n";
close YOUT;

### Call yorick, tex et al.
print STDERR "Calling yorick...\n" unless $quiet;
system("yorick -batch plot_contour2d.i");
exit($?<<8) if $?!=0;

if (-e "likeli1d_0.$output_format") {
    $cmd = "$path_bin/all_vs_all.pl -t \"$TITLE\" -b contour2d -e $output_format -l likeli1d > all_contour2d.tex";
} else {
    $cmd = "$path_bin/all_vs_all.pl -t \"$TITLE\" -b contour2d -e $output_format > all_contour2d.tex";
}
`$cmd`;
if (-e "density_0_1.$output_format") {
  $cmd = "$path_bin/all_vs_all.pl -t \"$TITLE\" -b density -e $output_format > all_density.tex";
  `$cmd`;
}

# Change line styles in eps files
if ($output_format eq "eps") {
  $i = 0;
  while ($i < $npar+$n_ded-1) {
    $j = $i + 1;
    while ($j < $npar+$n_ded) {
      $suff = "$i" . "_" . "$j";
      $name = "contour2d_$suff.eps";
      my $cmd = "$path_bin/change_lt_yor.pl $name > tmp.eps";
      `$cmd`;
      `mv tmp.eps $name`;
      $j++;
    }
    $i++;
  }
}


print STDERR "Calling latex...\n" unless $quiet;

if (-s "all_contour2d.tex") {
  $cmd = "$path_bin/ldp.sh all_contour2d.tex -q" if $output_format eq "eps";
  $cmd = "pdflatex all_contour2d.tex > /dev/null" if $output_format eq "pdf";
  `$cmd`;
  `rm -f all_contour2d.{log,dvi,aux}`;
}

if (-s "all_density.tex") {
  $cmd = "ldp.sh all_density.tex -q" if $output_format eq "eps";
  $cmd = "pdflatex all_density.tex > /dev/null" if $output_format eq "pdf";
  `$cmd`;
  `rm -f all_density.{log,dvi,aux}`;
}



##################

sub contour {
    my ($add_nosm) = @_;

    my $dd = $d+1;

    if ($add_nosm==0) {
	$w = $width;
	$first = 0;
    } else {
	$w = defined $smooth_fac_2d[$d]==0 ? $width : 1;
	$first = $dd==1? 1 : 0;
    }

    print YOUT "c2 = chi2(where(chi2>=0))\n";
    print YOUT "s = sort(c2)(::-1)\n";
    print YOUT "c2cum = c2(s)(cum)(2:)\n";
    print YOUT "n = sum(c2)\n";
    print YOUT "sig = [0.6827, 0.9545, 0.9973]\n";
    print YOUT "sig = sig(1:$sigma)\n";
    print YOUT "sigind = sig*n\n";
    print YOUT "levels1 = array(int, dimsof(sig)(2))\n";
    print YOUT "for (i=1; i<=dimsof(sig)(2); i++) { levels1(i) = where(c2cum>=sigind(i))(1); }\n";

    print YOUT "xlevels = array(double, dimsof(sig)(2))\n";
    #print YOUT "for (i=1; i<=dimsof(sig)(2); i++) { xlevels(i) = (sigind(i)-c2cum(levels1(i)-1))/(c2cum(levels1(i))-c2cum(levels1(i)-1)) + levels1(i)-1; }\n";
    print YOUT "y = xlevels - int(xlevels)\n";
    print YOUT "for (i=1; i<=dimsof(sig)(2); i++) { xlevels(i) = c2(s)(levels1(i)-1)*(1-y(i)) + c2(s)(levels1(i))*y(i); }\n";
    # Following line: no interpolation
    print YOUT "levels1 = xlevels\n";

    #print "*** $dd $key $key_str[$d]\n";
    $key_str[$d] = "" if ! defined $key_str[$d];
    printf YOUT "plot_chi2n, chi2, 1, 2, quiet=1, do_fma=%d, levels1=levels1, alllevels=0, width=$w, " .
      "type=%d, Col=[$color[$dd], $color[$dd], $color[$dd]]%s, key=%d, key_str=\"$key_str[$d]\", height=$title_height\n", 
	  $first, defined $options{S} ? 1 : $dd, $noshade==1 ? "" : ", shade=\"$color[$dd]\"", $key==0?0:$dd;

    # New
    if (defined $options{b}) {
        printf YOUT "ff = open(\"chi2_$i" . "_$j.block\", \"w\")\n";
        printf YOUT "write_matrix, file=ff, chi2\n";
        printf YOUT "close, ff\n";
    }
}

sub plot_1d {
    my ($add_nosm) = @_;

    my $dd = $d+1;

    if ($add_nosm==0) {
	$w = $width;
	$first = 0;
    } else {
	$w = $smooth_fac_1d[$d]==0 ? $width : 1;
	$first = $dd==1? 1 : 0;
    }

    $key_str[$d] = "" if ! defined $key_str[$d];
    printf YOUT "plot_chi2_1d, chi2, 1, meansig, pmax, mcmc=1, do_fma=%d, do_mean=$do_mean, do_sigma=$do_sigma, nice=1, " .
      "outputfile=\"likeli1d_$i.txt\", text=\"$parlist[$i]\", Type=%d, color=$color[$dd], width=$w, norm=$norm1d, " .
	"key=%d, key_str=\"$key_str[$d]\", height=$title_height\n",
	  $first, defined $options{S} ? 1 : $dd, $key==0?0:$dd;

    print YOUT "pmaxall($i+1,$d+1) = pmax\n";
}

sub usage {
    print STDERR "Usage: plot_contour2d.pl [OPTIONS] [DIR1 [DIR2 [...]]]\n";
    print STDERR "\nOPTIONS:\n";
#    print STDERR "  -f DO_COVAR    Plots covariance on top of contours\n";
#    print STDERR "                  0: No action (default)\n";
#    print STDERR "                  1: Final covariance (covar.fin)\n";
#    print STDERR "                  2: Fisher matrix (fisher)\n";
#    print STDERR "  -p DO_PROPOSAL Plots (PMC) proposal on top of contours\n";
#    print STDERR "                  0: No action (default)\n";
#    print STDERR "                  1: Weights (positions of mean)\n";
#    print STDERR "                  2: Directions (covariance)\n";
    print STDERR "  -i NITER       Number of iterations (needed if do_proposal=2)\n";
    print STDERR "  -c CONFIG_FILE Configuration file (default: in order config_mcmc, config_pmc)\n";
    print STDERR "  -t TITLE       Title string for each panel (default empty)\n";
    print STDERR "  -T TITLE       Title string for all_contour2d.{eps|pdf} (default empty)\n";
    print STDERR "  -n             No shade\n";
    print STDERR "  -w WIDTH       Line width WIDTH (default $width_def)\n";
    print STDERR "  -1 OPT         Add 1d posterior plots. OPT can contain the following letters:\n";
    print STDERR "                  m    Plot line at mean position\n";
    print STDERR "                  123  Plot line at 68%,95%,99.7 density\n";
    print STDERR "                  t    Write mean and 68% confidence intervals as text\n";
    print STDERR "                        (use with 'm' and '1'\n";
    print STDERR "                  n    None of the above\n";
    #print STDERR "  -d N_DED       Number of deduced parameters (default 0)\n";
    print STDERR "  -S             All contours with solid lines\n";
    print STDERR "  -s N           Outermost level is N sigma\n";
    print STDERR "  -r             Aspect ratio=1, changes plot limits such that dx=dy\n";
    print STDERR "  -g FACTOR      Gaussian smoothing of 2d-histograms with variance\n";
    print STDERR "                  box-width/|FACTOR|. If FACTOR is negativ, plots\n";
    print STDERR "                  unsmoothed histogram in addition (use with '-n').\n";
    print STDERR "                  Note: For multiple contours, use a list of values \"g1 g2 ...\"\n";
    print STDERR "  -G FACTOR      Gaussian smoothing of 1d-histograms (default: 2d factor)\n";
    print STDERR "  -C             Use covariance (file covar.fin) for Gaussian smoothing\n";
    print STDERR "  -N NORM        Normalisation of 1d posterior\n";
    print STDERR "                  'm'  Maximum = 1 (default)\n";
    print STDERR "                  'i'  Integral over posterior = 1\n";
    print STDERR "  -F NUM         Color scheme, NUM=0,1,2\n";
    print STDERR "  -k             Add key to plots\n";
    print STDERR "  -K \"KEY1 [KEY2 [...]]\" Key strings (default: directory names)\n";
    print STDERR "  -y FS          Font size FS (default 24)\n";
    print STDERR "  -o FORMAT      Output file format, FORMAT=eps|pdf (default: eps)\n";
    print STDERR "  -b             Writes the chi2 files in block format\n";
    print STDERR "  -m PAR         Plots a mark at position PAR (e.g. best-fit). PAR is white-space\n";
    print STDERR "                  separated list (use quotes or '\\ ', e.g. '0.3 0.8')\n";
    print STDERR "  -q             Run quietly, no verbose\n";
    print STDERR "  -h             This message\n";
    print STDERR "  DIR1 ...       List of directories containing histogram files (chi2_*_*)\n";
    print STDERR "                  Default: DIR1 = '.'\n";
    exit 1;
}
