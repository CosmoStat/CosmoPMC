#!/usr/bin/perl -w

# test_suite_cosmo_pmc.pl
# Martin Kilbinger 2010
# Performs test runs of Cosmo_PMC. The output can be compared
# between different versions or installations, to check for
# consistency.


use Fatal qw/ open close unlink chdir /;
use Getopt::Std;
use Cwd qw/ cwd abs_path /;



my %options = ();
getopts("rRn:cP:skvxh", \%options);

usage(0) if defined $options{h};

my $path_flag = "";
if (defined $options{P}) {
  my $path = abs_path($options{P});
  $path_flag = "-P $path";
}

my $verbose   = defined $options{v} ? 1 : 0;
my $verb_flag = defined $options{v} ? "-v" : "";
my $keep_flag = defined $options{k} ? "-k" : "";

my $cwd = cwd;
set_ENV_COSMOPMC();
die "Environment variable \$COSMOPMC = '$ENV{COSMOPMC}' not valid.\nSet \$COSMOPMC in the shell or use the '-P' option"
  unless -s "$ENV{COSMOPMC}/bin/cosmo_pmc.pl";

my $version = get_version();

my @subdir_pmc  = ("SN",
		   "Lensing/COSMOS-S10",
		   "BAO/distance_A",
		   "BAO/distance_d_z",
		   "WMAP_Distance_Priors",
		   "COSMOS-S10+SN+BAO",
		   "Lensing/CFHTLenS-K13",
		   "HOD/CFHTLS-T06");
my @options_pmc = ("", "-d", "-d", "-d", "-d", "", "", "-d -f '12.1 13.3 12.1 0.3 1.0'");
my @short_pmc   = (1, 0, 1, 1, 1, 0, 1, 0);

# Fiducial point for test_log_post
my @fid    = ("0.27 -1.0 19.31 1.6 -1.8",
	      "0.27 0.73 0.8 0.71 1.0",
	      "0.27 0.73",
	      "0.27 0.73",
	      "0.045 0.27 0.73 0.71",
	      "0.27 -1.0 0.8 0.7 1.0 19.31 1.6 -1.8",
	      "0.27 0.8 0.045 0.7 0.96",
	      "12.2 13.2 11.1 0.366 1.16");


####################
### Main program ###
####################

print STDERR "Running 'test_suite_cosmo_pmc.pl' in directory '$cwd'\n";

if (defined $options{x}) {
    print STDERR "Cleaning previous run\n";
    clean_all();
    exit 0;
}

unlink "suite.ps" if -e "suite.ps";

if (! defined $options{R}) {
  test_demo($cwd);
  $base = "numbers_logP";
  test_log_post($base, \@subdir_pmc, \@fid);
}

if (defined $options{r} or defined $options{R}) {
  my $ncpu = defined $options{n} ? $options{n} : 1;
  test_pmc($cwd, $ncpu);
}

my $title = "\"CosmoPMC v$version\"";
run_own_command("allps2tex.pl -t $title -f eps -w 12 Demo . > suite.tex");
run_own_command("ldp.sh suite.tex");
run_sys_command("pstopdf suite.ps");

my $cmd = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=suite_all.pdf suite.pdf";
foreach $s (@subdir_pmc) {
  $file = "Demo/MC_Demo/$s/essential.pdf";
  $cmd .= " $file" if -s $file;
}
print $cmd;
run_sys_command($cmd);
print "\n";

# Clean up
clean_parts();


###########
### end ###
###########


###################
### Subroutines ###
###################

# Performs tests on directory 'Demo'
sub test_demo {
  my ($cwd) = @_;

  chdir "Demo";
  `rm -rf *.eps`;


  ### Plots

  # Power spectrum, lensing
  run_sys_command("./lensingdemo");
  $plot{"logx"} = 1;
  $plot{"logy"} = 1;

  $plot{"ylabel"} = "P_{/Symbol d}(k, a)";
  plot_file("P_delta", %plot);

  $plot{"ylabel"} = "{/Symbol x}^{ij}_+({/Symbol q})";
  plot_file("xi_p", %plot);

  $plot{"ylabel"} = "{/Symbol x}^{ij}_-({/Symbol q})";
  plot_file("xi_m", %plot);

  # CMB
  if (defined $options{c} && check_CMB_support($cwd) ) {
    # The following command gets stuck in not used with 'system'
    system("./cmbdemo -c config_max_WMAP7");
    $plot{"logx"} = 1;
    $plot{"xlabel"} = "k[h/Mpc]";
    $plot{"ylabel"} = "P_{{/Symbol d}, lin}(k, a=0)";
    plot_file("Pdelta.dat", %plot);
    $plot{"logx"} = $plot{"logy"} = 0;
    $plot{"xlabel"} = "{/CMMI10 \140}";
    $plot{"ylabel"} = "{/CMMI10 \140}({/CMMI10 \140}+1)C_{/CMMI10 \140}/(2{/Symbol p})";
    plot_file("lensedCls.dat", %plot);
  }

  # SNIa
  run_sys_command("./sn1ademo");
  $plot{"logx"} = 0;
  $plot{"logy"} = 1;
  undef $plot{"xlabel"};
  $plot{"ylabel"} = "Luminosity distance [Mpc/h]";
  plot_file("D_Lum", %plot);


  ### Numbers

  # Power spectrum, lensing
  open(my $out_fh, ">numbers_0_Pdelta.txt");
  $txt{"func"} = "P_\\delta";
  print_number("P_delta", 109, $out_fh, %txt);
  print_number("P_delta", 230, $out_fh, %txt);
  close $out_fh;

  open($out_fh, ">numbers_1_xip.txt");
  $txt{"func"} = "\\xi_+";
  print_number("xi_p", 26, $out_fh, %txt);
  $txt{"func"} = "\\xi_-";
  print_number("xi_m", 62, $out_fh, %txt);
  close $out_fh;

  # CMB
  if (defined $options{c} && check_CMB_support($cwd)) {
    open($out_fh, ">numbers_2_Pdelta.dat.txt");
    $txt{"func"} = "P_{\\textrm{lin}}";
    print_number("Pdelta.dat", 200, $out_fh, %txt);
    print_number("Pdelta.dat", 450, $out_fh, %txt);
    close $out_fh;

    open($out_fh, ">numbers_3_Cl.dat.txt");
    $txt{"func"} = "\\ell(\\ell+1)C_\\ell/(2\\pi)";
    print_number("lensedCls.dat", 100, $out_fh, %txt);
    print_number("lensedCls.dat", 300, $out_fh, %txt);
    close $out_fh;
  }

  # SNIa
  open($out_fh, ">numbers_4_D_Lum.txt");
  $txt{"func"} = "D_{\\textrm{lum}}";
  print_number("D_Lum", 10, $out_fh, %txt);
  print_number("D_Lum", 80, $out_fh, %txt);

  txt2eps("numbers_0_Pdelta", "-m", "Cosmological quantities");
  txt2eps("numbers_1_xip");
  if (defined $options{c} && check_CMB_support($cwd)) {
     txt2eps("numbers_2_Pdelta.dat");
     txt2eps("numbers_3_Cl.dat");
   }
  txt2eps("numbers_4_D_Lum");

  # Clean up
  unlink <P_delta P_kappa xi_p xi_m mapsqr mapsqr_gauss gammasqr nz w_eos> unless defined $options{k};
  unlink <P_delta.gnu xi_p.gnu xi_m.gnu> unless defined $options{k};

  if (defined $options{c} && check_CMB_support($cwd)) {
    unlink <scalCls.dat lensedCls.dat Pdelta.dat> unless defined $options{k};
    unlink <lensedCls.dat.gnu Pdelta.dat.gnu> unless defined $options{k};
  }

  unlink <D_Lum D_Lum.gnu> unless defined $options{k};

  $base = "numbers_*_*";
  unlink <$base.log numbers_*_*.aux numbers_*_*.dvi>;
  unlink <numbers_*_*.txt numbers_*_*.tex numbers_*_*.tab> unless defined $options{k};

  chdir "..";

}

# Prints a line from an ascii data file
sub print_number {
  my ($name, $line, $out_fh, %txt) = @_;

  my $func = defined $txt{"func"} ? $txt{"func"} : "f";

  open(my $in_fh, "$name");
  my $nl = 0;
  my @header = ();
  while (<$in_fh>) {

    @F = split(" ", $_);
    if (/#/) {
      @header = @F;
      next;
    }

    if ($nl == $line) {
      printf {$out_fh} "$func(%g) =", $F[0];
      foreach $n (1 .. $#F) {
	printf {$out_fh} " %g ", $F[$n];
	printf {$out_fh} " (%s)", $header[$n+1] if $#header > -1;
	print {$out_fh} "," unless $n == $#F;
      }
      print {$out_fh} "\n";
      last;
    }

    $nl ++;
  }
  close $in_fh;
}

# Calls gnuplot to plot file with name $name. Creates output
# file 'fig_$name.eps'
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
  print {$gnu_fh} "set size square\n";
  print {$gnu_fh} "set output 'fig_$name.eps'\n";
  print {$gnu_fh} "set term post eps enhanced color 'Times-Roman' 20\n";
  print {$gnu_fh} "set logs x\n" if has_value($plot{"logx"}, 1);
  print {$gnu_fh} "set logs y\n" if has_value($plot{"logy"}, 1);
  if (defined $plot{"key"}) {
    $key = $plot{"key"};
    print {$gnu_fh} "set key '$key'\n";
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
    print {$gnu_fh} "u 1:(\$$n) w l t '$title'";
    print {$gnu_fh} ", " unless $n==$ncol;
  }
  print {$gnu_fh} "\n";
  print {$gnu_fh} "set output\n";
  close $gnu_fh;

  run_sys_command("gnuplot $name.gnu");
}

# Performs tests in MC_Demo sub-directories
sub test_log_post {
  my ($base, $subdir_ref, $fid_ref) = @_;

  # Copy arrays
  my $subdir = [@{$subdir_ref}];
  my $fid    = [@{$fid_ref}];

  push(@$subdir, "WMAP7");
  push(@fid, "0.045 0.27 0.087 0.96 2.41 0.71");

  # Note: 'WMAP7' has to be last item in the list (cmb check)

  open(my $out_fh, ">$base.txt");

  print {$out_fh} "{\\Large \\textbf{Log-posteriors}}\\bigskip\n\n";
  print {$out_fh} "\\begin{tabular}{lll}\n";

  my $nsubdir = $#{$subdir};
  foreach my $i (0 .. $nsubdir) {

    next if $i == $nsubdir && ! (defined $options{c} && check_CMB_support($cwd));

    my $cwd = cwd;
    my $dir = "Demo/MC_Demo/" . "$subdir->[$i]";
    chdir "$dir";
    print STDERR "*** Entering $dir ***\n";

    my $maxlogP = print_log_post($fid->[$i]);

    $subdir->[$i] =~ s/_/\\_/g;
    print {$out_fh} "$subdir->[$i]  fid & = ($fid->[$i])  log(posterior) & = $maxlogP \\\\\n";

    unlink <maxlogP log_max_post config_max_test_suite> unless defined $options{k};

    chdir "$cwd";
  }

  print {$out_fh} "\\end{tabular}\n";

  close $out_fh;

  run_own_command("txt2tex.pl $base.txt > $base.tex");
  run_own_command("lde.sh $base.tex");

  unlink <$base.log $base.aux $base.dvi>;
  unlink <$base.tex $base.txt> unless defined $options{k};
}


sub test_pmc {

  my ($cwd, $ncpu) = @_;

  foreach $i (0 .. $#subdir_pmc) {

    next if $short_pmc[$i] == 0 and defined $options{s};

    my $dir = "Demo/MC_Demo/$subdir_pmc[$i]";
    chdir "$dir";
    print STDERR "*** Entering $dir ***\n";
    run_own_command("cosmo_pmc.pl -n $ncpu -A n $options_pmc[$i]");

    # Look for config file with deduced parameters
    if (-s "config_pmc_ded") {
        mkdir "deduced" unless -d "deduced";
        my $last_iter = get_niter("config_pmc") - 1;
        run_own_command("add_deduced_halomodel -c config_pmc_ded -o deduced/pmcsim+ded iter_$last_iter/pmcsim");
        chdir "deduced";
        run_own_command("histograms_sample -c ../config_pmc_ded pmcsim+ded");
        run_own_command("plot_contour2d.pl -c ../config_pmc_ded -1 m1t -o pdf");
        chdir "..";
    }

    run_own_command("essential_cosmo_pmc_run.pl $path_flag $verb_flag $keep_flag");
    chdir "$cwd";
  }

}

# Returns niter, read from the config file
sub get_niter {
  my ($config) = @_;
  open(my $conf_fh, "$config");
  my $niter = -1;

  while (<$conf_fh>) {
    if (/^niter/) {
      my @F = split(" ", $_);
      $niter = $F[1];
    }
 }
 close $conf_fh;

 return $niter;
}


# Creates a eps file from the txt file '$base.txt'
sub txt2eps {
  my ($base, $mflag, $title) = @_;

  my $tflag = defined $title ? "-t '$title'" : "";
  $mflag = "-m" if ! defined $mflag;

  run_own_command("tab2tex.pl $mflag -l n $tflag $base.txt > $base.tab");
  run_own_command("txt2tex.pl $base.tab > $base.tex");
  run_own_command("lde.sh $base.tex");
}

# Runs 'max_post -m n'
sub print_log_post {
  my ($fid) = @_;

  run_own_command("config_pmc_to_max_and_fish.pl -M -f '$fid' > config_max_test_suite");
  #run_own_command("max_post -m n -c config_max_test_suite");
  system("$ENV{COSMOPMC}/bin/max_post -m n -c config_max_test_suite");

  my $maxlogP = get_max_logP();

  return $maxlogP;
}

# Prints max[log(posterior)] from the 'maxlogP' file
sub get_max_logP {

  open(my $max_fh, "maxlogP");
  my $line = <$max_fh>;
  close $max_fh;

  my ($maxlogP) = $line =~ /P\s*=\s*(\S*)\s*ftol/;
  return $maxlogP;
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

# Returns 1 if $var is defined and $var = $val
sub has_value {
  my ($var, $val) = @_;
  return 0 if ! defined $var;
  return $var == $val;
}

# Return the value of the variable 'DOWMAP' in '$cwd/Makefile.host'
sub check_CMB_support {
  my ($cwd) = @_;

  if (! -e "$cwd/Makefile.host") {
    print "No file 'Makefile.host' found\n";
    return;
  }

  open(my $mk_fh, "$cwd/Makefile.host");
  my $dowmap = 0;
  while (<$mk_fh>) {
    if (/DOWMAP\s*=\s*(\d+)/) {
      $dowmap = $1;
      last;
    }
  }
  close $mk_fh;

  return $dowmap;
}

# Returns the PMC version
sub get_version {
    open(my $fh, "$ENV{COSMOPMC}/tools/include/stdnames.h");
    my $version = -1;
    while (<$fh>) {
	if ($_ =~ /COSMO_PMC_V/) {
	    my @F = split(" ", $_);
	    $version = $F[2];
	    last;
	}
    }
    close $fh;

    return $version;
}

sub set_ENV_COSMOPMC {
  my ($cwd) = @_;

  # Copy from command argument
  if (defined $options{P}) {
    $ENV{COSMOPMC} = abs_path($options{P});
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

sub clean_parts {
    unlink <suite.log suite.aux suite.dvi>;
    unlink <suite.tex suite.ps> unless defined $options{k};
    unlink <$base.eps> if ! defined $options{R} and ! defined $options{k};
}

sub remove {
    my $dir  = $_[0];
    my @FILE = @{ $_[1] };
    my @EXT  = @{ $_[2] };

    foreach $file (@FILE) {
        foreach $ext (@EXT) {
            my $obj = $ext eq "NIX" ? "$dir/$file" : "$dir/$file.$ext";
            if (-e $obj) {
                print "Deleting $obj\n" if defined $options{v};
                unlink "$obj";
            }
        }
    }

}

sub clean_all {

    remove(".", ["suite"], [qw(log aux dvi ps pdf tex)]);
    remove(".", ["numbers_logP"], [qw(txt tex eps)]);
    remove(".", ["suite_all"], ["pdf"]);
    remove("Demo", [qw(xi_p xi_m D_Lum P_delta)], [qw(NIX  gnu)]);
    remove("Demo", [qw(fig_xi_p fig_xi_m fig_D_Lum fig_P_delta)], [qw(eps)]);
    remove("Demo", [qw(numbers_0_Pdelta numbers_1_xip numbers_4_D_Lum)], [qw(tex tab txt eps)]);
    remove("Demo", [qw(w_eos P_kappa mapsqr mapsqr_gauss gammasqr nz D_plus D_plus2)], ["NIX"]);
}

sub usage {
  my ($ex) = @_;

  print STDERR "Usage: test_suite_cosmo_pmc.pl [OPTIONS]\n";
  print STDERR "OPTIONS:\n";
  print STDERR "  -r             Do PMC test runs\n";
  print STDERR "  -R             Only do PMC test runs\n";
  print STDERR "  -n NCPU        Run PMC in parallel on NCPU cpus (default: 1)\n";
  print STDERR "  -c             Include CMB tests\n";
  print STDERR "  -P PATH        Use PATH as CosmoPMC root directory (default: environment\n";
  print STDERR "                  variable \$COSMOPMC)\n";
  print STDERR "  -s             Short, without time-taking PMC runs (e.g. Lensing/COSMOS-S10)\n";
  print STDERR "  -k             Keep temporary files\n";
  print STDERR "  -x             Clean previous run and exit\n";
  print STDERR "  -v             Verbose\n";
  print STDERR "  -h             This message\n";

  exit $ex if defined $ex;
}
