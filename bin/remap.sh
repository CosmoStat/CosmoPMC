#!/bin/csh

# remap.sh
# Martin Kilbinger 2008-2010
# Remaps the 1d and 2d histograms (files chi2_*) to another 
# parameter set.
# Also remaps mean and covariance

# Default parameters
set config      = "config_pmc"
set pmcsim      = "pmcsim"
set iter        = -1
set dir_in      = "."
set dir_out     = "remap"
set remap_file  = "remap.dat"
set npar        =
set n_ded       =

# Flags
set flag_config = 0
set flag_pmcsim  = 0
set flag_dir_in = 0

while ($#argv > 0)

    if ("$argv[1]" == "-c") then
	shift
	set config = $argv[1]
	set flag_config = 1
    endif
    if ("$argv[1]" == "-s") then
	shift
	set pmcsim = $argv[1]
	set flag_pmcsim = 1
    endif
    if ("$argv[1]" == "-i") then
	shift
        set dir_in = $argv[1]
	set flag_dir_in = 1
    endif
    if ("$argv[1]" == "-o") then
	shift
        set dir_out = $argv[1]
    endif
    if ("$argv[1]" == "-r") then
	shift
	set remap_file = $argv[1]
    endif
    if ("$argv[1]" == "-n") then
	shift
	set npar = $argv[1]
    endif
    if ("$argv[1]" == "-d") then
	shift
	set n_ded = $argv[1]
    endif
    if ("$argv[1]" == "-h") goto usage

    shift
end

if ("$dir_in" == ".") set dir_in = `pwd`

# Get absolute path of input directory
pushd $dir_in > /dev/null
#set hist_dir = $dir_in
set hist_dir = `pwd`
popd > /dev/null

if (! -s $remap_file) then
    echo "Remap file '$remap_file' does not exist or is empty, stopping script"
    exit 3
endif

set orig = (`head -n 1 $remap_file`)

if (! -e $dir_out) then
    echo "Creating directory for remapped run '$dir_out'"
    mkdir -p $dir_out
endif

if (-e $dir_out/chi2_0) rm -f $dir_out/chi2_*


@ i = 1
while ($i <= $#orig)
   #echo $i
   @ remapi = $i - 1
   ln -s $hist_dir/chi2_{$orig[$i]} $dir_out/chi2_{$remapi}

   @ j = $i + 1
   while ($j <= $#orig)
      @ remapj = $j - 1
      if ($orig[$i]>$orig[$j]) then
         echo "Swap $remapi $remapj"
	 # Swap columns and (re-)sort file
         cat $hist_dir/chi2_{$orig[$j]}_{$orig[$i]} | perl -ane 'print "$F[1] $F[0] $F[2]\n"' > $dir_out/tmp
	 sort -g $dir_out/tmp > $dir_out/chi2_{$remapi}_{$remapj}
      else
         echo "Link $remapi $remapj"
         ln -s $hist_dir/chi2_{$orig[$i]}_{$orig[$j]} $dir_out/chi2_{$remapi}_{$remapj}
      endif
      @ j ++
   end
   @ i ++
end

rm -f tmp

# Mean file
if (-e $dir_in/mean) then
   rm -f $dir_out/mean
   # Header line
   head -n 1 $dir_in/mean >> $dir_out/mean
   @ i = 1
   while ($i <= $#orig)
      @ line = $orig[$i] + 2
      head -n $line $dir_in/mean | tail -n 1 >> $dir_out/mean 
      @ i ++
   end
endif

if (-e $dir_in/covar+ded.fin) then
    remap_cov.pl $remap_file $dir_in/covar+ded.fin > $dir_out/covar+ded.fin
endif
if (-e $dir_in/covarinv+ded.fin) then
    remap_cov.pl $remap_file $dir_in/covarinv+ded.fin > $dir_out/covarinv_ded.fin
endif

remap_config.pl $remap_file $config $npar $n_ded > $dir_out/config_pmc

if ("$flag_pmcsim" == 1) then
    remap_sample.pl $remap_file $pmcsim > $dir_out/pmcsim
endif

exit 0


usage:
echo "Usage: remap.sh [OPTIONS]"
echo "OPTIONS:"
echo "   -c CONFIG            Input PMC configuration file (default './config_pmc')"
echo "   -i INPUT             Input directory INPUT (default '.')"
echo "   -s PMCSIM            Sample/PMC simulation file PMCSIM"
echo "   -o OUTPUT            Output directory OUTPUT (default './remap')"
echo "   -r REMAP             Remap file REMAP (default './remap.dat')"
echo "   -n NPAR              Number of parameters NPAR (default: read from remap file)"
echo "   -d N_DED             Number of deduced parameters N_DED (default: 0)"
echo "   -h                   This message"
exit 1
