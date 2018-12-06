#!/bin/csh

# Creates subdirectories for a MCMC run, similar to iter_? for PMC, each subdir iter_i corresponding to
# a subchain of the first i*nsamples points.


# Type of experiment (to copy the right links into subdirs)
set lens  =  1
set Nz    =  2
set Sn    =  4
set WMAP3 =  8
set WMAP5 = 16

set TYPE = $WMAP5

if ("$1" == "") then
    echo "Data types: Lens  $lens"
    echo "            Nz    $Nz"
    echo "            Sn    $Sn"
    echo "            WMAP3 $WMAP3"
    echo "            WMAP5 $WMAP5"

    echo -n "Data type? [$TYPE] "
    set type = $<
    if ("$type" == "") then
	set type = $TYPE
    endif
endif

# Number of sample points
set NSAMPLES = 10000
echo -n "number of sample points per step? [$NSAMPLES] "
set nsamples = $<
if ("$nsamples" == "") then
	set nsamples = $NSAMPLES
endif

set res = (`cut_mkmc.pl $nsamples`)

# Create subdirs and run mkmc
foreach i (iter_*)

    newdir_mkmc.sh $i <<EOF
y
$type
EOF

    cd $i
    ln -fs ../config_mcmc

    $mcsrc/mkmc

    cd ..

    set j = `perl -e '$ARGV[0] =~ s/iter_//; print "$ARGV[0]\n"' $i`
    echo "ln -sf $i/mean mean_$j"
    ln -sf $i/mean mean_$j

end

R < $COSMOPMC/R/plot_mean_var_iter.R $R_OPTIONS $res[1] $res[2] $res[3] 

# This is from plot_pmc_final.sh ... 
mkdir -p mean_var_iter
rm -f mean_var_iter/*.ps
mv mean_var_iter_*.ps mean_var_iter
cd mean_var_iter
allps2tex.pl > all_mean_var_iter.tex
ldp.sh all_mean_var_iter.tex -q
rm -f all_mean_var_iter.{log,dvi,aux}
cd ..

# end
