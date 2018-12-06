#!/bin/csh

# Martin Kilbinger 2008
# Creates all kinds of plots from a PMC simulation (iteration iter)

set version = 0.9

# For R
set Roptions = "-q  --slave --vanilla"
set Rpath    = $COSMOPMC/R/

### Command line and default parameters ###
set config = "config_pmc"
@ next = 0
while ($#argv)
    if ($next == 1) set config = $argv[1]

    if ("$argv[1]" ==  "-c") @ next = 1

    shift
end


###############
### Startup ###
###############

echo "plot_pmc_final.sh v$version"

# Parameters from config file
if (! -e $config) then
    echo "Erorr: config file '$config' not found, exiting."
    echo "       Use option -c CONFIG to change default name."
    exit 1
endif
echo Reading config file $config

set npar    = `cat $config | perl -ane 'print $F[1] if /npar/'`
set min     = (`cat $config | perl -ane 'print if /min/'`)
set max     = (`cat $config | perl -ane 'print if /max/'`)
set niter   = `cat $config | perl -ane 'print $F[1] if /niter/'`
set ncomp   = `cat $config | perl -ane 'print $F[1] if /ncomp/'`



############################
### PMC proposal weights ###
############################

echo "PMC proposal plots"
echo include, \"likeli.i\" > proposal.i
echo window, display=\"\", hcp=\"prop_weights\" >> proposal.i
echo mix_mvd = read_all_proposals\($niter\) >> proposal.i
echo plot_prop_weights_iter, mix_mvd >> proposal.i
echo hcps, \"prop_weights\" >> proposal.i
echo quit >> proposal.i
yorick -batch proposal.i


#########################
### Weight histograms ###
#########################

echo "Sample histograms"
@ i = 0
while ($i < $niter)
   histogram.pl -N 50 -L -n -q pmcsim_$i 0 > pmcsim_$i.hist
   if (-e iter_$i/sample_$i) then
      histogram.pl -L -N 20 -n -q iter_$i/sample_$i 0 > iter_$i/sample_$i.hist
   endif
   @ i ++ 
end

echo set nogrid > weight_hist.gnu
echo set logs y >> weight_hist.gnu
echo set term post eps enhanced color 24 >> weight_hist.gnu
echo set output \'weight_hist.ps\' >> weight_hist.gnu
echo set xlabel \'log\(weight\)\' >> weight_hist.gnu
echo set ylabel \'frequency\' >> weight_hist.gnu
echo set title \'PMC weights histograms\' >> weight_hist.gnu
echo set key left >> weight_hist.gnu
#echo -n "pl [-40:][0.001:] " >> weight_hist.gnu
echo -n "pl " >> weight_hist.gnu
@ i = 0
while ($i < $niter)
   echo -n \'pmcsim_$i.hist\' t \'iteration $i\' w l lw 3, >> weight_hist.gnu
   @ i += 3
end
echo 0 not >> weight_hist.gnu

echo set output \'sample_weight_hist.ps\' >> weight_hist.gnu
echo set xlabel \'weight\' >> weight_hist.gnu
echo set ylabel \'frequency\' >> weight_hist.gnu
echo set title \'Sample weights histograms\' >> weight_hist.gnu
echo -n "pl " >> weight_hist.gnu
@ i = 0
@ nsample = 0
while ($i < $niter - 1)
   empty.pl iter_$i/sample_$i.hist
   if ("$?" == "0") then
      echo -n \'iter_$i/sample_$i.hist\' t \'iteration $i\' w histeps, >> weight_hist.gnu
      @ nsample ++
   endif
   @ i ++
end
empty.pl iter_$i/sample_$i.hist
if ("$?" == "0") then
   echo \'iter_$i/sample_$i.hist\' t \'iteration $i\' w histeps >> weight_hist.gnu
   @ nsample ++
endif
if ($nsample == 0) then
   echo 1 >> weight_hist.gnu
endif
echo set output >> weight_hist.gnu

gnuplot weight_hist.gnu 


###########################################
### Mean, variance as fct. of iteration ###
###########################################

R $Roptions < $Rpath/plot_mean_var_iter.R

if ("$?" == "0") then
    mkdir -p mean_var_iter
    rm -f mean_var_iter/*.ps
    mv mean_var_iter_*.ps mean_var_iter
    cd mean_var_iter
    allps2tex.pl > all_mean_var_iter.tex
    ldp.sh all_mean_var_iter.tex -q
    rm -f all_mean_var_iter.{log,dvi,aux}
    cd ..
else
    echo Call to R with $Rpath/plot_mean_var_iter.R failed, skipping mean_var_iter
endif

################################
### PMC variance of the mean ###
################################

echo "PMC variance"
@ i = $niter - 1
R < $COSMOPMC/R/variance_pmc.R $Roptions  --args pmcsim_$i > variance_pmc

# end
