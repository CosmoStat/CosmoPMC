#!/bin/csh

# plot_cosmo_mcmc.sh
# Martin Kilbinger 2007, 2008
# Creates all kinds of plots from an MCMC chain (file chain.fin)

set version = 1.4

set pwd     = `pwd`
#set title = `basename $pwd`
set title = " "


### Tasks ###
set do_1dtime            = 0
set do_2dcloud           = 0
set do_1d2dlike_yor      = 1
set do_1d2dlike          = 0
set do_covar             = 0     # 1: covar.fin, 2: fisher
set do_proposal          = 0     # 1: Weights (mean), 2: Direction (covar)
set plot_proposal        = 1
set do_1dcorr_calc       = 0
set do_1dcorr_plot       = 0
set do_sample_histograms = 1
set do_variance_mcmc     = 1

set smooth_fac_2d = 30
set smooth_fac_1d = $smooth_fac_2d

set R_OPTIONS = "--slave --vanilla --args"

###############
### Startup ###
###############

echo "plot_cosmo_mcmc.sh v$version"


# parameters from config file
if (-e config_pmc) then
    set config  = config_pmc
    set pmc     = 1
else if (-e config_mcmc) then
    set config = config_mcmc
    set pmc    = 0
else if (-e config) then
    set config = config
    set pmc    = 1
else
    echo "No config file found, exiting."
    exit(2)
endif
echo Reading config file $config

set npar    = `cat $config | perl -ane 'print $F[1] if /npar/'`
set min     = (`cat $config | perl -ane 'print if /min/'`)
set max     = (`cat $config | perl -ane 'print if /max/'`)

if ($pmc == 0) then
    set nchain  = `cat $config | perl -ane 'print $F[1] if /nchain/'`
    set ncov    = `cat $config | perl -ane 'print $F[1] if /ncov/'`
    set fburnin = `cat $config | perl -ane 'print $F[1] if /fburnin/'`
    set ndecorr = `cat $config | perl -ane 'print $F[1] if /ndecorr/'`
    set niter   = -1
else
    set niter   = `cat $config | perl -ane 'print $F[1] if /niter/'`
    set ncomp   = `cat $config | perl -ane 'print $F[1] if /ncomp/'`
endif

# Some parameter combinations don't make sense
if ($pmc == 0) then
    set do_proposal = 0
    set plot_proposal = 0
    set do_sample_histograms = 0
endif
if ($pmc == 1) then
    set do_1dtime = 0
endif

# chain or sample file
if (-e chain.fin) then
    set sample = chain.fin
else
    @ nfinal = $niter - 1
    if (-e iter_$nfinal/sample_$nfinal) then
	set sample = iter_nfinal/sample_$nfinal
    else
        echo "No sample file found, exiting."
        exit(7)
    endif
endif


############################
### 1d time series plots ###
############################

if ($do_1dtime == 1 && -e chain.acc) then 

    echo "1d time series"
    echo unset logs > param1d.gnu
    echo unset grid >> param1d.gnu
    echo set size nosquare 1.5, 0.5 >> param1d.gnu
    unset grid >> param1d.gnu
    echo set xlabel \'i\' >> param1d.gnu
    echo set term post eps enhanced color 16 >> param1d.gnu
    @ i = 0
    while ($i < $npar)
	@ limi = $i + 2

        set ylabel = "`get_spar.pl -c $config -i $i gnuplot`"
	echo set ylabel \""$ylabel"\" >> param1d.gnu
	echo set title \"$title\" >> param1d.gnu
	echo set pointsize 0.5 >> param1d.gnu
	@ col = $i + 3
	set ii = `perl -e 'printf "%02d\n", '$i''`

	echo set output \'param1d_{$ii}.ps\' >> param1d.gnu
	echo pl "[][$min[$limi] : $max[$limi]]" \'chain.acc\' u 0:$col not w l >> param1d.gnu
	echo set output >> param1d.gnu

	set nfinal = `perl -e '$x = ($ARGV[0]-$ARGV[1]*$ARGV[2])/$ARGV[3]; print int($x), "\n"' $nchain $ncov $fburnin $ndecorr`
	@ start = $nfinal - 1000
	echo set output \'param1d_1k_{$ii}.ps\' >> param1d.gnu
	echo "pl [$start : $nfinal][$min[$limi] : $max[$limi]]" \'chain.acc\' u 0:$col not w l >> param1d.gnu
	echo set output >> param1d.gnu

	@ i++
    end

    gnuplot param1d.gnu

    mkdir tmptmp
    cp param1d*.ps tmptmp
    cd tmptmp
    allps2tex.pl 16 > all_param1d.tex
    ldp.sh all_param1d.tex -q
    mv all_param1d.{tex,ps} ..
    cd ..
    rm -rf tmptmp

endif


######################
### 2d cloud plots ###
######################

if ($do_2dcloud == 1) then

    echo "2d cloud plots"
    echo unset logs > param2d.gnu
    echo unset grid >> param2d.gnu
    echo set size squ >> param2d.gnu
    echo set term post eps enh color 16 >> param2d.gnu
    @ i = 0
    while ($i < $npar)
       @ j = $i + 1
       while ($j < $npar)
	   @ limi = $i + 2
	   @ limj = $j + 2
	   @ coli = $i + 3
	   @ colj = $j + 3

           set xlabel = "`get_spar.pl -c $config -i $i gnuplot`"
           set ylabel = "`get_spar.pl -c $config -i $j gnuplot`"

	   echo set xlabel \""$xlabel"\" >> param2d.gnu
	   echo set ylabel \""$ylabel"\" >> param2d.gnu
	   echo set title \"$title\" >> param2d.gnu
	   echo set output \'param2d_{$i}_{$j}.ps\' >> param2d.gnu
	   echo "pl [$min[$limi] : $max[$limi]][$min[$limj] : $max[$limj]]" \'$sample\' u $coli\:$colj not >> param2d.gnu
	   echo set output >> param2d.gnu
	   @ j++
       end
       @ i ++
    end

    gnuplot param2d.gnu

endif


#################################
### 2d contour plots (yorick) ###
#################################

if ($do_1d2dlike_yor == 1) then

    echo "2d contour plots (calling plot_contour2d.pl)"
    if ("$title" == " ") then
	set command = "plot_contour2d.pl -f $do_covar -p $do_proposal -1 m1t -g $smooth_fac_2d -G $smooth_fac_1d ."
    else
	set command = "plot_contour2d.pl -f $do_covar -p $do_proposal -t \"$title\" -1 -g $smooth_fac_2d -G $smooth_fac_1d ."
    endif
    echo $command > plot.sh
    `echo $command`

endif

after_contour2d:


#################################
### 1d, 2d marginal plots (R) ###
#################################

if ($do_1d2dlike == 1) then

    echo "1d and 2d likelihood plots"
    R < $COSMOPMC/R/plot_confidence.R $R_OPTIONS $sample
    all_vs_all.pl -b cont2d -e eps -l like1d -t $pwd > all_cont2d.tex
    ldp.sh all_cont2d.tex -q
    rm -f all_cont2d.{log,dvi,aux}

endif


#################
### Proposals ###
#################

if ($do_proposal & 2) then

    echo include, \"likeli.i\" > prop2d.i
    echo include, \"stuff.i\" >> prop2d.i

    @ i = 0
    while ($i < $npar)
	@ j = $i + 1
	while ($j < $npar)
	    @ coli = $i + 2
	    @ colj = $j + 2
	    @ i1   = $i + 1
	    @ j1   = $j + 1

	    @ p = 0
	    while ($p <= $niter)
		echo window, display=\"\", hcp=\"prop2d_{$p}__{$i}_{$j}\" >> prop2d.i
		echo mix_mvd = read_all_proposals\($p\) >> prop2d.i
		echo plot_mix_mvdens_covar, mix_mvd\($p\), [$i1, $j1] >> prop2d.i;
                set xlabel = "`get_spar.pl -c $config -i $i yorick`"
                set ylabel = "`get_spar.pl -c $config -i $j yorick`"
		echo xytitles, \""$xlabel"\", \""$ylabel"\", [-0.01, 0.015] >> prop2d.i
		echo pltitle, \"$title\" >> prop2d.i
		echo limits, $min[$coli], $max[$coli], $min[$colj], $max[$colj] >> prop2d.i
		echo hcps, \"prop2d_{$p}__{$i}_{$j}\" >> prop2d.i
		@ p ++
	    end

	    @ j ++
	end
    @ i ++
    end

    echo quit >> prop2d.i
    yorick -batch prop2d.i

    @ p = 0
    while ($p <= $niter)
	all_vs_all.pl -b prop2d_{$p}_ -e ps -l likeli1d -t $pwd> all_prop2d_{$p}.tex
	ldp.sh all_prop2d_{$p}.tex -q
	rm -f all_prop2d{$p}.{log,dvi,aux}
	@ p ++
    end

endif


############################
### PMC proposal weights ###
############################

if ($do_proposal & 1) then
    
    echo "PMC proposal plots"
    echo include, \"likeli.i\" > proposal.i
    echo window, display=\"\", hcp=\"prop_weights\" >> proposal.i
    echo mix_mvd = read_all_proposals\($niter\) >> proposal.i
    echo plot_prop_weights_iter, mix_mvd >> proposal.i
    echo hcps, \"prop_weights\" >> proposal.i
    echo quit >> proposal.i
    yorick -batch proposal.i

endif


#########################
### Weight histograms ###
#########################

if ($do_sample_histograms == 1) then

    echo "Sample histograms"
    @ i = 0
    while ($i < $niter)
	histogram.pl -N 100 -n -L -q pmcsim_$i 0 > pmcsim_$i.h 
	histogram.pl -L -N 20 -n -q sample_$i 0 > sample_$i.h 
	@ i ++ 
    end

    echo set nogrid > weight_hist.gnu
    echo set logs y >> weight_hist.gnu
    echo set term post eps enhanced color 16 >> weight_hist.gnu
    echo set output \'weight_hist.ps\' >> weight_hist.gnu
    echo set xlabel \'log\(weight\)\' >> weight_hist.gnu
    echo set ylabel \'frequency\' >> weight_hist.gnu
    echo set title \'PMC weights histograms\' >> weight_hist.gnu
    echo -n "pl [-100:] " >> weight_hist.gnu
    @ i = 0
    while ($i < $niter - 1)
	echo -n \'pmcsim_$i.h\' t \'iteration $i\' w l, >> weight_hist.gnu
	@ i ++
    end
    echo \'pmcsim_$i.h\' t \'iteration $i\' w l >> weight_hist.gnu

    echo set output \'sample_weight_hist.ps\' >> weight_hist.gnu
    echo set xlabel \'weight\' >> weight_hist.gnu
    echo set ylabel \'frequency\' >> weight_hist.gnu
    echo set title \'Sample weights histograms\' >> weight_hist.gnu
    echo -n "pl " >> weight_hist.gnu
    @ i = 0
    while ($i < $niter - 1)
	echo -n \'sample_$i.h\' t \'iteration $i\' w histeps, >> weight_hist.gnu
	@ i ++
    end
    echo \'sample_$i.h\' t \'iteration $i\' w histeps >> weight_hist.gnu
    echo set output >> weight_hist.gnu

    gnuplot weight_hist.gnu

endif


###############################
### 1d temporal correlation ###
###############################

if ($do_1dcorr_calc) then

    echo "1d temporal correlation"
    set imax = 300
    corr_temp.pl -m $imax $sample > corr

endif

if ($do_1dcorr_plot) then

    set imax = 300
    echo unset logs > corr1d.gnu
    echo unset grid >> corr1d.gnu
    echo set size 1.5, 0.5 >> corr1d.gnu
    echo set xlabel \'t\' >> corr1d.gnu
    echo set ylabel \'correlation\(t\)\' >> corr1d.gnu
    echo set title \"$title\" >> corr1d.gnu
    echo set term post eps enhanced color 16 >> corr1d.gnu
    echo set output \'corr1d.ps\' >> corr1d.gnu
    echo set multiplot >> corr1d.gnu
    set last = (`tail -n 1 corr`)
    @ xkey = $last[1] - 30
    set dy = `perl -e 'print 0.5/$ARGV[0], "\n"' $npar`

    @ i = 0
    while ($i < $npar)
       @ coli = $i + 2
       @ ip1  = $i + 1
       echo "ykey = 1 - ($i+1)*$dy" >> corr1d.gnu
       echo set key $xkey, ykey >> corr1d.gnu
       set title = "`get_spar.pl -c $config -i $i gnuplot`"
       echo pl [0 : $imax ][-0.2 : 1] \'corr\' u 1: $coli w l lt $ip1 t \""$title"\" >> corr1d.gnu
       @ i++
    end

    echo unset multiplot >> corr1d.gnu
    echo set output >> corr1d.gnu

    gnuplot corr1d.gnu

endif


#################################
### MCMC variance of the mean ###
#################################

if ($do_variance_mcmc) then

    echo "MCMC variance"
    R < $COSMOPMC/R/variance_mcmc.R $R_OPTIONS $sample > variance_mcmc

endif


echo "plot_cosmo_mcmc.sh finished"

# eof
