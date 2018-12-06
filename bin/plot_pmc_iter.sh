#!/bin/csh

# MK 5/2008
# Creates all kinds of plots from a PMC simulation (iteration iter)
# Usage: plot_pmc.sh

set version = 1.0

set pwd   = `pwd`
set title = " "

cd ..
set xxx = `pwd`
set TITLE = `basename $xxx`
cd $pwd


### Tasks ###
set do_sample_comp  = 0
set do_1d2dlike     = 0
set do_1d2dlike_yor = 1
set do_proposal     = 0     # Direction (covar)

### Command line and default parameters ###
set quiet = 0
set flag_quiet = ""
set config = "../config_pmc"
@ next = 0
while ($#argv)
    if ($next == 1) set config = $argv[1]

    if ("$argv[1]" == "-q") then
	set quiet = 1
	set flag_quiet = "-q"
    endif
    if ("$argv[1]" ==  "-c") @ next = 1

    shift
end


###############
### Startup ###
###############

if (! $quiet) echo "plot_pmc_iter.sh v$version $config"



# Parameters from config file
if (! -e $config) then
    echo "Erorr: config file '$config' not found, exiting."
    echo "       Use option -c CONFIG to change default name."
    exit 1
endif
if (! $quiet) echo Reading config file $config

set npar    = `cat $config | perl -ane 'print $F[1] if /^npar/'`
set min     = (`cat $config | perl -ane 'print if /^min/'`)
set max     = (`cat $config | perl -ane 'print if /^max/'`)
set ncomp   = `cat $config | perl -ane 'print $F[1] if /^ncomp/'`

set iter    = `perl -e '$ARGV[0] =~ /iter\_(\d?)$/; print "$1\n"' $pwd`

if (! $quiet) echo "Iteration $iter"


#####################################
### Sample points from components ###
#####################################

if ($do_sample_comp == 1) then

    if (! -e iter_$iter/sample) then
	if (! $quiet) echo "Warning: sample file iter_$iter/sample not found, skipping sample point plot"
    else

	if (! $quiet) echo  -n "Sample points from components"
	echo unset logs > samplec.gnu
	echo unset grid >> samplec.gnu
	echo set size squ >> samplec.gnu
	echo set term post eps enh color 20 >> samplec.gnu
	echo unset key >> samplec.gnu
	echo set pointsize 0.01 >> samplec.gnu
	echo ev = 5 >> samplec.gnu

	@ i = 0
	while ($i < $npar - 1)
	    @ j = $i + 1
		while ($j < $npar)
		
		    if (! $quiet) echo -n "."

		    @ limi = $i + 2
		    @ limj = $j + 2
		    @ coli = $i + 3
		    @ colj = $j + 3

		    set xlabel = "`$COSMOPMC/bin/get_spar.pl -c $config -i $i gnuplot`"
		    set ylabel = "`$COSMOPMC/bin/get_spar.pl -c $config -i $j gnuplot`"

		    echo set xlabel \""$xlabel"\" >> samplec.gnu
		    echo set ylabel \""$ylabel"\" >> samplec.gnu
		    echo set output \'samplec_{$i}_{$j}.ps\' >> samplec.gnu
		    echo -n "pl [$min[$limi] : $max[$limi]][$min[$limj] : $max[$limj]] " >> samplec.gnu

		    @ c = 0
		    while ($c < $ncomp)
			echo -n \'sample_{$iter}\' u $coli : \(\$2==-$c \? \$$colj : 1/0\) ev ev lt $c ps 0.5 >> samplec.gnu
			if ($c < $ncomp - 1) then
			    echo -n ', ' >> samplec.gnu
			else
			    echo >> samplec.gnu
			endif
			@ c ++
		    end
		    echo set output >> samplec.gnu
		@ j ++
		end
	    @ i ++
	end

	if (! $quiet) echo -n g
	gnuplot samplec.gnu


	if (! $quiet) echo -n l
	$COSMOPMC/bin/all_vs_all.pl -b samplec -e ps -t $pwd > all_samplec.tex
	$COSMOPMC/bin/ldp.sh all_samplec.tex -q
	rm -f all_samplec.{log,dvi,aux}
	if (! $quiet) echo

    endif

endif



#set x = $<
#################################
### 1d, 2d marginal plots (R) ###
#################################

if ($do_1d2dlike == 1) then

    if (! $quiet) echo "1d and 2d likelihood plots"
    $COSMOPMC/bin/plot_confidence.sh iter_$iter/sample

    $COSMOPMC/bin/all_vs_all.pl -b cont2d -e ps -l like1d -t $pwd > all_cont2d.tex
    $COSMOPMC/bin/ldp.sh all_cont2d.tex -q
    rm -f all_cont2d.{log,dvi,aux}


endif

if ($do_1d2dlike_yor == 1) then

    if (! $quiet) echo "1d and 2d likelihood plots (yorick)"
    #if (-e covar+ded.fin) then
       #set command = "$COSMOPMC/bin/plot_contour2d.pl -c $config -1 m1t -g -10 -C -n -T "'"'"$TITLE"'"'" $flag_quiet"
	echo "No smoothing of the contours!"
    #else
        set command = "$COSMOPMC/bin/plot_contour2d.pl -c $config -1 m1t -o pdf -g -25 -n -T "'"'"$TITLE"'"'" $flag_quiet"
    #endif

    #echo $command > plot.sh
    `echo $command`

endif


#################
### Proposals ###
#################

if ($do_proposal == 1) then

    if (-e proposal) then

       if (! $quiet) echo "2d plots of the proposal"
       echo include, \"$COSMOPMC/yorick/likeli.i\" > prop2d.i
       echo include, \"$COSMOPMC/yorick/stuff.i\" >> prop2d.i

       @ i = 0
       while ($i < $npar)
	   @ j = $i + 1
	   while ($j < $npar)
	       @ coli = $i + 2
	       @ colj = $j + 2
	       @ i1   = $i + 1
	       @ j1   = $j + 1

	       echo window, display=\"\", hcp=\"prop2d_{$i}_{$j}\" >> prop2d.i
	       echo mix_mvd = mix_mvdens_read\(\"proposal\"\) >> prop2d.i
	       echo plot_mix_mvdens_covar, mix_mvd, [$i1, $j1] >> prop2d.i;
               set xlabel = "`$COSMOPMC/bin/get_spar.pl -c $config -i $i yorick`"
               set ylabel = "`$COSMOPMC/bin/get_spar.pl -c $config -i $j yorick`"


#	       echo xytitles, \""$xlabel"\", \""$ylabel"\", [-0.01, 0.015] >> prop2d.i
	       echo "xytitles, " \"$xlabel\" ", " \"$ylabel\" ", [-0.01, 0.015]"  >> prop2d.i


	       echo pltitle, \"$title\" >> prop2d.i
	       echo limits, $min[$coli], $max[$coli], $min[$colj], $max[$colj] >> prop2d.i
	       echo hcps, \"prop2d_{$i}_{$j}\" >> prop2d.i

	       @ j ++
	   end
       @ i ++
       end

       echo quit >> prop2d.i
       yorick -batch prop2d.i

       $COSMOPMC/bin/all_vs_all.pl -b prop2d -e ps -t $pwd > all_prop2d.tex
       $COSMOPMC/bin/ldp.sh all_prop2d.tex -q
       rm -f all_prop2d.{log,dvi,aux}

    else

       if (! $quiet) echo "proposal does not exist, skipping plot_proposal"

    endif

endif

# end
