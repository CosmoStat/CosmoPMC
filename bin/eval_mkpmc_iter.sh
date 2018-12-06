#!/bin/tcsh

# MK 5/2008

set version = 0.2


set tsleep = 1

# Config file and info
if (-e config_pmc) then
    set config  = config_pmc
else
    echo "No config file found, exiting."
    exit(2)
endif

set niter = `cat $config | perl -ane 'print $F[1] if /niter/'`

@ iter = 0
while ($iter < $niter)

    set found = 0
    if (-e perplexity) then
       echo "perplexity :"
       tail perplexity
    endif
    @ ii = $iter + 1
    if (-e proposal_$ii) then
	echo "neff proposal_$ii :"
	neff_proposal.pl proposal_$ii
    endif
    echo "Waiting for iteration $iter..."
    while ($found == 0)
	if (-e iter_$iter) then
	    echo "mean_$iter :"
 	    cat mean_$iter
	    set found = 1
	endif
	sleep $tsleep
    end

    if (-e "iter_$iter/all_contour2d.ps") then
       echo "Skipping iteration $iter"
    else
       echo "Evaluating iteration $iter"
       cd iter_$iter
       plot_pmc_iter.sh -q
       cd ..
    endif

   @ iter++

end

plot_pmc_final.sh
proposal_mean.pl

echo "End eval_mkpmc_iter.sh"

# end
