#! /bin/sh

# -------------------------------------------------------------- #
# Example how to use cosmo_pmc to fit w(theta) with a HOD model 
# data are taken from Coupon et al. (2012) A&A, 542A, 5C
# -------------------------------------------------------------- #

# -------------------------------------------------------------- #
# Set the number of cpus 
# -------------------------------------------------------------- #

np=16

# -------------------------------------------------------------- #
# to clean up the directory: runme.sh clean
# -------------------------------------------------------------- #

if [ "$1" == "clean" ]; then
    echo "Cleaning directory..."
    rm -f maxlogP
    rm -f log_fish
    rm -f fisher
    rm -f config_fish
    rm -f rng_init_pmc
    rm -f proposal_0
    rm -f perplexity
    rm -f evidence_fisher
    rm -f evidence
    rm -f log_pmc
    rm -f enc
    rm -f cov_inv.dat
    rm -f config_max
    rm -f log_max_post
    #After pmc
    rm -rf iter_*
    rm -rf proposal_*
    rm -rf proposal_means
    rm -f  pmcsim*
    rm -f  config_pmc_ded
    rm -f  pmcsim_*
    rm -f  perplexity*
    rm -rf deduced    
    rm -f  diagnostics.*
    rm -f  essential.pdf
    rm -f  wtheta.model
    exit
fi


# -------------------------------------------------------------- #
# To compute the model
# runme.sh model
# -------------------------------------------------------------- #

if [ "$1" == "model" ]; then
    getHODModel halomodel.par -t wtheta -para 12.1795,13.3016,12.1779,0.328833,1.14355,1.0  -o wtheta.model
    exit
fi

# -------------------------------------------------------------- #
# MAX_POST: local maximum of likelihood
# -------------------------------------------------------------- #

# galaxy number density. For starting point use only 
Ngd=.0035188360

logMmin=`echo "scale=10; 10.0-0.85*l($Ngd)/l(10.0) " | bc -l`
logM1=`  echo "scale=10; $logMmin + 1.23  " | bc -l` 
logM0=`  echo "scale=10; $logMmin         " | bc -l` 

sstart=`echo $logMmin $logM1 $logM0 0.3 1.0`

config_pmc_to_max_and_fish.pl -M  -f "$sstart" > config_max
max_post

exit

# -------------------------------------------------------------- #
# fisher matrix (-d: diagonal fisher matrix)
# -------------------------------------------------------------- #

config_pmc_to_max_and_fish.pl -F -d -p maxlogP config_pmc > config_fish
mpirun -np $np go_fishing

# -------------------------------------------------------------- #
# cosmo_pmc
# -------------------------------------------------------------- #

cosmo_pmc.pl -n $np <<EOF
n
EOF
