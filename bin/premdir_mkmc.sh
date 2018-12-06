#!/bin/csh

# Creates a subdirectory ("prem" by default) of a MCMC run, copies chain.acc and modifies
# nchain in config_mcmc to the (premature) length of the chain.
# Run e.g. plot_mkmc.sh in the subdir.

if ("$1" == "") then
    set dir = "prem"
else
    set dir = "$1"
endif

newdir_mkmc.sh $dir

cd $dir
cp ../chain.acc .
set nchain = `cat chain.acc | perl -ane 'print if ! /#/' | wc -l`
echo $nchain
cat ../config_mcmc | perl -ane 's/nchain(\s*)(.*)/nchain\t\t'$nchain'/; print' > config_mcmc

###
