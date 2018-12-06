#!/bin/csh

# Martin Kilbinger 2008-2010
# Runs cosmo_mcmc, then plot_mcmc.sh.
# Deletes existing chains and Fisher matrix after confirmation.
# The following files are either read or recalculated by cosmo_mcmc:
#
# File      	Calculation step(s) if not existant
# -------------------------------------------------------
# fisher    	Maximum search and Fisher matrix
# chain.pre     - (If chain.pre exists, it is read and a
#               the MCMC run continued        
# chain.acc 	Runs MCM Chain with accepted points
# chain.fin 	From accepted chain removes burn-in phase
#		and decorrelate (thin-out). Generates
#		histograms
# -------------------------------------------------------
#
# For example, if after a first run of cosmo_mcmc a different burn-in
# phase is to be tried out, do the corresponding change in the 
# config file, start cosmo_mcmc.sh and confirm to delete chain.fin but
# not chain.acc and fisher.


# Check for existing files, delete if confirmed

if (-e chain.fin) then
    echo -n "Delete chain.fin? [n] " 
    set del = $<
    if ("$del" == "y") then
	rm -f chain.fin
    endif
endif

if (-e chain.acc) then
    echo -n "Delete chain.acc? [n] " 
    set del = $<
    if ("$del" == "y") then
	rm -f chain.acc
    endif
endif

if (-e fisher) then
    echo -n "Delete fisher? [n] " 
    set del = $<
    if ("$del" == "y") then
	rm -f fisher
    endif
endif

# Save old log file before overwriting
if (-e log_mcmc) then
    set name = `find . -maxdepth 1 -name log_mcmc -printf "%CY-%Cm-%Cd-%CT\n"`
    mv log_mcmc log_mcmc-$name
endif

# Copy command line arguments to string for cosmo_mcmc 
set str = ""
while (1 <= $#argv)
    set str = "$str $argv[1]"
    shift
end

# Run MCMC
$COSMOPMC/exec/cosmo_mcmc $str

# Preprocessing and plotting
if ($? == 0) then   
   plot_mcmc.sh
endif

# eof
