#!/bin/csh

# newdir_pmc.sh
# Martin Kilbinger 2008
# Creates a new directory and links necessary data files, according to the
# desired experiment(s). Different experiments are joined by bitwise
# adding the corresponding numbers.

if ("$1" == "-h") then
    echo "Usage: newdir_pmc.sh [DIR]"
    echo "Directory DIR (default: read on input) is created."
    echo 'Links are set to data files in \$COSMOPMC/data.'
    echo 'Parameter files are copied on request from \$COSMOPMC/par_files.'
    exit 1
endif

set bitcode = (1    2  4  8    16	            32  64    128 256)
set type    = (lens Nz Sn WMAP WMAP_Distance_Priors BAO Mring LSS halo)

# Default parameters
set DIR  = "new"
set TYPE = 8

if ("$1" == "") then
    echo -n "Directory name? [$DIR] "
    set dir = $<
    if ("$dir" == "") then
	set dir = $DIR
    endif
else
    set dir = $1
endif

mkdir -p $dir
cd $dir


@ i = 1
while ($i <= $#bitcode)
    perl -e 'printf "   %-15s %5d\n", '$type[$i]', '$bitcode[$i]';'
    @ i ++
end

echo -n "Data type? [$TYPE] "
set mytype = $<
if ("$mytype" == "") set mytype = $TYPE


# Copy data files
@ i = 1
while ($i <= $#bitcode)

    if ($mytype & $bitcode[$i]) then

	echo $bitcode[$i]

	# Choose subdirectory if existant
	if (-e $COSMOPMC/data/$type[$i]) then
	    set subdir = (`find $COSMOPMC/data/$type[$i] -mindepth 1 -maxdepth 1 -type d | grep -v \.svn`)
	    if ($#subdir > 0) then
		@ s = 1
		while ($s <= $#subdir)
		    set x = `basename $subdir[$s]`
		    perl -e 'printf "   %-15s %5d\n", $ARGV[0], $ARGV[1]' $x $s
		    @ s ++
		end
		echo -n "Data subtype? [1] "
		set subtype = $<
		if ("$subtype" == "") set subtype = 1
		set subdir = $subdir[$subtype]
	    else
		set subdir = $COSMOPMC/data/$type[$i]
	    endif

	    # Create links to files
	    foreach f ($subdir/*)
	    	echo Linking $f
		ln -sf $f
	    end

	endif

    endif

    @ i ++

end


# Type-specific additional links
if ($mytype & 8) then
   # WMAP
   ln -s $COSMOPMC/CMB/scamb
   if (-e $COSMOPMC/Makefile.host) then
      set wmap = (`grep MY_WMAP $COSMOPMC/Makefile.host | grep -v '#'`)
      ln -s $wmap[3] WMAP
   endif
endif

# Copy cosmo parameter files
set parfile = ()
if ($mytype & 1) set parfile = ($parfile cosmo.par cosmo_lens.par cosmo_3rd.par)
if ($mytype & 4) set parfile = ($parfile cosmo.par cosmo_SN.par)
if ($mytype & 8 || $mytype & 16 || $mytype & 32) set parfile = ($parfile cosmoDP.par)
if ($mytype & 64) set parfile = ($parfile cosmo.par cosmo_lens.par)
if ($mytype & 256) set parfile = ($parfile cosmo.par halomodel.par)

echo -n "Copy the following parameter files: $COSMOPMC/par_files/{$parfile}? [y] "
set ans = $<
if ("$ans" == "") then
    set ans = "y"
endif
if ("$ans" == "y") then
    foreach file ($parfile)
	cp $COSMOPMC/par_files/$file .
    end
endif

cd ..
echo $dir created

# end


