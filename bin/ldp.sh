#!/bin/csh

# Martin Kilbinger 2008
# Calls latex and dvips to create a ps from a latex file.
# Usage: ldp.sh file [-q] [-pdf]

set quiet  = 0
set dvi    = "dvips"
set suff   = "ps"
set flag_out = "-o"


while ($#argv > 0)
    if ("$argv[1]" == "-q") then
	set quiet = 1
    else if ("$argv[1]" == "-pdf") then
	set dvi  = "dvipdf"
	set suff = "pdf"
	set flag_out = ""
    else if ("$argv[1]" == "-h") then
	echo "Usage: lpd.sh FILE [-q] [-pdf]"
	exit 0
    else
	set name   = `basename $argv[1] .tex`
    endif
    shift
end

if ( "$2" == "-q" ) then
    latex $name > /dev/null && $dvi $name $flag_out tmptmp.$suff >& /dev/null && mv tmptmp.$suff $name.$suff
else
    latex $name && $dvi $name $flag_out tmptmp.$suff && mv tmptmp.$suff $name.$suff
endif
