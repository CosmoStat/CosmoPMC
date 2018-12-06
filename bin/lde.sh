#!/bin/tcsh

# Martin Kilbinger
# Usage: lde.sh file[.eps] [-q]

set name = `basename $1 .tex`

if ( "$2" == "-q" ) then
    latex $name > /dev/null && dvips -E $name -o $name.eps >& /dev/null
else
    latex $name && dvips -E $name -o $name.eps
endif

