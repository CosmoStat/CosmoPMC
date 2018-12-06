#!/bin/csh

# Martin Kilbinger 2008
# Prints the correlation coefficient of a matrix to stdout.

if ("$1" == "" || "$1" == "-h") then
   echo "Usage: corr_coeff filename [mvdens|block]"
   exit 1
endif

set format = 0
if ("$2" == "mvdens") then
    set format = 0
else if ("$2" == "block") then
    set format = 1
endif

echo include, \"stuff.i\" > cctmp.i
if ($format == 0) then
    echo mvd = mvdens_read\(\"$1\"\) >> cctmp.i
    echo "matrix = *(mvd.var)" >> cctmp.i
else if ($format == 1) then
    echo read_table, \"$1\", matrix >> cctmp.i
endif
# The following line does *not* write to stderr !?!
#echo 'write, stderr, format="%d x %d - matrix detected\\n", dimsof(matrix)(2), dimsof(matrix)(3)' >> cctmp.i
echo "r   = corr_coeff(matrix, zero=1e-20)" >> cctmp.i
echo write_matrix, r, format=\"% .3f \" >> cctmp.i
echo quit >> cctmp.i

yorick -batch cctmp.i

rm -rf cctmp.i

# end

