#!/bin/csh

# fit_s8Omalpha.sh
# Martin Kilbinger 2007
# Fits a power law sigma_8*Omega_m^alpha = const to the minimum region
# of the Omega_m-sigma_8 banana.
# Usage: fit_s8Omalpha.sh [chi2fname]

# Default parameters
set N         = 2
set CHI2THRES = 2.3
set OMEGAM0   = 0.27
set OUTFILE   = "out.txt"
set TITLE     = ""
set LIMIT     = (0.0 1.0 0.2 1.5)
set FILE      = "chi2_0_1"

if ("$1" != "") set FILE = $1

echo "fit_s8Omalpha.sh: file $FILE"
head -n 1 $FILE

echo -n "enter the number of parameters: [$N] "
set n = $<
if ("$n" == "") then
        set n = $N
endif

echo -n "enter chi^2 threshold: [$CHI2THRES] "
set chi2thres = $<
if ("$chi2thres" == "") then
    set chi2thres = $CHI2THRES
endif

echo -n "enter Omega_m_0: [$OMEGAM0] "
set Omegam0 = $<
if ("$Omegam0" == "") then
    set Omegam0 = $OMEGAM0
endif

echo -n "enter output file: [$OUTFILE] "
set outfile = $<
if ("$outfile" == "") then
    set outfile = $OUTFILE
endif

echo -n "enter title: [$TITLE] "
set title = $<
if ("$title" == "") then
    set title = $TITLE
endif
set title = "$title "

echo -n "limits (xl xu yl yu) [$LIMIT] "
set limit = ($<)
if ("$limit" == "") then
        set limit = ($LIMIT)
endif

set outps = `basename $outfile .txt`

rm -f fit_s8Oma.i
echo include, \"likeli.i\" > fit_s8Oma.i
echo include, \"yutils/lmfit.i\" >> fit_s8Oma.i
echo imin = read_chin\(\"$FILE\", $n, prior=0\) >> fit_s8Oma.i
echo "chi2 = array(double, nel(1), nel(2))" >> fit_s8Oma.i
#echo "for (i=k=1; i<=nel(1); i++) {" >> fit_s8Oma.i
#echo "  for (j=1; j<nel(2); j++,k++) {" >> fit_s8Oma.i
#echo "     chi2(i,j) = chinall(3,k); } }" >> fit_s8Oma.i
echo marginn, chinall, 1, 2, 1, 0 >> fit_s8Oma.i
echo fitpar = fit_power_law\(chi2, 1, 2, $chi2thres, Omegam0=$Omegam0, outname=\"$outfile\", mcmc=1\) >> fit_s8Oma.i;
#echo win, 0 >> fit_s8Oma.i
echo window, display=\"\", hcp=\"$outps\" >> fit_s8Oma.i


echo plot_chi2n, chi2, 1, 2, nice=1, chi2thres=$chi2thres, Omegam0=$Omegam0, fitpar=fitpar, alllevels=1 >> fit_s8Oma.i
echo limits, $limit[1], $limit[2], $limit[3], $limit[4] >> fit_s8Oma.i
echo pltitle, \"$title\" >> fit_s8Oma.i
echo eps, \"$outps\" >> fit_s8Oma.i
echo quit >> fit_s8Oma.i

yorick -batch fit_s8Oma.i

