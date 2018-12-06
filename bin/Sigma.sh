#!/bin/csh

if ("$1" == "") then
	set Omegam_fid = 0.27
else
	set Omegam_fid = $1
endif

fit_s8Omalpha.sh << EOF


$Omegam_fid 



EOF

set alpha = `cat out.txt | perl -ane 'print "$F[6]\n" if /Best/'`

if (! -s config_pmc) ln -s ../config_pmc
add_ded.pl pmcsim Sigma $alpha $Omegam_fid > pmcsim+ded
histograms_sample -c config_pmc+ded pmcsim+ded
meanvar_sample -c config_pmc+ded pmcsim+ded | tee mean+ded
