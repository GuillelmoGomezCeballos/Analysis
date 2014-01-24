#!/bin/tcsh

setenv MAINDIR /home/ceballos/releases/CMSSW_5_2_3_patch3/src/Ana/nt_scripts/limits_wh3l;

foreach mass (120 130 140 150 160 170 180 190 200)

setenv SECDIR sm_vhss;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_whss1_chan4_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/vhss1_cut.txt;
mv $MAINDIR/$SECDIR/histo_limits_whss2_chan4_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/vhss2_cut.txt;

setenv SECDIR fermiophobic_vhss;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_whss1_chan4_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/FF_vhss1_cut.txt;
mv $MAINDIR/$SECDIR/histo_limits_whss2_chan4_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/FF_vhss2_cut.txt;

end
