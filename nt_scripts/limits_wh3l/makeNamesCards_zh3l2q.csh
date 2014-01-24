#!/bin/tcsh

setenv MAINDIR /home/ceballos/releases/CMSSW_5_2_3_patch3/src/Ana/nt_scripts/limits_wh3l;

foreach mass (120 130 140 150 160 170 180 190 200)

setenv SECDIR sm_zh3l2q;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_zh3l2q_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/zh3l2q_cut.txt;

setenv SECDIR fermiophobic_zh3l2q;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_zh3l2q_mh${mass}_cut.txt   $MAINDIR/$SECDIR/$mass/FF_zh3l2q_cut.txt;

end
