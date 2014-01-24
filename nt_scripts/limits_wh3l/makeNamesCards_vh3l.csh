#!/bin/tcsh

setenv MAINDIR /home/ceballos/releases/CMSSW_5_2_8/src/Ana/nt_scripts/limits_wh3l;

foreach mass (110 115 120 125 130 135 140 150 160 170 180 190 200)

setenv SECDIR sm;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_vh3l1_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/vh3l1_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_vh3l2_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/vh3l2_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_vh3l1_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/vh3l1_cut_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_vh3l2_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/vh3l2_cut_8TeV.txt;

mv $MAINDIR/$SECDIR/histo_limits_vh3l1_mh${mass}_shape_7TeV.txt  $MAINDIR/$SECDIR/$mass/vh3l1_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/vh3l1_${mass}.input_7TeV.root		 $MAINDIR/$SECDIR/$mass/vh3l1_input_7TeV.root;
sed -ie s"/vh3l1_${mass}.input_7TeV.root/vh3l1_input_7TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l1_shape_7TeV.txt;

mv $MAINDIR/$SECDIR/histo_limits_vh3l2_mh${mass}_shape_7TeV.txt  $MAINDIR/$SECDIR/$mass/vh3l2_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/vh3l2_${mass}.input_7TeV.root		 $MAINDIR/$SECDIR/$mass/vh3l2_input_7TeV.root;
sed -ie s"/vh3l2_${mass}.input_7TeV.root/vh3l2_input_7TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l2_shape_7TeV.txt;

mv $MAINDIR/$SECDIR/histo_limits_vh3l1_mh${mass}_shape_8TeV.txt  $MAINDIR/$SECDIR/$mass/vh3l1_shape_8TeV.txt;
mv $MAINDIR/$SECDIR/vh3l1_${mass}.input_8TeV.root		 $MAINDIR/$SECDIR/$mass/vh3l1_input_8TeV.root;
sed -ie s"/vh3l1_${mass}.input_8TeV.root/vh3l1_input_8TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l1_shape_8TeV.txt;

mv $MAINDIR/$SECDIR/histo_limits_vh3l2_mh${mass}_shape_8TeV.txt  $MAINDIR/$SECDIR/$mass/vh3l2_shape_8TeV.txt;
mv $MAINDIR/$SECDIR/vh3l2_${mass}.input_8TeV.root		 $MAINDIR/$SECDIR/$mass/vh3l2_input_8TeV.root;
sed -ie s"/vh3l2_${mass}.input_8TeV.root/vh3l2_input_8TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l2_shape_8TeV.txt;

#setenv SECDIR sm4;
#mkdir -p $MAINDIR/$SECDIR/$mass;
#mv $MAINDIR/$SECDIR/histo_limits_vh3l_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/SM4_vh3l_cut_7TeV.txt;
#mv $MAINDIR/$SECDIR/histo_limits_vh3l_mh${mass}_shape_7TeV.txt $MAINDIR/$SECDIR/$mass/vh3l_shape_7TeV.txt;
#mv $MAINDIR/$SECDIR/vh3l${mass}.input_7TeV.root               $MAINDIR/$SECDIR/$mass/vh3l${mass}.input_7TeV.root;
#sed -ie s"/vh3l${mass}.input_7TeV.root/$SECDIR\/${mass}\/vh3l${mass}.input_7TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l_shape_7TeV.txt;

#setenv SECDIR fermiophobic;
#mkdir -p $MAINDIR/$SECDIR/$mass;
#mv $MAINDIR/$SECDIR/histo_limits_vh3l_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/FF_vh3l_cut_7TeV.txt;
#mv $MAINDIR/$SECDIR/histo_limits_vh3l_mh${mass}_shape_7TeV.txt $MAINDIR/$SECDIR/$mass/vh3l_shape_7TeV.txt;
#mv $MAINDIR/$SECDIR/vh3l${mass}.input_7TeV.root               $MAINDIR/$SECDIR/$mass/vh3l${mass}.input_7TeV.root;
#sed -ie s"/vh3l${mass}.input_7TeV.root/$SECDIR\/${mass}\/vh3l${mass}.input_7TeV.root/" $MAINDIR/$SECDIR/$mass/vh3l_shape_7TeV.txt;

end

rm -f $MAINDIR/*/*/vh3l*_shape_7TeV.txte;
rm -f $MAINDIR/*/*/vh3l*_shape_8TeV.txte;
