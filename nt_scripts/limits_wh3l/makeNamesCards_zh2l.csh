#!/bin/tcsh

setenv MAINDIR /home/ceballos/releases/CMSSW_5_2_8/src/Ana/nt_scripts/limits_wh3l;

foreach mass (105 115 125 135 145 175 200 300)

setenv SECDIR sm;
mkdir -p $MAINDIR/$SECDIR/$mass;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm0j_mh${mass}_shape_7TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee0j_mh${mass}_shape_7TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm1j_mh${mass}_shape_7TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee1j_mh${mass}_shape_7TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_shape_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm0j_mh${mass}_shape_8TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_shape_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee0j_mh${mass}_shape_8TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_shape_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm1j_mh${mass}_shape_8TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_shape_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee1j_mh${mass}_shape_8TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_shape_8TeV.txt;

mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm0j_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee0j_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm1j_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee1j_mh${mass}_cut_7TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_cut_7TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm0j_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_cut_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee0j_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_cut_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvmm1j_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_cut_8TeV.txt;
mv $MAINDIR/$SECDIR/histo_limits_zllhinvee1j_mh${mass}_cut_8TeV.txt   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_cut_8TeV.txt;

mv $MAINDIR/$SECDIR/zllhinvmm0j_${mass}_input_7TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvmm0j_input_7TeV.root;
mv $MAINDIR/$SECDIR/zllhinvee0j_${mass}_input_7TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvee0j_input_7TeV.root;
mv $MAINDIR/$SECDIR/zllhinvmm1j_${mass}_input_7TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvmm1j_input_7TeV.root;
mv $MAINDIR/$SECDIR/zllhinvee1j_${mass}_input_7TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvee1j_input_7TeV.root;
mv $MAINDIR/$SECDIR/zllhinvmm0j_${mass}_input_8TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvmm0j_input_8TeV.root;
mv $MAINDIR/$SECDIR/zllhinvee0j_${mass}_input_8TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvee0j_input_8TeV.root;
mv $MAINDIR/$SECDIR/zllhinvmm1j_${mass}_input_8TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvmm1j_input_8TeV.root;
mv $MAINDIR/$SECDIR/zllhinvee1j_${mass}_input_8TeV.root 		      $MAINDIR/$SECDIR/$mass/zllhinvee1j_input_8TeV.root;
sed -i s"/zllhinvmm0j_${mass}_input_7TeV.root/zllhinvmm0j_input_7TeV.root/"   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_shape_7TeV.txt;
sed -i s"/zllhinvee0j_${mass}_input_7TeV.root/zllhinvee0j_input_7TeV.root/"   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_shape_7TeV.txt;
sed -i s"/zllhinvmm1j_${mass}_input_7TeV.root/zllhinvmm1j_input_7TeV.root/"   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_shape_7TeV.txt;
sed -i s"/zllhinvee1j_${mass}_input_7TeV.root/zllhinvee1j_input_7TeV.root/"   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_shape_7TeV.txt;
sed -i s"/zllhinvmm0j_${mass}_input_8TeV.root/zllhinvmm0j_input_8TeV.root/"   $MAINDIR/$SECDIR/$mass/zmmhinv_nj0_shape_8TeV.txt;
sed -i s"/zllhinvee0j_${mass}_input_8TeV.root/zllhinvee0j_input_8TeV.root/"   $MAINDIR/$SECDIR/$mass/zeehinv_nj0_shape_8TeV.txt;
sed -i s"/zllhinvmm1j_${mass}_input_8TeV.root/zllhinvmm1j_input_8TeV.root/"   $MAINDIR/$SECDIR/$mass/zmmhinv_nj1_shape_8TeV.txt;
sed -i s"/zllhinvee1j_${mass}_input_8TeV.root/zllhinvee1j_input_8TeV.root/"   $MAINDIR/$SECDIR/$mass/zeehinv_nj1_shape_8TeV.txt;

end
