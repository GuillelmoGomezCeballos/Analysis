#/bin/sh

################### WW and H->WW
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim6.root\",6\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim6.root\",6\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundB.root\",\"backgroundB_skim6.root\",6\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_skim6.root\",6\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundC.root\",\"backgroundC_skim6.root\",6\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundD.root\",\"backgroundD_skim6.root\",6\)

################### bbWW
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim12.root\",12\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim12.root\",12\)

################### e-mu only, non-btagged
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim1.root\",1\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim1.root\",1\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_skim1.root\",1\)

################### trilepton
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_3l.root\",5\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundB.root\",\"backgroundB_3l.root\",5\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_3l.root\",5\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_3l.root\",5\)

################### Z(ll)H(invisible)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim10.root\",10\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim10.root\",10\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_skim10.root\",10\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"gamma.root\",\"gamma_skim10.root\",11\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundZH.root\",\"backgroundZH_skim10.root\",10\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"zh_syst.root\",\"zh_syst_skim10.root\",10\)

root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"dyll.root\",\"dyll_skim13.root\",13\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"dyllpt100.root\",\"dyllpt100_skim10.root\",10\)

if [ $1 == 1 ]; then
################### For both  same-sign dilepton and VBS
source ~/EVAL_SH65 5_3_14;
export INITDIR=$PWD;

cd ~/releases/CMSSW_5_3_14/src;
root -l -q -b Smurf/Analysis/Zll/WGammaStarScaleFactor.C+'("${INITDIR}/wwss_qed_2_qcd_99_ls0ls1.root","${INITDIR}/wwss_qcd_subtr_ls0ls1.root",-1.0)';
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qcd_subtr_ls0ls1.root","${INITDIR}/wwss_qcd_subtr2_ls0ls1.root",91)';
mv ${INITDIR}/wwss_qcd_subtr2_ls0ls1.root ${INITDIR}/wwss_qcd_subtr_ls0ls1.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qed_4_qcd_99_ls0ls1.root","${INITDIR}/wwss_qcdewk_ewkdstype_ls0ls1.root",91)';
cd $INITDIR;
hadd -f wwss_ewk_ewkdstype_ls0ls1.root wwss_qcdewk_ewkdstype_ls0ls1.root wwss_qcd_subtr_ls0ls1.root;

cd ~/releases/CMSSW_5_3_14/src;
root -l -q -b Smurf/Analysis/Zll/WGammaStarScaleFactor.C+'("${INITDIR}/wwss_qed_2_qcd_99_lt012.root","${INITDIR}/wwss_qcd_subtr_lt012.root",-1.0)';
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qcd_subtr_lt012.root","${INITDIR}/wwss_qcd_subtr2_lt012.root",91)';
mv ${INITDIR}/wwss_qcd_subtr2_lt012.root ${INITDIR}/wwss_qcd_subtr_lt012.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qed_4_qcd_99_lt012.root","${INITDIR}/wwss_qcdewk_ewkdstype_lt012.root",91)';
cd $INITDIR;
hadd -f wwss_ewk_ewkdstype_lt012.root wwss_qcdewk_ewkdstype_lt012.root wwss_qcd_subtr_lt012.root;

cd ~/releases/CMSSW_5_3_14/src;
root -l -q -b Smurf/Analysis/Zll/WGammaStarScaleFactor.C+'("${INITDIR}/wwss_qed_2_qcd_99_lm0123.root","${INITDIR}/wwss_qcd_subtr_lm0123.root",-1.0)';
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qcd_subtr_lm0123.root","${INITDIR}/wwss_qcd_subtr2_lm0123.root",91)';
mv ${INITDIR}/wwss_qcd_subtr2_lm0123.root ${INITDIR}/wwss_qcd_subtr_lm0123.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${INITDIR}/wwss_qed_4_qcd_99_lm0123.root","${INITDIR}/wwss_qcdewk_ewkdstype_lm0123.root",91)';
cd $INITDIR;
hadd -f wwss_ewk_ewkdstype_lm0123.root wwss_qcdewk_ewkdstype_lm0123.root wwss_qcd_subtr_lm0123.root;

################### same-sign dilepton
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim8.root\",8\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim8.root\",8\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_skim8.root\",8\)

hadd -f aux1.root backgroundA_skim8.root wwss_qcdewk_ewkdstype_ls0ls1.root wwss_qcd_subtr_ls0ls1.root wwss_qed_2_qcd_99_ls0ls1.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux1.root\",\"backgroundA_skim8_ls0ls1.root\",8\)
rm aux1.root

hadd -f aux1.root backgroundA_skim8.root wwss_qcdewk_ewkdstype_lt012.root wwss_qcd_subtr_lt012.root wwss_qed_2_qcd_99_lt012.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux1.root\",\"backgroundA_skim8_lt012.root\",8\)
rm aux1.root

hadd -f aux1.root backgroundA_skim8.root wwss_qcdewk_ewkdstype_lm0123.root wwss_qcd_subtr_lm0123.root wwss_qed_2_qcd_99_lm0123.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux1.root\",\"backgroundA_skim8_lm0123.root\",8\)
rm aux1.root

################### VBS
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"backgroundA.root\",\"backgroundA_skim14.root\",14\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"data.root\",\"data_skim14.root\",14\)
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"hww_syst.root\",\"hww_syst_skim14.root\",14\)

hadd -f aux2.root backgroundA_skim14.root wwss_qcdewk_ewkdstype_ls0ls1.root wwss_qcd_subtr_ls0ls1.root wwss_qed_2_qcd_99_ls0ls1.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux2.root\",\"backgroundA_skim14_ls0ls1.root\",14\)
rm aux2.root

hadd -f aux2.root backgroundA_skim14.root wwss_qcdewk_ewkdstype_lt012.root wwss_qcd_subtr_lt012.root wwss_qed_2_qcd_99_lt012.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux2.root\",\"backgroundA_skim14_lt012.root\",14\)
rm aux2.root

hadd -f aux2.root backgroundA_skim14.root wwss_qcdewk_ewkdstype_lm0123.root wwss_qcd_subtr_lm0123.root wwss_qed_2_qcd_99_lm0123.root ww_dps.root
root -l -q -b ~/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/skimFastNtuple.C+\(\"aux2.root\",\"backgroundA_skim14_lm0123.root\",14\)
rm aux2.root

fi
