#!/bin/tcsh

source $HOME/EVAL65 5_3_14;

### NSEL
### 1 ==> WW selection
### 2 ==> Full WH->3l selection
### 3 ==> WW within H->WW signal region
### 4 ==> Z->ll selection
### 5 ==> ZZ->llnn selection
### 6 ==> Full ZH->3l+2jets selection
### 7 ==> WW same-sign

setenv dirB /home/ceballos/condor/old_53x/histo_s12-wwj-v7a_all_noskim.root;
setenv NSEL 1;
setenv NJETS $2;

if      ($1 == '1') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-wwj-v7a_all_noskim.root;
  setenv NSEL 1;
else if ($1 == '2') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ww-mcatnlo-v7a_all_noskim.root;
  setenv NSEL 1;
else if ($1 == '3') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ww-powheg-v7c_all_noskim.root;
  setenv NSEL 1;
else if ($1 == '4') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ggww2l-v7a_all_noskim.root;
  setenv NSEL 1;
else if ($1 == '5'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-h125ww2l-gf-v7a_all_noskim.root;
  setenv NSEL 1;

else if ($1 == '6') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-wwj-v7a_all_noskim.root;
  setenv NSEL 3;
else if ($1 == '7') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ww-mcatnlo-v7a_all_noskim.root;
  setenv NSEL 3;
else if ($1 == '8') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ww-powheg-v7c_all_noskim.root;
  setenv NSEL 3;
else if ($1 == '9') then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-ggww2l-v7a_all_noskim.root;
  setenv NSEL 3;
else if ($1 == '10'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-h125ww2l-gf-v7a_all_noskim.root;
  setenv NSEL 3;

else if ($1 == '11'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-zz2l2n-v7a_all_noskim.root;
  setenv NSEL 5;
else if ($1 == '12'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-wz3ln-v7a_all_noskim.root;
  setenv NSEL 5;
else if ($1 == '13'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-h125inv-zh-v7a_all_noskim.root;
  setenv NSEL 5

else if ($1 == '14'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_wwss_qed_4_qcd_99_lt012_all_noskim.root;
  setenv NSEL 7
else if ($1 == '15'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-wz3ln-v7a_all_noskim.root;
  setenv NSEL 7

else if ($1 == '16'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-zz4l-v7a_all_noskim.root;
  setenv NSEL 8;
else if ($1 == '17'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-wz3ln-v7a_all_noskim.root;
  setenv NSEL 8;
else if ($1 == '18'  ) then
  setenv dirB /home/ceballos/condor/old_53x/histo_s12-h95invg-zh-v7c_all_noskim.root;
  setenv NSEL 8

else if ($1 != '0'  ) then
  exit;
endif

rm -f pdf_cteq66.txt;
touch pdf_cteq66.txt;
setenv MAXSETS 52;
@ count00 = 0;
while ($count00 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count00},\"CT10.LHgrid\",${NSEL},${NJETS}\);
  cat compare.txt >> pdf_cteq66.txt;
  @ count00++;
end

rm -f pdf_nnpdf.txt;
touch pdf_nnpdf.txt;
setenv MAXSETS 100;
@ count01 = 0;
while ($count01 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count01},\"NNPDF21_100.LHgrid\",${NSEL},${NJETS}\);
  cat compare.txt >> pdf_nnpdf.txt;
  @ count01++;
end

rm -f pdf_mstw.txt;
touch pdf_mstw.txt;
setenv MAXSETS 40;
@ count02 = 0;
while ($count02 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count02},\"MSTW2008nlo68cl.LHgrid\",${NSEL},${NJETS}\);
  cat compare.txt >> pdf_mstw.txt;
  @ count02++;
end

rm -f pdf_cteq66_alphas.txt;
touch pdf_cteq66_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"CT10as.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_cteq66_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",10,\"CT10as.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_cteq66_alphas.txt;

rm -f pdf_mstw_alphas.txt;
touch pdf_mstw_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz+68cl.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_mstw_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz-68cl.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_mstw_alphas.txt;

#setenv LHAPATH /afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets
rm -f pdf_nnpdf_alphas.txt;
touch pdf_nnpdf_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF21_as_0117_100.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_nnpdf_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF21_as_0121_100.LHgrid\",${NSEL},${NJETS}\);
cat compare.txt >> pdf_nnpdf_alphas.txt;

rm -f compare.txt;
root -l -q -b final_pdf_error.C+;

exit 0;
