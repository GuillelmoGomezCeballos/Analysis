#!/bin/tcsh

source $HOME/EVAL65 5_2_3_patch3;

setenv dirB /home/ceballos/condor/old_53x/histo_s12-wwj-v7a_all_noskim.root;

rm -f pdf_cteq66.txt;
touch pdf_cteq66.txt;
setenv MAXSETS 44;
@ count00 = 0;
while ($count00 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count00},\"cteq66.LHgrid\",2\);
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count00},\"cteq66.LHgrid\",3\);
  root -l -q -b wcuts.C+\(1\);
  cat compare.txt >> pdf_cteq66.txt;
  @ count00++;
end

rm -f pdf_nnpdf.txt;
touch pdf_nnpdf.txt;
setenv MAXSETS 100;
@ count01 = 0;
while ($count01 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count01},\"NNPDF20_100.LHgrid\",2\);
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count01},\"NNPDF20_100.LHgrid\",3\);
  root -l -q -b wcuts.C+\(1\);
  cat compare.txt >> pdf_nnpdf.txt;
  @ count01++;
end

rm -f pdf_mstw.txt;
touch pdf_mstw.txt;
setenv MAXSETS 40;
@ count02 = 0;
while ($count02 <= $MAXSETS)
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count02},\"MSTW2008nlo68cl.LHgrid\",2\);
  root -l -q -b runSinglePDF.C+\(\"${dirB}\",${count02},\"MSTW2008nlo68cl.LHgrid\",3\);
  root -l -q -b wcuts.C+\(1\);
  cat compare.txt >> pdf_mstw.txt;
  @ count02++;
end

rm -f pdf_cteq66_alphas.txt;
touch pdf_cteq66_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"cteq66alphas.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"cteq66alphas.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_cteq66_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",4,\"cteq66alphas.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",4,\"cteq66alphas.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_cteq66_alphas.txt;

rm -f pdf_mstw_alphas.txt;
touch pdf_mstw_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz+68cl.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz+68cl.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_mstw_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz-68cl.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"MSTW2008nlo68cl_asmz-68cl.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_mstw_alphas.txt;

setenv LHAPATH /afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets
rm -f pdf_nnpdf_alphas.txt;
touch pdf_nnpdf_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF20_as_0117_100.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF20_as_0117_100.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_nnpdf_alphas.txt;
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF20_as_0121_100.LHgrid\",2\);
root -l -q -b runSinglePDF.C+\(\"${dirB}\",0,\"NNPDF20_as_0121_100.LHgrid\",3\);
root -l -q -b wcuts.C+\(1\);
cat compare.txt >> pdf_nnpdf_alphas.txt;

rm compare.txt newfile_reco.root newfile_gen.root;

root -l -q -b final_pdf_error.C+;

exit 0;
