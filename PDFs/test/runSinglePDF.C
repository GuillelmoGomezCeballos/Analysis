#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Ana/PDFs/interface/pdf.h"
void runSinglePDF(std::string iName,int iPDF,std::string iPDFName, int nsel, unsigned int nJetsType = 0) {
  mithep::pdf a(iName, iPDF, iPDFName, nsel, nJetsType);
}
// iName     == input file
// iPDF      == number of the PDF to run
// iPDFName  == PDF name (cteq66.LHgrid, NNPDF20_100.LHgrid, MSTW2008nlo68cl.LHgrid)
// nsel      == 0 (after ww cuts), 1 (no cuts), 2 (after wh->3l cuts)
// nJetsType == number of reconstructed jets (used in some cases)
