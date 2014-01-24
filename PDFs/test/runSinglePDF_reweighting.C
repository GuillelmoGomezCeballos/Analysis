#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Ana/PDFs/interface/pdf_reweighting.h"
void runSinglePDF_reweighting(std::string iName, int iPDF, std::string iPDFName) {
  mithep::pdf_reweighting a(iName, iPDF, iPDFName);
}
// iName    == input file
// iPDF     == number of the PDF to run
// iPDFName == PDF name (cteq66.LHgrid, NNPDF20_100.LHgrid, MSTW2008nlo68cl.LHgrid)
// nsel     == 0
