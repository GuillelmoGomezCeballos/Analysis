#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/PDFs/interface/pdfFiller.h"
void runPDFFiller(std::string iName) {
  mithep::pdfFiller a(iName);
}
// iName     == input file
