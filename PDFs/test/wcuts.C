#include <iostream>
#include <string>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "Ana/PDFs/interface/MitPDFNtupleEvent.h"

void wcuts(int nsel = 0) {
  if(nsel == 0 || nsel == 1){ // 0 == RECO vs. GEN, 1 == RECO1 vs. RECO2
    double sumRECO = 0.0;
    double sumGEN = 0.0;

    SmurfTree bgdEvent;
    bgdEvent.LoadTree("newfile_reco.root",0);
    bgdEvent.InitTree(0);
    for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
      if(i0 % 100000 == 0) std::cout << "=== RECO Processed ===> " << i0 << std::endl;
      bgdEvent.tree_->GetEntry(i0);
      sumRECO = sumRECO + bgdEvent.scale1fb_;
    }

    TFile *iFilePDF = new TFile("newfile_gen.root");
    if     (nsel == 0){
      TTree *iTreePDF = (TTree*) iFilePDF->FindObjectAny("PDFTree");
      MitPDFNtupleEvent lMitPDFNtupleEvent(iTreePDF);
      for(int i0 = 0; i0 < iTreePDF->GetEntries(); i0++) {
	if(i0 % 100000 == 0) std::cout << "=== GEN Processed ===> " << i0 << std::endl;
	iTreePDF->GetEntry(i0);
	sumGEN = sumGEN + lMitPDFNtupleEvent.H_weight;
      }
    }
    else if(nsel == 1){
      SmurfTree genEvent;
      genEvent.LoadTree("newfile_gen.root",0);
      genEvent.InitTree(0);
      for(int i0 = 0; i0 < genEvent.tree_->GetEntries(); i0++) {
	if(i0 % 100000 == 0) std::cout << "=== GEN-RECO Processed ===> " << i0 << std::endl;
	genEvent.tree_->GetEntry(i0);
	sumGEN = sumGEN + genEvent.scale1fb_;
      }
    }

    ofstream  oFile("compare.txt");
    std::cout << setprecision (25) << sumRECO/sumGEN << " " << setprecision (25) << sumRECO << std::endl;
    oFile     << setprecision (25) << sumRECO/sumGEN << " " << setprecision (25) << sumRECO << std::endl;
    oFile.close();
  }
}
