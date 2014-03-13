#include <iostream>
#include <string>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"

void wcuts() {
  double sumRECO = 0.0;
  double sumGEN = 0.0;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree("newfile_reco.root",0);
  bgdEvent.InitTree(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== RECO Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);
    sumRECO = sumRECO + bgdEvent.scale1fb_;
  }

  SmurfTree genEvent;
  genEvent.LoadTree("newfile_gen.root",0);
  genEvent.InitTree(0);
  for(int i0 = 0; i0 < genEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== GEN Processed ===> " << i0 << std::endl;
    genEvent.tree_->GetEntry(i0);
    sumGEN = sumGEN + genEvent.scale1fb_;
  }

  ofstream  oFile("compare.txt");
  std::cout << std::setprecision(25) << sumRECO/sumGEN << " " << std::setprecision(25) << sumRECO << std::endl;
  oFile     << std::setprecision(25) << sumRECO/sumGEN << " " << std::setprecision(25) << sumRECO << std::endl;
  oFile.close();
}
