#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "LHAPDF/LHAPDF.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "Ana/PDFs/interface/pdf_reweighting.h"

using namespace mithep;
ClassImp(pdf_reweighting)

pdf_reweighting::pdf_reweighting(std::string iName, int iPDF, std::string iPDFName) { 
  std::cout << "====> " << iPDFName << std::endl;
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(iPDFName.c_str());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(iPDF);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(iName.c_str(),0);
  bgdEvent.InitTree(0);

  TFile *lFRECO = new TFile("newfile_pdf.root","recreate");
  TTree *newtree = bgdEvent.tree_->CloneTree(0);

  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 100000 == 0) std::cout << "=== RECO Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    double lx1  = bgdEvent.x1_;
    double lx2  = bgdEvent.x2_;
    double lx1_new  = bgdEvent.x1_*3.5/4.0;
    double lx2_new  = bgdEvent.x2_*3.5/4.0;
    int lid1 = bgdEvent.id1_;
    int lid2 = bgdEvent.id2_;
    double lq   = bgdEvent.Q_;
    double lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    double lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    double lxf1_new = LHAPDF::xfx(lx1_new,lq,lid1)/lx1_new;
    double lxf2_new = LHAPDF::xfx(lx2_new,lq,lid2)/lx2_new;
    
    double weight = (lxf1_new*lxf2_new)/(lxf1*lxf2);
    bgdEvent.scale1fb_ = bgdEvent.scale1fb_ * weight;
    newtree->Fill();
  }
  newtree->Write();
  lFRECO->Close();
}
