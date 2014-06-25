#include <iostream>
#include <string>
#include <fstream>
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "LHAPDF/LHAPDF.h"
#include "/afs/cern.ch/user/c/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "Analysis/PDFs/interface/pdfFiller.h"

using namespace mithep;
ClassImp(pdfFiller)

pdfFiller::pdfFiller(std::string iName) {
  TString pdfA = "CT10.LHgrid";            const int NpdfA = 52;
  TString pdfB = "NNPDF21_100.LHgrid";     const int NpdfB = 100;
  TString pdfC = "MSTW2008nlo68cl.LHgrid"; const int NpdfC = 40;
  TString pdfAlphaA0 = "CT10as.LHgrid";
  TString pdfAlphaA1 = "CT10as.LHgrid";
  TString pdfAlphaB0 = "NNPDF21_as_0117_100.LHgrid";
  TString pdfAlphaB1 = "NNPDF21_as_0121_100.LHgrid";
  TString pdfAlphaC0 = "MSTW2008nlo68cl_asmz+68cl.LHgrid";
  TString pdfAlphaC1 = "MSTW2008nlo68cl_asmz-68cl.LHgrid";
  std::map <int , std::vector<double> > pdfAWeights_;
  std::map <int , std::vector<double> > pdfBWeights_;
  std::map <int , std::vector<double> > pdfCWeights_;
  std::map <int , std::vector<double> > pdfDWeights_;

  std::string pName; std::string pFile; 
  float lq = 0; float lx1 = 0; float lx2 = 0; 
  int lid1  = 0; int  lid2 = 0;; 

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(iName.c_str(),-1);
  bgdEvent.InitTree(0);
  double lxf1,lxf2;

  // CT10
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfA.Data());
  LHAPDF::getDescription();
  for(int j=0; j<=NpdfA; j++){
    LHAPDF::usePDFMember(j);
    for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
      if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
      bgdEvent.tree_->GetEntry(i0);

      lx1  = bgdEvent.x1_;
      lx2  = bgdEvent.x2_;
      lid1 = bgdEvent.id1_;
      lid2 = bgdEvent.id2_;
      lq   = bgdEvent.Q_;
      lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
      lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
      
      pdfAWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
    }
  }

  // NNPDF
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfB.Data());
  LHAPDF::getDescription();
  for(int j=0; j<=NpdfB; j++){
    LHAPDF::usePDFMember(j);
    for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
      if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
      bgdEvent.tree_->GetEntry(i0);

      lx1  = bgdEvent.x1_;
      lx2  = bgdEvent.x2_;
      lid1 = bgdEvent.id1_;
      lid2 = bgdEvent.id2_;
      lq   = bgdEvent.Q_;
      lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
      lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
      
      pdfBWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
    }
  }

  // MSTW
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfC.Data());
  LHAPDF::getDescription();
  for(int j=0; j<=NpdfC; j++){
    LHAPDF::usePDFMember(j);
    for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
      if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
      bgdEvent.tree_->GetEntry(i0);

      lx1  = bgdEvent.x1_;
      lx2  = bgdEvent.x2_;
      lid1 = bgdEvent.id1_;
      lid2 = bgdEvent.id2_;
      lq   = bgdEvent.Q_;
      lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
      lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
      
      pdfCWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
    }
  }

  // CT10-Alpha0
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaA0.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }
  
  // CT10-Alpha1
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaA1.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(10);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }

  // NNPDF-Alpha0
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaB0.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }
  
  // NNPDF-Alpha1
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaB1.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }
  
  // MSTW-Alpha0
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaC0.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }
  
  // MSTW-Alpha1
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(pdfAlphaC1.Data());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(0);
  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    
    pdfDWeights_[i0].push_back(lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0);
  }

  std::vector<double> pdfWeightsA;
  std::vector<double> pdfWeightsB;
  std::vector<double> pdfWeightsC;
  std::vector<double> pdfWeightsD;

  TString ofn(iName);
  ofn.ReplaceAll(".root","_pdf.root");
  TFile *out = TFile::Open( Form("%s" ,ofn.Data() ) ,"RECREATE" );
  out->cd();
  TTree *clone = bgdEvent.tree_->CloneTree(0);
  clone->Branch("pdfWeightsA",    "std::vector<double>",   &pdfWeightsA);
  clone->Branch("pdfWeightsB",    "std::vector<double>",   &pdfWeightsC);
  clone->Branch("pdfWeightsC",    "std::vector<double>",   &pdfWeightsB);
  clone->Branch("pdfWeightsD",    "std::vector<double>",   &pdfWeightsD);

  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 100000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);
    pdfWeightsA.clear();
    pdfWeightsB.clear();
    pdfWeightsC.clear();
    pdfWeightsD.clear();
    for(unsigned int j=1; j<pdfAWeights_[i0].size(); j++) pdfWeightsA.push_back(pdfAWeights_[i0][j]/pdfAWeights_[i0][0]);
    for(unsigned int j=1; j<pdfBWeights_[i0].size(); j++) pdfWeightsB.push_back(pdfBWeights_[i0][j]/pdfBWeights_[i0][0]);
    for(unsigned int j=1; j<pdfCWeights_[i0].size(); j++) pdfWeightsC.push_back(pdfCWeights_[i0][j]/pdfCWeights_[i0][0]);

    pdfWeightsD.push_back(pdfDWeights_[i0][0]/pdfAWeights_[i0][0]);
    pdfWeightsD.push_back(pdfDWeights_[i0][1]/pdfAWeights_[i0][0]);

    pdfWeightsD.push_back(pdfDWeights_[i0][2]/pdfBWeights_[i0][0]);
    pdfWeightsD.push_back(pdfDWeights_[i0][3]/pdfBWeights_[i0][0]);

    pdfWeightsD.push_back(pdfDWeights_[i0][4]/pdfCWeights_[i0][0]);
    pdfWeightsD.push_back(pdfDWeights_[i0][5]/pdfCWeights_[i0][0]);
    clone->Fill();
  }
  
  clone->Write(); 
  out->Close();
}
