#include <iostream>
#include <string>
#include <fstream>
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "LHAPDF/LHAPDF.h"
#include "/afs/cern.ch/user/c/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/afs/cern.ch/user/c/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/afs/cern.ch/user/c/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "Analysis/PDFs/interface/pdf.h"

using namespace mithep;
ClassImp(pdf)

// 1 ==> WW selection
// 2 ==> Full WH->3l selection
// 3 ==> WW within H->WW signal region
// 4 ==> Z->ll selection
// 5 ==> ZZ->llnn selection
// 6 ==> Full ZH->3l+2jets selection
// 7 ==> WW same-sign

pdf::pdf(std::string iName,int iPDF,std::string iPDFName, int nsel, unsigned int nJetsType) { 
  std::cout << "====> " << iPDFName << std::endl;
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(iPDFName.c_str());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(iPDF);
  std::string pName; std::string pFile; 
  float lq = 0; float lx1 = 0; float lx2 = 0; 
  int lid1  = 0; int  lid2 = 0;; 

  double sumRECO = 0.0;
  double sumGEN = 0.0;

  bool useDYMVA = true;
  bool WWXSSel = false; double ptLepMin = 10.0;
  if(nsel == 1) WWXSSel = true;
  if(WWXSSel == true) ptLepMin = 20.;
  unsigned int patternTopVeto = SmurfTree::TopVeto;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(iName.c_str(),0);
  bgdEvent.InitTree(0);

  if(nsel <= 0 || nsel >= 8) assert(0);

  for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
    if(i0 % 1000000 == 0) std::cout << "=== Processed ===> " << i0 << std::endl;
    bgdEvent.tree_->GetEntry(i0);

    lx1  = bgdEvent.x1_;
    lx2  = bgdEvent.x2_;
    lid1 = bgdEvent.id1_;
    lid2 = bgdEvent.id2_;
    lq   = bgdEvent.Q_;
    double lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
    double lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
    sumGEN = sumGEN + lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0;

    if(bgdEvent.lep1_.Pt() < 1.0) continue;

    bool pass = false;
    if     (nsel == 1 || nsel == 3) {
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
      double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
      bool   passMET = usedMet > 20.0;
      if(useDYMVA == false){
	if     (bgdEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
	else if(bgdEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
	else                          passMET = passMET && (bgdEvent.met_ > 45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      } else {
	if     (bgdEvent.njets_ == 0) passMET = passMET && (bgdEvent.dymva_ >  0.88 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
	else if(bgdEvent.njets_ == 1) passMET = passMET && (bgdEvent.dymva_ >  0.84 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
	else                          passMET = passMET && (bgdEvent.met_   >  45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      }
      bool dPhiDiLepJetCut = true;
      if(useDYMVA == false){
	if(bgdEvent.njets_ <= 1) dPhiDiLepJetCut = bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
      	                                           bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
	else                     dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
    	                                           bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
      }
      if(bgdEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                           bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;

      double ptLLCut = 30.0; if(bgdEvent.type_ == SmurfTree::ee || bgdEvent.type_ == SmurfTree::mm) ptLLCut = 45.0;
      bool passNjets         = bgdEvent.njets_ == nJetsType; if(nJetsType == 3) passNjets = bgdEvent.njets_ <= 1;
      bool preselCuts        = bgdEvent.lep1_.Pt() > 20.0 && bgdEvent.lep2_.Pt() > ptLepMin && bgdEvent.dilep_.Pt() > ptLLCut;
      bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
      bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
      bool passMass          = bgdEvent.dilep_.M() > 12 && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);

      pass = 
  	charge == 0 &&
  	passNjets && preselCuts && passBtagVeto && pass3rLVeto && passMass && dPhiDiLepJetCut && passMET;

      if     (nsel == 3){ // WW within H->WW signal region
        pass = pass &&
          bgdEvent.dilep_.M() < 50.0 &&
	  bgdEvent.mt_ < 125.0;
      }
    }
    else if(nsel == 2) {
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      double massZMin = 999.0; double massMin = 999.0; double dRMin = 999.0;
      if(bgdEvent.lid3_ != 0){
  	massZMin = trilepton_info(0,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
  				    bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
        			    bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
        			    bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
  	massMin  = trilepton_info(1,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
  				    bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
        			    bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
        			    bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
  	dRMin = trilepton_info(2,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
  				 bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
        			 bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
        			 bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      }
      pass = 
  	bgdEvent.lid3_ != 0 &&
  	abs(charge) == 1 &&
  	bgdEvent.lep1_.Pt() > 20. &&
  	bgdEvent.lep2_.Pt() > 10. &&
  	bgdEvent.lep3_.Pt() > 10. &&
  	TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > 40. &&
  	(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
  	bgdEvent.jet1_.Pt() <= 40. &&
        massMin > 12 && massMin < 100 &&
        dRMin < 2.0 &&
	massZMin > 25;
    }
    else if(nsel == 4) {
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
      pass = 
      (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
       charge == 0 &&
       bgdEvent.lep1_.Pt() > 20. &&
       bgdEvent.lep2_.Pt() > 20. &&
      (fabs(bgdEvent.dilep_.M()-91.1876) < 15.);
    }
    else if(nsel == 5) {
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
      pass = 
      (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
      (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
       charge == 0 &&
       bgdEvent.lep1_.Pt() > 20. &&
       bgdEvent.lep2_.Pt() > 20. &&
      (fabs(bgdEvent.dilep_.M()-91.1876) < 15.) &&
      bgdEvent.met_ > 120.0 &&
      DeltaPhi(bgdEvent.dilep_.Phi() ,bgdEvent.metPhi_) > 2.7 &&
      fabs(bgdEvent.met_-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < 0.25 &&
      bgdEvent.njets_ == nJetsType;
    }
    else if(nsel == 6) {
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      double massZMin = 999.0; double massMin = 999.0;
      if(bgdEvent.lid3_ != 0){
  	massZMin = trilepton_info(0,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
  				    bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
        			    bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
        			    bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
  	massMin  = trilepton_info(1,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
  				    bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
        			    bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
        			    bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      }
      pass = 
  	bgdEvent.lid3_ != 0 &&
  	abs(charge) == 1 &&
  	bgdEvent.lep1_.Pt() > 20. &&
  	bgdEvent.lep2_.Pt() > 10. &&
  	bgdEvent.lep3_.Pt() > 10. &&
  	bgdEvent.njets_ >= 2 &&
  	massZMin < 15 &&
        massMin > 12;
    }
    else if(nsel == 7) {
      int centrality = 0;
      if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
  	  (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
  	 ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
  	  (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
      double metMin = 30.0; if(bgdEvent.type_ == SmurfTree::ee) metMin = 40.0;
      int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);

      int newId=int(bgdEvent.jet1McId_);
      int qDisAgree=int((newId%1000-newId%100)/100);
      int hasZCand=int(newId/1000);
      int trackSel[4] = {int((bgdEvent.jet2McId_%100-bgdEvent.jet2McId_%10)/10),int((bgdEvent.jet2McId_%1000-bgdEvent.jet2McId_%100)/10),int((bgdEvent.jet2McId_%10000-bgdEvent.jet2McId_%1000)/10),int(bgdEvent.jet2McId_/10000)};

      bool passNjets         = bgdEvent.njets_ >= 2;
      bool passMET           = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > metMin;
      bool preselCuts        = bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 20.;
      bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
      bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && hasZCand == 0 && trackSel[0]+trackSel[1]+trackSel[2] == 0;
      bool passMass          = bgdEvent.dilep_.M() > 50.0 && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ != SmurfTree::ee);
      bool passVBFSel        = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && centrality == 1;

      pass = 
       charge != 0 &&
       passNjets  == true && passMET == true &&
       preselCuts == true && bgdEvent.dilep_.M() > 15.0 && qDisAgree == 0 &&
       passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto;
    }

    if(pass == false) continue;

    sumRECO = sumRECO + lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0;
  }

  ofstream  oFile("compare.txt");
  std::cout << std::setprecision(25) << sumRECO/sumGEN << " " << std::setprecision(25) << sumRECO << std::endl;
  oFile     << std::setprecision(25) << sumRECO/sumGEN << " " << std::setprecision(25) << sumRECO << std::endl;
  oFile.close();

}
