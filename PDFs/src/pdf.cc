#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "LHAPDF/LHAPDF.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Ana/nt_scripts/trilepton.h"
#include "Ana/PDFs/interface/MitPDFNtupleEvent.h"
#include "Ana/PDFs/interface/pdf.h"

using namespace mithep;
ClassImp(pdf)

// 0 ==> num - WW selection
// 1 ==> den - no cuts
// 2 ==> num - Full WH->3l selection
// 3 ==> den - WH->3l selection, reverting Z veto
// 4 ==> num - WW within H->WW signal region
// 5 ==> den - WW control region
// 6 ==> num - Z->ll selection
// 7 ==> num - ZZ->llnn selection
// 8 ==> num - Full ZH->3l+2jets selection
// 9 ==> num - WW same-sign

pdf::pdf(std::string iName,int iPDF,std::string iPDFName, int nsel, unsigned int nJetsType) { 
  std::cout << "====> " << iPDFName << std::endl;
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(iPDFName.c_str());
  LHAPDF::getDescription();
  LHAPDF::usePDFMember(iPDF);
  std::string pName; std::string pFile; 
  float lq = 0; float lx1 = 0; float lx2 = 0; 
  int lid1  = 0; int  lid2 = 0;; 

  TFile *iFile = new TFile(iName.c_str());
  bool WWXSSel = true;
  double ptLepMin = 10.0;
  if(WWXSSel == true) ptLepMin = 20.;
  unsigned int patternTopVeto         = SmurfTree::TopVeto;

  if(nsel == 0 || nsel == 2 || nsel == 3 || nsel == 4 || nsel == 5 || nsel == 6 || nsel == 7 || nsel == 8 || nsel == 9){
    SmurfTree bgdEvent;
    bgdEvent.LoadTree(iName.c_str(),0);
    bgdEvent.InitTree(0);

    SmurfTree         fSmurfTree;
    TFile *lFRECO;
    if     (nsel == 0 || nsel == 2 || 
            nsel == 4 || nsel == 6 || 
	    nsel == 7 || nsel == 8 || 
	    nsel == 9)   lFRECO = new TFile("newfile_reco.root","recreate");
    else if(nsel == 3 || nsel == 5)   lFRECO = new TFile("newfile_gen.root","recreate");
    fSmurfTree.CreateTree(0);

    for(int i0 = 0; i0 < bgdEvent.tree_->GetEntries(); i0++) {
      if(i0 % 100000 == 0) std::cout << "=== RECO Processed ===> " << i0 << std::endl;
      bgdEvent.tree_->GetEntry(i0);

      bool pass = false;
      if     (nsel == 6) {
	int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
        pass = 
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
        (fabs(bgdEvent.dilep_.M()-91.1876) < 15.);
      }
      else if(nsel == 9) {
        int centrality = 0;
        if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
            (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
           ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
            (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
        double metMin = 30.0; if(bgdEvent.type_ == SmurfTree::ee) metMin = 40.0;
	int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
        bool passNjets    = bgdEvent.njets_ >= 2;
        bool passMET      = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > metMin;
        bool preselCuts   = bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 20.;
        bool passBtagVeto = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
        bool pass3rLVeto  = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
        bool passMass     = bgdEvent.dilep_.M() > 200.0;
        bool passVBFSel   = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && centrality == 1;
        pass = 
         charge != 0 &&
	 passNjets    == true && passMET     == true && preselCuts == true && 
	 passBtagVeto == true && pass3rLVeto == true && passMass   == true && 
	 passVBFSel   == true;
      }
      else if(nsel == 0 || nsel == 4 || nsel == 5) {
        double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
        bool   passMET = usedMet > 20.;
        if     (bgdEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
        else if(bgdEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
        else                          passMET = passMET && (bgdEvent.met_ > 45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);

	bool dPhiDiLepJetCut = true;
	if(bgdEvent.njets_ <= 1) dPhiDiLepJetCut = bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
   						   bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
	else		         dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
        					   bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;

	unsigned int Njet3 = bgdEvent.njets_;
        if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
          if(bgdEvent.jet3_.pt() <= 30)					                                         Njet3 = 2;
          else if(bgdEvent.jet3_.pt() > 30 && (
    	    (bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() < 0) ||
    	    (bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() < 0)))   Njet3 = 0;
          else							                                                 Njet3 = 2;
          if(bgdEvent.njets_ < 2 || bgdEvent.njets_ > 3)                                                         Njet3 = 0;
        }

	int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
        pass = 
	  bgdEvent.dilep_.M()   > 12 &&
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          Njet3 == nJetsType &&
          bgdEvent.lep1_.Pt() > 20. &&
          bgdEvent.lep2_.Pt() > ptLepMin &&
          passMET == true &&
	  bgdEvent.dilep_.Pt() > 45 &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dPhiDiLepJetCut == true;
	if     (nsel == 4){ // WW within H->WW signal region
	  pass = pass &&
	    bgdEvent.dilep_.M() < 50.0;
        }
	else if(nsel == 5){ // WW control region
	  pass = pass &&
	    bgdEvent.dilep_.M() > 100.0;
	}
      }
      else if(nsel == 2 || nsel == 3) {
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
	  dRMin < 2.0;
	if     (nsel == 2){ // Full WH->3l selection
	  pass = pass &&
	    massZMin > 25;
        }
	else if(nsel == 3){ // WH->3l selection, reverting Z veto
	  pass = pass &&
	    massZMin <= 25;
	}
      }
      else if(nsel == 7) {
	int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
        pass = 
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
	(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
        (fabs(bgdEvent.dilep_.M()-91.1876) < 15.) &&
	bgdEvent.met_ > 120.0 &&
	DeltaPhi(bgdEvent.dilep_.Phi() ,bgdEvent.metPhi_)*180.0/TMath::Pi() > 160.0 &&
	fabs(bgdEvent.met_-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < 0.25 &&
	bgdEvent.njets_ == nJetsType;
      }
      else if(nsel == 8) {
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

      if(pass == false) continue;

      lx1  = bgdEvent.x1_;
      lx2  = bgdEvent.x2_;
      lid1 = bgdEvent.id1_;
      lid2 = bgdEvent.id2_;
      lq   = bgdEvent.Q_;
      double lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
      double lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
      fSmurfTree.scale1fb_ = bgdEvent.scale1fb_ *
                           lxf1*lxf2/bgdEvent.pdf1_/bgdEvent.pdf2_/1000.0;
      fSmurfTree.tree_->Fill();
    }
    fSmurfTree.tree_->Write();
    lFRECO->Close();
  }
  else if(nsel == 1){
    // No cuts part
    TTree *iTreeGEN = (TTree*) iFile->FindObjectAny("PDFTree");
    MitPDFNtupleEvent lMitPDFNtupleEvent(iTreeGEN);

    TFile *lFGEN = new TFile("newfile_gen.root","recreate");
    //lFGEN->cd();
    TTree *lFTreeGEN = (TTree*) iTreeGEN->CloneTree(0);
    for(int i0 = 0; i0 < iTreeGEN->GetEntries(); i0++) {
      if(i0 % 100000 == 0) std::cout << "=== GEN Processed ===> " << i0 << std::endl;
      iTreeGEN->GetEntry(i0);
      lx1  = lMitPDFNtupleEvent.H_lx1;
      lx2  = lMitPDFNtupleEvent.H_lx2;
      lid1 = lMitPDFNtupleEvent.H_lid1;
      lid2 = lMitPDFNtupleEvent.H_lid2;
      lq   = lMitPDFNtupleEvent.H_lq;
      double lxf1 = LHAPDF::xfx(lx1,lq,lid1)/lx1;
      double lxf2 = LHAPDF::xfx(lx2,lq,lid2)/lx2;
      lMitPDFNtupleEvent.H_weight = lMitPDFNtupleEvent.H_weight *
                                    lxf1*lxf2/lMitPDFNtupleEvent.H_lpdf1/lMitPDFNtupleEvent.H_lpdf2;
      lFTreeGEN->Fill();
    }
    lFTreeGEN->Write();
    lFGEN->Close();
  }
}
