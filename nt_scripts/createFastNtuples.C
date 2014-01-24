#include "TFile.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void createFastNtuples(TString inputFile = "input.root",TString outputFile = "output.root"){

TFile *_file0 = TFile::Open(inputFile);
_file0->cd();
TTree* tree = (TTree*)_file0->Get("tree");
TTree* oldtree = tree;
cout << oldtree->GetEntries() << endl;

// Additional variables
float	       auxVar0_;
auxVar0_       = -999;

TFile *_file2 = TFile::Open(outputFile,"recreate");
_file2->cd();
TTree *newtree2 = oldtree->CopyTree("");
newtree2->SetBranchStatus("*", 1);
TBranch* br_test20 = newtree2->Branch("auxVar0",	   &auxVar0_	  ,   "auxVar0/F");
 for (Long64_t ievt=0; ievt<newtree2->GetEntries();ievt++) {
   newtree2->GetEntry(ievt);
   if (ievt%10000 == 0) cout << ievt<<endl;
   br_test20->Fill();
 }
cout << newtree2->GetEntries() << endl;
_file2->Write();
_file2->Close();

/*
TFile *_file1 = TFile::Open("data_lfake.root","recreate");
_file1->cd();
TTree *newtree1 = oldtree->CopyTree("((cuts & 4) != 4) || ((cuts & 512) != 512)");
  // Additional variables
  LorentzVector  lep3_;
  int            lq3_;
  int            lid3_;
  LorentzVector  jet3_;
  float          jet3Btag_;
  int            lep3McId_;
  int            lep3MotherMcId_;
  int            jet3McId_;
  float          dPhiLep3Jet1_;
  float          dRLep3Jet1_;
  float          dPhiLep3MET_;
  float          mt3_;
  float          jetLowBtag_;
  unsigned int   nSoftMuons_;
  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          higgsPt_;
  float          hPtWeight_;
  LorentzVector* lepPtr3_;
  lq3_  	 = 0;
  lid3_ 	 = 0;
  LorentzVector* jetPtr3_;
  jet3Btag_	 = -999.;
  lep3McId_	 = 0;
  lep3MotherMcId_= 0;
  jet3McId_	 = 0;
  dPhiLep3Jet1_  = -999.;
  dRLep3Jet1_	 = -999.;
  dPhiLep3MET_   = -999.;
  mt3_  	 = -999.;
  jetLowBtag_	 = -999.;
  nSoftMuons_	 = 0;
  Q_		 = -999.;
  id1_  	 = -999.;
  x1_		 = -999.;
  pdf1_ 	 = -999.;  
  id2_  	 = -999.;  
  x2_		 = -999.;
  pdf2_ 	 = -999.;  
  processId_	 = 0;
  higgsPt_	 = -999;
  hPtWeight_	 = -999;

  newtree1->Branch("lep3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
  newtree1->Branch("lq3",	    &lq3_,	    "lq3/I");
  newtree1->Branch("lid3",	    &lid3_,	    "lid3/I");
  newtree1->Branch("jet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
  newtree1->Branch("jet3Btag",      &jet3Btag_,      "jet3Btag/F");
  newtree1->Branch("lep3McId",      &lep3McId_,      "lep3McId/I");
  newtree1->Branch("lep3MotherMcId",&lep3MotherMcId_,"lep3MotherMcId/I");
  newtree1->Branch("jet3McId",      &jet3McId_,      "jet3McId/I");
  newtree1->Branch("dPhiLep3Jet1",  &dPhiLep3Jet1_,  "dPhiLep3Jet1/F");
  newtree1->Branch("dRLep3Jet1",    &dRLep3Jet1_,    "dRLep3Jet1/F");
  newtree1->Branch("dPhiLep3MET",   &dPhiLep3MET_,   "dPhiLep3MET/F");
  newtree1->Branch("mt3",	    &mt3_,	     "mt3/F");
  newtree1->Branch("jetLowBtag",    &jetLowBtag_,    "jetLowBtag/F");
  newtree1->Branch("nSoftMuons",    &nSoftMuons_,    "nSoftMuons/i");
  newtree1->Branch("Q", 	    &Q_     ,	  "Q/F");
  newtree1->Branch("id1",	    &id1_   ,	  "id1/F");
  newtree1->Branch("x1",	    &x1_    ,	  "x1/F");
  newtree1->Branch("pdf1",	    &pdf1_  ,	  "pdf1/F");
  newtree1->Branch("id2",	    &id2_   ,	  "id2/F");
  newtree1->Branch("x2",	    &x2_    ,	  "x2/F");
  newtree1->Branch("pdf2",	    &pdf2_  ,	  "pdf2/F");
cout << newtree1->GetEntries() << endl;
_file1->Write();
_file1->Close();
TFile *_file2 = TFile::Open("data_2l.root","recreate");
_file2->cd();
//TTree *newtree2 = oldtree->CopyTree("((cuts & 4) == 4) && ((cuts & 512) == 512) && met > 20&&(type==1||type==2||TMath::Abs(dilep->mass()-91.1876)>15)");
TTree *newtree2 = oldtree->CopyTree("((cuts & 4) == 4) && ((cuts & 512) == 512)");
  newtree2->SetBranchStatus("*", 1);
  //TTree *clone = t->CloneTree(-1, "fast");
  TBranch* br_test00 = newtree2->Branch("lep3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
  TBranch* br_test01 = newtree2->Branch("lq3",	    &lq3_,	    "lq3/I");
  TBranch* br_test02 = newtree2->Branch("lid3",	    &lid3_,	    "lid3/I");
  TBranch* br_test03 = newtree2->Branch("jet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
  TBranch* br_test04 = newtree2->Branch("jet3Btag",      &jet3Btag_,      "jet3Btag/F");
  TBranch* br_test05 = newtree2->Branch("lep3McId",      &lep3McId_,      "lep3McId/I");
  TBranch* br_test06 = newtree2->Branch("lep3MotherMcId",&lep3MotherMcId_,"lep3MotherMcId/I");
  TBranch* br_test07 = newtree2->Branch("jet3McId",      &jet3McId_,      "jet3McId/I");
  TBranch* br_test08 = newtree2->Branch("dPhiLep3Jet1",  &dPhiLep3Jet1_,  "dPhiLep3Jet1/F");
  TBranch* br_test09 = newtree2->Branch("dRLep3Jet1",    &dRLep3Jet1_,    "dRLep3Jet1/F");
  TBranch* br_test10 = newtree2->Branch("dPhiLep3MET",   &dPhiLep3MET_,   "dPhiLep3MET/F");
  TBranch* br_test11 = newtree2->Branch("mt3",	    &mt3_,	     "mt3/F");
  TBranch* br_test12 = newtree2->Branch("jetLowBtag",    &jetLowBtag_,    "jetLowBtag/F");
  TBranch* br_test13 = newtree2->Branch("nSoftMuons",    &nSoftMuons_,    "nSoftMuons/i");
  TBranch* br_test14 = newtree2->Branch("Q", 	    &Q_     ,	  "Q/F");
  TBranch* br_test15 = newtree2->Branch("id1",	    &id1_   ,	  "id1/F");
  TBranch* br_test16 = newtree2->Branch("x1",	    &x1_    ,	  "x1/F");
  TBranch* br_test17 = newtree2->Branch("pdf1",	   &pdf1_  ,	 "pdf1/F");
  TBranch* br_test18 = newtree2->Branch("id2",	   &id2_   ,	 "id2/F");
  TBranch* br_test19 = newtree2->Branch("x2", 	   &x2_    ,	 "x2/F");
  TBranch* br_test20 = newtree2->Branch("pdf2",	   &pdf2_  ,	 "pdf2/F");
    for (Long64_t ievt=0; ievt<newtree2->GetEntries();ievt++) {
      newtree2->GetEntry(ievt);
      if (ievt%10000 == 0) cout << ievt<<endl;
      br_test00->Fill();
      br_test01->Fill();
      br_test02->Fill();
      br_test03->Fill();
      br_test04->Fill();
      br_test05->Fill();
      br_test06->Fill();
      br_test07->Fill();
      br_test08->Fill();
      br_test09->Fill();
      br_test10->Fill();
      br_test11->Fill();
      br_test12->Fill();
      br_test13->Fill();
      br_test14->Fill();
      br_test15->Fill();
      br_test16->Fill();
      br_test17->Fill();
      br_test18->Fill();
      br_test19->Fill();
      br_test20->Fill();
    }
cout << newtree2->GetEntries() << endl;
_file2->Write();
_file2->Close();
*/
}
