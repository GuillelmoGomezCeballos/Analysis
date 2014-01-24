#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Ana/nt_scripts/trilepton.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"

const int verboseLevel =   1;

//------------------------------------------------------------------------------
// dataEstimations_EPS
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void dataEstimations_EPS
(
 int     mH  	 = 5,
 double mhAna = 4,
 TString bgdInputFile    = "ntuples_42x_EPS/background42x.root",
 TString signalInputFile = "ntuples_42x_EPS/hww160.root",
 TString dataInputFile   = "ntuples_42x_EPS/data_2l.root"
 )
{
  unsigned int nJetsType = (int)mhAna/10;
  int lDecay    = (int)mhAna%10;
  printf("nJetsType: %d, lDecay = %d\n",nJetsType,lDecay);
  double lumi = 1.092;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  //TFile *fLeptonFRFileM = TFile::Open("/data/smurf/data/EPS/auxiliar/FakeRates_SmurfV6.root");
  //TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  TFile *fLeptonFRFileM = TFile::Open("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_mu_fr.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("mu_fr_m2_15"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  //TFile *fLeptonFRFileE = TFile::Open("/data/smurf/data/EPS/auxiliar/FakeRates_SmurfV6.root");
  //TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  TFile *fLeptonFRFileE = TFile::Open("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_el_fr.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("el_fr_v4_35"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");

  int channel = mH;
  int binc = -1;
  if     (mH == 115) binc = 0;
  else if(mH == 120) binc = 1;
  else if(mH == 130) binc = 2;
  else if(mH == 140) binc = 3;
  else if(mH == 150) binc = 4;
  else if(mH == 160) binc = 5;
  else if(mH == 170) binc = 6;
  else if(mH == 180) binc = 7;
  else if(mH == 190) binc = 8;
  else if(mH == 200) binc = 9;
  else if(mH == 210) binc = 10;
  else if(mH == 220) binc = 11;
  else if(mH == 230) binc = 12;
  else if(mH == 250) binc = 13;
  else if(mH == 300) binc = 14;
  else if(mH == 350) binc = 15;
  else if(mH == 400) binc = 16;
  else if(mH == 450) binc = 17;
  else if(mH == 500) binc = 18;
  else if(mH == 550) binc = 19;
  else if(mH == 600) binc = 20;
  else               mH   = 115;

  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/data/EPS/auxiliar/ggHWW_KFactors_PowhegToHQT.root");
  TH1D *HiggsPtKFactor;
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", mH);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;

  double cutMassLow[21]       = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  //double cutMetLow[21]        = { 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35};
  //double cutMetLowEM[21]      = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
  double RFixedMM[22]; RFixedMM[21] = 0.32;
  double RFixedEE[22]; RFixedEE[21] = 0.26;
  double RFixedLL[22]; RFixedLL[21] = 0.29;
  if     (nJetsType == 1){
    RFixedMM[21] = 0.27;
    RFixedEE[21] = 0.16;
    RFixedLL[21] = 0.21;
  }
  else if(nJetsType == 1){
    RFixedMM[21] = 0.34;
    RFixedEE[21] = 0.24;
    RFixedLL[21] = 0.29;
  }
  const int nbins= 5;
  const float bins [nbins+1] = {20, 23, 26, 31, 40, 50};
  char sb[200];
  char proc[200];
  TH1D* hDMET[20];
  for(int i=0; i<20; i++){
    if     (i%8==0) sprintf(proc,"mc_mm");
    else if(i%8==1) sprintf(proc,"data_mm");
    else if(i%8==2) sprintf(proc,"mc_ee");
    else if(i%8==3) sprintf(proc,"data_ee");
    else if(i%8==4) sprintf(proc,"mc_em");
    else if(i%8==5) sprintf(proc,"data_em");
    else if(i%8==6) sprintf(proc,"mc_ll");
    else if(i%8==7) sprintf(proc,"data_ll");
    sprintf(sb,"hDMET_%s_%d",proc,i); hDMET[i] = new TH1D(sb,sb, nbins, bins);
    hDMET[i]->Sumw2();
  }

  double btag_lowpt_1j_den[50],btag_lowpt_1j_num[50],btag_lowpt_0j_den[50],btag_lowpt_0j_num[50];
  double btag_highestpt_2j_den[50],btag_highestpt_2j_num[50],btag_highestpt_1j_den[50],btag_highestpt_1j_num[50];
  double zll_ee[50],zll_mm[50];
  for(int i=0; i<50; i++){
    btag_lowpt_1j_den[i] = 0.0;
    btag_lowpt_1j_num[i] = 0.0;
    btag_lowpt_0j_den[i] = 0.0;
    btag_lowpt_0j_num[i] = 0.0;
    btag_highestpt_2j_den[i] = 0.0;
    btag_highestpt_2j_num[i] = 0.0;
    btag_highestpt_1j_den[i] = 0.0;
    btag_highestpt_1j_num[i] = 0.0;;
    zll_ee[i] = 0.0;
    zll_mm[i] = 0.0;
  }
  double nZ_METCut[20]; for(int i=0; i<20; i++) nZ_METCut[i] = 0.0;
  double nZ_METMTCut[20]; for(int i=0; i<20; i++) nZ_METMTCut[i] = 0.0;

  unsigned int patternTopVeto         = SmurfTree::TopVeto;
  unsigned int patternTopTagNotInJets = SmurfTree::TopTagNotInJets;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%10000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);
    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::tw    ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz    ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data  ) fDecay =  1;
    else                                          {fDecay = 0;cout << bgdEvent.dstype_ << endl;assert(0);}
    if(fDecay == -1 || fDecay > 100) fDecay = 0;//44;
    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
    if((fDecay == 27 || fDecay == 28) && channel != 128) {
      if(bgdEvent.lep1MotherMcId_ == 23 && bgdEvent.lep2MotherMcId_ == 23) {
        fDecay = 31;
      }
    }
    int Njet3 = 0;
    if(bgdEvent.jet3_.Pt() <= 30)									 Njet3 = 0;
    else if(bgdEvent.jet3_.Pt() > 30 && (
      (bgdEvent.jet1_.Eta()-bgdEvent.jet3_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.jet3_.Eta() < 0) ||
      (bgdEvent.jet2_.Eta()-bgdEvent.jet3_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.jet3_.Eta() < 0))) Njet3 = 2;
    else												 Njet3 = 1;
    int centrality = 0;
    if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
       ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    double newMet      = bgdEvent.met_;
    double newTrackMet = bgdEvent.trackMet_;
    bool applyCor = false;
    if(applyCor == true && bgdEvent.dstype_ != SmurfTree::data){
      double metx=0.0;double mety=0.0;double trkmetx=0.0;double trkmety=0.0;
      if      (bgdEvent.njets_ == 0){
        metx	= newMet*cos(bgdEvent.metPhi_)+gRandom->Gaus(0.0,4.8);
        mety	= newMet*sin(bgdEvent.metPhi_)+gRandom->Gaus(0.0,4.8);
        trkmetx = newTrackMet*cos(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,1.4);
        trkmety = newTrackMet*sin(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,1.4);
      }
      else if(bgdEvent.njets_ == 1){
        metx	= newMet*cos(bgdEvent.metPhi_)+gRandom->Gaus(0.0,4.9);
        mety	= newMet*sin(bgdEvent.metPhi_)+gRandom->Gaus(0.0,4.9);
        trkmetx = newTrackMet*cos(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.4);
        trkmety = newTrackMet*sin(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.4);
      }
      else if(bgdEvent.njets_ >= 2){
        metx	= newMet*cos(bgdEvent.metPhi_)+gRandom->Gaus(0.0,5.0);
        mety	= newMet*sin(bgdEvent.metPhi_)+gRandom->Gaus(0.0,5.0);
        trkmetx = newTrackMet*cos(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.8);
        trkmety = newTrackMet*sin(bgdEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.8);
      }
      newMet      = sqrt(metx*metx+mety*mety);
      newTrackMet = sqrt(trkmetx*trkmetx+trkmety*trkmety);
    }
    double deltaPhiA[3] = {TMath::Abs(bgdEvent.lep1_.Phi()-bgdEvent.metPhi_),TMath::Abs(bgdEvent.lep2_.Phi()-bgdEvent.metPhi_),0.0};
    while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
    while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
    deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
    double pmetA = newMet;
    if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

    double deltaPhiB[3] = {TMath::Abs(bgdEvent.lep1_.Phi()-bgdEvent.trackMetPhi_),TMath::Abs(bgdEvent.lep2_.Phi()-bgdEvent.trackMetPhi_),0.0};
    while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
    while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
    deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
    double pmetB = newTrackMet;
    if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

    double deltaPhiC[3] = {TMath::Abs(bgdEvent.jet1_.Phi()-bgdEvent.metPhi_),TMath::Abs(bgdEvent.jet2_.Phi()-bgdEvent.metPhi_),0.0};
    while(deltaPhiC[0]>TMath::Pi()) deltaPhiC[0] = TMath::Abs(deltaPhiC[0] - 2*TMath::Pi());
    while(deltaPhiC[1]>TMath::Pi()) deltaPhiC[1] = TMath::Abs(deltaPhiC[1] - 2*TMath::Pi());
    double pmetC = newMet;

    double deltaPhiD[3] = {TMath::Abs(bgdEvent.jet1_.Phi()-bgdEvent.trackMetPhi_),TMath::Abs(bgdEvent.jet2_.Phi()-bgdEvent.trackMetPhi_),0.0};
    while(deltaPhiD[0]>TMath::Pi()) deltaPhiD[0] = TMath::Abs(deltaPhiD[0] - 2*TMath::Pi());
    while(deltaPhiD[1]>TMath::Pi()) deltaPhiD[1] = TMath::Abs(deltaPhiD[1] - 2*TMath::Pi());
    double pmetD = newTrackMet;

    if     (bgdEvent.njets_ == 0){
      pmetC = pmetA;
      pmetD = pmetB;
    }
    else if(bgdEvent.njets_ == 1){
      deltaPhiC[2] = deltaPhiC[0];
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = deltaPhiD[0];
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }
    else if(bgdEvent.njets_ >= 2){
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiC[1]);
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiD[1]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }

    double usedMet = 0.0;
    bool passMET = false;
    usedMet = TMath::Min(pmetA,pmetB);
    if     (bgdEvent.njets_ == 0){
       passMET = usedMet > 20. &&
                (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }
    else if(bgdEvent.njets_ == 1){
       passMET = usedMet > 20. &&
                (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }
    else if(bgdEvent.njets_ >= 2){
       passMET = usedMet > 20. &&
                (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }
    Float_t mTWMin   = bgdEvent.mt1_;
    Float_t mTWMax   = bgdEvent.mt1_;
    if(mTWMin > bgdEvent.mt2_                        ) mTWMin = bgdEvent.mt2_;
    if(mTWMin > bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMin = bgdEvent.mt3_;
    if(mTWMax < bgdEvent.mt2_	 		     ) mTWMax = bgdEvent.mt2_;
    if(mTWMax < bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMax = bgdEvent.mt3_; 
    double pxzll = bgdEvent.dilep_.Px() + newMet*cos( bgdEvent.metPhi_);
    double pyzll = bgdEvent.dilep_.Py() + newMet*sin( bgdEvent.metPhi_);
    double mtHZZ = TMath::Power(sqrt(bgdEvent.dilep_.Pt()*bgdEvent.dilep_.Pt()+91.1876*91.1876)+
                                sqrt(newMet       *newMet       +91.1876*91.1876),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;

    double mtHWW = TMath::Power(sqrt(bgdEvent.dilep_.Pt()*bgdEvent.dilep_.Pt()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+
                                sqrt(newMet       *newMet       +bgdEvent.dilep_.M()*bgdEvent.dilep_.M()),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHWW >0) mtHWW = sqrt(mtHWW); else mtHWW = 0.0;

    double pznunu = 91.1876*91.1876 - newMet*newMet;
    if(pznunu > 0) pznunu = sqrt(pznunu); else pznunu = 0.0;
    double mHZZ = (sqrt(bgdEvent.dilep_.P()*bgdEvent.dilep_.P()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+91.1876)*
                  (sqrt(bgdEvent.dilep_.P()*bgdEvent.dilep_.P()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+91.1876)-
                 -pxzll*pxzll-pyzll*pyzll-(bgdEvent.dilep_.Pz()+pznunu)*(bgdEvent.dilep_.Pz()+pznunu);
    if(mHZZ > 0) mHZZ = sqrt(mHZZ); else mHZZ = 0.001;

    double massZMin = 999.0;
    double massMin = 999.0;
    double dRMin = 999.0; double mTW3 = 0.0; double type3l = 0.0;
    if(bgdEvent.lid3_ != 0){
      usedMet = TMath::Min(newMet,newTrackMet);
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
      mTW3  = trilepton_info(3,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      type3l= trilepton_info(4,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
    }    

    // begin computing weights
    double theWeight = 0.0;
    double add       = 1.0;
    int nFake = 0;
    if(bgdEvent.dstype_ == SmurfTree::data ){
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
    } else {
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV1)  == SmurfTree::Lep1LooseMuV1)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV1)  == SmurfTree::Lep2LooseMuV1)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV1)  == SmurfTree::Lep3LooseMuV1)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
    }
    if(nFake > 1){
      theWeight = 0.0;
    }
    else if(nFake == 1){
      if(bgdEvent.dstype_ == SmurfTree::data){
    	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        										(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        										(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        theWeight	       = add;
      }
      else if(TMath::Abs(bgdEvent.lep1McId_)*TMath::Abs(bgdEvent.lep2McId_) > 0 || bgdEvent.dstype_ == SmurfTree::wgamma){
    	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV1)  == SmurfTree::Lep1LooseMuV1  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        										(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV1)  == SmurfTree::Lep2LooseMuV1  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        										(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
    	add = add*scaleFactor(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.nvtx_, 2);
        fDecay  	       = 3;
        theWeight	       = -1.0 * bgdEvent.scale1fb_*lumi*add;
      }
    }
    else if(bgdEvent.dstype_ != SmurfTree::data){
      add = add*scaleFactor(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.nvtx_, 2);
      if(fDecay == 27 || fDecay == 28 || fDecay == 19 || fDecay == 29 || fDecay == 30 || fDecay == 31){
    	if((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
        if((bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
          add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
    	double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
        							 fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        							 TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
    	add = add*trigEff;
      }
      if(channel == 5){
        if((fDecay == 9) && bgdEvent.type_ == SmurfTree::mm) {
          if(bgdEvent.njets_ == 0) add=add*4.215;
          if(bgdEvent.njets_ == 1) add=add*2.570;
          if(bgdEvent.njets_ >= 2) add=add*3.855;
        }
        if((fDecay == 9) && bgdEvent.type_ == SmurfTree::ee) {
          if(bgdEvent.njets_ == 0) add=add*4.215;
          if(bgdEvent.njets_ == 1) add=add*2.570;
          if(bgdEvent.njets_ >= 2) add=add*3.855;
        }
      }
      if(channel == 6 || (channel >= 100 && channel <= 800) || channel == 900 || channel == 901){
        if((fDecay == 5 ||fDecay == 13)) {
          if(bgdEvent.njets_ == 0) add=add*1.74;
          if(bgdEvent.njets_ == 1) add=add*1.38; 
          if(bgdEvent.njets_ >= 2) add=add*1.00; 
        }
      }
      if(fDecay == 3) {
        add=add*1.95; 
      }
      // CAREFUL HERE, no data-driven corrections
      //add = 1.0;
      theWeight = bgdEvent.scale1fb_*lumi*add;
    }
    if(channel == 5){ // btag study
      if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
         bgdEvent.dstype_ != SmurfTree::data) &&
         bgdEvent.dilep_.M()   > 12 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
        (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
	(bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
	 1 == 1
	){
	int classType = 0;
	if(fDecay == 5 ) classType = 1;
	if(fDecay == 13) classType = 2;
	if(bgdEvent.jet1Btag_ >= 2.1 && bgdEvent.njets_ == 1){
	  btag_lowpt_1j_den[classType+10*4]		 += theWeight;
	  btag_lowpt_1j_den[classType+10*bgdEvent.type_] += theWeight;
	  if((bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
	    btag_lowpt_1j_num[classType+10*4]		   += theWeight;
	    btag_lowpt_1j_num[classType+10*bgdEvent.type_] += theWeight;
	  }
        }
	if(bgdEvent.njets_ == 0){
	  btag_lowpt_0j_den[classType+10*4]		 += theWeight;
	  btag_lowpt_0j_den[classType+10*bgdEvent.type_] += theWeight;
	  if((bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
	    btag_lowpt_0j_num[classType+10*4]		   += theWeight;
	    btag_lowpt_0j_num[classType+10*bgdEvent.type_] += theWeight;
	  }
        }
	if(bgdEvent.jet2Btag_ >= 2.1 && bgdEvent.njets_ == 2){
	  btag_highestpt_2j_den[classType+10*4] 	     += theWeight;
	  btag_highestpt_2j_den[classType+10*bgdEvent.type_] += theWeight;
	  if(bgdEvent.jet1Btag_ >= 2.1){
	    btag_highestpt_2j_num[classType+10*4]	       += theWeight;
	    btag_highestpt_2j_num[classType+10*bgdEvent.type_] += theWeight;
	  }
        }
	if((bgdEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && bgdEvent.njets_ == 1){
	  btag_highestpt_1j_den[classType+10*4] 	     += theWeight;
	  btag_highestpt_1j_den[classType+10*bgdEvent.type_] += theWeight;
	  if(bgdEvent.jet1Btag_ >= 2.1){
	    btag_highestpt_1j_num[classType+10*4]	       += theWeight;
	    btag_highestpt_1j_num[classType+10*bgdEvent.type_] += theWeight;
	  }
        }
      }
    } // btag study
    // Z counting
    if(
     (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
        bgdEvent.dstype_ != SmurfTree::data) &&
      (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
       bgdEvent.lid3_ == 0 &&
       charge == 0 &&
       bgdEvent.lep1_.Pt() > 20. &&
       bgdEvent.lep2_.Pt() > 10. &&
       bgdEvent.njets_ == nJetsType &&
      (fabs(bgdEvent.dilep_.M()-91.1876) < 15.) &&
     !(bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
       usedMet > 20 &&
       1 == 1
      ){
      int classType = 0;
      if(fDecay == 9) classType = 1;
      if     (bgdEvent.type_ == SmurfTree::ee) zll_ee[classType] += theWeight;
      else if(bgdEvent.type_ == SmurfTree::mm) zll_mm[classType] += theWeight;
    }
    if(channel == 6 || (channel >= 100 && channel <= 800)){ // DY study
      if(
      (((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
       ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
         bgdEvent.dstype_ != SmurfTree::data) &&
         bgdEvent.dilep_.M()   > 12 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
	 bgdEvent.njets_ == nJetsType &&
        (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         usedMet > 20. &&
	(bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165.) &&
	(fDecay == 9 || fDecay == 31) &&
	 1 == 1
	){
        int classType = 0;
        if     (bgdEvent.type_ == SmurfTree::mm) classType += 0;
        else if(bgdEvent.type_ == SmurfTree::ee) classType += 2;
        else                                     classType += 4;
	if(fabs(bgdEvent.dilep_.M()-91.1876) < 15.) classType += 10;
	if(channel == 6) {
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
        }
	if((channel >= 100 && channel <= 800) &&
          bgdEvent.dilep_.M() > cutMassLow[binc] &&
          bgdEvent.dilep_.M() < cutMassHigh[binc] &&
	  fabs(bgdEvent.dilep_.M()-91.1876) >= 15. &&
          bgdEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
          bgdEvent.lep2_.Pt() > cutPtMinLow[binc] &&
          bgdEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
	  1 == 1){
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
	  if(usedMet > 40.0 &&
             bgdEvent.mt_  > cutMTLow[binc] &&
             bgdEvent.mt_  < cutMTHigh[binc]) nZ_METMTCut[classType] += theWeight;
        }
        else if((channel >= 100 && channel <= 800) &&
	  fabs(bgdEvent.dilep_.M()-91.1876) < 15. &&
          bgdEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
          bgdEvent.lep2_.Pt() > cutPtMinLow[binc] &&
          bgdEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
	  1 == 1){
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
	  if(usedMet > 40.0 &&
             bgdEvent.mt_  > cutMTLow[binc] &&
             bgdEvent.mt_  < cutMTHigh[binc]) nZ_METMTCut[classType] += theWeight;
        }
      }
    } // DY study
  } // end background loop

  int nData=dataEvent.tree_->GetEntries();
  for (int i=0; i<nData; ++i) {

    if (i%10000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);
    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;
    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    int Njet3 = 0;
    if(dataEvent.jet3_.Pt() <= 30)									     Njet3 = 0;
    else if(dataEvent.jet3_.Pt() > 30 && (
      (dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta() < 0) ||
      (dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta() < 0))) Njet3 = 2;
    else												     Njet3 = 1;
    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    double deltaPhiA[3] = {TMath::Abs(dataEvent.lep1_.Phi()-dataEvent.metPhi_),TMath::Abs(dataEvent.lep2_.Phi()-dataEvent.metPhi_),0.0};
    while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
    while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
    deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
    double pmetA = dataEvent.met_;
    if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

    double deltaPhiB[3] = {TMath::Abs(dataEvent.lep1_.Phi()-dataEvent.trackMetPhi_),TMath::Abs(dataEvent.lep2_.Phi()-dataEvent.trackMetPhi_),0.0};
    while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
    while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
    deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
    double pmetB = dataEvent.trackMet_;
    if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

    double deltaPhiC[3] = {TMath::Abs(dataEvent.jet1_.Phi()-dataEvent.metPhi_),TMath::Abs(dataEvent.jet2_.Phi()-dataEvent.metPhi_),0.0};
    while(deltaPhiC[0]>TMath::Pi()) deltaPhiC[0] = TMath::Abs(deltaPhiC[0] - 2*TMath::Pi());
    while(deltaPhiC[1]>TMath::Pi()) deltaPhiC[1] = TMath::Abs(deltaPhiC[1] - 2*TMath::Pi());
    double pmetC = dataEvent.met_;

    double deltaPhiD[3] = {TMath::Abs(dataEvent.jet1_.Phi()-dataEvent.trackMetPhi_),TMath::Abs(dataEvent.jet2_.Phi()-dataEvent.trackMetPhi_),0.0};
    while(deltaPhiD[0]>TMath::Pi()) deltaPhiD[0] = TMath::Abs(deltaPhiD[0] - 2*TMath::Pi());
    while(deltaPhiD[1]>TMath::Pi()) deltaPhiD[1] = TMath::Abs(deltaPhiD[1] - 2*TMath::Pi());
    double pmetD = dataEvent.trackMet_;

    if     (dataEvent.njets_ == 0){
      pmetC = pmetA;
      pmetD = pmetB;
    }
    else if(dataEvent.njets_ == 1){
      deltaPhiC[2] = deltaPhiC[0];
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = deltaPhiD[0];
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }
    else if(dataEvent.njets_ >= 2){
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiC[1]);
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiD[1]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }

    double usedMet = 0.0;
    bool passMET = false;
    usedMet = TMath::Min(pmetA,pmetB);
    if     (dataEvent.njets_ == 0){
       passMET = usedMet > 20. &&
                (usedMet > 40. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    else if(dataEvent.njets_ == 1){
       passMET = usedMet > 20. &&
                (usedMet > 40. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    else if(dataEvent.njets_ >= 2){
       passMET = usedMet > 20. &&
                (usedMet > 40. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    double massZMin = 999.0;
    double massMin = 999.0;
    double dRMin = 999.0; double mTW3 = 0.0; double type3l = 0.0;
    Float_t mTWMin   = dataEvent.mt1_;
    Float_t mTWMax   = dataEvent.mt1_;
    if(mTWMin > dataEvent.mt2_  		       ) mTWMin = dataEvent.mt2_;
    if(mTWMin > dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMin = dataEvent.mt3_;
    if(mTWMax < dataEvent.mt2_  		       ) mTWMax = dataEvent.mt2_;
    if(mTWMax < dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMax = dataEvent.mt3_;
    double pxzll = dataEvent.dilep_.Px() + dataEvent.met_*cos( dataEvent.metPhi_);
    double pyzll = dataEvent.dilep_.Py() + dataEvent.met_*sin( dataEvent.metPhi_);
    double mtHZZ = TMath::Power(sqrt(dataEvent.dilep_.Pt()*dataEvent.dilep_.Pt()+91.1876*91.1876)+
                              sqrt(dataEvent.met_       *dataEvent.met_       +91.1876*91.1876),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;
    double mtHWW = TMath::Power(sqrt(dataEvent.dilep_.Pt()*dataEvent.dilep_.Pt()+dataEvent.dilep_.M()*dataEvent.dilep_.M())+
                                sqrt(dataEvent.met_       *dataEvent.met_       +dataEvent.dilep_.M()*dataEvent.dilep_.M()),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHWW >0) mtHWW = sqrt(mtHWW); else mtHWW = 0.0;
    double pznunu = 91.1876*91.1876 - dataEvent.met_*dataEvent.met_;
    if(pznunu > 0) pznunu = sqrt(pznunu); else pznunu = 0.0;
    double mHZZ = (sqrt(dataEvent.dilep_.P()*dataEvent.dilep_.P()+dataEvent.dilep_.M()*dataEvent.dilep_.M())+91.1876)*
                  (sqrt(dataEvent.dilep_.P()*dataEvent.dilep_.P()+dataEvent.dilep_.M()*dataEvent.dilep_.M())+91.1876)-
                 -pxzll*pxzll-pyzll*pyzll-(dataEvent.dilep_.Pz()+pznunu)*(dataEvent.dilep_.Pz()+pznunu);
    if(mHZZ > 0) mHZZ = sqrt(mHZZ); else mHZZ = 0.001;

    if(dataEvent.lid3_ != 0){
      usedMet = TMath::Min(dataEvent.met_,dataEvent.trackMet_);
      massZMin = trilepton_info(0,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      massMin  = trilepton_info(1,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      dRMin = trilepton_info(2,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      mTW3  = trilepton_info(3,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      type3l= trilepton_info(4,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
    }
    double theWeight = 1.0;
    if(channel == 5){ // btag study
      if(
	(dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
	(dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
         dataEvent.dilep_.M()   > 12 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         passMET == true &&
        (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
	(dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
	 1 == 1
	){
	int classType = 3;
	if(dataEvent.jet1Btag_ >= 2.1 && dataEvent.njets_ == 1){
	  btag_lowpt_1j_den[classType+10*4]		  += theWeight;
	  btag_lowpt_1j_den[classType+10*dataEvent.type_] += theWeight;
	  if((dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
	    btag_lowpt_1j_num[classType+10*4]		    += theWeight;
	    btag_lowpt_1j_num[classType+10*dataEvent.type_] += theWeight;
	  }
        }
	if(dataEvent.njets_ == 0){
	  btag_lowpt_0j_den[classType+10*4]		  += theWeight;
	  btag_lowpt_0j_den[classType+10*dataEvent.type_] += theWeight;
	  if((dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets){
	    btag_lowpt_0j_num[classType+10*4]		    += theWeight;
	    btag_lowpt_0j_num[classType+10*dataEvent.type_] += theWeight;
	  }
        }
	if(dataEvent.jet2Btag_ >= 2.1 && dataEvent.njets_ == 2){
	  btag_highestpt_2j_den[classType+10*4] 	      += theWeight;
	  btag_highestpt_2j_den[classType+10*dataEvent.type_] += theWeight;
	  if(dataEvent.jet1Btag_ >= 2.1){
	    btag_highestpt_2j_num[classType+10*4]	        += theWeight;
	    btag_highestpt_2j_num[classType+10*dataEvent.type_] += theWeight;
	  }
        }
	if((dataEvent.cuts_ & patternTopTagNotInJets) != patternTopTagNotInJets && dataEvent.njets_ == 1){
	  btag_highestpt_1j_den[classType+10*4] 	      += theWeight;
	  btag_highestpt_1j_den[classType+10*dataEvent.type_] += theWeight;
	  if(dataEvent.jet1Btag_ >= 2.1){
	    btag_highestpt_1j_num[classType+10*4]	        += theWeight;
	    btag_highestpt_1j_num[classType+10*dataEvent.type_] += theWeight;
	  }
        }
      }
    } // btag study
    // Z counting
    if(
      (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
      (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
      (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
       dataEvent.lid3_ == 0 &&
       charge == 0 &&
       dataEvent.lep1_.Pt() > 20. &&
       dataEvent.lep2_.Pt() > 10. &&
       dataEvent.njets_ == nJetsType &&
      (fabs(dataEvent.dilep_.M()-91.1876) < 15.) &&
     !(dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
       usedMet > 20 &&
       1 == 1
      ){
      int classType = 2;
      if     (dataEvent.type_ == SmurfTree::ee) zll_ee[classType] += theWeight;
      else if(dataEvent.type_ == SmurfTree::mm) zll_mm[classType] += theWeight;
    }
    if(channel == 6 || (channel >= 100 && channel <= 800)){ // DY study
      if(
        (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
        (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
         dataEvent.dilep_.M()   > 12 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
	 dataEvent.njets_ == nJetsType &&
        (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         usedMet > 20. &&
	(dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165.) &&
	 1 == 1
	){
        int classType = 1;
        if     (dataEvent.type_ == SmurfTree::mm) classType += 0;
        else if(dataEvent.type_ == SmurfTree::ee) classType += 2;
        else                                      classType += 4;
	if(fabs(dataEvent.dilep_.M()-91.1876) < 15.) classType += 10;
	if(channel == 6) {
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
        }
        if((channel >= 100 && channel <= 800) &&
          dataEvent.dilep_.M() > cutMassLow[binc] &&
          dataEvent.dilep_.M() < cutMassHigh[binc] &&
	  fabs(dataEvent.dilep_.M()-91.1876) >= 15. &&
          dataEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
          dataEvent.lep2_.Pt() > cutPtMinLow[binc] &&
          dataEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
	  1 == 1){
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
	  if(usedMet > 40.0 &&
             dataEvent.mt_  > cutMTLow[binc] &&
             dataEvent.mt_  < cutMTHigh[binc]) nZ_METMTCut[classType] += theWeight;
        }
        else if((channel >= 100 && channel <= 800) &&
	  fabs(dataEvent.dilep_.M()-91.1876) < 15. &&
          dataEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
          dataEvent.lep2_.Pt() > cutPtMinLow[binc] &&
          dataEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
	  1 == 1){
	  hDMET[classType]->Fill(TMath::Min(usedMet,49.999), theWeight);
	  if(usedMet > 40.0) nZ_METCut[classType] += theWeight;
	  if(usedMet > 40.0 &&
             dataEvent.mt_  > cutMTLow[binc] &&
             dataEvent.mt_  < cutMTHigh[binc]) nZ_METMTCut[classType] += theWeight;
        }
      }
    } // DY study
  } // End loop data
  
  double k_em[2] = {0.0, 0.0};
  if(zll_mm[1] != 0){
    k_em[0] = zll_ee[1]/zll_mm[1];
    k_em[1] = (zll_ee[2]-zll_ee[0])/(zll_mm[2]-zll_mm[0]);
    for(int i=0; i<2; i++) if(k_em[i] > 0) k_em[i] = sqrt(k_em[i]); else k_em[i] = 0.0;
    printf("**********k ee/mm**********\n");
    printf("k ee/mm ratio: %6.3f - %6.3f\n",k_em[0],k_em[1]);
  }
  if(channel == 5){
    printf("**********eff highest pt jet 2-j**********\n");
    double effttMC_btag_highestpt_2j[5],effttDA_btag_highestpt_2j[5],effttDA_btag_highestpt_2j_Err[5];
    for(int i=0; i<5; i++) effttMC_btag_highestpt_2j[i] = (btag_highestpt_2j_num[1+i*10]+btag_highestpt_2j_num[2+i*10])/(btag_highestpt_2j_den[1+i*10]+btag_highestpt_2j_den[2+i*10]);
    for(int i=0; i<5; i++) effttDA_btag_highestpt_2j[i] = (btag_highestpt_2j_num[3+i*10]-btag_highestpt_2j_num[0+i*10])/(btag_highestpt_2j_den[3+i*10]-btag_highestpt_2j_den[0+i*10]);
    for(int i=0; i<5; i++) effttDA_btag_highestpt_2j_Err[i] = sqrt((1-effttDA_btag_highestpt_2j[i])*effttDA_btag_highestpt_2j[i]/btag_highestpt_2j_den[3+i*10]);
    for(int i=0; i<5; i++) printf("numerator(%d)   --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_highestpt_2j_num[3+i*10],btag_highestpt_2j_num[0+i*10],(btag_highestpt_2j_num[1+i*10]+btag_highestpt_2j_num[2+i*10]));
    for(int i=0; i<5; i++) printf("denominator(%d) --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_highestpt_2j_den[3+i*10],btag_highestpt_2j_den[0+i*10],(btag_highestpt_2j_den[1+i*10]+btag_highestpt_2j_den[2+i*10]));
    for(int i=0; i<5; i++) printf("eff (%d): %6.3f --> %6.3f +/- %6.3f\n",i,effttMC_btag_highestpt_2j[i],effttDA_btag_highestpt_2j[i],effttDA_btag_highestpt_2j_Err[i]);

    printf("**********eff highest pt jet 1-j**********\n");
    double effttMC_btag_highestpt_1j[5],effttDA_btag_highestpt_1j[5],effttDA_btag_highestpt_1j_Err[5];
    for(int i=0; i<5; i++) effttMC_btag_highestpt_1j[i] = (btag_highestpt_1j_num[1+i*10]+btag_highestpt_1j_num[2+i*10])/(btag_highestpt_1j_den[1+i*10]+btag_highestpt_1j_den[2+i*10]);
    for(int i=0; i<5; i++) effttDA_btag_highestpt_1j[i] = (btag_highestpt_1j_num[3+i*10]-btag_highestpt_1j_num[0+i*10])/(btag_highestpt_1j_den[3+i*10]-btag_highestpt_1j_den[0+i*10]);
    for(int i=0; i<5; i++) effttDA_btag_highestpt_1j_Err[i] = sqrt((1-effttDA_btag_highestpt_1j[i])*effttDA_btag_highestpt_1j[i]/btag_highestpt_1j_den[3+i*10]);
    for(int i=0; i<5; i++) printf("numerator(%d)   --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_highestpt_1j_num[3+i*10],btag_highestpt_1j_num[0+i*10],(btag_highestpt_1j_num[1+i*10]+btag_highestpt_1j_num[2+i*10]));
    for(int i=0; i<5; i++) printf("denominator(%d) --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_highestpt_1j_den[3+i*10],btag_highestpt_1j_den[0+i*10],(btag_highestpt_1j_den[1+i*10]+btag_highestpt_1j_den[2+i*10]));
    for(int i=0; i<5; i++) printf("eff (%d): %6.3f - %6.3f +/- %6.3f (indicative)\n",i,effttMC_btag_highestpt_1j[i],effttDA_btag_highestpt_1j[i],effttDA_btag_highestpt_1j_Err[i]);

    // we use the combined efficiency obtained in the 2-j bin, instead of the obtained final state by final state
    double estimationMC_btag_highestpt_1j[5]; for(int i=0; i<5; i++) estimationMC_btag_highestpt_1j[i] = (1-effttMC_btag_highestpt_2j[4])/effttMC_btag_highestpt_2j[4]*(btag_highestpt_1j_num[1+i*10]+btag_highestpt_1j_num[2+i*10]);
    double estimationDA_btag_highestpt_1j[5]; for(int i=0; i<5; i++) estimationDA_btag_highestpt_1j[i] = (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*(btag_highestpt_1j_num[3+i*10]-btag_highestpt_1j_num[0+i*10]);
    double estimationDA_btag_highestpt_1j_Err[5]; for(int i=0; i<5; i++) estimationDA_btag_highestpt_1j_Err[i] = sqrt(((btag_highestpt_1j_num[3+i*10]-btag_highestpt_1j_num[0+i*10])*effttDA_btag_highestpt_2j_Err[4]/effttDA_btag_highestpt_2j[4]/effttDA_btag_highestpt_2j[4])*
                                                                                                                      ((btag_highestpt_1j_num[3+i*10]-btag_highestpt_1j_num[0+i*10])*effttDA_btag_highestpt_2j_Err[4]/effttDA_btag_highestpt_2j[4]/effttDA_btag_highestpt_2j[4])+
										                                      (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*
														      (1-effttDA_btag_highestpt_2j[4])/effttDA_btag_highestpt_2j[4]*btag_highestpt_1j_num[3+i*10]);
    for(int i=0; i<5; i++) printf("top 1-jet(%d): %6.3f  %6.3f  %6.3f +/- %6.3f\n",i,(btag_highestpt_1j_den[1+i*10]+btag_highestpt_1j_den[2+i*10])-(btag_highestpt_1j_num[1+i*10]+btag_highestpt_1j_num[2+i*10]),
                                                                                   estimationMC_btag_highestpt_1j[i],estimationDA_btag_highestpt_1j[i],estimationDA_btag_highestpt_1j_Err[i]);

    printf("**********eff low pt jet 1-j**********\n");
    double effttMC_btag_lowpt_1j[5],effttDA_btag_lowpt_1j[5],effttDA_btag_lowpt_1j_Err[5];
    for(int i=0; i<5; i++) effttMC_btag_lowpt_1j[i]     = (btag_lowpt_1j_num[1+10*i]+btag_lowpt_1j_num[2+10*i])/ (btag_lowpt_1j_den[1+10*i]+btag_lowpt_1j_den[2+10*i]);
    for(int i=0; i<5; i++) effttDA_btag_lowpt_1j[i]     = (btag_lowpt_1j_num[3+10*i]-btag_lowpt_1j_num[0+10*i])/ (btag_lowpt_1j_den[3+10*i]-btag_lowpt_1j_den[0+10*i]);
    for(int i=0; i<5; i++) effttDA_btag_lowpt_1j_Err[i] = sqrt((1-effttDA_btag_lowpt_1j[i])*effttDA_btag_lowpt_1j[i]/btag_lowpt_1j_den[3+10*i]);
    for(int i=0; i<5; i++) printf("numerator(%d)   --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_lowpt_1j_num[3+i*10],btag_lowpt_1j_num[0+i*10],(btag_lowpt_1j_num[1+i*10]+btag_lowpt_1j_num[2+i*10]));
    for(int i=0; i<5; i++) printf("denominator(%d) --> data: %4.0f, background: %7.3f, signal: %7.3f\n",i,btag_lowpt_1j_den[3+i*10],btag_lowpt_1j_den[0+i*10],(btag_lowpt_1j_den[1+i*10]+btag_lowpt_1j_den[2+i*10]));
    for(int i=0; i<5; i++) printf("eff (%d): %6.3f --> %6.3f +/- %6.3f\n",i,effttMC_btag_lowpt_1j[i],effttDA_btag_lowpt_1j[i],effttDA_btag_lowpt_1j_Err[i]);

    printf("**********eff low pt jet 0-j**********\n");
    double N_top_expected_0j[5]; for(int i=0; i<5; i++) N_top_expected_0j[i] = (btag_lowpt_0j_den[1+10*i]+btag_lowpt_0j_den[2+10*i])-(btag_lowpt_0j_num[1+10*i]+btag_lowpt_0j_num[2+10*i]);
    double fttbar[5]; for(int i=0; i<5; i++) fttbar[i] = btag_lowpt_0j_den[1+10*i]/(btag_lowpt_0j_den[1+10*i]+btag_lowpt_0j_den[2+10*i]);
   
    double effMC_btag_lowpt_0j[5]; for(int i=0; i<5; i++) effMC_btag_lowpt_0j[i] = (fttbar[i]*(1-(1-effttMC_btag_lowpt_1j[i])*(1-effttMC_btag_lowpt_1j[i]))+(1-fttbar[i])*effttMC_btag_lowpt_1j[i]);
    double effDA_btag_lowpt_0j[5]; for(int i=0; i<5; i++) effDA_btag_lowpt_0j[i] = (fttbar[i]*(1-(1-effttDA_btag_lowpt_1j[i])*(1-effttDA_btag_lowpt_1j[i]))+(1-fttbar[i])*effttDA_btag_lowpt_1j[i]);
    double sigma_ftop[2]={0.00,0.17};
    double effttMC_btag_lowpt_1j_Err[5] = {0.01,0.01,0.01,0.01,0.01};
    double effMC_btag_lowpt_0j_Err[5];
    double effDA_btag_lowpt_0j_Err[5];
    for(int i=0; i<5; i++){
      effMC_btag_lowpt_0j_Err[i] = sqrt(sigma_ftop[0]*sigma_ftop[0]*(effMC_btag_lowpt_0j[i]*(1-effMC_btag_lowpt_0j[i]))*(effMC_btag_lowpt_0j[i]*(1-effMC_btag_lowpt_0j[i]))+
    					effttMC_btag_lowpt_1j_Err[i]*effttMC_btag_lowpt_1j_Err[i]*(fttbar[i]*(1-2*effMC_btag_lowpt_0j[i])+1)*(fttbar[i]*(1-2*effMC_btag_lowpt_0j[i])+1));
      effDA_btag_lowpt_0j_Err[i] = sqrt(sigma_ftop[1]*sigma_ftop[1]*(effDA_btag_lowpt_0j[i]*(1-effDA_btag_lowpt_0j[i]))*(effDA_btag_lowpt_0j[i]*(1-effDA_btag_lowpt_0j[i]))+
    					effttDA_btag_lowpt_1j_Err[i]*effttDA_btag_lowpt_1j_Err[i]*(fttbar[i]*(1-2*effDA_btag_lowpt_0j[i])+1)*(fttbar[i]*(1-2*effDA_btag_lowpt_0j[i])+1));
    }
    for(int i=0; i<5; i++) printf("fttbar= %5.3f, eff btag lowpt 0j(%d): %6.3f, %6.3f +/- %6.3f --> %6.3f +/- %6.3f\n",fttbar[i],i
                                                                                                       ,(btag_lowpt_0j_num[1+10*i]+btag_lowpt_0j_num[2+10*i])/(btag_lowpt_0j_den[1+10*i]+btag_lowpt_0j_den[2+10*i])
                                                                                                       ,effMC_btag_lowpt_0j[i],effMC_btag_lowpt_0j_Err[i],effDA_btag_lowpt_0j[i],effDA_btag_lowpt_0j_Err[i]);

    double sigma_0f_bck = 0.20;
    // we use the combined efficiency obtained in the 1-j bin, instead of the obtained final state by final state
    double estimationMC_btag_lowpt_0j[5]; for(int i=0; i<5; i++) estimationMC_btag_lowpt_0j[i] = (1-effMC_btag_lowpt_0j[4])/effMC_btag_lowpt_0j[4]*(btag_lowpt_0j_num[1+i*10]+btag_lowpt_0j_num[2+i*10]);
    double estimationDA_btag_lowpt_0j[5]; for(int i=0; i<5; i++) estimationDA_btag_lowpt_0j[i] = (1-effDA_btag_lowpt_0j[4])/effDA_btag_lowpt_0j[4]*(btag_lowpt_0j_num[3+i*10]-btag_lowpt_0j_num[0+i*10]);
    
    double estimationMC_btag_lowpt_0j_Err[5]; for(int i=0; i<5; i++) estimationMC_btag_lowpt_0j_Err[i] = sqrt(((btag_lowpt_0j_num[1+i*10]+btag_lowpt_0j_num[2+i*10])*effMC_btag_lowpt_0j_Err[4]/effMC_btag_lowpt_0j[4]/effMC_btag_lowpt_0j[4])*
                                                                                                              ((btag_lowpt_0j_num[1+i*10]+btag_lowpt_0j_num[2+i*10])*effMC_btag_lowpt_0j_Err[4]/effMC_btag_lowpt_0j[4]/effMC_btag_lowpt_0j[4]));
    double estimationDA_btag_lowpt_0j_Err[5]; for(int i=0; i<5; i++) estimationDA_btag_lowpt_0j_Err[i] = sqrt(((btag_lowpt_0j_num[3+i*10]-btag_lowpt_0j_num[0+i*10])*effDA_btag_lowpt_0j_Err[4]/effDA_btag_lowpt_0j[4]/effDA_btag_lowpt_0j[4])*
                                                                                                              ((btag_lowpt_0j_num[3+i*10]-btag_lowpt_0j_num[0+i*10])*effDA_btag_lowpt_0j_Err[4]/effDA_btag_lowpt_0j[4]/effDA_btag_lowpt_0j[4])+
										                              (1-estimationDA_btag_lowpt_0j[i])/estimationDA_btag_lowpt_0j[i]*
													      (1-estimationDA_btag_lowpt_0j[i])/estimationDA_btag_lowpt_0j[i]*btag_lowpt_0j_num[3+i*10]+
													      TMath::Power(sigma_0f_bck*btag_lowpt_0j_num[0+i*10]*(1-effMC_btag_lowpt_0j[4])/effMC_btag_lowpt_0j[4],2));
    for(int i=0; i<5; i++) printf("top 0-jet(%d), data(%2d), bck(%6.3f): %6.3f - %6.3f +/- %6.3f -  %6.3f +/- %6.3f\n",i,
                                 (int)btag_lowpt_0j_num[3+i*10],btag_lowpt_0j_num[0+i*10],N_top_expected_0j[i],
			          estimationMC_btag_lowpt_0j[i],estimationMC_btag_lowpt_0j_Err[i],
			          estimationDA_btag_lowpt_0j[i],estimationDA_btag_lowpt_0j_Err[i]);
  }
  if(channel == 6 || (channel >= 100 && channel <= 800)){
    hDMET[ 6]->Add( hDMET[ 0],1.0);
    hDMET[ 6]->Add( hDMET[ 2],1.0);
    hDMET[ 6]->Add( hDMET[ 4],-1*0.5*(1.0/k_em[0]+k_em[0]));
    hDMET[16]->Add( hDMET[10],1.0);
    hDMET[16]->Add( hDMET[12],1.0);
    hDMET[16]->Add( hDMET[14],-1*0.5*(1.0/k_em[0]+k_em[0]));
    hDMET[ 7]->Add( hDMET[ 1],1.0);
    hDMET[ 7]->Add( hDMET[ 3],1.0);
    hDMET[ 7]->Add( hDMET[ 5],-1*0.5*(1.0/k_em[1]+k_em[1]));
    hDMET[17]->Add( hDMET[11],1.0);
    hDMET[17]->Add( hDMET[13],1.0);
    hDMET[17]->Add( hDMET[15],-1*0.5*(1.0/k_em[1]+k_em[1]));
    for(int i=1; i<=hDMET[ 0]->GetNbinsX(); i++){
      printf("Nout(%d) mm/ee/em/ll mc - data: %7.2f  %7.2f  %7.2f  %7.2f - %7.2f  %7.2f  %7.2f  %7.2f\n",i,
                 hDMET[ 0]->GetBinContent(i),hDMET[ 2]->GetBinContent(i),hDMET[ 4]->GetBinContent(i),hDMET[ 6]->GetBinContent(i),
                 hDMET[ 1]->GetBinContent(i),hDMET[ 3]->GetBinContent(i),hDMET[ 5]->GetBinContent(i),hDMET[ 7]->GetBinContent(i));
    }
    for(int i=1; i<=hDMET[ 0]->GetNbinsX(); i++){
      printf("Nin-(%d) mm/ee/em/ll mc - data: %7.2f  %7.2f  %7.2f  %7.2f - %7.2f  %7.2f  %7.2f  %7.2f\n",i,
                 hDMET[10]->GetBinContent(i),hDMET[12]->GetBinContent(i),hDMET[14]->GetBinContent(i),hDMET[16]->GetBinContent(i),
                 hDMET[11]->GetBinContent(i),hDMET[13]->GetBinContent(i),hDMET[15]->GetBinContent(i),hDMET[17]->GetBinContent(i));
    }
    hDMET[ 6]->Divide(hDMET[16]);
    hDMET[ 7]->Divide(hDMET[17]);
    hDMET[ 6]->SetMinimum(0.0); hDMET[ 6]->SetMaximum(2.0);
    hDMET[ 7]->SetMinimum(0.0); hDMET[ 7]->SetMaximum(2.0);

    hDMET[ 0]->Add( hDMET[ 4],-1*0.5/k_em[0]);
    hDMET[ 2]->Add( hDMET[ 4],-1*0.5*k_em[0]);
    hDMET[10]->Add( hDMET[14],-1*0.5/k_em[0]);
    hDMET[12]->Add( hDMET[14],-1*0.5*k_em[0]);
    hDMET[ 1]->Add( hDMET[ 5],-1*0.5/k_em[1]);
    hDMET[ 3]->Add( hDMET[ 5],-1*0.5*k_em[1]);
    hDMET[11]->Add( hDMET[15],-1*0.5/k_em[1]);
    hDMET[13]->Add( hDMET[15],-1*0.5*k_em[1]);
    hDMET[ 0]->Divide(hDMET[10]);
    hDMET[ 2]->Divide(hDMET[12]);
    hDMET[ 1]->Divide(hDMET[11]);
    hDMET[ 3]->Divide(hDMET[13]);
    hDMET[ 0]->SetMinimum(0.0); hDMET[ 0]->SetMaximum(2.0);
    hDMET[ 2]->SetMinimum(0.0); hDMET[ 2]->SetMaximum(2.0);
    hDMET[ 1]->SetMinimum(0.0); hDMET[ 1]->SetMaximum(2.0);
    hDMET[ 3]->SetMinimum(0.0); hDMET[ 3]->SetMaximum(2.0);
    TCanvas* c1 = new TCanvas("c1","c1",100,100,700,500);
    c1->Divide(3,2);
    c1->cd(1);
    hDMET[ 0]->Draw("e");
    c1->cd(2);
    hDMET[ 2]->Draw("e");
    c1->cd(3);
    hDMET[ 6]->Draw("e");
    c1->cd(4);
    hDMET[ 1]->Draw("e");
    c1->cd(5);
    hDMET[ 3]->Draw("e");
    c1->cd(6);
    hDMET[ 7]->Draw("e");
    for(int i=1; i<=hDMET[ 0]->GetNbinsX(); i++){
      printf("R(%d) mm/ee/ll mc - data: %9.3f  %9.3f  %9.3f - %9.3f  %9.3f  %9.3f\n",i,
                 hDMET[ 0]->GetBinContent(i),hDMET[ 2]->GetBinContent(i),hDMET[ 6]->GetBinContent(i),
                 hDMET[ 1]->GetBinContent(i),hDMET[ 3]->GetBinContent(i),hDMET[ 7]->GetBinContent(i));
    }
    printf("MC-out-mm/ee/em, MET>40: %7.3f  %7.3f  %7.3f\n",nZ_METCut[ 0],nZ_METCut[ 2],nZ_METCut[ 4]);
    printf("MC- in-mm/ee/em, MET>40: %7.3f  %7.3f  %7.3f\n",nZ_METCut[10],nZ_METCut[12],nZ_METCut[14]);
    printf("DA-out-mm/ee/em, MET>40: %7.3f  %7.3f  %7.3f\n",nZ_METCut[ 1],nZ_METCut[ 3],nZ_METCut[ 5]);
    printf("DA- in-mm/ee/em, MET>40: %7.3f  %7.3f  %7.3f\n",nZ_METCut[11],nZ_METCut[13],nZ_METCut[15]);    
    printf("MC-out-mm/ee/em, MET-MT: %7.3f  %7.3f  %7.3f\n",nZ_METMTCut[ 0],nZ_METMTCut[ 2],nZ_METMTCut[ 4]);
    printf("MC- in-mm/ee/em, MET-MT: %7.3f  %7.3f  %7.3f\n",nZ_METMTCut[10],nZ_METMTCut[12],nZ_METMTCut[14]);
    printf("DA-out-mm/ee/em, MET-MT: %7.3f  %7.3f  %7.3f\n",nZ_METMTCut[ 1],nZ_METMTCut[ 3],nZ_METMTCut[ 5]);
    printf("DA- in-mm/ee/em, MET-MT: %7.3f  %7.3f  %7.3f\n",nZ_METMTCut[11],nZ_METMTCut[13],nZ_METMTCut[15]);
    printf("DA-EST-mm/ee/em, MET>40: %7.3f  %7.3f  %7.3f\n",RFixedMM[21]*(nZ_METCut[11]-0.5*nZ_METCut[15]/k_em[1]),
                                                            RFixedEE[21]*(nZ_METCut[13]-0.5*nZ_METCut[15]*k_em[1]),
							    RFixedLL[21]*(nZ_METCut[11]+nZ_METCut[13]-0.5*nZ_METCut[15]*(1.0/k_em[1]+k_em[1]))); 
  }
  return;
}
