#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_4_2_2/src/Analysis/nt_scripts/trilepton.h"
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
// optimalCuts_41x
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCuts_41x
(
 int     mH  	 = 29,
 Char_t xTitle[]="myX", Char_t yTitle[]="Fraction",
 int thePlot = 0,
 TString bgdInputFile    = "ntuples_41x/background.root",
 TString signalInputFile = "ntuples_41x/hww160.root",
 TString dataInputFile   = "ntuples_41x/data_2l.root",
 bool fillInfoNote = false,
 double mhAna = 4
 )
{
  unsigned int nJetsType = (int)mhAna/10;
  int lDecay    = (int)mhAna%10;
  printf("nJetsType: %d, lDecay = %d\n",nJetsType,lDecay);
  bool makeGoodPlots = true;
  double lumi = 0.1876;

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/Run2011_Spring11_SmurfV6/lepton_eff/efficiency_results_v6.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/FakeRates_SmurfV6.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open("/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/FakeRates_SmurfV6.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup("/data/smurf/data/Run2011_Spring11_SmurfV6/lepton_eff/efficiency_results_v6.root");

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

  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ggHWW_KFactors_PowhegToHQT.root");
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

  double wwScaleFactor0jCut[21]    = {1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,1};

  double cutMassLow[21]       = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  double cutMetLow[21]        = { 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35};
  double cutMetLowEM[21]      = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};

  int nBin    = 100;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 4; xminPlot = -0.5; xmaxPlot = 3.5;}
  else if(thePlot >= 20 && thePlot <= 22) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 23 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 30 && thePlot <= 33) {nBinPlot = 200; xminPlot = -5.0; xmaxPlot = 15.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  1000.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 38) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  100.0;}
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  5.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  nBin = nBinPlot;

  TH1D* hDHSig[3];
  hDHSig[0] = new TH1D("hDHSig_0", "hDHSig_0", 1, -0.5, 1.5);
  hDHSig[1] = new TH1D("hDHSig_1", "hDHSig_1", 1, -0.5, 1.5);
  hDHSig[2] = new TH1D("hDHSig_2", "hDHSig_2", 1, -0.5, 1.5);
  hDHSig[0]->Sumw2();
  hDHSig[1]->Sumw2();
  hDHSig[2]->Sumw2();
  hDHSig[0]->Scale(0.0);
  hDHSig[1]->Scale(0.0);
  hDHSig[2]->Scale(0.0);
  double S0[nBin],S1[nBin];
  double B0[nBin],B1[nBin];
  for(int i=0; i<nBin; i++){
    S0[i] = 0.0; S1[i] = 0.0;
    B0[i] = 0.0; B1[i] = 0.0;
  }
  TH1D* hDSignif[2];
  hDSignif[0] = new TH1D("hDSignif_0", "hDSignif_0", nBin, -0.5, nBin-0.5);
  hDSignif[1] = new TH1D("hDSignif_1", "hDSignif_1", nBin, -0.5, nBin-0.5);
  TH1D* hDSigOpt = new TH1D("hDSigOpt", "hDSigOpt", nBin, xmin, xmax);
  TH1D* hDBckOpt = new TH1D("hDBckOpt", "hDBckOpt", nBin, xmin, xmax);
  hDSigOpt->Sumw2();
  hDBckOpt->Sumw2();

  TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
  histos->Sumw2();
  TH1D* histo0 = (TH1D*) histos->Clone("histo0");
  TH1D* histo1 = (TH1D*) histos->Clone("histo1");
  TH1D* histo2 = (TH1D*) histos->Clone("histo2");
  TH1D* histo3 = (TH1D*) histos->Clone("histo3");
  TH1D* histo4 = (TH1D*) histos->Clone("histo4");
  TH1D* histo5 = (TH1D*) histos->Clone("histo5");
  histos->Scale(0.0);
  histo0->Scale(0.0);
  histo1->Scale(0.0);
  histo2->Scale(0.0);
  histo3->Scale(0.0);
  histo4->Scale(0.0);
  histo5->Scale(0.0);

  char sb[200];
  TH1D *hDnjets[70];
  TH1D *hDpttag[60], *hDetatag[60];
  if(channel == 5 && thePlot == 101){
    for(int nj=0; nj<=64; nj++){
      sprintf(sb,"hDnjets_%d",nj); hDnjets[nj] = new TH1D(sb,sb, 6, -0.5, 5.5);
      hDnjets[nj]->Sumw2();
    }
    for(int nj=65; nj<=67; nj++){
      sprintf(sb,"hDnjets_%d",nj); hDnjets[nj] = new TH1D(sb,sb, 4, -0.5, 3.5);
      hDnjets[nj]->Sumw2();
    }
    for(int nj=0; nj<60; nj++){
      sprintf(sb,"hDpttag_%d",nj);   hDpttag[nj]   = new TH1D(sb, sb,   10,    0,  200);
      sprintf(sb,"hDetatag_%d",nj);  hDetatag[nj]  = new TH1D(sb, sb,   20,    0,  5.0);
      hDpttag[nj] ->Sumw2();
      hDetatag[nj]->Sumw2();
    }
  }

  double bgdDecayFake[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double bgdDecayFakeE[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  double bgdDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double weiDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  unsigned int patternTopVeto         = SmurfTree::TopVeto;
  unsigned int patternTopTag          = SmurfTree::TopTag;
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
    else                                           fDecay = 0;
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
       //usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }
    else if(bgdEvent.njets_ == 1){
       //usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }
    else if(bgdEvent.njets_ >= 2){
       //usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
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
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 29){ // WW selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
	 bgdEvent.njets_ == nJetsType &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         //usedMet > 35. &&
         //fabs(bgdEvent.dilep_.M()-91.1876) < 15. && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         //(bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets &&
	 //(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.nSoftMuons_ == 0) &&
	 //bgdEvent.jet1Btag_ >= 2.1 &&
	 //bgdEvent.jet2Btag_ >= 2.1 &&
	 (bgdEvent.njets_ != 1 ||
	 (bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) &&
	 //bgdEvent.dilep_.M()   > 100 &&
	 //(bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
	 //bgdEvent.type_ == SmurfTree::mm &&
	 (bgdEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 29 || fDecay == 30 || fDecay == 14) isSignalDecay = true;	     
      }
    } // WW selection

    if(channel == 6){ // Z selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0. &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
        (fabs(bgdEvent.dilep_.M()-91.1876) < 15.) && 
         //bgdEvent.njets_ == 0 &&
	 //fabs(bgdEvent.dilep_.Rapidity()) < 1 &&
	 !(bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
	 //bgdEvent.type_ == SmurfTree::mm &&
	 //bgdEvent.type_ == SmurfTree::ee &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 9) isSignalDecay = true;
      }
    } // Z selection

    if(channel == 5){ // ttbar selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > 20. &&
        (TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > 35. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
	 (bgdEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 5 || fDecay == 13) isSignalDecay = true;
      }
    } // ttbar selection

    if(channel == 28 || channel == 128){ // ZZllnn selection
      if(
        // bgdEvent.lid3_ == 0. &&
        // charge == 0 &&
        // fabs(bgdEvent.dilep_.M()-91.1876) < 15. && 
	//(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee) &&
        // bgdEvent.lep1_.Pt() > 20. &&
        // bgdEvent.lep2_.Pt() > 20. &&
        // bgdEvent.dilep_.Pt() > 55. &&
        //(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
        // newMet > mhAna/3 &&
	// mtHZZ > 4*mhAna/5.+  25. &&
	// mtHZZ < 5*mhAna/5.+  10. &&
        //(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.lid3_ == 0. &&
         charge == 0 &&
         fabs(bgdEvent.dilep_.M()-91.1876) < 15. && 
	(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee) &&
         bgdEvent.lep1_.Pt() > 30. &&
         bgdEvent.lep2_.Pt() > 20. &&
         bgdEvent.dilep_.Pt() > 55. &&
	(bgdEvent.jet1Btag_ < 2 || bgdEvent.jet1_.Pt() < 30) &&
	(bgdEvent.jet2Btag_ < 2 || bgdEvent.jet2_.Pt() < 30) &&
	(bgdEvent.jet3Btag_ < 2 || bgdEvent.jet3_.Pt() < 30) &&
         bgdEvent.met_ > hzz2l_cuts(mhAna,0) &&
         bgdEvent.met_ < hzz2l_cuts(mhAna,1) &&
	 mtHZZ  > hzz2l_cuts(mhAna,2) &&
	 mtHZZ  < hzz2l_cuts(mhAna,3) &&
	 1 == 1
	){
	//
	passCuts = true;
	if(fDecay == 28 && channel != 128) isSignalDecay = true;
      }
    } // ZZllnn selection

    if(channel == 27){ // WZ selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
        (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.jet1_.Pt() <= 40. &&
	 massZMin < 25 &&
	 1 == 1
	){
	passCuts = true;
	if(bgdEvent.dstype_ == SmurfTree::wz) isSignalDecay = true;
      }
    } // WZ selection

    if(channel == 9){ // Ztautau selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.dilep_.M()   < 80 &&
         bgdEvent.lid3_ == 0. &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
	 TMath::Min(bgdEvent.mt1_,bgdEvent.mt2_) < 40. &&
	 (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 9 || fDecay == 8) isSignalDecay = true;
      }
    } // Ztautau selection

    if(channel == 900 || channel == 901){ // qqH selection
      if(
        (bgdEvent.dilep_.M()   < 100 || channel == 901) &&
         bgdEvent.dilep_.M()   > 12  &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
	 bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 && Njet3 != 2 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
        (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
        (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         (bgdEvent.jet1_+bgdEvent.jet2_).M() > 450. &&
         //TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 &&
         bgdEvent.jet1_.Eta()*bgdEvent.jet2_.Eta() < 0 &&
	 centrality == 1 &&
	 (bgdEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // qqH selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection

      double metCut = cutMetLow[binc];
      if(bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) metCut = cutMetLowEM[binc];
      if(
         bgdEvent.mt_	   > cutMTLow[binc] &&
         bgdEvent.mt_	   < cutMTHigh[binc] &&
         bgdEvent.dilep_.M() > cutMassLow[binc] &&
         bgdEvent.dilep_.M() < cutMassHigh[binc] &&
         bgdEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
         bgdEvent.lep2_.Pt() > cutPtMinLow[binc] &&
         bgdEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
	 bgdEvent.njets_ == nJetsType &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (bgdEvent.njets_ != 1 ||
	 (bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) &&
	 //(bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
	 //bgdEvent.type_ == SmurfTree::mm &&
	 (bgdEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // H->WW selection

    if(channel == 1000){ // HW->2l selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0 &&
         charge != 0 &&
	 bgdEvent.njets_ >= 2 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (bgdEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->2l selection

    if(channel == 1100){ // HW->3l selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.jet1_.Pt() <= 40. &&
	 massZMin > 25 &&
	 massMin < 100 &&
	 mTWMax > 100 &&
	 dRMin < 2.0 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection
    
    if(passCuts == true){
      double theWeight = 0.0;
      double add       = 1.0;
      if(bgdEvent.dstype_ == SmurfTree::data || 
      ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1) || ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1) || ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV1) == SmurfTree::Lep3LooseMuV1) || ((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) ||
      ((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2) == SmurfTree::Lep3LooseMuV2) ){
      add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV1) == SmurfTree::Lep1LooseMuV1), ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4));
      add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV1) == SmurfTree::Lep2LooseMuV1), ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4));
      add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2), 0);
      add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2), 0);
        if(bgdEvent.dstype_ == SmurfTree::data){
	  theWeight              = add;
	  bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += theWeight*theWeight;
	}
	else if(TMath::Abs(bgdEvent.lep1McId_)*TMath::Abs(bgdEvent.lep2McId_) > 0 || bgdEvent.dstype_ == SmurfTree::wgamma){
          add = add*scaleFactor(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.nvtx_, 2);
	  fDecay                 = 3;
	  isSignalDecay          = false;
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add;
          bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += TMath::Abs(theWeight)*TMath::Abs(theWeight);
	}
	else {
	  theWeight = 0.0;
	}
      }
      else {
	add = add*scaleFactor(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.nvtx_, 2);
	if(fDecay == 27 || fDecay == 28 || fDecay == 19){
          if((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
	    add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  if((bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
	    add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	                                                           fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
								   TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
          add = add*trigEff;
        }
  	if((channel == 29 || (channel >= 100 && channel <= 800) || channel == 900 || channel == 901 || channel == 5) &&
	   channel != 128){
	  /*
	  if(fDecay == 9) {				    	
  	    if(bgdEvent.njets_ == 1) add=add*0.95; 
  	    if(bgdEvent.njets_ == 2) add=add*1.60; 
  	    if(bgdEvent.njets_ >= 3) add=add*1.90; 
  	  }
	  */
  	  if((fDecay == 9 || fDecay == 31) && bgdEvent.type_ == SmurfTree::mm) {
  	    if(bgdEvent.njets_ == 0) add=add*5.003;
  	    if(bgdEvent.njets_ == 1) add=add*2.429;
  	    if(bgdEvent.njets_ >= 2) add=add*3.477;
  	  }
  	  if((fDecay == 9 || fDecay == 31) && bgdEvent.type_ == SmurfTree::ee) {
  	    if(bgdEvent.njets_ == 0) add=add*3.018;
  	    if(bgdEvent.njets_ == 1) add=add*1.592;
  	    if(bgdEvent.njets_ >= 2) add=add*1.981;
  	  }
  	  if((fDecay == 5 ||fDecay == 13)) {
  	    if(bgdEvent.njets_ == 0) add=add*2.50;
  	    if(bgdEvent.njets_ == 1) add=add*1.16; 
  	    if(bgdEvent.njets_ >= 2) add=add*1.00; 
  	  }
  	  if(fDecay == 3) {
  	    add=add*1.95; 
          }
	  if((fDecay == 29 || fDecay == 30 || fDecay == 14) &&
	     (channel >= 100 && channel <= 800) && bgdEvent.njets_ == 0.){     
  	    add=add*wwScaleFactor0jCut[binc]; 
          }
	}
	// CAREFUL HERE, no data-driven corrections, just Higgs k-factors
	//add = 1.0;
        if (bgdEvent.processId_ == 10010) {
          add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(bgdEvent.higgsPt_));
        }
	theWeight              = bgdEvent.scale1fb_*lumi*add;
        bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += theWeight*theWeight;
      }
      double myVar = newMet;
      if     (thePlot == 1) myVar = bgdEvent.lep1_.Pt();
      else if(thePlot == 2) myVar = bgdEvent.lep2_.Pt();
      else if(thePlot == 3) myVar = bgdEvent.lep3_.Pt();
      else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = bgdEvent.dilep_.M();
      else if(thePlot == 8) myVar = bgdEvent.mt_;
      else if(thePlot == 9) myVar = bgdEvent.mt1_;
      else if(thePlot ==10) myVar = bgdEvent.mt2_;
      else if(thePlot ==11) myVar = bgdEvent.mt3_;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = newMet/bgdEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = bgdEvent.njets_;
      else if(thePlot ==18) myVar = bgdEvent.nvtx_;
      else if(thePlot ==19) myVar = bgdEvent.nSoftMuons_;
      else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(bgdEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(bgdEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(bgdEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(bgdEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(bgdEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(bgdEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = bgdEvent.jet1Btag_;
      else if(thePlot ==31) myVar = bgdEvent.jet2Btag_;
      else if(thePlot ==32) myVar = bgdEvent.jet3Btag_;
      else if(thePlot ==33) myVar = bgdEvent.jetLowBtag_;
      else if(thePlot ==34) myVar = (bgdEvent.jet1_+bgdEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = dRMin;
      else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = mtHWW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = newMet*cos(bgdEvent.metPhi_);
      else if(thePlot ==50) myVar = newMet*sin(bgdEvent.metPhi_);
      else if(thePlot ==51) myVar = newTrackMet*cos(bgdEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = newTrackMet*sin(bgdEvent.trackMetPhi_);

      if(channel == 5 && thePlot == 101){
	                                                                        hDnjets[ 0+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jetLowBtag_ >= 2.1)                    	    		hDnjets[ 1+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.nSoftMuons_ > 0)                       	    		hDnjets[ 2+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.nSoftMuons_ > 0)             hDnjets[ 3+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if((bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets) hDnjets[ 4+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if((bgdEvent.cuts_ & patternTopTag) == patternTopTag)         	        hDnjets[ 5+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jet1_.Pt() > 30.0 && bgdEvent.jet1Btag_ >= 2.1)             hDnjets[ 6+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jet1_.Pt() > 30.0 && bgdEvent.jet1Btag_ >= 2.1 &&
	   (bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets) hDnjets[ 7+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.nSoftMuons_ == 0)             hDnjets[ 8+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.nSoftMuons_ == 0 &&
	   bgdEvent.jet1_.Pt() > 30.0 && bgdEvent.jet1Btag_ >= 2.1)             hDnjets[ 9+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if((bgdEvent.jet1_.Pt() <= 30.0 || bgdEvent.jet1Btag_ < 2.1) && 
	   (bgdEvent.jet2_.Pt() <= 30.0 || bgdEvent.jet2Btag_ < 2.1))           hDnjets[10+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if((bgdEvent.jet1_.Pt() <= 30.0 || bgdEvent.jet1Btag_ < 2.1) && 
	   (bgdEvent.jet2_.Pt() <= 30.0 || bgdEvent.jet2Btag_ < 2.1) &&
	   (bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.nSoftMuons_ > 0))           hDnjets[11+20*(int)isSignalDecay]->Fill( bgdEvent.njets_,theWeight);
	if(bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.nSoftMuons_ > 0){
	  if     (fabs(bgdEvent.jet1_.Eta()) <= fabs(bgdEvent.jet2_.Eta()) && bgdEvent.jet2Btag_ < 2.1){
	    hDetatag[0+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    if     (bgdEvent.jet1_.Pt()<60)
	      hDetatag[1+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    else if(bgdEvent.jet1_.Pt()<90)
	      hDetatag[2+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    else
	      hDetatag[3+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    if(bgdEvent.jet1Btag_ > 2.1){
	      hDetatag[4+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      if     (bgdEvent.jet1_.Pt()<60)
		hDetatag[5+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      else if(bgdEvent.jet1_.Pt()<90)
		hDetatag[6+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      else
		hDetatag[7+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    }
	  }
	  else if(fabs(bgdEvent.jet2_.Eta()) <= fabs(bgdEvent.jet1_.Eta()) && bgdEvent.jet1Btag_ < 2.1){
	    hDetatag[0+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    if     (bgdEvent.jet2_.Pt()<60)
	      hDetatag[1+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    else if(bgdEvent.jet2_.Pt()<90)
	      hDetatag[2+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    else
	      hDetatag[3+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    if(bgdEvent.jet2Btag_ > 2.1){
	      hDetatag[4+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      if     (bgdEvent.jet2_.Pt()<60)
		hDetatag[5+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      else if(bgdEvent.jet2_.Pt()<90)
		hDetatag[6+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      else
		hDetatag[7+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    }
	  }
	}
	int jtype = 0;
	if     (fabs(bgdEvent.jet1McId_) == 5 && fabs(bgdEvent.jet2McId_) == 5) jtype = 1;
	else if(fabs(bgdEvent.jet1McId_) == 5 || fabs(bgdEvent.jet2McId_) == 5) jtype = 2;
	else									jtype = 3;
	if(bgdEvent.njets_ == 2){
	  hDnjets[65]->Fill(jtype,theWeight);
	}
	if(isSignalDecay == true &&
	   TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 &&
	   bgdEvent.jet1_.Eta()*bgdEvent.jet2_.Eta() < 0 &&
	   (bgdEvent.jet1_+bgdEvent.jet2_).M() > 450 &&
	   bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 && Njet3 != 2){
	                                                                  hDnjets[60]->Fill(bgdEvent.njets_,theWeight);
	  if(bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.nSoftMuons_ > 0)     hDnjets[61]->Fill( bgdEvent.njets_,theWeight);
	  if((bgdEvent.jet1_.Pt() <= 30.0 || bgdEvent.jet1Btag_ < 2.1) && 
	     (bgdEvent.jet2_.Pt() <= 30.0 || bgdEvent.jet2Btag_ < 2.1))   hDnjets[62]->Fill( bgdEvent.njets_,theWeight);
	  if((bgdEvent.jet1_.Pt() <= 30.0 || bgdEvent.jet1Btag_ < 2.1) && 
	     (bgdEvent.jet2_.Pt() <= 30.0 || bgdEvent.jet2Btag_ < 2.1) &&
	     (bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.nSoftMuons_ > 0))   hDnjets[63]->Fill( bgdEvent.njets_,theWeight);
	  if((bgdEvent.cuts_ & patternTopTag) == patternTopTag)           hDnjets[64]->Fill(bgdEvent.njets_,theWeight);

	  hDnjets[66]->Fill(jtype,theWeight);
          if((bgdEvent.jet1_.Pt() <= 30.0 || bgdEvent.jet1Btag_ < 2.1) && 
	     (bgdEvent.jet2_.Pt() <= 30.0 || bgdEvent.jet2Btag_ < 2.1)) 
	    hDnjets[67]->Fill(jtype,theWeight);

	  if(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.nSoftMuons_ == 0){
	    if     (fabs(bgdEvent.jet1_.Eta()) <= fabs(bgdEvent.jet2_.Eta()) && bgdEvent.jet2Btag_ < 2.1){
	      hDetatag[8+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      if     (bgdEvent.jet1_.Pt()<60)
		hDetatag[9+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      else if(bgdEvent.jet1_.Pt()<90)
		hDetatag[10+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      else
		hDetatag[11+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      if(bgdEvent.jet1Btag_ > 2.1){
		hDetatag[12+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
		if     (bgdEvent.jet1_.Pt()<60)
		  hDetatag[13+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
		else if(bgdEvent.jet1_.Pt()<90)
		  hDetatag[14+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
		else
		  hDetatag[15+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      }
	    }
	    else if(fabs(bgdEvent.jet2_.Eta()) <= fabs(bgdEvent.jet1_.Eta()) && bgdEvent.jet1Btag_ < 2.1){
	      hDetatag[8+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      if     (bgdEvent.jet2_.Pt()<60)
		hDetatag[9+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      else if(bgdEvent.jet2_.Pt()<90)
		hDetatag[10+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      else
		hDetatag[11+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      if(bgdEvent.jet2Btag_ > 2.1){
		hDetatag[12+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
		if     (bgdEvent.jet2_.Pt()<60)
		  hDetatag[13+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
		else if(bgdEvent.jet2_.Pt()<90)
		  hDetatag[14+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
		else
		  hDetatag[15+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      }
	    }
	  }
	  if     (bgdEvent.jet1Btag_ >= 2.1 || bgdEvent.jet2Btag_ >= 2.1){
	    if(fabs(bgdEvent.jet1_.Eta()) <= fabs(bgdEvent.jet2_.Eta())){
	      hDpttag[16+20*(int)isSignalDecay] ->Fill(bgdEvent.jet1_.Pt(),theWeight);
	      hDetatag[16+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      hDpttag[17+20*(int)isSignalDecay] ->Fill(bgdEvent.jet2_.Pt(),theWeight);
	      hDetatag[17+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    }
	    else {
	      hDpttag[16+20*(int)isSignalDecay] ->Fill(bgdEvent.jet2_.Pt(),theWeight);
	      hDetatag[16+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      hDpttag[17+20*(int)isSignalDecay] ->Fill(bgdEvent.jet1_.Pt(),theWeight);
	      hDetatag[17+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    }
	  }
	  else if(bgdEvent.jet1Btag_ < 2.1 && bgdEvent.jet2Btag_ < 2.1){
	    if(fabs(bgdEvent.jet1_.Eta()) <= fabs(bgdEvent.jet2_.Eta())){
	      hDpttag[18+20*(int)isSignalDecay] ->Fill(bgdEvent.jet1_.Pt(),theWeight);
	      hDetatag[18+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	      hDpttag[19+20*(int)isSignalDecay] ->Fill(bgdEvent.jet2_.Pt(),theWeight);
	      hDetatag[19+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	    }
	    else {
	      hDpttag[18+20*(int)isSignalDecay] ->Fill(bgdEvent.jet2_.Pt(),theWeight);
	      hDetatag[18+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet2_.Eta()),theWeight);
	      hDpttag[19+20*(int)isSignalDecay] ->Fill(bgdEvent.jet1_.Pt(),theWeight);
	      hDetatag[19+20*(int)isSignalDecay]->Fill(fabs(bgdEvent.jet1_.Eta()),theWeight);
	    }
	  }
        }
      }

      if     (fDecay == 14 || fDecay == 29 || fDecay == 30){
        histo0->Fill(myVar,theWeight);
      }
      else if(fDecay == 0 || fDecay == 1 || fDecay == 2 || fDecay == 3 || fDecay == 19){
        histo4->Fill(myVar,theWeight);
      }
      else if(fDecay == 6 || fDecay == 7 || fDecay == 8 || fDecay == 9 || fDecay == 10 || fDecay == 31){
        histo1->Fill(myVar,theWeight);
      }
      else if(fDecay == 27 || fDecay == 28){
        histo3->Fill(myVar,theWeight);
      }
      else if(fDecay == 4 || fDecay == 5 || fDecay == 11 || fDecay == 12 || fDecay == 13){
        histo2->Fill(myVar,theWeight);
      }
      else {
        printf("NOOOOOOOOOOOOOOOOOOOO\n");
      }
      myVar = myVar / xmaxPlot;
      if     (myVar <= 0) myVar = 0.001;
      else if(myVar >= 1) myVar = 0.999;
      for(int n1=0; n1<nBin; n1++){
        if(myVar > 1.0*n1/nBin){
          if(isSignalDecay == true ) S0[n1] = S0[n1] + theWeight;
          if(isSignalDecay == false) B0[n1] = B0[n1] + theWeight;
        }
        if(myVar < 1.0*n1/nBin){
          if(isSignalDecay == true ) S1[n1] = S1[n1] + theWeight;
          if(isSignalDecay == false) B1[n1] = B1[n1] + theWeight;
        }
      }
      if(isSignalDecay == true ) hDSigOpt->Fill(myVar,theWeight);
      if(isSignalDecay == false) hDBckOpt->Fill(myVar,theWeight);
    }
  } // end background loop

  if((channel >= 0 && channel <= 8000)){
    int nSig=sigEvent.tree_->GetEntries();
    for (int i=0; i<nSig; ++i) {

      if (i%10000 == 0 && verboseLevel > 0)
	printf("--- reading Signal event %5d of %5d\n",i,nSig);
      sigEvent.tree_->GetEntry(i);
    int charge = (int)(sigEvent.lq1_ + sigEvent.lq2_);

    int Njet3 = 0;
    if(sigEvent.jet3_.Pt() <= 30)									 Njet3 = 0;
    else if(sigEvent.jet3_.Pt() > 30 && (
      (sigEvent.jet1_.Eta()-sigEvent.jet3_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.jet3_.Eta() < 0) ||
      (sigEvent.jet2_.Eta()-sigEvent.jet3_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.jet3_.Eta() < 0))) Njet3 = 2;
    else												 Njet3 = 1;
    int centrality = 0;
    if(((sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() < 0)) &&
       ((sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() < 0))) centrality = 1; 
    
    double newMet      = sigEvent.met_;
    double newTrackMet = sigEvent.trackMet_;
    double deltaPhiA[3] = {TMath::Abs(sigEvent.lep1_.Phi()-sigEvent.metPhi_),TMath::Abs(sigEvent.lep2_.Phi()-sigEvent.metPhi_),0.0};
    while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
    while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
    deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
    double pmetA = newMet;
    if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

    double deltaPhiB[3] = {TMath::Abs(sigEvent.lep1_.Phi()-sigEvent.trackMetPhi_),TMath::Abs(sigEvent.lep2_.Phi()-sigEvent.trackMetPhi_),0.0};
    while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
    while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
    deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
    double pmetB = newTrackMet;
    if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

    double deltaPhiC[3] = {TMath::Abs(sigEvent.jet1_.Phi()-sigEvent.metPhi_),TMath::Abs(sigEvent.jet2_.Phi()-sigEvent.metPhi_),0.0};
    while(deltaPhiC[0]>TMath::Pi()) deltaPhiC[0] = TMath::Abs(deltaPhiC[0] - 2*TMath::Pi());
    while(deltaPhiC[1]>TMath::Pi()) deltaPhiC[1] = TMath::Abs(deltaPhiC[1] - 2*TMath::Pi());
    double pmetC = newMet;

    double deltaPhiD[3] = {TMath::Abs(sigEvent.jet1_.Phi()-sigEvent.trackMetPhi_),TMath::Abs(sigEvent.jet2_.Phi()-sigEvent.trackMetPhi_),0.0};
    while(deltaPhiD[0]>TMath::Pi()) deltaPhiD[0] = TMath::Abs(deltaPhiD[0] - 2*TMath::Pi());
    while(deltaPhiD[1]>TMath::Pi()) deltaPhiD[1] = TMath::Abs(deltaPhiD[1] - 2*TMath::Pi());
    double pmetD = newTrackMet;

    if     (sigEvent.njets_ == 0){
      pmetC = pmetA;
      pmetD = pmetB;
    }
    else if(sigEvent.njets_ == 1){
      deltaPhiC[2] = deltaPhiC[0];
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = deltaPhiD[0];
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }
    else if(sigEvent.njets_ >= 2){
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiC[1]);
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiD[1]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }

    double usedMet = 0.0;
    bool passMET = false;
    if     (sigEvent.njets_ == 0){
       usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
    }
    else if(sigEvent.njets_ == 1){
       usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
    }
    else if(sigEvent.njets_ >= 2){
       usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
    }
    double massZMin = 999.0;
    double massMin = 999.0;
    double dRMin = 999.0; double mTW3 = 0.0; double type3l = 0.0;
    Float_t mTWMin   = sigEvent.mt1_;
    Float_t mTWMax   = sigEvent.mt1_;
    if(mTWMin > sigEvent.mt2_			     ) mTWMin = sigEvent.mt2_;
    if(mTWMin > sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMin = sigEvent.mt3_;
    if(mTWMax < sigEvent.mt2_			     ) mTWMax = sigEvent.mt2_;
    if(mTWMax < sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMax = sigEvent.mt3_;
    double pxzll = sigEvent.dilep_.Px() + newMet*cos( sigEvent.metPhi_);
    double pyzll = sigEvent.dilep_.Py() + newMet*sin( sigEvent.metPhi_);
    double mtHZZ = TMath::Power(sqrt(sigEvent.dilep_.Pt()*sigEvent.dilep_.Pt()+91.1876*91.1876)+
                              sqrt(newMet       *newMet       +91.1876*91.1876),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;
    double mtHWW = TMath::Power(sqrt(sigEvent.dilep_.Pt()*sigEvent.dilep_.Pt()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+
                                sqrt(newMet       *newMet       +sigEvent.dilep_.M()*sigEvent.dilep_.M()),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHWW >0) mtHWW = sqrt(mtHWW); else mtHWW = 0.0;

    double pznunu = 91.1876*91.1876 - newMet*newMet;
    if(pznunu > 0) pznunu = sqrt(pznunu); else pznunu = 0.0;
    double mHZZ = (sqrt(sigEvent.dilep_.P()*sigEvent.dilep_.P()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+91.1876)*
                  (sqrt(sigEvent.dilep_.P()*sigEvent.dilep_.P()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+91.1876)-
                 -pxzll*pxzll-pyzll*pyzll-(sigEvent.dilep_.Pz()+pznunu)*(sigEvent.dilep_.Pz()+pznunu);
    if(mHZZ > 0) mHZZ = sqrt(mHZZ); else mHZZ = 0.001;

    if(sigEvent.lid3_ != 0){
      usedMet = TMath::Min(newMet,newTrackMet);
      massZMin = trilepton_info(0,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
      massMin  = trilepton_info(1,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
      dRMin = trilepton_info(2,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
      mTW3  = trilepton_info(3,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
      type3l= trilepton_info(4,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
    }    
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 29){ // WW selection
      if(
         sigEvent.dilep_.M()   > 12 &&
         sigEvent.lid3_ == 0 &&
         charge == 0 &&
	 sigEvent.njets_ == nJetsType &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.njets_ != 1 ||
	 (sigEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) &&
	 //sigEvent.dilep_.M()   > 100 &&
	 //(sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) &&
	 //sigEvent.type_ == SmurfTree::mm &&
	 (sigEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WW selection

    if(channel == 6){ // Z selection
      if(
         sigEvent.dilep_.M()   > 12 &&
         sigEvent.lid3_ == 0. &&
         charge == 0 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
        (sigEvent.lep2_.Pt() > 15. || sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::em) &&
         (fabs(sigEvent.dilep_.M()-91.1876) < 15.) && 
         //sigEvent.njets_ == 1 &&
	 //fabs(sigEvent.dilep_.Rapidity()) < 1 &&
	 //(sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) &&
	 //sigEvent.type_ == SmurfTree::mm &&
	 //sigEvent.type_ == SmurfTree::ee &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Z selection

    if(channel == 5){ // ttbar selection
      if(
         sigEvent.dilep_.M()   > 12 &&
         sigEvent.lid3_ == 0 &&
         charge == 0 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_) > 20. &&
        (TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_) > 35. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
	 (sigEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // ttbar selection

    if(channel == 28 || channel == 128){ // ZZllnn selection
      if(
        //(sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         sigEvent.lid3_ == 0. &&
         charge == 0 &&
         fabs(sigEvent.dilep_.M()-91.1876) < 15. && 
	(sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee) &&
         sigEvent.lep1_.Pt() > 30. &&
         sigEvent.lep2_.Pt() > 20. &&
         sigEvent.dilep_.Pt() > 55. &&
	(sigEvent.jet1Btag_ < 2 || sigEvent.jet1_.Pt() < 30) &&
	(sigEvent.jet2Btag_ < 2 || sigEvent.jet2_.Pt() < 30) &&
	(sigEvent.jet3Btag_ < 2 || sigEvent.jet3_.Pt() < 30) &&
         sigEvent.met_ > hzz2l_cuts(mhAna,0) &&
         sigEvent.met_ < hzz2l_cuts(mhAna,1) &&
	 mtHZZ  > hzz2l_cuts(mhAna,2) &&
	 mtHZZ  < hzz2l_cuts(mhAna,3) &&
	 1 == 1
	){
	//
	passCuts = true;
        isSignalDecay = true;
      }
    } // ZZllnn selection

    if(channel == 27){ // WZ selection
      charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
      if(
         sigEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         sigEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
        (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         sigEvent.jet1_.Pt() <= 40. &&
	 massZMin < 25 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WZ selection

    if(channel == 9){ // Ztautau selection
      if(
         sigEvent.dilep_.M()   > 12 &&
         sigEvent.dilep_.M()   < 80 &&
         sigEvent.lid3_ == 0. &&
         charge == 0 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
	 TMath::Min(sigEvent.mt1_,sigEvent.mt2_) < 40. &&
	 (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Ztautau selection

    if(channel == 900 || channel == 901){ // qqH selection
      if(
        (sigEvent.dilep_.M()   < 100 || channel == 901) &&
         sigEvent.dilep_.M()   > 12  &&
         sigEvent.lid3_ == 0 &&
         charge == 0 &&
	 sigEvent.njets_ >= 2 && sigEvent.njets_ <= 3 && Njet3 != 2 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
        (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
        (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         (sigEvent.jet1_+sigEvent.jet2_).M() > 450. &&
          //TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta()) > 3.5 &&
         sigEvent.jet1_.Eta()*sigEvent.jet2_.Eta() < 0 &&
	 centrality == 1 &&
	 (sigEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // qqH selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection
      double metCut = cutMetLow[binc];
      if(sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) metCut = cutMetLowEM[binc];
      if(
         sigEvent.mt_	     > cutMTLow[binc] &&
         sigEvent.mt_	     < cutMTHigh[binc] &&
         sigEvent.dilep_.M() > cutMassLow[binc] &&
         sigEvent.dilep_.M() < cutMassHigh[binc] &&
         sigEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
         sigEvent.lep2_.Pt() > cutPtMinLow[binc] &&
         sigEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
         sigEvent.lid3_ == 0 &&
         charge == 0 &&
	 sigEvent.njets_ == nJetsType &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.njets_ != 1 ||
	 (sigEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) &&
	 //(sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) &&
	 //sigEvent.type_ == SmurfTree::mm &&
	 (sigEvent.type_ == lDecay || lDecay == 4) &&
	  1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // H->WW selection

    if(channel == 1000){ // HW->2l selection
      if(
         sigEvent.dilep_.M()   > 12 &&
         sigEvent.lid3_ == 0 &&
         charge != 0 &&
	 sigEvent.njets_ >= 2 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 10. || sigEvent.type_ != SmurfTree::ee) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // HW->2l selection

    if(channel == 1100){ // HW->3l selection
      charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
      if(
         sigEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         sigEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         sigEvent.jet1_.Pt() <= 40. &&
	 massZMin > 25 &&
	 massMin < 100 &&
	 mTWMax > 100 &&
	 dRMin < 2.0 &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // HW->3l selection
    
    //if(passCuts == true && sigEvent.processId_==10001){
    if(passCuts == true){
      double add = 1.;
      add = scaleFactor(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), sigEvent.nvtx_, 2);
      if((sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
        add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
      if((sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
        add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
     							       fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
	   						       TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
      add = add*trigEff;

      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      //add = 1.0;
      if (sigEvent.processId_ == 10010) {
        add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(sigEvent.higgsPt_));
	//add = add * enhancementFactor(channel);
      }

      double myWeight = sigEvent.scale1fb_*lumi*add;
      double myVar = newMet;
      if     (thePlot == 1) myVar = sigEvent.lep1_.Pt();
      else if(thePlot == 2) myVar = sigEvent.lep2_.Pt();
      else if(thePlot == 3) myVar = sigEvent.lep3_.Pt();
      else if(thePlot == 4) myVar = sigEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = sigEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = sigEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = sigEvent.dilep_.M();
      else if(thePlot == 8) myVar = sigEvent.mt_;
      else if(thePlot == 9) myVar = sigEvent.mt1_;
      else if(thePlot ==10) myVar = sigEvent.mt2_;
      else if(thePlot ==11) myVar = sigEvent.mt3_;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = newMet/sigEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = sigEvent.njets_;
      else if(thePlot ==18) myVar = sigEvent.nvtx_;
      else if(thePlot ==19) myVar = sigEvent.nSoftMuons_;
      else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = sigEvent.dPhiDiLepJet1_*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(sigEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(sigEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(sigEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(sigEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(sigEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(sigEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = sigEvent.jet1Btag_;
      else if(thePlot ==31) myVar = sigEvent.jet2Btag_;
      else if(thePlot ==32) myVar = sigEvent.jet3Btag_;
      else if(thePlot ==33) myVar = sigEvent.jetLowBtag_;
      else if(thePlot ==34) myVar = (sigEvent.jet1_+sigEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = dRMin;
      else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,sigEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = mtHWW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = newMet*cos(sigEvent.metPhi_);
      else if(thePlot ==50) myVar = newMet*sin(sigEvent.metPhi_);
      else if(thePlot ==51) myVar = newTrackMet*cos(sigEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = newTrackMet*sin(sigEvent.trackMetPhi_);
      histos->Fill(myVar,myWeight);
      if     (sigEvent.processId_==10010) hDHSig[0]->Fill(1.0,myWeight);
      else if(sigEvent.processId_==10001) hDHSig[1]->Fill(1.0,myWeight);
      else if(sigEvent.processId_==26  ||
              sigEvent.processId_==24  ||
              sigEvent.processId_==121 ||
              sigEvent.processId_==122)   hDHSig[2]->Fill(1.0,myWeight);

      myVar = myVar / xmaxPlot;

      if     (myVar <= 0) myVar = 0.001;
      else if(myVar >= 1) myVar = 0.999;
      for(int n1=0; n1<nBin; n1++){
        if(myVar > 1.0*n1/nBin){
          if(isSignalDecay == true ) S0[n1] = S0[n1] + myWeight;
        }
        if(myVar < 1.0*n1/nBin){
          if(isSignalDecay == true ) S1[n1] = S1[n1] + myWeight;
        }
      }
      if(isSignalDecay == true ) hDSigOpt->Fill(myVar,myWeight);
    }
    } // Loop over signal
  }

  int nData=dataEvent.tree_->GetEntries();
  int nSelectedData = 0;
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
       //usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    else if(dataEvent.njets_ == 1){
       //usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    else if(dataEvent.njets_ >= 2){
       //usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
       passMET = usedMet > 20. &&
                (usedMet > 35. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
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
    bool passCuts = false;
    if(channel == 29){ // WW selection
      if(
         dataEvent.dilep_.M()   > 12 &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
	 dataEvent.njets_ == nJetsType &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
	(dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
	(dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
         passMET == true &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
         //usedMet > 35. &&
         //fabs(dataEvent.dilep_.M()-91.1876) < 15. && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         //(dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets &&
	 //(dataEvent.jetLowBtag_ < 2.1 && dataEvent.nSoftMuons_ == 0) &&
	 //dataEvent.jet1Btag_ >= 2.1 &&
	 //dataEvent.jet2Btag_ >= 2.1 &&
	 (dataEvent.njets_ != 1 ||
	 (dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me)) &&
	 //dataEvent.dilep_.M()   > 100 &&
	 //(dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
	 //dataEvent.type_ == SmurfTree::mm &&
	 (dataEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
	//cout << "DATA "  << dataEvent.type_ << " " << dataEvent.run_ << " " << dataEvent.event_ << " " << dataEvent.lumi_ << endl;
      }
    } // WW selection
    if(channel == 6){ // Z selection
      if(
         dataEvent.dilep_.M()   > 12 &&
         dataEvent.lid3_ == 0. &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
        (fabs(dataEvent.dilep_.M()-91.1876) < 15.) && 
         //dataEvent.njets_ == 0 &&
	 //fabs(dataEvent.dilep_.Rapidity()) < 1 &&
	 !(dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
	 //dataEvent.type_ == SmurfTree::mm &&
	 //dataEvent.type_ == SmurfTree::ee &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Z selection

    if(channel == 5){ // ttbar selection
      if(
         dataEvent.dilep_.M()   > 12 &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_) > 20. &&
        (TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_) > 35. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
	 (dataEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // ttbar selection

    if(channel == 28 || channel == 128){ // ZZllnn selection
      if(
        //(dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.lid3_ == 0. &&
         charge == 0 &&
         fabs(dataEvent.dilep_.M()-91.1876) < 15. && 
	(dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee) &&
         dataEvent.lep1_.Pt() > 30. &&
         dataEvent.lep2_.Pt() > 20. &&
         dataEvent.dilep_.Pt() > 55. &&
	(dataEvent.jet1Btag_ < 2 || dataEvent.jet1_.Pt() < 30) &&
	(dataEvent.jet2Btag_ < 2 || dataEvent.jet2_.Pt() < 30) &&
	(dataEvent.jet3Btag_ < 2 || dataEvent.jet3_.Pt() < 30) &&
         dataEvent.met_ > hzz2l_cuts(mhAna,0) &&
         dataEvent.met_ < hzz2l_cuts(mhAna,1) &&
	 mtHZZ  > hzz2l_cuts(mhAna,2) &&
	 mtHZZ  < hzz2l_cuts(mhAna,3) &&
	 1 == 1
	){
	//
	passCuts = true;
      }
    } // ZZllnn selection

    if(channel == 27){ // WZ selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
        (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.jet1_.Pt() <= 40. &&
	 massZMin < 25 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WZ selection

    if(channel == 9){ // Ztautau selection
      if(
         dataEvent.dilep_.M()   > 12 &&
         dataEvent.dilep_.M()   < 80 &&
         dataEvent.lid3_ == 0. &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
	 TMath::Min(dataEvent.mt1_,dataEvent.mt2_) < 40. &&
	 (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Ztautau selection

    if(channel == 900 || channel == 901){ // qqH selection
      if(
        (dataEvent.dilep_.M()   < 100 || channel == 901) &&
         dataEvent.dilep_.M()   > 12  &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
	 dataEvent.njets_ >= 2 && dataEvent.njets_ <= 3 && Njet3 != 2 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         passMET == true &&
        (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
        (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         (dataEvent.jet1_+dataEvent.jet2_).M() > 450. &&
         // TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 &&
         dataEvent.jet1_.Eta()*dataEvent.jet2_.Eta() < 0 &&
	 centrality == 1 &&
	 (dataEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
	//cout << dataEvent.lep1_.Pt() << " " << dataEvent.lep2_.Pt() << " " << dataEvent.mt_ << " " << dataEvent.dPhi_*180.0/TMath::Pi()<< " " << dataEvent.type_<<endl;
      }
    } // qqH selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection
      double metCut = cutMetLow[binc];
      if(dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) metCut = cutMetLowEM[binc];
      if(
         dataEvent.mt_	     > cutMTLow[binc] &&
         dataEvent.mt_	     < cutMTHigh[binc] &&
         dataEvent.dilep_.M() > cutMassLow[binc] &&
         dataEvent.dilep_.M() < cutMassHigh[binc] &&
         dataEvent.lep1_.Pt() > cutPtMaxLow[binc] &&
         dataEvent.lep2_.Pt() > cutPtMinLow[binc] &&
         dataEvent.dPhi_*180.0/TMath::Pi() < cutDeltaphilHigh[binc] &&
         dataEvent.lid3_ == 0 &&
         charge == 0 &&
	 dataEvent.njets_ == nJetsType &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         passMET == true &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (dataEvent.njets_ != 1 ||
	 (dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me)) &&
	 //(dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) &&
	 //dataEvent.type_ == SmurfTree::mm &&
	 (dataEvent.type_ == lDecay || lDecay == 4) &&
	  1 == 1
	){
	passCuts = true;
      }
    } // H->WW selection

    if(channel == 1000){ // HW->2l selection
      if(
         dataEvent.dilep_.M()   > 12 &&
         dataEvent.lid3_ == 0 &&
         charge != 0 &&
	 dataEvent.njets_ >= 2 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
        (dataEvent.lep2_.Pt() > 15. || dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::em) &&
         passMET == true &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (dataEvent.type_ == lDecay || lDecay == 4) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->2l selection

    if(channel == 1100){ // HW->3l selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 20. &&
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.jet1_.Pt() <= 40. &&
	 massZMin > 25 &&
	 massMin < 100 &&
	 mTWMax > 100 &&
	 dRMin < 2.0 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection
    
    if(passCuts == true){
      double myVar = dataEvent.met_;
      if     (thePlot == 1) myVar = dataEvent.lep1_.Pt();
      else if(thePlot == 2) myVar = dataEvent.lep2_.Pt();
      else if(thePlot == 3) myVar = dataEvent.lep3_.Pt();
      else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = dataEvent.dilep_.M();
      else if(thePlot == 8) myVar = dataEvent.mt_;
      else if(thePlot == 9) myVar = dataEvent.mt1_;
      else if(thePlot ==10) myVar = dataEvent.mt2_;
      else if(thePlot ==11) myVar = dataEvent.mt3_;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = dataEvent.met_/dataEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = dataEvent.njets_;
      else if(thePlot ==18) myVar = dataEvent.nvtx_;
      else if(thePlot ==19) myVar = dataEvent.nSoftMuons_;
      else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(dataEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(dataEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(dataEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(dataEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(dataEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(dataEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = dataEvent.jet1Btag_;
      else if(thePlot ==31) myVar = dataEvent.jet2Btag_;
      else if(thePlot ==32) myVar = dataEvent.jet3Btag_;
      else if(thePlot ==33) myVar = dataEvent.jetLowBtag_;
      else if(thePlot ==34) myVar = (dataEvent.jet1_+dataEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = dRMin;
      else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = mtHWW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = dataEvent.met_*cos(dataEvent.metPhi_);
      else if(thePlot ==50) myVar = dataEvent.met_*sin(dataEvent.metPhi_);
      else if(thePlot ==51) myVar = dataEvent.trackMet_*cos(dataEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = dataEvent.trackMet_*sin(dataEvent.trackMetPhi_);

      if(channel == 5 && thePlot == 101){
	                                                                          hDnjets[ 0+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jetLowBtag_ >= 2.1)                    	    		  hDnjets[ 1+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.nSoftMuons_ > 0)                       	    		  hDnjets[ 2+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jetLowBtag_ >= 2.1 || dataEvent.nSoftMuons_ > 0)             hDnjets[ 3+40]->Fill( dataEvent.njets_,1.0);
	if((dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets)  hDnjets[ 4+40]->Fill( dataEvent.njets_,1.0);
	if((dataEvent.cuts_ & patternTopTag) == patternTopTag)         	          hDnjets[ 5+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jet1_.Pt() > 30.0 && dataEvent.jet1Btag_ >= 2.1)             hDnjets[ 6+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jet1_.Pt() > 30.0 && dataEvent.jet1Btag_ >= 2.1 &&
	   (dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets)  hDnjets[ 7+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jetLowBtag_ < 2.1 && dataEvent.nSoftMuons_ == 0)             hDnjets[ 8+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jetLowBtag_ < 2.1 && dataEvent.nSoftMuons_ == 0 &&
	   dataEvent.jet1_.Pt() > 30.0 && dataEvent.jet1Btag_ >= 2.1)             hDnjets[ 9+40]->Fill( dataEvent.njets_,1.0);
	if((dataEvent.jet1_.Pt() <= 30.0 || dataEvent.jet1Btag_ < 2.1) && 
	   (dataEvent.jet2_.Pt() <= 30.0 || dataEvent.jet2Btag_ < 2.1))           hDnjets[10+40]->Fill( dataEvent.njets_,1.0);
	if((dataEvent.jet1_.Pt() <= 30.0 || dataEvent.jet1Btag_ < 2.1) && 
	   (dataEvent.jet2_.Pt() <= 30.0 || dataEvent.jet2Btag_ < 2.1) &&
	   (dataEvent.jetLowBtag_ >= 2.1 || dataEvent.nSoftMuons_ > 0))           hDnjets[11+40]->Fill( dataEvent.njets_,1.0);
	if(dataEvent.jetLowBtag_ >= 2.1 || dataEvent.nSoftMuons_ > 0){
	  if     (fabs(dataEvent.jet1_.Eta()) <= fabs(dataEvent.jet2_.Eta()) && dataEvent.jet2Btag_ < 2.1){
	    hDetatag[0+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	    if     (dataEvent.jet1_.Pt()<60)
	      hDetatag[1+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	    else if(dataEvent.jet1_.Pt()<90)
	      hDetatag[2+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	    else
	      hDetatag[3+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	    if(dataEvent.jet1Btag_ > 2.1){
	      hDetatag[4+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	      if     (dataEvent.jet1_.Pt()<60)
		hDetatag[5+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	      else if(dataEvent.jet1_.Pt()<90)
		hDetatag[6+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	      else
		hDetatag[7+40]->Fill(fabs(dataEvent.jet1_.Eta()),1.0);
	    }
	  }
	  else if(fabs(dataEvent.jet2_.Eta()) <= fabs(dataEvent.jet1_.Eta()) && dataEvent.jet1Btag_ < 2.1){
	    hDetatag[0+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	    if     (dataEvent.jet2_.Pt()<60)
	      hDetatag[1+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	    else if(dataEvent.jet2_.Pt()<90)
	      hDetatag[2+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	    else
	      hDetatag[3+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	    if(dataEvent.jet2Btag_ > 2.1){
	      hDetatag[4+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	      if     (dataEvent.jet2_.Pt()<60)
		hDetatag[5+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	      else if(dataEvent.jet2_.Pt()<90)
		hDetatag[6+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	      else
		hDetatag[7+40]->Fill(fabs(dataEvent.jet2_.Eta()),1.0);
	    }
	  }
	}
      }

      nSelectedData = nSelectedData + 1;
      histo5->Fill(myVar,1.0);
    }
  } // End loop data
  printf("data: %d\n",nSelectedData);

  if(channel == 5 && thePlot == 101){
    char output[200];
    sprintf(output,"histo_njets.root");     
    TFile* outFileNjets = new TFile(output,"recreate");
    outFileNjets->cd();
    for(int nj=0; nj<68; nj++) hDnjets[nj]->Write();
    for(int nj=0; nj<60; nj++) hDpttag[nj]->Write();
    for(int nj=0; nj<60; nj++) hDetatag[nj]->Write();
    outFileNjets->Close();
  }

  if(makeGoodPlots){
    char output[200];
    sprintf(output,"histo_nice.root");     
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
      histos->Write();
      histo0->Write();
      histo1->Write();
      histo2->Write();
      histo3->Write();
      histo4->Write();
      histo5->Write();
    outFilePlotsNote->Close();
  }

  if(channel == 1000 || channel == 1100 ||
     (channel >= 100 && channel <= 800)){
    for(int n1=0; n1<nBin; n1++){
      double signif0 = 0.0;
      if(S0[n1] != 0 && B0[n1] != 0) signif0 = scpFast(S0[n1],B0[n1],0.35*B0[n1]);
      double signif1 = 0.0;
      if(S1[n1] != 0 && B1[n1] != 0) signif1 = scpFast(S1[n1],B1[n1],0.35*B1[n1]);
      hDSignif[0]->SetBinContent(n1, signif0);
      hDSignif[1]->SetBinContent(n1, signif1);
    }
    printf("S: %f, B: %f, Signif: %f, S/B: %f\n",S0[0],B0[0],
    	   scpFast(S0[0],B0[0],0.35*B0[0]),S0[0]/B0[0]);
  }
  else if(B0[0] != 0){
    for(int n1=0; n1<nBin; n1++){
      double signif0 = 0.0;
      if(B0[n1] != 0) signif0 = S0[n1]/sqrt(S0[n1]+B0[n1]+0.35*0.35*B0[n1]*B0[n1]);
      double signif1 = 0.0;
      if(B1[n1] != 0) signif1 = S1[n1]/sqrt(S1[n1]+B1[n1]+0.35*0.35*B1[n1]*B1[n1]);
      hDSignif[0]->SetBinContent(n1, signif0);
      hDSignif[1]->SetBinContent(n1, signif1);
    }
    printf("S: %f, B: %f, Signif: %f, S/B: %f\n",S0[0],B0[0],
    	   S0[0]/sqrt(S0[0]+B0[0]+0.35*0.35*B0[0]*B0[0]),S0[0]/B0[0]);
  }
  
  bgdDecay[29] = bgdDecay[29] + bgdDecay[14];
  weiDecay[29] = weiDecay[29] + weiDecay[14];
  bgdDecay[14] = 0.0;
  weiDecay[14] = 0.0;
  // W/Z/tX/Vgamma
  double bgdCombined[4] = {bgdDecay[0]+bgdDecay[1]+bgdDecay[2]+bgdDecay[3],
                           bgdDecay[6]+bgdDecay[7]+bgdDecay[8]+bgdDecay[9]+bgdDecay[10]+bgdDecay[31],
			   bgdDecay[4]+bgdDecay[5]+bgdDecay[11]+bgdDecay[12]+bgdDecay[13],
			   bgdDecay[17]+bgdDecay[18]+bgdDecay[19]};
  double bgdCombinedE[4] = {sqrt(weiDecay[0]+weiDecay[1]+weiDecay[2]+weiDecay[3]),
                            sqrt(weiDecay[6]+weiDecay[7]+weiDecay[8]+weiDecay[9] +weiDecay[10]+weiDecay[31]),
                            sqrt(bgdDecay[4]*weiDecay[4]+weiDecay[5]+weiDecay[11]+weiDecay[12]+weiDecay[13]),
                            sqrt(weiDecay[17]+weiDecay[18]+weiDecay[19])};
  for(int i=0; i<45; i++){
    if(bgdDecay[i] != 0) printf("bdg(%2d) = %f +/- %f\n",i,bgdDecay[i],sqrt(weiDecay[i]));
  }
  for(int i=0; i<16; i++){
    if(bgdDecayFake[i] != 0) printf("bgdDecayFake(%2d) = %4.1f +/- %4.1f\n",i,bgdDecayFake[i],sqrt(bgdDecayFakeE[i]));
  }
  if(bgdCombined[0] != 0) printf("bdg( W) = %f +/- %f\n",bgdCombined[0],bgdCombinedE[0]);
  if(bgdCombined[1] != 0) printf("bdg( Z) = %f +/- %f\n",bgdCombined[1],bgdCombinedE[1]);
  if(bgdCombined[2] != 0) printf("bdg(tX) = %f +/- %f\n",bgdCombined[2],bgdCombinedE[2]);
  if(bgdCombined[3] != 0) printf("bdg(Vg) = %f +/- %f\n",bgdCombined[3],bgdCombinedE[3]);
  // ggWW, WW, WZ, ZZ, W, Z, tX, Vg
  double Syst[8] = {0.00, 0.00, 0.18, 0.18, 0.50, 1.00, 1.00, 0.35};
  //double Syst[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if(channel >=100 && channel <= 200) {Syst[0] = 0.50; Syst[1] = 0.55;}
  if(channel > 200 && channel <= 800) {Syst[0] = 0.50; Syst[1] = 0.50;}
  double total_stat = sqrt(weiDecay[30] +
                           weiDecay[29] +
                           weiDecay[28] +
                           weiDecay[27] +
                           bgdCombinedE[0]*bgdCombinedE[0]+
		           bgdCombinedE[1]*bgdCombinedE[1]+
		           bgdCombinedE[2]*bgdCombinedE[2]+
		           bgdCombinedE[3]*bgdCombinedE[3]);
  double totalS = sqrt(weiDecay[30]+bgdDecay[30]*Syst[0]*bgdDecay[30]*Syst[0] +
                       weiDecay[29]+bgdDecay[29]*Syst[1]*bgdDecay[29]*Syst[1] +
                       weiDecay[28]+bgdDecay[28]*Syst[2]*bgdDecay[28]*Syst[2] +
                       weiDecay[27]+bgdDecay[27]*Syst[3]*bgdDecay[27]*Syst[3] +
                       bgdCombinedE[0]*bgdCombinedE[0]+bgdCombined[0]*Syst[4]*bgdCombined[0]*Syst[4]+
		       bgdCombinedE[1]*bgdCombinedE[1]+bgdCombined[1]*Syst[5]*bgdCombined[1]*Syst[5]+
		       bgdCombinedE[2]*bgdCombinedE[2]+bgdCombined[2]*Syst[6]*bgdCombined[2]*Syst[6]+
		       bgdCombinedE[3]*bgdCombinedE[3]+bgdCombined[3]*Syst[7]*bgdCombined[3]*Syst[7]);
  printf("totalstat: %f\n",total_stat);
  printf("totalS   : %f\n",totalS);
  printf("totalE   : %f\n",sqrt(hDSigOpt->GetSumOfWeights()+hDBckOpt->GetSumOfWeights()));
  hDSignif[0]->SetDirectory(0);
  hDSignif[1]->SetDirectory(0);
  hDSigOpt->SetDirectory(0);
  hDBckOpt->SetDirectory(0);
  //hDSigOpt->Divide(hDBckOpt);
  //hDSigOpt->Draw();
  if(fillInfoNote == true){
    double jeteff_E = 1.02;
    double topXS_E  = 1.37;
    double  wwXS_E  = 1.30;
    double XS_QCDscale_ggH[3];
    XS_QCDscale_ggH[0] = 1.110;
    XS_QCDscale_ggH[1] = 0.910;
    XS_QCDscale_ggH[2] = 1.000;
    if     (nJetsType == 1) {
      jeteff_E = 1.05;
      topXS_E  = 1.22;
      wwXS_E   = 1.30;
      XS_QCDscale_ggH[0] = 1.000;
      XS_QCDscale_ggH[1] = 1.300;
      XS_QCDscale_ggH[2] = 0.890;
    }
    else if(nJetsType == 2) {
      jeteff_E = 1.05;
      topXS_E  = 2.00;
      wwXS_E   = 1.30;
      XS_QCDscale_ggH[0] = 1.000;
      XS_QCDscale_ggH[1] = 1.000;
      XS_QCDscale_ggH[2] = 1.620;
    }

    for(int i=0; i<45; i++) if(bgdDecay[i] < 0) bgdDecay[i] = 0.0;
    bgdDecay[27] = bgdDecay[27] + bgdDecay[28];
    weiDecay[27] = sqrt(weiDecay[27] + weiDecay[28]);
    if(bgdDecay[27] > 0) weiDecay[27] = weiDecay[27]/bgdDecay[27];
    if(bgdDecay[29] > 0) weiDecay[29] = weiDecay[29]/bgdDecay[29];
    if(bgdDecay[30] > 0) weiDecay[30] = weiDecay[30]/bgdDecay[30];
    if(bgdCombined[0] > 0) bgdCombinedE[0] =  bgdCombinedE[0] / bgdCombined[0];
    if(bgdCombined[1] > 0) bgdCombinedE[1] =  bgdCombinedE[1] / bgdCombined[1];
    if(bgdCombined[2] > 0) bgdCombinedE[2] =  bgdCombinedE[2] / bgdCombined[2];
    if(bgdCombined[3] > 0) bgdCombinedE[3] =  bgdCombinedE[3] / bgdCombined[3];
    
    char outputLimitsCut[200];
    sprintf(outputLimitsCut,"histo_limits_%dj_chan%d_mh%d_cut.txt",nJetsType,lDecay,mH);     
    ofstream newcardCut;
    newcardCut.open(outputLimitsCut);
    newcardCut << Form("imax 1 number of channels\n");
    newcardCut << Form("jmax 7 number of background\n");
    newcardCut << Form("kmax * number of nuisance parameters\n");
    newcardCut << Form("Observation %d\n",nSelectedData);
    newcardCut << Form("bin 1 1 1 1 1 1 1 1 1 1\n");
    newcardCut << Form("process VH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma\n");
    newcardCut << Form("process -2 -1 0 1 2 3 4 5 6 7\n");
    newcardCut << Form("rate %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",TMath::Max(hDHSig[2]->GetSumOfWeights(),0.001),TMath::Max(hDHSig[1]->GetSumOfWeights(),0.001),TMath::Max(hDHSig[0]->GetSumOfWeights(),0.001),bgdDecay[29],bgdDecay[30],bgdDecay[27],bgdCombined[2],bgdCombined[1],bgdCombined[0],bgdCombined[3]);
    newcardCut << Form("lumi                  lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.000 1.000 1.000 1.040\n");			   
    newcardCut << Form("CMS_eff_m             lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n");			   
    newcardCut << Form("CMS_trigger_m         lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n");			   
    newcardCut << Form("CMS_eff_e             lnN 1.025 1.025 1.025 1.025 1.025 1.025 1.000 1.000 1.000 1.025\n");			   
    newcardCut << Form("CMS_trigger_e         lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n");			   
    newcardCut << Form("CMS_p_scale_m         lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n");			   
    newcardCut << Form("CMS_p_scale_e         lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n");			   
    newcardCut << Form("CMS_met               lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.000 1.000 1.000 1.020\n");			   
    newcardCut << Form("CMS_jes               lnN %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f 1.000 1.000 1.000 %5.3f\n",jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E,jeteff_E);			   
    newcardCut << Form("CMS_fake              lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.400 1.000\n");
    newcardCut << Form("UE_PS                 lnN 1.030 1.030 1.030 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("pdf_gg                lnN 1.000 1.000 1.100 1.000 1.040 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("pdf_qqbar             lnN 1.050 1.050 1.000 1.040 1.000 1.040 1.000 1.000 1.000 1.040\n");
    newcardCut << Form("QCDscale_ggH          lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[0]);  
    newcardCut << Form("QCDscale_ggH1in       lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[1]);  
    newcardCut << Form("QCDscale_ggH2in       lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",XS_QCDscale_ggH[2]);  
    newcardCut << Form("QCDscale_qqH          lnN 1.000 1.010 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_VH           lnN 1.020 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");			   
    newcardCut << Form("QCDscale_VV           lnN 1.000 1.000 1.000 1.000 1.000 1.040 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_ggVV         lnN 1.000 1.000 1.000 1.000 1.500 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_ggH_EXTRAP   lnN 1.000 1.000 1.100 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_qqH_EXTRAP   lnN 1.000 1.100 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("QCDscale_VH_EXTRAP    lnN 1.100 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n");
    newcardCut << Form("CMS_ttbar             lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",topXS_E);		  
    newcardCut << Form("CMS_Z                 lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.500 1.000 1.000\n");			  
    newcardCut << Form("CMS_WW                lnN 1.000 1.000 1.000 %5.3f %5.3f 1.000 1.000 1.000 1.000 1.000\n",wwXS_E,wwXS_E);			  
    newcardCut << Form("CMS_stat_VH_%1dj      lnN %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,0.03+1.0);	     
    newcardCut << Form("CMS_stat_qqH_%1dj     lnN 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,0.03+1.0);	     
    newcardCut << Form("CMS_stat_ggH_%1dj     lnN 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,0.03+1.0);	     
    newcardCut << Form("CMS_stat_WW_%1dj      lnN 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000 1.000\n",nJetsType,weiDecay[29]+1.0);	 
    newcardCut << Form("CMS_stat_ggWW_%1dj    lnN 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000 1.000\n",nJetsType,weiDecay[30]+1.0);	 
    newcardCut << Form("CMS_stat_VV_%1dj      lnN 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000 1.000\n",nJetsType,weiDecay[27]+1.0);	 
    newcardCut << Form("CMS_stat_ttbar_%1dj   lnN 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000 1.000\n",nJetsType,bgdCombinedE[2]+1.0); 
    newcardCut << Form("CMS_stat_Z_%1dj       lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000 1.000\n",nJetsType,bgdCombinedE[1]+1.0); 
    newcardCut << Form("CMS_stat_Wjets_%1dj   lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f 1.000\n",nJetsType,bgdCombinedE[0]+1.0); 
    newcardCut << Form("CMS_stat_Wgamma_%1dj  lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %5.3f\n",nJetsType,bgdCombinedE[3]+1.0); 
  }
  return;
  bool doNicePlot = false;
  if(doNicePlot == false){
    TCanvas* c1 = new TCanvas("c1","c1",100,100,700,800);
    c1->Divide(2,2);
    c1->cd(1);
    hDSignif[0]->Draw();
    c1->cd(2);
    hDSignif[1]->Draw();
    c1->cd(3);
    hDSigOpt->Draw();
    c1->cd(4);
    hDBckOpt->Draw();
    char output[200];
    sprintf(output,"histo_tmva_%d.root",mH);
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
      hDSigOpt->Write();
      hDBckOpt->Write();
    outFilePlotsNote->Close();
  }
  if(doNicePlot == true){
    TCanvas* newc = new TCanvas("newc","newc",100,100,700,800);
    newc->Divide(1);
    newc->cd(1);
    hDSigOpt->SetTitle("");
    hDSigOpt->SetLineColor(2);
    hDSigOpt->SetLineStyle(0);
    hDSigOpt->SetLineWidth(4);
    hDSigOpt->GetXaxis()->SetTitle(xTitle);
    hDSigOpt->GetXaxis()->SetTitleOffset(1.0);
    hDSigOpt->GetYaxis()->SetTitle(yTitle);
    hDSigOpt->GetYaxis()->SetTitleOffset(1.3);

    hDBckOpt->SetTitle("");
    hDBckOpt->SetFillColor(4);
    //hDBckOpt->SetFillStyle(1001);
    hDBckOpt->SetLineStyle(0);
    hDBckOpt->SetLineWidth(4);
    hDBckOpt->GetXaxis()->SetTitle(xTitle);
    hDBckOpt->GetXaxis()->SetTitleOffset(1.0);
    hDBckOpt->GetYaxis()->SetTitle(yTitle);
    hDBckOpt->GetYaxis()->SetTitleOffset(1.3);

    //hDBckOpt->Add(hDSigOpt);
    hDSigOpt->SetMinimum(0.001);
    hDBckOpt->SetMinimum(0.001);

    hDSigOpt->Scale(1./hDSigOpt->GetSumOfWeights());
    hDBckOpt->Scale(1./hDBckOpt->GetSumOfWeights());

    hDSigOpt->Draw("hist");
    hDBckOpt->Draw("same,hist,e");
    //hDBckOpt->Draw("hist,e");
    //hDSigOpt->Draw("same,hist");

    TLegend* leg = new TLegend(0.5,0.7,0.8,0.8);    	        				
    leg->SetFillColor(10);			    	        				
    leg->SetTextSize(0.03);			    	        			 
    leg->AddEntry(hDBckOpt,"Background","l");    	   
    //leg->AddEntry(hDSigOpt,"WZ #rightarrow l#null","l");	    	        	       
    //leg->AddEntry(hDSigOpt,"WW #rightarrow l#nul#nu","l");	    	        	       
    leg->AddEntry(hDSigOpt,"H(160) #rightarrow WW","l");	    	        	       
    //leg->AddEntry(hDSigOpt,"Signal","l");	    	        	       
    leg->Draw("same");
    TPaveText* labelcms  = new TPaveText(0.7,0.95,0.97,0.97,"NDCBR");
    labelcms->SetTextAlign(12);
    labelcms->SetTextSize(0.05);
    labelcms->SetFillColor(0);
    labelcms->AddText("L = 1 fb^{-1}");
    //labelcms->AddText("L = 100 pb^{-1}");
    labelcms->SetBorderSize(0);
    //labelcms->Draw("same");
  }

  return;

}
