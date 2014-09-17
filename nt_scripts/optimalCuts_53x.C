#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
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
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/HiggsWWStarMassBoundNoROOT.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/HWWCuts.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/HiggsQCDScaleSystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/PSUESystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/PDFgHHSystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/InterfgHHSystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/HWWKinematics.cc"

double caloMet(double met, double metPhi, LorentzVector l1,LorentzVector l2);

const int verboseLevel =   1;
const bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC

// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCuts_53x
(
 int     mH  	 = 29,
 int thePlot = 7,
 TString bgdInputFile    = "ntuples_53x/backgroundA_skim6.root",
 TString signalInputFile = "ntuples_53x/hww125.root",
 TString dataInputFile   = "ntuples_53x/data_skim6.root",
 bool fillInfoNote = true,
 double mhAna = 16,
 int period = 3
 //,double VarAux0 = 0.0, double VarAux1 = 0.0, int VarAux2 = 0
 )
{
  int category = 0;
  if(period > 9) {period = period - 10; category = 1;}
  unsigned int nJetsType = (int)mhAna/10;
  int lDecay    = (int)mhAna%10;
  printf("nJetsType: %d, lDecay = %d\n",nJetsType,lDecay);
  bool makeGoodPlots = true;
  double lumi = 1.0;
  bool isFermioPhobic = false;
  bool isSM4          = false;
 
  bool fCheckProblem = true;

  bool WWXSSel = false;
  double ptLepMin = 10.0;
  if(WWXSSel == true) ptLepMin = 20.;

  SmurfTree bgdEvent;
  if(thePlot != 63 && thePlot != 64) bgdEvent.LoadTree(bgdInputFile,-1);
  else                               bgdEvent.LoadTree(bgdInputFile,0);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  if	 (period == 0){ // Full2012-Summer12-V9-3500ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    lumi     = 3.553;minRun =      0;maxRun = 999999;
  }
  else if(period == -1){ // Full2012-Summer12-V9 NoJetId
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-04_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-04_V1/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3500ipb.root";
    lumi     = 3.553;minRun =      0;maxRun = 999999;
  }
  else if(period == 1){ // Full2012-Summer12-V9-5000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-06_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_5000ipb_71mb.root";
    lumi     = 5.098;minRun =      0;maxRun = 999999;
    //effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/efficiency_results_HWWICHEP2012.root";
    //fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/fakerate_results_HWWICHEP2012.root";
  }
  else if(period == 2){ // Full2012-Summer12-V9-12000ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-07_HCP_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_12p1ifb.root";
    lumi     = 12.1;minRun =      0;maxRun = 999999;
  }
  else if(period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;minRun =      0;maxRun = 999999;
    //effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/summary.root";
    //fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_HCP.root";
    //puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_12p1ifb.root";
    //lumi     = 12.1;minRun =      0;maxRun = 203002;
    //effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-09/summary_PostHCP_V0.root";
    //fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_PostHCP.root";
    //puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_postHCP.root";
    //lumi     = 7.367;minRun =	  203003;maxRun = 999999;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }
  
  double special[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};special[0] = -9.0;
  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open(fakePath.Data());
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  LeptonScaleLookup trigLookup(effPath.Data());

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int channel = mH;
  if(channel > 1000 && channel < 2000) mH = channel-1000;
  int binc = HiggsMassIndex(mH)-1;
  if(binc == -1) mH = 125;
  bool useDYMVA = true;
  if(channel > 1000 && channel < 2000) useDYMVA = false;

  int newMH = mH;
  if(newMH == 110) newMH = 115; // there is no correction for mh=110!
  if(newMH  > 600) newMH = 600; // there is no correction for mh>600!
  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/ggHWW_KFactors_PowhegToHQT_WithAdditionalMassPoints.root");
  TH1D *HiggsPtKFactor;
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", newMH);
  printf("KFactor_PowhegToHQT: %s\n",kfactorHistName);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;

  double theCutMassLow       = 12;
  double theCutMassHigh	     = cutMassHigh (mH);
  double theCutPtMaxLow	     = cutPtMaxLow (mH);
  double theCutDeltaphilHigh = cutDeltaphiHigh (mH);
  double theCutMTLow         = cutMTLow (mH);
  double theCutMTHigh        = cutMTHigh (mH);
  
  if(nJetsType >= 2) theCutMTLow = 30.0;

  int nBin    = 100;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;

  if     (thePlot >=  7 && thePlot <=  7) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot = 200.0;}
  else if(thePlot >=  8 && thePlot <=  8) {nBinPlot =  30; xminPlot = 60.0; xmaxPlot = 280.0;}
  else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot = 500.0;}
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 30; xminPlot = -0.5; xmaxPlot = 29.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 20 && thePlot <= 22) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 26 && thePlot <= 26) {nBinPlot = 100; xminPlot =-2.5; xmaxPlot = 2.5;}
  else if(thePlot >= 23 && thePlot <= 28) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 29 && thePlot <= 29) {nBinPlot = 100; xminPlot = -1.0; xmaxPlot = 1.0;}
  else if(thePlot >= 30 && thePlot <= 30) {nBinPlot = 200; xminPlot = -5.0; xmaxPlot = 15.0;}
  else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 300; xminPlot = 0.0; xmaxPlot = 600.0;}
  else if(thePlot >= 33 && thePlot <= 33) {nBinPlot = 90; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  1000.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =  2000.0;}
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  1.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 44 && thePlot <= 44) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 49) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot = 0.5;}
  else if(thePlot >= 50 && thePlot <= 50) {nBinPlot = 40; xminPlot = 140.0; xmaxPlot = 180.0;}
  else if(thePlot >= 51 && thePlot <= 51) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot = 50.0;}
  else if(thePlot >= 52 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  else if(thePlot >= 53 && thePlot <= 55) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 56 && thePlot <= 56) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 57 && thePlot <= 57) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}
  else if(thePlot >= 58 && thePlot <= 59) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 60 && thePlot <= 60) {nBinPlot = 40; xminPlot = 0.0; xmaxPlot = 4.0;}
  else if(thePlot >= 61 && thePlot <= 61) {nBinPlot = 180; xminPlot =  0.0; xmaxPlot =  180.0;}
  else if(thePlot >= 62 && thePlot <= 62) {nBinPlot = 72; xminPlot = -1.0; xmaxPlot =  1.0;}
  else if(thePlot >= 63 && thePlot <= 64) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot =  400.0;}
  nBin = nBinPlot;

  double nSigCut[6]  = {0,0,0,0,0,0};
  double nSigECut[6] = {0,0,0,0,0,0};
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

  TH1D* histoggH = (TH1D*) histos->Clone("histoggH");
  TH1D* histoqqH = (TH1D*) histos->Clone("histoqqH");
  TH1D* histoVH  = (TH1D*) histos->Clone("histoVH");

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

  float ewkMVA = -999.;
  bgdEvent.tree_->SetBranchAddress("ewkMVA", &ewkMVA );

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);
    
    if(thePlot == 63 || thePlot == 64){
      //if     (thePlot == 63) histo0->Fill(TMath::Min((bgdEvent.genlep1_+bgdEvent.genlep2_+bgdEvent.genmet_).Pt(),399.999),1.0);
      //else if(thePlot == 64) histo0->Fill(TMath::Min((bgdEvent.genPart1_+bgdEvent.genPart2_+bgdEvent.genPart3_+bgdEvent.genPart4_).Pt(),399.999),1.0);
      continue;
    }

    if(!(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  bgdEvent.dstype_ != SmurfTree::data)) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;if(special[0]<0)assert(0);}
    if(fDecay == -1 || fDecay > 100) fDecay = 0;//44;
    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
    if((fDecay == 27 || fDecay == 28) && channel != 128 && WWXSSel == false && 
       (channel <= 1000 || channel >= 2000)) {
      if(bgdEvent.lep1MotherMcId_ == 23 && bgdEvent.lep2MotherMcId_ == 23) {
        fDecay = 31;
      }
    }
    int centrality = 0;
    if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
       ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    unsigned int Njet3 = bgdEvent.njets_;bool VBFCuts = true;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(bgdEvent.jet3_.pt() <= 30)					                                     Njet3 = 2;
      else if(bgdEvent.jet3_.pt() > 30 && (
    	(bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() < 0) ||
    	(bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() < 0)))   Njet3 = 0;
      else							                                             Njet3 = 2;
      if(bgdEvent.njets_ < 2 || bgdEvent.njets_ > 3)                                                         Njet3 = 0;
      if((channel >= 100 && channel <= 800) || channel == 30 || channel == 31) VBFCuts = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500. && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    }
    double deltaPhiQQL[3] = {DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep1_.Phi()),DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep2_.Phi()),0.0};
    if(bgdEvent.njets_ == 1){deltaPhiQQL[0]=DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep1_.Phi());deltaPhiQQL[1]=DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep2_.Phi());}
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double recoilx   = bgdEvent.met_*cos(bgdEvent.metPhi_)-bgdEvent.lep1_.Px()-bgdEvent.lep2_.Px();
    double recoily   = bgdEvent.met_*sin(bgdEvent.metPhi_)-bgdEvent.lep1_.Py()-bgdEvent.lep2_.Py();
    double recoilPhi = TMath::ATan2(recoily,recoilx);
    double dPhiDilepRecoil = TMath::Abs(bgdEvent.dilep_.Phi()-recoilPhi);
    while(dPhiDilepRecoil>TMath::Pi()) dPhiDilepRecoil = TMath::Abs(dPhiDilepRecoil - 2*TMath::Pi());
    double ptww[3] = {bgdEvent.met_*cos(bgdEvent.metPhi_)+bgdEvent.lep1_.Px()+bgdEvent.lep2_.Px(),
                      bgdEvent.met_*sin(bgdEvent.metPhi_)+bgdEvent.lep1_.Py()+bgdEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);

    bool passNewCuts = true;
    if(bgdEvent.dilep_.Pt() <= 30) passNewCuts = false;
    if(bgdEvent.dilep_.Pt() <= 45 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) passNewCuts = false;
    if(bgdEvent.dilep_.Pt() <= 45 && !(channel == 29 || channel == 30 || channel == 31)) passNewCuts = false;

    double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool   passMET = usedMet > 20.;
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
    if(channel > 1000 && channel < 2000) {passMET = usedMet > 30.; if(bgdEvent.type_ == SmurfTree::ee) passMET = usedMet > 45.;}

    Float_t mTWMin   = bgdEvent.mt1_;
    Float_t mTWMax   = bgdEvent.mt1_;
    if(mTWMin > bgdEvent.mt2_                        ) mTWMin = bgdEvent.mt2_;
    if(mTWMin > bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMin = bgdEvent.mt3_;
    if(mTWMax < bgdEvent.mt2_	 		     ) mTWMax = bgdEvent.mt2_;
    if(mTWMax < bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMax = bgdEvent.mt3_; 
    double pxzll = bgdEvent.lep1_.Px() + bgdEvent.lep2_.Px() + bgdEvent.met_*cos( bgdEvent.metPhi_) + bgdEvent.jet1_.Px() + bgdEvent.jet2_.Px();
    double pyzll = bgdEvent.lep1_.Py() + bgdEvent.lep2_.Py() + bgdEvent.met_*cos( bgdEvent.metPhi_) + bgdEvent.jet1_.Py() + bgdEvent.jet2_.Py();
    double mtHZZ = TMath::Power(TMath::Sqrt((bgdEvent.lep1_.Pt() + bgdEvent.lep2_.Pt() + bgdEvent.met_)*(bgdEvent.lep1_.Pt() + bgdEvent.lep2_.Pt() + bgdEvent.met_)+80.4*80.4) + 
                   		TMath::Sqrt((bgdEvent.jet1_.Pt() + bgdEvent.jet2_.Pt())*(bgdEvent.jet1_.Pt() + bgdEvent.jet2_.Pt())+80.4*80.4),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;

    double mtHWW = TMath::Power(sqrt(bgdEvent.dilep_.Pt()*bgdEvent.dilep_.Pt()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+
                                sqrt(bgdEvent.met_       *bgdEvent.met_       +bgdEvent.dilep_.M()*bgdEvent.dilep_.M()),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHWW >0) mtHWW = sqrt(mtHWW); else mtHWW = 0.0;

    double pznunu = 91.1876*91.1876 - bgdEvent.met_*bgdEvent.met_;
    if(pznunu > 0) pznunu = sqrt(pznunu); else pznunu = 0.0;
    double mHZZ = (sqrt(bgdEvent.dilep_.P()*bgdEvent.dilep_.P()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+91.1876)*
                  (sqrt(bgdEvent.dilep_.P()*bgdEvent.dilep_.P()+bgdEvent.dilep_.M()*bgdEvent.dilep_.M())+91.1876)-
                 -pxzll*pxzll-pyzll*pyzll-(bgdEvent.dilep_.Pz()+pznunu)*(bgdEvent.dilep_.Pz()+pznunu);
    if(mHZZ > 0) mHZZ = sqrt(mHZZ); else mHZZ = 0.001;

    TLorentzVector Lep1(bgdEvent.lep1_.Px(),bgdEvent.lep1_.Py(),bgdEvent.lep1_.Pz(),bgdEvent.lep1_.P());
    TLorentzVector Lep2(bgdEvent.lep2_.Px(),bgdEvent.lep2_.Py(),bgdEvent.lep2_.Pz(),bgdEvent.lep2_.P());
    TVector3 theMet( bgdEvent.met_*cos( bgdEvent.metPhi_),bgdEvent.met_*sin( bgdEvent.metPhi_),0);
    HWWKinematics HWWKin(Lep1,Lep2,theMet);
    TVector3 theWW(bgdEvent.lep1_.Px()+bgdEvent.lep2_.Px()+bgdEvent.met_*cos(bgdEvent.metPhi_),bgdEvent.lep1_.Py()+bgdEvent.lep2_.Py()+bgdEvent.met_*sin(bgdEvent.metPhi_),bgdEvent.lep1_.Pz()+bgdEvent.lep2_.Pz());
    TVector3 theHH(theWW.Px()+(bgdEvent.jet1_+bgdEvent.jet2_).Px(),theWW.Py()+(bgdEvent.jet1_+bgdEvent.jet2_).Py(),theWW.Pz()+(bgdEvent.jet1_+bgdEvent.jet2_).Pz());
    double ptFracHH = TMath::Abs((bgdEvent.jet1_+bgdEvent.jet2_).Pt()-theWW.Pt())/theWW.Pt();
    double dPhiHH = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),theWW.Phi())*180.0/TMath::Pi();
    double dRMinJL = 999.0;
    double dRJL[4] = {DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta()),
                      DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta()),
                      DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta()),
                      DeltaR(bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta())};
    for(int nj=0; nj<4; nj++) {if(dRJL[nj] < dRMinJL) dRMinJL = dRJL[nj];}

    double Mjj = (bgdEvent.jet1_+bgdEvent.jet2_).M(); double qqM = (bgdEvent.jet3_+bgdEvent.jet4_).M(); double qqDeltaEta = TMath::Abs(bgdEvent.jet3_.Eta()-bgdEvent.jet4_.Eta());
    if(bgdEvent.njets_ >= 4){
      mHZZ = (bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.lep1_+bgdEvent.lep2_).M();
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((bgdEvent.jet1_+bgdEvent.jet3_).M()-90.1876)) {Mjj = (bgdEvent.jet1_+bgdEvent.jet3_).M();qqM = (bgdEvent.jet2_+bgdEvent.jet4_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet4_.Eta());mHZZ = (bgdEvent.jet1_+bgdEvent.jet3_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((bgdEvent.jet1_+bgdEvent.jet4_).M()-90.1876)) {Mjj = (bgdEvent.jet1_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet2_+bgdEvent.jet3_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet3_.Eta());mHZZ = (bgdEvent.jet1_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((bgdEvent.jet2_+bgdEvent.jet3_).M()-90.1876)) {Mjj = (bgdEvent.jet2_+bgdEvent.jet3_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet4_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet4_.Eta());mHZZ = (bgdEvent.jet2_+bgdEvent.jet3_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((bgdEvent.jet2_+bgdEvent.jet4_).M()-90.1876)) {Mjj = (bgdEvent.jet2_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet3_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet3_.Eta());mHZZ = (bgdEvent.jet2_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((bgdEvent.jet3_+bgdEvent.jet4_).M()-90.1876)) {Mjj = (bgdEvent.jet3_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet2_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());mHZZ = (bgdEvent.jet3_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet4_.Eta())) {Mjj = (bgdEvent.jet1_+bgdEvent.jet3_).M();qqM = (bgdEvent.jet2_+bgdEvent.jet4_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet4_.Eta());mHZZ = (bgdEvent.jet1_+bgdEvent.jet3_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet3_.Eta())) {Mjj = (bgdEvent.jet1_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet2_+bgdEvent.jet3_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet2_.Eta()-bgdEvent.jet3_.Eta());mHZZ = (bgdEvent.jet1_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet4_.Eta())) {Mjj = (bgdEvent.jet2_+bgdEvent.jet3_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet4_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet4_.Eta());mHZZ = (bgdEvent.jet2_+bgdEvent.jet3_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet3_.Eta())) {Mjj = (bgdEvent.jet2_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet3_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet3_.Eta());mHZZ = (bgdEvent.jet2_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta())) {Mjj = (bgdEvent.jet3_+bgdEvent.jet4_).M();qqM = (bgdEvent.jet1_+bgdEvent.jet2_).M();qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());mHZZ = (bgdEvent.jet3_+bgdEvent.jet4_+bgdEvent.lep1_+bgdEvent.lep2_).M();}
    } else {qqDeltaEta = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());}
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 29 || channel == 30 || channel == 31){ // WW selection
      if(
	  bgdEvent.dilep_.M()   > 12 &&
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          Njet3 == nJetsType &&
          bgdEvent.lep1_.Pt() > 20. &&
          bgdEvent.lep2_.Pt() > ptLepMin &&
          passMET == true &&
	  passNewCuts == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         //(bgdEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets &&
	 //(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.nSoftMuons_ == 0) &&
	 //bgdEvent.jet1Btag_ < 2.1 &&
	 //bgdEvent.jet2Btag_ >= 2.1 &&
	 dPhiDiLepJetCut == true &&
	(channel != 30 || (bgdEvent.mt_ > 30 && (bgdEvent.mt_ > 60 || nJetsType == 2) && bgdEvent.mt_ < 280 && bgdEvent.dilep_.M() < 200 && bgdEvent.lep1_.Pt() > 20. && VBFCuts == true)) &&
	(channel != 31 || (bgdEvent.mt_ > 30 && (bgdEvent.mt_ > 80 || nJetsType == 2) && bgdEvent.mt_ < 600 && bgdEvent.dilep_.M() < 600 && bgdEvent.lep1_.Pt() > 50. && VBFCuts == true)) &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 29 || fDecay == 30 || fDecay == 14) isSignalDecay = true;	     
      }
    } // WW selection

    if(channel == 6){ // Z selection
      //bool passBtagVeto = bgdEvent.nSoftMuons_ == 0 && (bgdEvent.jet1_.Pt() <= 20 || bgdEvent.jet1ProbBtag_ < 0.244) && (bgdEvent.jet2_.Pt() <= 20 || bgdEvent.jet2ProbBtag_ < 0.244) && (bgdEvent.jet3_.Pt() <= 20 || bgdEvent.jet3ProbBtag_ < 0.244) && (bgdEvent.jet4_.Pt() <= 20 || bgdEvent.jet4ProbBtag_ < 0.244);
      bool passBtagVeto = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
      if(
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
        passBtagVeto &&
       //((charge == 0 && VarAux2 == 0) || (charge != 0 && VarAux2 == 1)) &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
	 int((bgdEvent.jet1McId_%1000-bgdEvent.jet1McId_%100)/100) == 0 &&
         Njet3 == nJetsType &&
       //((abs(bgdEvent.lep1_.Eta())>=VarAux0&&abs(bgdEvent.lep1_.Eta())<VarAux0+0.5&&abs(bgdEvent.lep2_.Eta())>=VarAux1&&abs(bgdEvent.lep2_.Eta())<VarAux1+0.5)||
       // (abs(bgdEvent.lep1_.Eta())>=VarAux1&&abs(bgdEvent.lep1_.Eta())<VarAux1+0.5&&abs(bgdEvent.lep2_.Eta())>=VarAux0&&abs(bgdEvent.lep2_.Eta())<VarAux0+0.5)) &&
        (fabs(bgdEvent.dilep_.M()-91.1876) < 15.) && 
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 19) passCuts = false;
	if(fDecay == 9) isSignalDecay = true;
      }
    } // Z selection
    if(channel == 5){ // bbmm selection
      int nb = 0;
      if(bgdEvent.jet1_.Pt() > 30 && bgdEvent.jet1ProbBtag_ > 0.783) nb++;
      if(bgdEvent.jet2_.Pt() > 30 && bgdEvent.jet2ProbBtag_ > 0.783) nb++;
      if(bgdEvent.jet3_.Pt() > 30 && bgdEvent.jet3ProbBtag_ > 0.783) nb++;
      if(bgdEvent.jet4_.Pt() > 30 && bgdEvent.jet4ProbBtag_ > 0.783) nb++;
      if(
	  bgdEvent.dilep_.M() > 12 && bgdEvent.dilep_.M() < 70 && 
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          bgdEvent.njets_ < 2 &&
          bgdEvent.lep1_.Pt() > 25. &&
          bgdEvent.lep2_.Pt() > 25. &&
	  DeltaR(bgdEvent.lep1_.Phi(),bgdEvent.lep1_.Eta(),bgdEvent.lep2_.Phi(),bgdEvent.lep2_.Eta()) < 1.1 &&
          nb >= 1 &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	//
	passCuts = true;
      }
    } // bbmm selection

    if(channel==15||channel==16 ||channel==17||channel==18){ // HHbbww selection
      bool passBtag[3] = {false,false,false};
      if     ( bgdEvent.jet1ProbBtag_ > 0.679 && bgdEvent.jet2ProbBtag_ > 0.679) passBtag[0] = true;
      else if((bgdEvent.jet1ProbBtag_ > 0.679 && bgdEvent.jet2ProbBtag_ > 0.244)||
              (bgdEvent.jet1ProbBtag_ > 0.244 && bgdEvent.jet2ProbBtag_ > 0.679)) passBtag[1] = true;
      else if( bgdEvent.jet1ProbBtag_ > 0.679 || bgdEvent.jet2ProbBtag_ > 0.679)  passBtag[2] = true;
      if(
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          bgdEvent.lep1_.Pt() > 20. &&
          bgdEvent.lep2_.Pt() > 20. &&
          usedMet > 20. &&
	  bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
        ((passBtag[0]&&channel==15)||(passBtag[1]&&channel==16)||(passBtag[2]&&channel==17)||
	((passBtag[0]||passBtag[1]||passBtag[2])&&channel==18)) &&
	  bgdEvent.dilep_.M() > 12 && bgdEvent.dilep_.M() < 200 && bgdEvent.mt_ < 200 && (bgdEvent.jet1_+bgdEvent.jet2_).M() < 200.0 &&
	  bgdEvent.dilep_.M() <  50 && bgdEvent.mt_ > 40 && bgdEvent.mt_ < 125 && fabs((bgdEvent.jet1_+bgdEvent.jet2_).M()-125.0) < 20.0 &&
         (bgdEvent.jet1_+bgdEvent.jet2_).Pt() > 100 && theWW.Pt() > 100 &&
	  dPhiHH > 160.0 && dRMinJL > 1.8 &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	//
	passCuts = true;
      }
    } // HHbbww selection

    if(channel == 9){ // Ztautau selection
      if(
         bgdEvent.dilep_.M() > 12 && bgdEvent.dilep_.M() < 100 &&
         Njet3 == nJetsType &&
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
	 bgdEvent.dilep_.Pt() < 30. &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	if(fDecay == 9 || fDecay == 8) isSignalDecay = true;
      }
    } // Ztautau selection

    if(channel == 902){ // VH,qqlnln selection
      if(
         bgdEvent.dilep_.M()   < 70  &&
         bgdEvent.dilep_.M()   > 12  &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 bgdEvent.njets_ == 2 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
        (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
	 dPhiDiLepJetCut == true &&
        (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 bgdEvent.jet1Btag_ < 1.6 && bgdEvent.jet2Btag_ < 1.6 && 
        (bgdEvent.jet1_+bgdEvent.jet2_).M() > 60. && (bgdEvent.jet1_+bgdEvent.jet2_).M() < 110. &&
         TMath::Abs(bgdEvent.jet1_.Eta()) < 2.5 &&
         TMath::Abs(bgdEvent.jet2_.Eta()) < 2.5 &&
         TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) < 2.1 &&
         bgdEvent.dR_   < 1.3 &&
         bgdEvent.mt_	> 50 &&
         bgdEvent.mt_	< 160 &&
	 DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180/TMath::Pi() > 75.0 &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // VH,qqlnln selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection
      double theCutPtMinLow = cutPtMinLow (mH, bgdEvent.type_);
      if(
         bgdEvent.mt_  > theCutMTLow &&
         bgdEvent.mt_  < theCutMTHigh &&
         bgdEvent.dilep_.M() > theCutMassLow &&
         bgdEvent.dilep_.M() < theCutMassHigh &&
         bgdEvent.lep1_.Pt() > theCutPtMaxLow &&
         bgdEvent.lep2_.Pt() > theCutPtMinLow &&
         bgdEvent.dPhi_*180.0/TMath::Pi() < theCutDeltaphilHigh &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 Njet3 == nJetsType &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dPhiDiLepJetCut == true &&
	 VBFCuts == true &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // H->WW selection

    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
         bgdEvent.dilep_.M() > 12 && bgdEvent.dilep_.M() < 200 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 bgdEvent.njets_ <= 2. && bgdEvent.njets_ >= 0. &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         (bgdEvent.jet1_.Pt() <= 30 || TMath::Abs(bgdEvent.jet1_.Eta()) < 2.5) && (bgdEvent.jet2_.Pt() <= 30 || TMath::Abs(bgdEvent.jet2_.Eta()) < 2.5) && (bgdEvent.jet3_.Pt() <= 30 || TMath::Abs(bgdEvent.jet3_.Eta()) < 2.5) &&
         passMET == true &&
         bgdEvent.mt_ > 0. && bgdEvent.mt_ < 200. && 
         (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
         (dPhiDiLepJetCut == true || bgdEvent.type_ != SmurfTree::ee) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	
	passCuts = true;
	bool isRealLepton = false;
        if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
           (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;
	if(isRealLepton == false &&
	   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
	    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
	    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // HW->2l selection
 
    if(channel == 2000){ // H->ZZ->2l2q selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 qqM > 200 &&
	 qqDeltaEta > 2.5 &&
	 Mjj > 50 && Mjj < 130 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
	 bgdEvent.njets_ >= 4 &&
         abs(bgdEvent.jet1_.Eta()) < 5.0 && abs(bgdEvent.jet2_.Eta()) < 5.0 && abs(bgdEvent.jet3_.Eta()) < 5.0 && abs(bgdEvent.jet4_.Eta()) < 5.0 &&
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // H->ZZ->2l2q selection


    if(passCuts == true && category == 1){
      passCuts = (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection ||
                 (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection;
    }
    if(passCuts == true){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(category == 1) {
        if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV1)  == SmurfTree::Lep1LooseEleV1)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake--;
        if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV1)  == SmurfTree::Lep2LooseEleV1)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake--;
      }
      if(nFake < 0) assert(0);
 
      bool isRealLepton = false;
      if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;
      if(nFake > 1){
	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	fDecay = 22;
	if(category == 1) add = add*1.10; // HACK!!!!
	theWeight	       = -1.0*add;
	bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += theWeight*theWeight;
      }
      else if(nFake == 1){
        if(bgdEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001 && category == 0)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  if(category == 1) add = add*1.10; // HACK!!!!
	  // new category, W+jetsM
	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	  bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += theWeight*theWeight;
	  //cout << "SSS " <<bgdEvent.type_ << " " << add << " " <<bgdEvent.run_ << " " << bgdEvent.event_ << " " <<bgdEvent.lumi_<<endl;
	}
	else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,bgdEvent.npu_);
          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt(), 
								   fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						   TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
          add = add*trigEff;
	  if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001 && category == 0)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  fDecay                 = 3;
	  isSignalDecay          = false;
	  if(category == 1) add = add*1.10; // HACK!!!!
	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add;
          bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += TMath::Abs(theWeight)*TMath::Abs(theWeight);
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					        fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							        TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        double sf_eff = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
        	        leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        theWeight = ZttScaleFactor(period,bgdEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
        bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += TMath::Abs(theWeight)*TMath::Abs(theWeight);
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);
        double add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add2 = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
								 fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						 TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        add = add1*add2*trigEff;
        if(channel > 1000 && channel < 2000 && isRealLepton != true && bgdEvent.dstype_ != SmurfTree::wgamma) add = 0.0;
  	if(special[0] > -1.0 && (fDecay == 29 || fDecay == 30)) {
	  special[0]+=bgdEvent.scale1fb_*lumi*add1;
	  special[1]+=bgdEvent.scale1fb_*lumi*add2;
	  special[2]+=bgdEvent.scale1fb_*lumi*trigEff;
	  special[3]+=bgdEvent.scale1fb_*lumi*add1*add2*trigEff;
	  special[4]+=bgdEvent.scale1fb_*lumi;
	  special[5]+=1.0;
	}

  	//if(bgdEvent.dstype_ == SmurfTree::www) add = bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_;

	//double myWeight = bgdEvent.scale1fb_*lumi; // MC study
	//double myWeight = bgdEvent.scale1fb_*lumi*(bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
        if(fCheckProblem == true && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001 && category == 0)
	  printf("PROBLEMC(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",bgdEvent.event_,add1,add2,trigEff,add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
/*
  	if(channel > 1000 && channel < 2000) {
	  if(bgdEvent.dstype_ == SmurfTree::qqww)  add = add*WrongChargeScaleFactor_VHqqll();
	  if(bgdEvent.dstype_ == SmurfTree::ggww)  add = add*WrongChargeScaleFactor_VHqqll();
	  if(fDecay == 9)                          add = add*WrongChargeScaleFactor_VHqqll()*DYMCScaleFactor_VHqqll(nJetsType);
	  if(bgdEvent.dstype_ == SmurfTree::ttbar) add = add*TopMCScaleFactor_VHqqll(nJetsType);
	  if(bgdEvent.dstype_ == SmurfTree::tw)    add = add*TopMCScaleFactor_VHqqll(nJetsType);
	}
*/
	if((channel == 29 || channel == 30 || channel == 31 || (channel >= 100 && channel <= 800)) &&
	   (channel != 128)){
  	  if(fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) {
  	    if(bgdEvent.njets_ == 0) add=add*DYBkgScaleFactor(0,0);
  	    if(bgdEvent.njets_ == 1) add=add*DYBkgScaleFactor(0,1);
  	    if(bgdEvent.njets_ >= 2) add=add*DYBkgScaleFactor(0,2);
  	  }
  	  if((fDecay == 5 ||fDecay == 13)) {
  	    if     (bgdEvent.njets_ == 0) add=add*TopBkgScaleFactor(0);
  	    else if(bgdEvent.njets_ == 1) add=add*TopBkgScaleFactor(1); 
  	    else if(bgdEvent.njets_ >= 2) add=add*TopBkgScaleFactor(2);
	    if(bgdEvent.njets_ >= 2 && ((channel >= 100 && channel <= 800) || channel == 30 || channel == 31)) add=add*TopVBFBkgScaleFactor(0)/TopBkgScaleFactor(2);
  	  }
  	  if(fDecay == 3) add=add*WJetsMCScaleFactor();

	  if((fDecay == 29 || fDecay == 30 || fDecay == 14) &&
	     ((channel >= 100 && channel <= 800) || WWXSSel == true ||
	       channel == 30 || channel == 31)){
	     if(WWXSSel == true) {
	       //if(bgdEvent.njets_ == 0) add=add*1.21;
	       //else                     add=add*0.90;
	     }
	     else if(channel == 30 || channel == 31){
	       if(bgdEvent.njets_ == 0) add=add*WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mH,115),200),0); 
	       else                     add=add*WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mH,115),200),1);
	     }
	     else {
	       if(bgdEvent.njets_ == 0) add=add*WWBkgScaleFactorCutBased(TMath::Min(TMath::Max((int)mH,115),200),0); 
	       else                     add=add*WWBkgScaleFactorCutBased(TMath::Min(TMath::Max((int)mH,115),200),1);
	     }
          }
	  if((channel >= 100 && channel <= 800) &&
	      fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)){
	    if(mH <= 300 && (lDecay == 4 || lDecay == 5)){
  	      //add=add*DYBkgScaleFactor(TMath::Max((int)mH,115),(int)bgdEvent.njets_)/DYBkgScaleFactor(0,(int)bgdEvent.njets_);
  	      add=0.0; // we will add it later
            }
	  }
	  if((channel == 30 || channel == 31) &&
	      fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)){
	    if(mH <= 300 && (lDecay == 4 || lDecay == 5)){
  	      //add=add*DYBkgScaleFactor(TMath::Max((int)mH,115),(int)bgdEvent.njets_)/DYBkgScaleFactor(0,(int)bgdEvent.njets_);
  	      add=0.0; // we will add it later
            }
	  }
	}

	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,bgdEvent.met_);

        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add;
        bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += theWeight*theWeight;
      }
      double myVar = bgdEvent.met_;
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
      else if(thePlot ==11) myVar = qqM;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = bgdEvent.met_/bgdEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = bgdEvent.njets_;
      else if(thePlot ==18) myVar = bgdEvent.nvtx_;
      else if(thePlot ==19) myVar = bgdEvent.pTrackMet_;
      else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(bgdEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(bgdEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(bgdEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(bgdEvent.lid1_)==13) myVar = bgdEvent.lep1_.Eta();else if(fabs(bgdEvent.lid2_)==13) myVar = bgdEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(TMath::Min((double)ewkMVA,0.999),-0.999);
      else if(thePlot ==30) myVar = TMath::Max(bgdEvent.jet1Btag_,bgdEvent.jet2Btag_);
      else if(thePlot ==31) myVar = 2*HWWKin.CalcMR();
      else if(thePlot ==32) myVar = 2*HWWKin.CalcMRNEW();
      else if(thePlot ==33) myVar = HWWKin.CalcDeltaPhiRFRAME()*180.0/TMath::Pi();
      else if(thePlot ==34) myVar = Mjj;
      else if(thePlot ==35) myVar = qqDeltaEta;
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = TMath::Min((bgdEvent.jet1_+bgdEvent.jet2_).M(),1999.999);
      else if(thePlot ==38) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
      else if(thePlot ==39) myVar = TMath::Min(TMath::Max((1.3491*log((bgdEvent.jet1_+bgdEvent.jet2_).M())-0.01163*(bgdEvent.jet1_.Eta()*bgdEvent.jet2_.Eta())+1.433*TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta())-7.647)/15.00,0.001),0.999);
      else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = (bgdEvent.jet1_+bgdEvent.jet2_).Pt();
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = theWW.Pt();
      else if(thePlot ==48) myVar = bgdEvent.type_;
      else if(thePlot ==49) myVar = ptFracHH;
      else if(thePlot ==50) myVar = dPhiHH;
      else if(thePlot ==51) myVar = theHH.Pt();
      else if(thePlot ==52) myVar = bgdEvent.trackMet_*sin(bgdEvent.trackMetPhi_);
      else if(thePlot ==53) myVar = DeltaPhi(bgdEvent.jet3_.Phi(),bgdEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==54) myVar = dPhiDilepRecoil*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = bgdEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==56) myVar = poorManMetSyst(bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                                   bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,bgdEvent.met_,bgdEvent.metPhi_,
                                    		   bgdEvent.jet1_,bgdEvent.jet2_,bgdEvent.jet3_,bgdEvent.jet4_,0);
      else if(thePlot ==57) myVar = bgdEvent.dR_;
      else if(thePlot ==58) myVar = Lester::higgsWWStarMassBoundNoROOT(bgdEvent.lep1_.e(),bgdEvent.lep1_.px(),bgdEvent.lep1_.py(),bgdEvent.lep1_.pz(),
                                                                       bgdEvent.lep2_.e(),bgdEvent.lep2_.px(),bgdEvent.lep2_.py(),bgdEvent.lep2_.pz(),
								       bgdEvent.met_*cos(bgdEvent.metPhi_),bgdEvent.met_*sin(bgdEvent.metPhi_),80.41);
      else if(thePlot ==59) myVar = mt_atlas(bgdEvent.dilep_,bgdEvent.met_,bgdEvent.metPhi_);
      else if(thePlot ==60) myVar = dRMinJL;
      else if(thePlot ==61) myVar = TMath::Min(DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep1_.Phi()),TMath::Min(DeltaPhi(bgdEvent.jet1_.Phi(),bgdEvent.lep2_.Phi()),TMath::Min(DeltaPhi(bgdEvent.jet2_.Phi(),bgdEvent.lep1_.Phi()),DeltaPhi(bgdEvent.jet2_.Phi(),bgdEvent.lep2_.Phi()))))*180.0/TMath::Pi();

      if     (fDecay == 14 || fDecay == 29 || fDecay == 30){
        histo0->Fill(myVar,theWeight);
      }
      else if(fDecay == 0 || fDecay == 1 || fDecay == 2 || fDecay == 3 || fDecay == 19 || fDecay == 20 || fDecay == 23){
        histo4->Fill(myVar,theWeight);
      }
      else if(fDecay == 6 || fDecay == 7 || fDecay == 8 || fDecay == 9 || fDecay == 10 || fDecay == 31){
        histo1->Fill(myVar,theWeight);
      }
      else if(fDecay == 27 || fDecay == 28 || fDecay == 21){
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

  if(channel >= 100 && channel <= 800){
    if(mH <= 300 && (lDecay == 4 || lDecay == 5)){
      bgdDecay[9] += DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType);
      weiDecay[9] += (DYBkgScaleFactorKappa(TMath::Max((int)mH,115),nJetsType)-1)*DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType)
                    *(DYBkgScaleFactorKappa(TMath::Max((int)mH,115),nJetsType)-1)*DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType);
      B0[0] += DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType);
    }
  }
  if(channel == 30 || channel == 31){
    if(mH <= 300 && (lDecay == 4 || lDecay == 5)){
      bgdDecay[9] += DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType);
      weiDecay[9] += (DYBkgScaleFactorKappa(TMath::Max((int)mH,115),nJetsType)-1)*DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType)
                    *(DYBkgScaleFactorKappa(TMath::Max((int)mH,115),nJetsType)-1)*DYBkgScaleFactor(TMath::Max((int)mH,115),nJetsType);
      B0[0] += DYBkgScaleFactorBDT(TMath::Max((int)mH,115),nJetsType);
    }
  }

  if(special[0] > -1.0){
    printf("%d %11.3f %8.6f %8.6f %8.6f %8.6f %11.3f\n",(int)special[5],special[4],special[0]/special[4],special[1]/special[4],special[2]/special[4],special[3]/special[4],special[3]);
  }

  if((channel >= 0 && channel <= 8000)){
    sigEvent.tree_->SetBranchAddress("ewkMVA", &ewkMVA );
    int nSig=sigEvent.tree_->GetEntries();
    for (int i=0; i<nSig; ++i) {

     if (i%100000 == 0 && verboseLevel > 0)
	printf("--- reading Signal event %5d of %5d\n",i,nSig);
     sigEvent.tree_->GetEntry(i);

     bool lId = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
     if(category == 1) lId = ((sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (sigEvent.cuts_ & SmurfTree::Lep2LooseEleV1)    == SmurfTree::Lep2LooseEleV1   ) ||
                             ((sigEvent.cuts_ & SmurfTree::Lep1LooseEleV1)    == SmurfTree::Lep1LooseEleV1    && (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
     if(!lId) continue;

    int charge = (int)(sigEvent.lq1_ + sigEvent.lq2_);

    int centrality = 0;
    if(((sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() < 0)) &&
       ((sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() < 0))) centrality = 1; 
    unsigned int Njet3 = sigEvent.njets_;bool VBFCuts = true;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(sigEvent.jet3_.pt() <= 30)					                                     Njet3 = 2;
      else if(sigEvent.jet3_.pt() > 30 && (
    	(sigEvent.jet1_.eta()-sigEvent.jet3_.eta() > 0 && sigEvent.jet2_.eta()-sigEvent.jet3_.eta() < 0) ||
    	(sigEvent.jet2_.eta()-sigEvent.jet3_.eta() > 0 && sigEvent.jet1_.eta()-sigEvent.jet3_.eta() < 0)))   Njet3 = 0;
      else							                                             Njet3 = 2;
      if(sigEvent.njets_ < 2 || sigEvent.njets_ > 3)                                                         Njet3 = 0;
      if((channel >= 100 && channel <= 800) || channel == 30 || channel == 31) VBFCuts = (sigEvent.jet1_+sigEvent.jet2_).M() > 500. && TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    }
    double deltaPhiQQL[3] = {DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep1_.Phi()),DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep2_.Phi()),0.0};
    if(sigEvent.njets_ == 1){deltaPhiQQL[0]=DeltaPhi(sigEvent.jet1_.Phi(),sigEvent.lep1_.Phi());deltaPhiQQL[1]=DeltaPhi(sigEvent.jet1_.Phi(),sigEvent.lep2_.Phi());}
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double recoilx   = sigEvent.met_*cos(sigEvent.metPhi_)-sigEvent.lep1_.Px()-sigEvent.lep2_.Px();
    double recoily   = sigEvent.met_*sin(sigEvent.metPhi_)-sigEvent.lep1_.Py()-sigEvent.lep2_.Py();
    double recoilPhi = TMath::ATan2(recoily,recoilx);
    double dPhiDilepRecoil = TMath::Abs(sigEvent.dilep_.Phi()-recoilPhi);
    while(dPhiDilepRecoil>TMath::Pi()) dPhiDilepRecoil = TMath::Abs(dPhiDilepRecoil - 2*TMath::Pi());
    double ptww[3] = {sigEvent.met_*cos(sigEvent.metPhi_)+sigEvent.lep1_.Px()+sigEvent.lep2_.Px(),
                      sigEvent.met_*sin(sigEvent.metPhi_)+sigEvent.lep1_.Py()+sigEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);

    bool passNewCuts = true;
    if(sigEvent.dilep_.Pt() <= 30) passNewCuts = false;
    if(sigEvent.dilep_.Pt() <= 45 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) passNewCuts = false;
    if(sigEvent.dilep_.Pt() <= 45 && !(channel == 29 || channel == 30 ||channel == 31)) passNewCuts = false;

    double usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
    bool   passMET = usedMet > 20.;
    if(useDYMVA == false){
      if     (sigEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
      else if(sigEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (sigEvent.met_ > 45.0 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
    } else {
      if     (sigEvent.njets_ == 0) passMET = passMET && (sigEvent.dymva_ >  0.88 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
      else if(sigEvent.njets_ == 1) passMET = passMET && (sigEvent.dymva_ >  0.84 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (sigEvent.met_   >  45.0 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);
    }
    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(sigEvent.njets_ <= 1) dPhiDiLepJetCut = sigEvent.jet1_.Pt() <= 15. || sigEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
      	                                         sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
    	                                         sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
    }
    if(sigEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
    if(channel > 1000 && channel < 2000) {passMET = usedMet > 30.; if(sigEvent.type_ == SmurfTree::ee) passMET = usedMet > 45.;}

    TLorentzVector Lep1(sigEvent.lep1_.Px(),sigEvent.lep1_.Py(),sigEvent.lep1_.Pz(),sigEvent.lep1_.P());
    TLorentzVector Lep2(sigEvent.lep2_.Px(),sigEvent.lep2_.Py(),sigEvent.lep2_.Pz(),sigEvent.lep2_.P());
    TVector3 theMet( sigEvent.met_*cos( sigEvent.metPhi_),sigEvent.met_*sin( sigEvent.metPhi_),0);
    HWWKinematics HWWKin(Lep1,Lep2,theMet);
    TVector3 theWW(sigEvent.lep1_.Px()+sigEvent.lep2_.Px()+sigEvent.met_*cos(sigEvent.metPhi_),sigEvent.lep1_.Py()+sigEvent.lep2_.Py()+sigEvent.met_*sin(sigEvent.metPhi_),sigEvent.lep1_.Pz()+sigEvent.lep2_.Pz());
    TVector3 theHH(theWW.Px()+(sigEvent.jet1_+sigEvent.jet2_).Px(),theWW.Py()+(sigEvent.jet1_+sigEvent.jet2_).Py(),theWW.Pz()+(sigEvent.jet1_+sigEvent.jet2_).Pz());
    double ptFracHH = TMath::Abs((sigEvent.jet1_+sigEvent.jet2_).Pt()-theWW.Pt())/theWW.Pt();
    double dPhiHH = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),theWW.Phi())*180.0/TMath::Pi();
    double dRMinJL = 999.0;
    double dRJL[4] = {DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.lep1_.Phi(),sigEvent.lep1_.Eta()),
                      DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.lep2_.Phi(),sigEvent.lep2_.Eta()),
                      DeltaR(sigEvent.jet2_.Phi(),sigEvent.jet2_.Eta(),sigEvent.lep1_.Phi(),sigEvent.lep1_.Eta()),
                      DeltaR(sigEvent.jet2_.Phi(),sigEvent.jet2_.Eta(),sigEvent.lep2_.Phi(),sigEvent.lep2_.Eta())};
    for(int nj=0; nj<4; nj++) {if(dRJL[nj] < dRMinJL) dRMinJL = dRJL[nj];}

    Float_t mTWMin   = sigEvent.mt1_;
    Float_t mTWMax   = sigEvent.mt1_;
    if(mTWMin > sigEvent.mt2_			     ) mTWMin = sigEvent.mt2_;
    if(mTWMin > sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMin = sigEvent.mt3_;
    if(mTWMax < sigEvent.mt2_			     ) mTWMax = sigEvent.mt2_;
    if(mTWMax < sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMax = sigEvent.mt3_;
    double pxzll = sigEvent.dilep_.Px() + sigEvent.met_*cos( sigEvent.metPhi_);
    double pyzll = sigEvent.dilep_.Py() + sigEvent.met_*sin( sigEvent.metPhi_);
    double mtHZZ = TMath::Power(TMath::Sqrt((sigEvent.lep1_.Pt() + sigEvent.lep2_.Pt() + sigEvent.met_)*(sigEvent.lep1_.Pt() + sigEvent.lep2_.Pt() + sigEvent.met_)+80.4*80.4) + 
                   		TMath::Sqrt((sigEvent.jet1_.Pt() + sigEvent.jet2_.Pt())*(sigEvent.jet1_.Pt() + sigEvent.jet2_.Pt())+80.4*80.4),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHZZ >0) mtHZZ = sqrt(mtHZZ); else mtHZZ = 0.0;
    double mtHWW = TMath::Power(sqrt(sigEvent.dilep_.Pt()*sigEvent.dilep_.Pt()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+
                                sqrt(sigEvent.met_       *sigEvent.met_       +sigEvent.dilep_.M()*sigEvent.dilep_.M()),2)
		 -pxzll*pxzll-pyzll*pyzll;
    if(mtHWW >0) mtHWW = sqrt(mtHWW); else mtHWW = 0.0;

    double pznunu = 91.1876*91.1876 - sigEvent.met_*sigEvent.met_;
    if(pznunu > 0) pznunu = sqrt(pznunu); else pznunu = 0.0;
    double mHZZ = (sqrt(sigEvent.dilep_.P()*sigEvent.dilep_.P()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+91.1876)*
                  (sqrt(sigEvent.dilep_.P()*sigEvent.dilep_.P()+sigEvent.dilep_.M()*sigEvent.dilep_.M())+91.1876)-
                 -pxzll*pxzll-pyzll*pyzll-(sigEvent.dilep_.Pz()+pznunu)*(sigEvent.dilep_.Pz()+pznunu);
    if(mHZZ > 0) mHZZ = sqrt(mHZZ); else mHZZ = 0.001;

    double Mjj = (sigEvent.jet1_+sigEvent.jet2_).M(); double qqM = (sigEvent.jet3_+sigEvent.jet4_).M(); double qqDeltaEta = TMath::Abs(sigEvent.jet3_.Eta()-sigEvent.jet4_.Eta());
    if(sigEvent.njets_ >= 4){
      mHZZ = (sigEvent.jet1_+sigEvent.jet2_+sigEvent.lep1_+sigEvent.lep2_).M();
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((sigEvent.jet1_+sigEvent.jet3_).M()-90.1876)) {Mjj = (sigEvent.jet1_+sigEvent.jet3_).M();qqM = (sigEvent.jet2_+sigEvent.jet4_).M();qqDeltaEta = TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet4_.Eta());mHZZ = (sigEvent.jet1_+sigEvent.jet3_+sigEvent.lep1_+sigEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((sigEvent.jet1_+sigEvent.jet4_).M()-90.1876)) {Mjj = (sigEvent.jet1_+sigEvent.jet4_).M();qqM = (sigEvent.jet2_+sigEvent.jet3_).M();qqDeltaEta = TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet3_.Eta());mHZZ = (sigEvent.jet1_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((sigEvent.jet2_+sigEvent.jet3_).M()-90.1876)) {Mjj = (sigEvent.jet2_+sigEvent.jet3_).M();qqM = (sigEvent.jet1_+sigEvent.jet4_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet4_.Eta());mHZZ = (sigEvent.jet2_+sigEvent.jet3_+sigEvent.lep1_+sigEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((sigEvent.jet2_+sigEvent.jet4_).M()-90.1876)) {Mjj = (sigEvent.jet2_+sigEvent.jet4_).M();qqM = (sigEvent.jet1_+sigEvent.jet3_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet3_.Eta());mHZZ = (sigEvent.jet2_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((sigEvent.jet3_+sigEvent.jet4_).M()-90.1876)) {Mjj = (sigEvent.jet3_+sigEvent.jet4_).M();qqM = (sigEvent.jet1_+sigEvent.jet2_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());mHZZ = (sigEvent.jet3_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet4_.Eta())) {Mjj = (sigEvent.jet1_+sigEvent.jet3_).M();qqM = (sigEvent.jet2_+sigEvent.jet4_).M();qqDeltaEta = TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet4_.Eta());mHZZ = (sigEvent.jet1_+sigEvent.jet3_+sigEvent.lep1_+sigEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet3_.Eta())) {Mjj = (sigEvent.jet1_+sigEvent.jet4_).M();qqM = (sigEvent.jet2_+sigEvent.jet3_).M();qqDeltaEta = TMath::Abs(sigEvent.jet2_.Eta()-sigEvent.jet3_.Eta());mHZZ = (sigEvent.jet1_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet4_.Eta())) {Mjj = (sigEvent.jet2_+sigEvent.jet3_).M();qqM = (sigEvent.jet1_+sigEvent.jet4_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet4_.Eta());mHZZ = (sigEvent.jet2_+sigEvent.jet3_+sigEvent.lep1_+sigEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet3_.Eta())) {Mjj = (sigEvent.jet2_+sigEvent.jet4_).M();qqM = (sigEvent.jet1_+sigEvent.jet3_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet3_.Eta());mHZZ = (sigEvent.jet2_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta())) {Mjj = (sigEvent.jet3_+sigEvent.jet4_).M();qqM = (sigEvent.jet1_+sigEvent.jet2_).M();qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());mHZZ = (sigEvent.jet3_+sigEvent.jet4_+sigEvent.lep1_+sigEvent.lep2_).M();}
    } else {qqDeltaEta = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());}
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 29 || channel == 30 || channel == 31){ // WW selection
      if(
	  sigEvent.dilep_.M()   > 12 &&
         (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          Njet3 == nJetsType &&
          sigEvent.lep1_.Pt() > 20. &&
          sigEvent.lep2_.Pt() > ptLepMin &&
          passMET == true &&
	  passNewCuts == true &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         //(sigEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets &&
	 //(sigEvent.jetLowBtag_ < 2.1 && sigEvent.nSoftMuons_ == 0) &&
	 //sigEvent.jet1Btag_ < 2.1 &&
	 //sigEvent.jet2Btag_ >= 2.1 &&
	 dPhiDiLepJetCut == true &&
	(channel != 30 || (sigEvent.mt_ > 30 && (sigEvent.mt_ > 60 || nJetsType == 2) && sigEvent.mt_ < 280 && sigEvent.dilep_.M() < 200 && sigEvent.lep1_.Pt() > 20. && VBFCuts == true)) &&
	(channel != 31 || (sigEvent.mt_ > 30 && (sigEvent.mt_ > 80 || nJetsType == 2) && sigEvent.mt_ < 600 && sigEvent.dilep_.M() < 600 && sigEvent.lep1_.Pt() > 50. && VBFCuts == true)) &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WW selection

    if(channel==15||channel==16 ||channel==17||channel==18){ // HHbbww selection
      bool passBtag[3] = {false,false,false};
      if     ( sigEvent.jet1ProbBtag_ > 0.679 && sigEvent.jet2ProbBtag_ > 0.679) passBtag[0] = true;
      else if((sigEvent.jet1ProbBtag_ > 0.679 && sigEvent.jet2ProbBtag_ > 0.244)||
              (sigEvent.jet1ProbBtag_ > 0.244 && sigEvent.jet2ProbBtag_ > 0.679)) passBtag[1] = true;
      else if( sigEvent.jet1ProbBtag_ > 0.679 || sigEvent.jet2ProbBtag_ > 0.679)  passBtag[2] = true;
      if(
         (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          sigEvent.lep1_.Pt() > 20. &&
          sigEvent.lep2_.Pt() > 20. &&
          usedMet > 20. &&
	  sigEvent.njets_ >= 2 && sigEvent.njets_ <= 3 &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
        ((passBtag[0]&&channel==15)||(passBtag[1]&&channel==16)||(passBtag[2]&&channel==17)||
	((passBtag[0]||passBtag[1]||passBtag[2])&&channel==18)) &&
	  sigEvent.dilep_.M() > 12 && sigEvent.dilep_.M() < 200 && sigEvent.mt_ < 200 && (sigEvent.jet1_+sigEvent.jet2_).M() < 200.0 &&
	  sigEvent.dilep_.M() <  50 && sigEvent.mt_ > 40 && sigEvent.mt_ < 125 && fabs((sigEvent.jet1_+sigEvent.jet2_).M()-125.0) < 25.0 &&
         (sigEvent.jet1_+sigEvent.jet2_).Pt() > 100 && theWW.Pt() > 100 &&
	  dPhiHH > 160.0 && dRMinJL > 1.8 &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	//
	passCuts = true;
        isSignalDecay = true;
      }
    } // HHbbww selection

    if(channel == 9){ // Ztautau selection
      if(
         sigEvent.dilep_.M() > 12 && sigEvent.dilep_.M() < 100 &&
         Njet3 == nJetsType &&
         (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
	 sigEvent.dilep_.Pt() < 30. &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Ztautau selection

    if(channel == 902){ // VH,qqlnln selection
      if(
         sigEvent.dilep_.M()   < 70  &&
         sigEvent.dilep_.M()   > 12  &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 sigEvent.njets_ == 2 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
        (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
	 dPhiDiLepJetCut == true &&
        (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 sigEvent.jet1Btag_ < 1.6 && sigEvent.jet2Btag_ < 1.6 && 
        (sigEvent.jet1_+sigEvent.jet2_).M() > 60. && (sigEvent.jet1_+sigEvent.jet2_).M() < 110. &&
         TMath::Abs(sigEvent.jet1_.Eta()) < 2.5 &&
         TMath::Abs(sigEvent.jet2_.Eta()) < 2.5 &&
         TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta()) < 2.1 &&
         sigEvent.dR_   < 1.3 &&
         sigEvent.mt_	> 50 &&
         sigEvent.mt_	< 160 &&
	 DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180/TMath::Pi() > 75.0 &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // VH,qqlnln selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection
      double theCutPtMinLow = cutPtMinLow (mH, sigEvent.type_);
      if(
         sigEvent.mt_  > theCutMTLow &&
         sigEvent.mt_  < theCutMTHigh &&
         sigEvent.dilep_.M() > theCutMassLow &&
         sigEvent.dilep_.M() < theCutMassHigh &&
         sigEvent.lep1_.Pt() > theCutPtMaxLow &&
         sigEvent.lep2_.Pt() > theCutPtMinLow &&
         sigEvent.dPhi_*180.0/TMath::Pi() < theCutDeltaphilHigh &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 Njet3 == nJetsType &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 15. || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dPhiDiLepJetCut == true &&
	 VBFCuts == true &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	  1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // H->WW selection

    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
         sigEvent.dilep_.M() > 12 && sigEvent.dilep_.M() < 200 &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 sigEvent.njets_ <= 2. && sigEvent.njets_ >= 0. &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         (sigEvent.jet1_.Pt() <= 30 || TMath::Abs(sigEvent.jet1_.Eta()) < 2.5) && (sigEvent.jet2_.Pt() <= 30 || TMath::Abs(sigEvent.jet2_.Eta()) < 2.5) && (sigEvent.jet3_.Pt() <= 30 || TMath::Abs(sigEvent.jet3_.Eta()) < 2.5) &&
         passMET == true &&
         sigEvent.mt_ > 0. && sigEvent.mt_ < 200. && 
         (fabs(sigEvent.dilep_.M()-91.1876) > 10. || sigEvent.type_ != SmurfTree::ee) && 
         (dPhiDiLepJetCut == true || sigEvent.type_ != SmurfTree::ee) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
	bool isRealLepton = false;
        if((TMath::Abs(sigEvent.lep1McId_) == 11 || TMath::Abs(sigEvent.lep1McId_) == 13) &&
           (TMath::Abs(sigEvent.lep2McId_) == 11 || TMath::Abs(sigEvent.lep2McId_) == 13)) isRealLepton = true;
	if(isRealLepton == false) {passCuts = false; isSignalDecay = false;}
      }
    } // HW->2l selection
    
    if(channel == 2000){ // H->ZZ->2l2q selection
      if(
         sigEvent.dilep_.M()   > 12 &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 qqM > 200 &&
	 qqDeltaEta > 2.5 &&
	 Mjj > 50 && Mjj < 130 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
	 sigEvent.njets_ >= 4 &&
         abs(sigEvent.jet1_.Eta()) < 5.0 && abs(sigEvent.jet2_.Eta()) < 5.0 && abs(sigEvent.jet3_.Eta()) < 5.0 && abs(sigEvent.jet4_.Eta()) < 5.0 &&
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // H->ZZ->2l2q selection

    //if(passCuts == true && sigEvent.processId_==10001){
    if(passCuts == true){
      double add = 1.;
      double addPU = 1.;
      addPU = nPUScaleFactor2012(fhDPU,sigEvent.npu_);
      add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
      add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
     							       fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
	   						       TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
      add = add*trigEff*addPU;

      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;

      if(isFermioPhobic == true) {
        add = add * enhancementFactor(mH,2); // FF BR(H->WW) enhancement factor
	if(sigEvent.processId_==121 || sigEvent.processId_==122) add = 0.0;
      }

      if(isSM4 == true) add = add * enhancementFactor(mH,1); // SM4 BR(H->WW) enhancement factor

      double addIntf = 1.0;
      if (sigEvent.processId_ == 10010) {
        // new interference term
        addIntf = InterfgHHSystematics(mH);

        //add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(sigEvent.higgsPt_));
	if(isFermioPhobic == true) add = 0.0;
	if(isSM4 == true) add = add * enhancementFactor(mH,0); // SM4 ggH enhancement factor
      }

      double myWeight = sigEvent.scale1fb_*lumi*add*addIntf;
      //double myWeight = sigEvent.scale1fb_*lumi*(sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_*sigEvent.sfWeightHPt_);
      if(fCheckProblem == true && TMath::Abs((sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_)-add)/add>0.0001 && category == 0) {
        printf("PROBLEM: %f - %f %f %f %f %f = %f\n",add,sigEvent.sfWeightFR_,sigEvent.sfWeightPU_,sigEvent.sfWeightEff_,sigEvent.sfWeightTrig_,sigEvent.sfWeightHPt_,sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_*sigEvent.sfWeightHPt_);
        printf("		  %f %f %f %f %f\n",1.0,addPU,leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_)*
              leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_), trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
        																		     fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        																		     TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_)),
               HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(sigEvent.higgsPt_)));
      }
      //if(sigEvent.event_==89870){
      //  printf("SSS %d %f %f %f\n",sigEvent.event_,2*HWWKin.CalcMR(),2*HWWKin.CalcMRNEW(),TMath::Abs(HWWKin.CalcDeltaPhiRFRAME()));
      //}
      double myVar = sigEvent.met_;
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
      else if(thePlot ==11) myVar = qqM;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = sigEvent.met_/sigEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = sigEvent.njets_;
      else if(thePlot ==18) myVar = sigEvent.nvtx_;
      else if(thePlot ==19) myVar = sigEvent.pTrackMet_;
      else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(sigEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(sigEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(sigEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(sigEvent.lid1_)==13) myVar = sigEvent.lep1_.Eta();else if(fabs(sigEvent.lid2_)==13) myVar = sigEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(TMath::Min((double)ewkMVA,0.999),-0.999);
      else if(thePlot ==30) myVar = TMath::Max(sigEvent.jet1Btag_,sigEvent.jet2Btag_);
      else if(thePlot ==31) myVar = 2*HWWKin.CalcMR();
      else if(thePlot ==32) myVar = 2*HWWKin.CalcMRNEW();
      else if(thePlot ==33) myVar = HWWKin.CalcDeltaPhiRFRAME()*180.0/TMath::Pi();
      else if(thePlot ==34) myVar = Mjj;
      else if(thePlot ==35) myVar = qqDeltaEta;
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = TMath::Min((sigEvent.jet1_+sigEvent.jet2_).M(),1999.999);
      else if(thePlot ==38) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
      else if(thePlot ==39) myVar = TMath::Min(TMath::Max((1.3491*log((sigEvent.jet1_+sigEvent.jet2_).M())-0.01163*(sigEvent.jet1_.Eta()*sigEvent.jet2_.Eta())+1.433*TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta())-7.647)/15.00,0.001),0.999);
      else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,sigEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = (sigEvent.jet1_+sigEvent.jet2_).Pt();
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = theWW.Pt();
      else if(thePlot ==48) myVar = sigEvent.type_;
      else if(thePlot ==49) myVar = ptFracHH;
      else if(thePlot ==50) myVar = dPhiHH;
      else if(thePlot ==51) myVar = theHH.Pt();
      else if(thePlot ==52) myVar = sigEvent.trackMet_*sin(sigEvent.trackMetPhi_);
      else if(thePlot ==53) myVar = DeltaPhi(sigEvent.jet3_.Phi(),sigEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==54) myVar = dPhiDilepRecoil*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = sigEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==56) myVar = poorManMetSyst(sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                                   sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,sigEvent.met_,sigEvent.metPhi_,
                                    		   sigEvent.jet1_,sigEvent.jet2_,sigEvent.jet3_,sigEvent.jet4_,0);
      else if(thePlot ==57) myVar = sigEvent.dR_;
      else if(thePlot ==58) myVar = Lester::higgsWWStarMassBoundNoROOT(sigEvent.lep1_.e(),sigEvent.lep1_.px(),sigEvent.lep1_.py(),sigEvent.lep1_.pz(),
                                                                       sigEvent.lep2_.e(),sigEvent.lep2_.px(),sigEvent.lep2_.py(),sigEvent.lep2_.pz(),
								       sigEvent.met_*cos(sigEvent.metPhi_),sigEvent.met_*sin(sigEvent.metPhi_),80.41);
      else if(thePlot ==59) myVar = mt_atlas(sigEvent.dilep_,sigEvent.met_,sigEvent.metPhi_);
      else if(thePlot ==60) myVar = dRMinJL;
      else if(thePlot ==61) myVar = TMath::Min(DeltaPhi(sigEvent.jet1_.Phi(),sigEvent.lep1_.Phi()),TMath::Min(DeltaPhi(sigEvent.jet1_.Phi(),sigEvent.lep2_.Phi()),TMath::Min(DeltaPhi(sigEvent.jet2_.Phi(),sigEvent.lep1_.Phi()),DeltaPhi(sigEvent.jet2_.Phi(),sigEvent.lep2_.Phi()))))*180.0/TMath::Pi();

      int nSigBin = -1;
      if     (sigEvent.processId_==121 ||
    	      sigEvent.processId_==122)   nSigBin = 1;
      else if(sigEvent.processId_==24)    nSigBin = 2;
      else if(sigEvent.processId_==26)    nSigBin = 3;
      else if(sigEvent.processId_==10001) nSigBin = 4;
      else if(sigEvent.processId_==10010) nSigBin = 5;
      else  {if(channel != 1200) return;}
      nSigCut[0]  = nSigCut[0]  + myWeight;
      nSigECut[0] = nSigECut[0] + myWeight*myWeight;
      nSigCut[nSigBin]  = nSigCut[nSigBin]  + myWeight;
      nSigECut[nSigBin] = nSigECut[nSigBin] + myWeight*myWeight;

      histos->Fill(myVar,myWeight);

      if     (sigEvent.processId_==10010) {hDHSig[0]->Fill(1.0,myWeight); histoggH->Fill(myVar,myWeight);}
      else if(sigEvent.processId_==10001) {hDHSig[1]->Fill(1.0,myWeight); histoqqH->Fill(myVar,myWeight);}
      else if(sigEvent.processId_==26  ||
              sigEvent.processId_==24  ||
              sigEvent.processId_==121 ||
              sigEvent.processId_==122)   {hDHSig[2]->Fill(1.0,myWeight); histoVH ->Fill(myVar,myWeight);}

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
  dataEvent.tree_->SetBranchAddress("ewkMVA", &ewkMVA );
  int nData=dataEvent.tree_->GetEntries();
  double nSelectedData = 0;
  for (int i=0; i<nData; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if(category == 1) lId = ((dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2LooseEleV1)    == SmurfTree::Lep2LooseEleV1   ) ||
    			    ((dataEvent.cuts_ & SmurfTree::Lep1LooseEleV1)    == SmurfTree::Lep1LooseEleV1    && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;
    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    unsigned int Njet3 = dataEvent.njets_;bool VBFCuts = true;
    if(nJetsType == 2){ // nJetsType = 0/1/2-jet selection
      if(dataEvent.jet3_.pt() <= 30)					                                     Njet3 = 2;
      else if(dataEvent.jet3_.pt() > 30 && (
    	(dataEvent.jet1_.eta()-dataEvent.jet3_.eta() > 0 && dataEvent.jet2_.eta()-dataEvent.jet3_.eta() < 0) ||
    	(dataEvent.jet2_.eta()-dataEvent.jet3_.eta() > 0 && dataEvent.jet1_.eta()-dataEvent.jet3_.eta() < 0)))   Njet3 = 0;
      else							                                             Njet3 = 2;
      if(dataEvent.njets_ < 2 || dataEvent.njets_ > 3)                                                         Njet3 = 0;
      if((channel >= 100 && channel <= 800) || channel == 30 || channel == 31) VBFCuts = (dataEvent.jet1_+dataEvent.jet2_).M() > 500. && TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 && centrality == 1;
    }
    double deltaPhiQQL[3] = {DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep1_.Phi()),DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep2_.Phi()),0.0};
    if(dataEvent.njets_ == 1){deltaPhiQQL[0]=DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep1_.Phi());deltaPhiQQL[1]=DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep2_.Phi());}
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double recoilx   = dataEvent.met_*cos(dataEvent.metPhi_)-dataEvent.lep1_.Px()-dataEvent.lep2_.Px();
    double recoily   = dataEvent.met_*sin(dataEvent.metPhi_)-dataEvent.lep1_.Py()-dataEvent.lep2_.Py();
    double recoilPhi = TMath::ATan2(recoily,recoilx);
    double dPhiDilepRecoil = TMath::Abs(dataEvent.dilep_.Phi()-recoilPhi);
    while(dPhiDilepRecoil>TMath::Pi()) dPhiDilepRecoil = TMath::Abs(dPhiDilepRecoil - 2*TMath::Pi());
    double ptww[3] = {dataEvent.met_*cos(dataEvent.metPhi_)+dataEvent.lep1_.Px()+dataEvent.lep2_.Px(),
                      dataEvent.met_*sin(dataEvent.metPhi_)+dataEvent.lep1_.Py()+dataEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);

    bool passNewCuts = true;
    if(dataEvent.dilep_.Pt() <= 30) passNewCuts = false;
    if(dataEvent.dilep_.Pt() <= 45 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) passNewCuts = false;
    if(dataEvent.dilep_.Pt() <= 45 && !(channel == 29 || channel == 30 || channel == 31)) passNewCuts = false;

    double usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool   passMET = usedMet > 20.;
    if(useDYMVA == false){
      if     (dataEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                           passMET = passMET && (dataEvent.met_ > 45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    } else {
      if     (dataEvent.njets_ == 0) passMET = passMET && (dataEvent.dymva_ >  0.88 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (dataEvent.dymva_ >  0.84 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                           passMET = passMET && (dataEvent.met_   >  45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(dataEvent.njets_ <= 1) dPhiDiLepJetCut = dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
      	                                          dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
      else                      dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
    	                                          dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    }
    if(dataEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    if(channel > 1000 && channel < 2000) {passMET = usedMet > 30.; if(dataEvent.type_ == SmurfTree::ee) passMET = usedMet > 45.;}

    TLorentzVector Lep1(dataEvent.lep1_.Px(),dataEvent.lep1_.Py(),dataEvent.lep1_.Pz(),dataEvent.lep1_.P());
    TLorentzVector Lep2(dataEvent.lep2_.Px(),dataEvent.lep2_.Py(),dataEvent.lep2_.Pz(),dataEvent.lep2_.P());
    TVector3 theMet( dataEvent.met_*cos( dataEvent.metPhi_),dataEvent.met_*sin( dataEvent.metPhi_),0);
    HWWKinematics HWWKin(Lep1,Lep2,theMet);
    TVector3 theWW(dataEvent.lep1_.Px()+dataEvent.lep2_.Px()+dataEvent.met_*cos(dataEvent.metPhi_),dataEvent.lep1_.Py()+dataEvent.lep2_.Py()+dataEvent.met_*sin(dataEvent.metPhi_),dataEvent.lep1_.Pz()+dataEvent.lep2_.Pz());
    TVector3 theHH(theWW.Px()+(dataEvent.jet1_+dataEvent.jet2_).Px(),theWW.Py()+(dataEvent.jet1_+dataEvent.jet2_).Py(),theWW.Pz()+(dataEvent.jet1_+dataEvent.jet2_).Pz());
    double ptFracHH = TMath::Abs((dataEvent.jet1_+dataEvent.jet2_).Pt()-theWW.Pt())/theWW.Pt();
    double dPhiHH = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),theWW.Phi())*180.0/TMath::Pi();
    double dRMinJL = 999.0;
    double dRJL[4] = {DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta()),
                      DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta()),
                      DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta()),
                      DeltaR(dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta())};
    for(int nj=0; nj<4; nj++) {if(dRJL[nj] < dRMinJL) dRMinJL = dRJL[nj];}

    Float_t mTWMin   = dataEvent.mt1_;
    Float_t mTWMax   = dataEvent.mt1_;
    if(mTWMin > dataEvent.mt2_  		       ) mTWMin = dataEvent.mt2_;
    if(mTWMin > dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMin = dataEvent.mt3_;
    if(mTWMax < dataEvent.mt2_  		       ) mTWMax = dataEvent.mt2_;
    if(mTWMax < dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMax = dataEvent.mt3_;
    double pxzll = dataEvent.dilep_.Px() + dataEvent.met_*cos( dataEvent.metPhi_);
    double pyzll = dataEvent.dilep_.Py() + dataEvent.met_*sin( dataEvent.metPhi_);
    double mtHZZ = TMath::Power(TMath::Sqrt((dataEvent.lep1_.Pt() + dataEvent.lep2_.Pt() + dataEvent.met_)*(dataEvent.lep1_.Pt() + dataEvent.lep2_.Pt() + dataEvent.met_)+80.4*80.4) + 
                   		TMath::Sqrt((dataEvent.jet1_.Pt() + dataEvent.jet2_.Pt())*(dataEvent.jet1_.Pt() + dataEvent.jet2_.Pt())+80.4*80.4),2)
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

    double Mjj = (dataEvent.jet1_+dataEvent.jet2_).M(); double qqM = (dataEvent.jet3_+dataEvent.jet4_).M(); double qqDeltaEta = TMath::Abs(dataEvent.jet3_.Eta()-dataEvent.jet4_.Eta());
    if(dataEvent.njets_ >= 4){
      mHZZ = (dataEvent.jet1_+dataEvent.jet2_+dataEvent.lep1_+dataEvent.lep2_).M();
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((dataEvent.jet1_+dataEvent.jet3_).M()-90.1876)) {Mjj = (dataEvent.jet1_+dataEvent.jet3_).M();qqM = (dataEvent.jet2_+dataEvent.jet4_).M();qqDeltaEta = TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet4_.Eta());mHZZ = (dataEvent.jet1_+dataEvent.jet3_+dataEvent.lep1_+dataEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((dataEvent.jet1_+dataEvent.jet4_).M()-90.1876)) {Mjj = (dataEvent.jet1_+dataEvent.jet4_).M();qqM = (dataEvent.jet2_+dataEvent.jet3_).M();qqDeltaEta = TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta());mHZZ = (dataEvent.jet1_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((dataEvent.jet2_+dataEvent.jet3_).M()-90.1876)) {Mjj = (dataEvent.jet2_+dataEvent.jet3_).M();qqM = (dataEvent.jet1_+dataEvent.jet4_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet4_.Eta());mHZZ = (dataEvent.jet2_+dataEvent.jet3_+dataEvent.lep1_+dataEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((dataEvent.jet2_+dataEvent.jet4_).M()-90.1876)) {Mjj = (dataEvent.jet2_+dataEvent.jet4_).M();qqM = (dataEvent.jet1_+dataEvent.jet3_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta());mHZZ = (dataEvent.jet2_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
      //if(TMath::Abs(Mjj-90.1876) > TMath::Abs((dataEvent.jet3_+dataEvent.jet4_).M()-90.1876)) {Mjj = (dataEvent.jet3_+dataEvent.jet4_).M();qqM = (dataEvent.jet1_+dataEvent.jet2_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());mHZZ = (dataEvent.jet3_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet4_.Eta())) {Mjj = (dataEvent.jet1_+dataEvent.jet3_).M();qqM = (dataEvent.jet2_+dataEvent.jet4_).M();qqDeltaEta = TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet4_.Eta());mHZZ = (dataEvent.jet1_+dataEvent.jet3_+dataEvent.lep1_+dataEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta())) {Mjj = (dataEvent.jet1_+dataEvent.jet4_).M();qqM = (dataEvent.jet2_+dataEvent.jet3_).M();qqDeltaEta = TMath::Abs(dataEvent.jet2_.Eta()-dataEvent.jet3_.Eta());mHZZ = (dataEvent.jet1_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet4_.Eta())) {Mjj = (dataEvent.jet2_+dataEvent.jet3_).M();qqM = (dataEvent.jet1_+dataEvent.jet4_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet4_.Eta());mHZZ = (dataEvent.jet2_+dataEvent.jet3_+dataEvent.lep1_+dataEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta())) {Mjj = (dataEvent.jet2_+dataEvent.jet4_).M();qqM = (dataEvent.jet1_+dataEvent.jet3_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet3_.Eta());mHZZ = (dataEvent.jet2_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
      if(qqDeltaEta < TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta())) {Mjj = (dataEvent.jet3_+dataEvent.jet4_).M();qqM = (dataEvent.jet1_+dataEvent.jet2_).M();qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());mHZZ = (dataEvent.jet3_+dataEvent.jet4_+dataEvent.lep1_+dataEvent.lep2_).M();}
    } else {qqDeltaEta = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());}
    bool passCuts = false;
    if(channel == 29 || channel == 30 || channel == 31){ // WW selection
      if(
	  dataEvent.dilep_.M()   > 12 &&
         (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          Njet3 == nJetsType &&
          dataEvent.lep1_.Pt() > 20. &&
          dataEvent.lep2_.Pt() > ptLepMin &&
          passMET == true &&
	  passNewCuts == true &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         //(dataEvent.cuts_ & patternTopTagNotInJets) == patternTopTagNotInJets &&
	 //(dataEvent.jetLowBtag_ < 2.1 && dataEvent.nSoftMuons_ == 0) &&
	 //dataEvent.jet1Btag_ < 2.1 &&
	 //dataEvent.jet2Btag_ >= 2.1 &&
	 dPhiDiLepJetCut == true &&
	(channel != 30 || (dataEvent.mt_ > 30 && (dataEvent.mt_ > 60 || nJetsType == 2) && dataEvent.mt_ < 280 && dataEvent.dilep_.M() < 200 && dataEvent.lep1_.Pt() > 20. && VBFCuts == true)) &&
	(channel != 31 || (dataEvent.mt_ > 30 && (dataEvent.mt_ > 80 || nJetsType == 2) && dataEvent.mt_ < 600 && dataEvent.dilep_.M() < 600 && dataEvent.lep1_.Pt() > 50. && VBFCuts == true)) &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	//cout << "DATA "  << dataEvent.njets_ << " " << dataEvent.type_ << " " << dataEvent.run_ << " " << dataEvent.event_ << " " << dataEvent.lumi_ << endl;
      }
    } // WW selection
    if(channel == 6){ // Z selection
      //bool passBtagVeto = dataEvent.nSoftMuons_ == 0 && (dataEvent.jet1_.Pt() <= 20 || dataEvent.jet1ProbBtag_ < 0.244) && (dataEvent.jet2_.Pt() <= 20 || dataEvent.jet2ProbBtag_ < 0.244) && (dataEvent.jet3_.Pt() <= 20 || dataEvent.jet3ProbBtag_ < 0.244) && (dataEvent.jet4_.Pt() <= 20 || dataEvent.jet4ProbBtag_ < 0.244);
      bool passBtagVeto = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
      if(
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
        passBtagVeto &&
       //((charge == 0 && VarAux2 == 0) || (charge != 0 && VarAux2 == 1)) &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20. &&
	 int((dataEvent.jet1McId_%1000-dataEvent.jet1McId_%100)/100) == 0 &&
         Njet3 == nJetsType &&
       //((abs(dataEvent.lep1_.Eta())>=VarAux0&&abs(dataEvent.lep1_.Eta())<VarAux0+0.5&&abs(dataEvent.lep2_.Eta())>=VarAux1&&abs(dataEvent.lep2_.Eta())<VarAux1+0.5)||
       // (abs(dataEvent.lep1_.Eta())>=VarAux1&&abs(dataEvent.lep1_.Eta())<VarAux1+0.5&&abs(dataEvent.lep2_.Eta())>=VarAux0&&abs(dataEvent.lep2_.Eta())<VarAux0+0.5)) &&
        (fabs(dataEvent.dilep_.M()-91.1876) < 15.) && 
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	//cout << "DATA "  << dataEvent.njets_ << " " << dataEvent.type_ << " " << dataEvent.run_ << " " << dataEvent.event_ << " " << dataEvent.lumi_ << endl;
      }
    } // Z selection

    if(channel == 5){ // bbmm selection
      int nb = 0;
      if(dataEvent.jet1_.Pt() > 30 && dataEvent.jet1ProbBtag_ > 0.783) nb++;
      if(dataEvent.jet2_.Pt() > 30 && dataEvent.jet2ProbBtag_ > 0.783) nb++;
      if(dataEvent.jet3_.Pt() > 30 && dataEvent.jet3ProbBtag_ > 0.783) nb++;
      if(dataEvent.jet4_.Pt() > 30 && dataEvent.jet4ProbBtag_ > 0.783) nb++;
      if(
	  dataEvent.dilep_.M() > 12 && dataEvent.dilep_.M() < 70 && 
         (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          dataEvent.njets_ < 2 &&
          dataEvent.lep1_.Pt() > 25. &&
          dataEvent.lep2_.Pt() > 25. &&
	  DeltaR(dataEvent.lep1_.Phi(),dataEvent.lep1_.Eta(),dataEvent.lep2_.Phi(),dataEvent.lep2_.Eta()) < 1.1 &&
          nb >= 1 &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	//
	passCuts = true;
      }
    } // bbmm selection

    if(channel==15||channel==16 ||channel==17||channel==18){ // HHbbww selection
      bool passBtag[3] = {false,false,false};
      if     ( dataEvent.jet1ProbBtag_ > 0.679 && dataEvent.jet2ProbBtag_ > 0.679)  passBtag[0] = true;
      else if((dataEvent.jet1ProbBtag_ > 0.679 && dataEvent.jet2ProbBtag_ > 0.244)||
              (dataEvent.jet1ProbBtag_ > 0.244 && dataEvent.jet2ProbBtag_ > 0.679)) passBtag[1] = true;
      else if( dataEvent.jet1ProbBtag_ > 0.679 || dataEvent.jet2ProbBtag_ > 0.679)  passBtag[2] = true;
      if(
         (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
          charge == 0 &&
          dataEvent.lep1_.Pt() > 20. &&
          dataEvent.lep2_.Pt() > 20. &&
          usedMet > 20. &&
	  dataEvent.njets_ >= 2 && dataEvent.njets_ <= 3 &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
        ((passBtag[0]&&channel==15)||(passBtag[1]&&channel==16)||(passBtag[2]&&channel==17)||
	((passBtag[0]||passBtag[1]||passBtag[2])&&channel==18)) &&
	  dataEvent.dilep_.M() > 12 && 
	  dataEvent.dilep_.M() < 200 && dataEvent.mt_ < 200 && (dataEvent.jet1_+dataEvent.jet2_).M() < 200.0 &&
	  dataEvent.dilep_.M() <  50 && dataEvent.mt_ > 40 && dataEvent.mt_ < 125 && fabs((dataEvent.jet1_+dataEvent.jet2_).M()-125.0) < 20.0 &&
         (dataEvent.jet1_+dataEvent.jet2_).Pt() > 100 && theWW.Pt() > 100 &&
          dPhiHH > 160.0 && dRMinJL > 1.8 &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	//
	passCuts = true;
      }
    } // HHbbww selection

    if(channel == 9){ // Ztautau selection
      if(
         dataEvent.dilep_.M() > 12 && dataEvent.dilep_.M() < 100 &&
         Njet3 == nJetsType &&
         (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
	 dataEvent.dilep_.Pt() < 30. &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Ztautau selection

    if(channel == 902){ // VH,qqlnln selection
      if(
         dataEvent.dilep_.M()   < 70  &&
         dataEvent.dilep_.M()   > 12  &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 dataEvent.njets_ == 2 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
        (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
	 dPhiDiLepJetCut == true &&
        (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dataEvent.jet1Btag_ < 1.6 && dataEvent.jet2Btag_ < 1.6 && 
        (dataEvent.jet1_+dataEvent.jet2_).M() > 60. && (dataEvent.jet1_+dataEvent.jet2_).M() < 110. &&
         TMath::Abs(dataEvent.jet1_.Eta()) < 2.5 &&
         TMath::Abs(dataEvent.jet2_.Eta()) < 2.5 &&
         TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) < 2.1 &&
         dataEvent.dR_   < 1.3 &&
         dataEvent.mt_	> 50 &&
         dataEvent.mt_	< 160 &&
	 DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180/TMath::Pi() > 75.0 &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // VH,qqlnln selection

    if((channel >= 100 && channel <= 800) && channel != 128){ // H->WW selection
      double theCutPtMinLow = cutPtMinLow (mH, dataEvent.type_);
      if(
         dataEvent.mt_  > theCutMTLow &&
         dataEvent.mt_  < theCutMTHigh &&
         dataEvent.dilep_.M() > theCutMassLow &&
         dataEvent.dilep_.M() < theCutMassHigh &&
         dataEvent.lep1_.Pt() > theCutPtMaxLow &&
         dataEvent.lep2_.Pt() > theCutPtMinLow &&
         dataEvent.dPhi_*180.0/TMath::Pi() < theCutDeltaphilHigh &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 Njet3 == nJetsType &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         passMET == true &&
	 passNewCuts == true &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 15. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dPhiDiLepJetCut == true &&
	 VBFCuts == true &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
         1 == 1
	){
	passCuts = true;
	//printf("%6d %12d %4d - %d %5.1f %5.1f %5.1f %5.1f - %5.1f %5.1f %5.1f - %5.1f %5.1f %5.1f - %5.1f %5.1f - %5.1f %5.1f %5.1f %5.1f\n",
	//dataEvent.run_,dataEvent.event_,dataEvent.lumi_,
	//dataEvent.type_,dataEvent.met_,dataEvent.pmet_,dataEvent.trackMet_,dataEvent.pTrackMet_,
	//dataEvent.lep1_.Pt(),dataEvent.lep2_.Pt(),dataEvent.jet1_.Pt(),
	//dataEvent.lep1_.Phi(),dataEvent.lep2_.Phi(),dataEvent.jet1_.Phi(),
	//dataEvent.metPhi_,dataEvent.trackMetPhi_,usedMet,dataEvent.dilep_.M(),dataEvent.dPhi_*180.0/TMath::Pi(),dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi());
      }
    } // H->WW selection

    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
         dataEvent.dilep_.M() > 12 && dataEvent.dilep_.M() < 200 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 dataEvent.njets_ <= 2. && dataEvent.njets_ >= 0. &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         (dataEvent.jet1_.Pt() <= 30 || TMath::Abs(dataEvent.jet1_.Eta()) < 2.5) && (dataEvent.jet2_.Pt() <= 30 || TMath::Abs(dataEvent.jet2_.Eta()) < 2.5) && (dataEvent.jet3_.Pt() <= 30 || TMath::Abs(dataEvent.jet3_.Eta()) < 2.5) &&
         passMET == true &&
         dataEvent.mt_ > 0. && dataEvent.mt_ < 200. && 
         (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
         (dPhiDiLepJetCut == true || dataEvent.type_ != SmurfTree::ee) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->2l selection

    if(channel == 2000){ // H->ZZ->2l2q selection
      if(
         dataEvent.dilep_.M()   > 12 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge == 0 &&
	 qqM > 200 &&
	 qqDeltaEta > 2.5 &&
	 Mjj > 50 && Mjj < 130 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
	 dataEvent.njets_ >= 4 &&
         abs(dataEvent.jet1_.Eta()) < 5.0 && abs(dataEvent.jet2_.Eta()) < 5.0 && abs(dataEvent.jet3_.Eta()) < 5.0 && abs(dataEvent.jet4_.Eta()) < 5.0 &&
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
      }
    } // H->ZZ->2l2q selection

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
      else if(thePlot ==11) myVar = qqM;
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = dataEvent.met_/dataEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = dataEvent.njets_;
      else if(thePlot ==18) myVar = dataEvent.nvtx_;
      else if(thePlot ==19) myVar = dataEvent.pTrackMet_;
      else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(dataEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(dataEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(dataEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(dataEvent.lid1_)==13) myVar = dataEvent.lep1_.Eta();else if(fabs(dataEvent.lid2_)==13) myVar = dataEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(TMath::Min((double)ewkMVA,0.999),-0.999);
      else if(thePlot ==30) myVar = TMath::Max(dataEvent.jet1Btag_,dataEvent.jet2Btag_);
      else if(thePlot ==31) myVar = 2*HWWKin.CalcMR();
      else if(thePlot ==32) myVar = 2*HWWKin.CalcMRNEW();
      else if(thePlot ==33) myVar = HWWKin.CalcDeltaPhiRFRAME()*180.0/TMath::Pi();
      else if(thePlot ==34) myVar = Mjj;
      else if(thePlot ==35) myVar = qqDeltaEta;
      else if(thePlot ==36) myVar = Njet3;
      else if(thePlot ==37) myVar = TMath::Min((dataEvent.jet1_+dataEvent.jet2_).M(),1999.999);
      else if(thePlot ==38) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
      else if(thePlot ==39) myVar = TMath::Min(TMath::Max((1.3491*log((dataEvent.jet1_+dataEvent.jet2_).M())-0.01163*(dataEvent.jet1_.Eta()*dataEvent.jet2_.Eta())+1.433*TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta())-7.647)/15.00,0.001),0.999);
      else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = (dataEvent.jet1_+dataEvent.jet2_).Pt();
      else if(thePlot ==45) myVar = mtHZZ;
      else if(thePlot ==46) myVar = mHZZ;
      else if(thePlot ==47) myVar = theWW.Pt();
      else if(thePlot ==48) myVar = dataEvent.type_;
      else if(thePlot ==49) myVar = ptFracHH;
      else if(thePlot ==50) myVar = dPhiHH;
      else if(thePlot ==51) myVar = theHH.Pt();
      else if(thePlot ==52) myVar = dataEvent.trackMet_*sin(dataEvent.trackMetPhi_);
      else if(thePlot ==53) myVar = DeltaPhi(dataEvent.jet3_.Phi(),dataEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==54) myVar = dPhiDilepRecoil*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = dataEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==56) myVar = poorManMetSyst(dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                                   dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,dataEvent.met_,dataEvent.metPhi_,
                                    		   dataEvent.jet1_,dataEvent.jet2_,dataEvent.jet3_,dataEvent.jet4_,0);
      else if(thePlot ==57) myVar = dataEvent.dR_;
      else if(thePlot ==58) myVar = Lester::higgsWWStarMassBoundNoROOT(dataEvent.lep1_.e(),dataEvent.lep1_.px(),dataEvent.lep1_.py(),dataEvent.lep1_.pz(),
                                                                       dataEvent.lep2_.e(),dataEvent.lep2_.px(),dataEvent.lep2_.py(),dataEvent.lep2_.pz(),
								       dataEvent.met_*cos(dataEvent.metPhi_),dataEvent.met_*sin(dataEvent.metPhi_),80.41);
      else if(thePlot ==59) myVar = mt_atlas(dataEvent.dilep_,dataEvent.met_,dataEvent.metPhi_);
      else if(thePlot ==60) myVar = dRMinJL;
      else if(thePlot ==61) myVar = TMath::Min(DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep1_.Phi()),TMath::Min(DeltaPhi(dataEvent.jet1_.Phi(),dataEvent.lep2_.Phi()),TMath::Min(DeltaPhi(dataEvent.jet2_.Phi(),dataEvent.lep1_.Phi()),DeltaPhi(dataEvent.jet2_.Phi(),dataEvent.lep2_.Phi()))))*180.0/TMath::Pi();

      double wei = 1.0;
      nSelectedData = nSelectedData + wei;
      histo5->Fill(myVar,1.0*wei);
    }
  } // End loop data
  printf("data: %f\n",nSelectedData);
  //printf("%d %d %d %d\n",(int)histo5->GetBinContent(1),(int)histo5->GetBinContent(2),(int)histo5->GetBinContent(3),(int)histo5->GetBinContent(4));

  //histos->Draw();
  if(makeGoodPlots){
    char output[200];
    sprintf(output,"histo_nice8.root");     
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
      histos->Write();
      histo0->Write();
      histo1->Write();
      histo2->Write();
      histo3->Write();
      histo4->Write();
      histo5->Write();
      histoggH->Write();
      histoqqH->Write();
      histoVH->Write();
    outFilePlotsNote->Close();
  }

  if((channel > 1000 && channel < 2000) ||
     (channel >= 100 && channel <= 800)){
    for(int n1=0; n1<nBin; n1++){
      double signif0 = 0.0;
      if(S0[n1] != 0 && B0[n1] != 0) signif0 = scpFast(S0[n1],B0[n1],0.35*B0[n1]);
      double signif1 = 0.0;
      if(S1[n1] != 0 && B1[n1] != 0) signif1 = scpFast(S1[n1],B1[n1],0.35*B1[n1]);
      hDSignif[0]->SetBinContent(n1, signif0);
      hDSignif[1]->SetBinContent(n1, signif1);
    }
    printf("S: %f  B: %f  Signif: %f  S/B: %f\n",S0[0],B0[0],
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
    printf("S: %f  B: %f  Signif: %f  S/B: %f\n",S0[0],B0[0],
    	   S0[0]/sqrt(S0[0]+B0[0]+0.35*0.35*B0[0]*B0[0]),S0[0]/B0[0]);
  }
  
  bgdDecay[29] = bgdDecay[29] + bgdDecay[14];
  weiDecay[29] = weiDecay[29] + weiDecay[14];
  bgdDecay[14] = 0.0;
  weiDecay[14] = 0.0;
  // W/Z/tX/Vgamma
  double bgdCombined[4] = {bgdDecay[0]+bgdDecay[1]+bgdDecay[2]+bgdDecay[3]+bgdDecay[23],
                           bgdDecay[6]+bgdDecay[7]+bgdDecay[8]+bgdDecay[9]+bgdDecay[10]+bgdDecay[31],
			   bgdDecay[4]+bgdDecay[5]+bgdDecay[11]+bgdDecay[12]+bgdDecay[13],
			   bgdDecay[17]+bgdDecay[18]+bgdDecay[19]+bgdDecay[20]};
  double bgdCombinedE[4] = {sqrt(weiDecay[0]+weiDecay[1]+weiDecay[2]+weiDecay[3]+weiDecay[23]),
                            sqrt(weiDecay[6]+weiDecay[7]+weiDecay[8]+weiDecay[9]+weiDecay[10]+weiDecay[31]),
                            sqrt(bgdDecay[4]*weiDecay[4]+weiDecay[5]+weiDecay[11]+weiDecay[12]+weiDecay[13]),
                            sqrt(weiDecay[17]+weiDecay[18]+weiDecay[19]+weiDecay[20])};
  for(int i=0; i<45; i++){
    if(bgdDecay[i] != 0 ||
       (i==27||i==28||i==29||i==30)) printf("bdg(%2d) = %8.3f +/- %6.3f\n",i,bgdDecay[i],sqrt(weiDecay[i]));
  }
  for(int i=0; i<16; i++){
    if(bgdDecayFake[i] != 0) printf("bgdDecayFake(%2d) = %4.1f +/- %4.1f\n",i,bgdDecayFake[i],sqrt(bgdDecayFakeE[i]));
  }
  printf("bdg(xW) = %8.3f +/- %6.3f\n",bgdCombined[0],bgdCombinedE[0]);
  printf("bdg(xZ) = %8.3f +/- %6.3f\n",bgdCombined[1],bgdCombinedE[1]);
  printf("bdg(tX) = %8.3f +/- %6.3f\n",bgdCombined[2],bgdCombinedE[2]);
  printf("bdg(Vg) = %8.3f +/- %6.3f\n",bgdCombined[3],bgdCombinedE[3]);
  // ggWW, WW, WZ, ZZ, W, Z, tX, Vg, Ztt
  double Syst[9] = {0.00, 0.00, 0.18, 0.18, 0.50, 1.00, 1.00, 0.35, 0.10};
  //double Syst[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if(channel >=100 && channel <= 200) {Syst[0] = 0.50; Syst[1] = 0.55;}
  if(channel > 200 && channel <= 800) {Syst[0] = 0.50; Syst[1] = 0.50;}
  double total_stat = sqrt(weiDecay[21] +
                           weiDecay[30] +
                           weiDecay[29] +
                           weiDecay[28] +
                           weiDecay[27] +
                           bgdCombinedE[0]*bgdCombinedE[0]+
		           bgdCombinedE[1]*bgdCombinedE[1]+
		           bgdCombinedE[2]*bgdCombinedE[2]+
		           bgdCombinedE[3]*bgdCombinedE[3]+
			   weiDecay[10]);
  double totalS = sqrt(weiDecay[21]+bgdDecay[21]*Syst[0]*bgdDecay[21]*Syst[0] +
                       weiDecay[30]+bgdDecay[30]*Syst[0]*bgdDecay[30]*Syst[0] +
                       weiDecay[29]+bgdDecay[29]*Syst[1]*bgdDecay[29]*Syst[1] +
                       weiDecay[28]+bgdDecay[28]*Syst[2]*bgdDecay[28]*Syst[2] +
                       weiDecay[27]+bgdDecay[27]*Syst[3]*bgdDecay[27]*Syst[3] +
                       bgdCombinedE[0]*bgdCombinedE[0]+bgdCombined[0]*Syst[4]*bgdCombined[0]*Syst[4]+
		       bgdCombinedE[1]*bgdCombinedE[1]+bgdCombined[1]*Syst[5]*bgdCombined[1]*Syst[5]+
		       bgdCombinedE[2]*bgdCombinedE[2]+bgdCombined[2]*Syst[6]*bgdCombined[2]*Syst[6]+
		       bgdCombinedE[3]*bgdCombinedE[3]+bgdCombined[3]*Syst[7]*bgdCombined[3]*Syst[7]+
		       weiDecay[10]+bgdDecay[10]);
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
    printf("histo -> s: %8.2f d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f\n",histos->GetSumOfWeights(),histo5->GetSumOfWeights(),histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights());
  
    for(int i=0; i<6; i++) if(nSigCut[i] > 0) nSigECut[i] = sqrt(nSigECut[i]);
    for(int i=0; i<6; i++) printf("nSig(%1d) = %6.3f +/- %6.3f\n",i,nSigCut[i],nSigECut[i]);
    if(WWXSSel == false){
      double totalB   =      bgdDecay[21]+bgdDecay[27]+bgdDecay[28]+bgdDecay[29]+bgdDecay[30]+                 bgdCombined[0]+                 bgdCombined[1]+                 bgdCombined[2]+                 bgdCombined[3];
      double totalB_E = sqrt(weiDecay[21]+weiDecay[27]+weiDecay[28]+weiDecay[29]+weiDecay[30]+bgdCombinedE[0]*bgdCombinedE[0]+bgdCombinedE[1]*bgdCombinedE[1]+bgdCombinedE[2]*bgdCombinedE[2]+bgdCombinedE[3]*bgdCombinedE[3]);
      printf("Higgs     &    data   & SumBkg      &   qqWW   &   ggWW   & Top   & Wjets  & VV & VVV &  Zjets & Vgamma & Vgamma* \\\\\n"); 
      printf("%6.3f \\pm %5.3f & %3d & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f  \\pm %5.2f & %6.2f  \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f \\\\\n",
             nSigCut[0],nSigECut[0],(int)nSelectedData,totalB,totalB_E,bgdDecay[29],sqrt(weiDecay[29]),bgdDecay[30],sqrt(weiDecay[30]),bgdCombined[2],bgdCombinedE[2],bgdCombined[0],bgdCombinedE[0],
	     bgdDecay[27]+bgdDecay[28],sqrt(weiDecay[27]+weiDecay[28]),bgdDecay[21],sqrt(weiDecay[21]),bgdCombined[1],bgdCombinedE[1],bgdDecay[19],sqrt(weiDecay[19]),bgdDecay[20],sqrt(weiDecay[20]));
    } else {
      double totalB   =      bgdDecay[21]+bgdDecay[27]+bgdDecay[28]+ nSigCut[0]+                 bgdCombined[0]+                 bgdCombined[1]+                 bgdCombined[2]+                 bgdCombined[3];
      double totalB_E = sqrt(weiDecay[21]+weiDecay[27]+weiDecay[28]+nSigECut[0]+bgdCombinedE[0]*bgdCombinedE[0]+bgdCombinedE[1]*bgdCombinedE[1]+bgdCombinedE[2]*bgdCombinedE[2]+bgdCombinedE[3]*bgdCombinedE[3]);
      printf("data  &   qqWW   &   ggWW & SumBkg     &   Higgs      & Top   & Wjets  & VV  & VVV & Zjets & Vgamma & Vgamma* \\\\\n"); 
      printf("%3d & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f  \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f  & %6.2f \\pm %5.2f \\\\\n",
             (int)nSelectedData,bgdDecay[29],sqrt(weiDecay[29]),bgdDecay[30],sqrt(weiDecay[30]),totalB,totalB_E,nSigCut[0],nSigECut[0],
	     bgdCombined[2],bgdCombinedE[2],bgdCombined[0],bgdCombinedE[0],
	     bgdDecay[27]+bgdDecay[28],sqrt(weiDecay[27]+weiDecay[28]),bgdDecay[21],sqrt(weiDecay[21]),bgdCombined[1],bgdCombinedE[1],
	     bgdDecay[19],sqrt(weiDecay[19]),bgdDecay[20],sqrt(weiDecay[20]));
    }
  }

  return;

}

double caloMet(double met, double metPhi, LorentzVector l1,LorentzVector l2){
double metx = met * cos(metPhi) - l1.Px() - l2.Px();
double mety = met * sin(metPhi) - l1.Py() - l2.Py();

return sqrt(metx*metx+mety*mety);
}
