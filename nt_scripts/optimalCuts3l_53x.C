#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
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
const bool UseDyttDataDriven = false; // if true, then remove em events in dyll MC
SmurfTree ttEvent;
SmurfTree systEvent;

//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCuts3l_53x
(
 int     mH  	 = 1100,
 int thePlot = 0,
 TString bgdInputFile     = "ntuples_53x/backgroundA_3l.root",
 TString signalInputFile  = "ntuples_53x/hww125.root",
 TString dataInputFile    = "ntuples_53x/data_3l.root",
 TString signalInputFile2 = "ntuples_53x/vhtt125.root",
 TString systInputFile    = "ntuples_53x/hww_syst_3l.root",
 bool fillInfoNote = true,
 double mhAna = 125,
 int period = 13
 )
{
  bool applyFRCorr   = true;
  bool makeGoodPlots = true;
  double lumi = 1.545;
  bool isFermioPhobic = false;
  bool isSM4          = false;

  int nsel = (int)period/10;
  period   = (int)period%10;
  printf("nsel: %d, period: %d\n",nsel,period);

  const bool useTemplates = true; // otherwise, if false, will crash
  bool useAlternativeStatTemplates = false;
  int rebinMVAHist = 1;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  if(signalInputFile2 != ""){
    ttEvent.LoadTree(signalInputFile2,-1);
    ttEvent.InitTree(0);
  }

  if(systInputFile != "" && useTemplates == true){
    systEvent.LoadTree(systInputFile,-1);
    systEvent.InitTree(0);
  }

  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  if	 (period == 0){ // Full2012-Summer12-V9-3500ipb
  }
  else if(period == -1){ // Full2012-Summer12-V9 NoJetId
  }
  else if(period == 1){ // Full2012-Summer12-V9-5000ipb
  }
  else if(period == 2){ // Full2012-Summer12-V9-12100ipb
    effPath  = "/data/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/summary.root";
    fakePath = "/data/smurf/dlevans/FakeRates/V00-02-07_HCP_V0/summary.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_12p1ifb.root";
    lumi     = 12.1;
  }
  else if(period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

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

  //Fake rate systematics
  TFile *fLeptonFRFileSyst = TFile::Open(fakePath.Data());
  TH2D *fhDFRMuSyst = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  TH2D *fhDFRElSyst = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold20_PtEta"));
  assert(fhDFRMuSyst);
  assert(fhDFRElSyst);
  fhDFRMuSyst->SetDirectory(0);
  fhDFRElSyst->SetDirectory(0);
  fLeptonFRFileSyst->Close();
  delete fLeptonFRFileSyst;

  TFile *fPUS4File = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUS4File;

  int channel = mH;
  if(channel == 1200 || channel == 1099 || channel == 1301) applyFRCorr = false;
  if(channel == 1102) applyFRCorr = false;
  int newMH = (int)mhAna;
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

  int nBin    = 100;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  double maxHis = 1.0; double minHis = -1.0;
  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 5; xminPlot = -1.0; xmaxPlot = 1.0;}
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 60; xminPlot = 60.0; xmaxPlot = 120.0;}
  else if(thePlot >= 20 && thePlot <= 22) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 23 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 30 && thePlot <= 31) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot =  400.0;}
  else if(thePlot >= 32 && thePlot <= 33) {nBinPlot = 120; xminPlot = 0; xmaxPlot = 12.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  1000.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 100; xminPlot = 50; xmaxPlot =  800;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  100.0;}
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  200.0;}
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 40; xminPlot = 0.0; xmaxPlot =  4.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 45) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 8; xminPlot = 0.5; xmaxPlot = 8.5;}
  else if(thePlot >= 49 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  else if(thePlot >= 53 && thePlot <= 53) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 54 && thePlot <= 54) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 61 && thePlot <= 61) {nBinPlot = 20;  xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 62 && thePlot <= 62) {nBinPlot = 20;  xminPlot = 0.0; xmaxPlot = 800.0;}
  nBin = nBinPlot;

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
  TH1D* histo6 = (TH1D*) histos->Clone("histo6");
  TH1D* histoWW= (TH1D*) histos->Clone("histoWW");
  TH1D* histoTT= (TH1D*) histos->Clone("histoTT");
  TH1D *histo_WZ     = (TH1D*) histos->Clone("histo_WZ");
  TH1D *histo_ZZ     = (TH1D*) histos->Clone("histo_ZZ");
  TH1D *histo_Wjets  = (TH1D*) histos->Clone("histo_Wjets");
  TH1D *histo_Wgamma = (TH1D*) histos->Clone("histo_Wgamma");
  TH1D *histo_VVV    = (TH1D*) histos->Clone("histo_VVV");
  TH1D *histo_WH_htt_SM  = (TH1D*) histos->Clone("histo_WH_htt_SM");
  TH1D *histo_WH_hww_SM  = (TH1D*) histos->Clone("histo_WH_hww_SM");
  histos ->Scale(0.0);
  histo0 ->Scale(0.0);
  histo1 ->Scale(0.0);
  histo2 ->Scale(0.0);
  histo3 ->Scale(0.0);
  histo4 ->Scale(0.0);
  histo5 ->Scale(0.0);
  histo_WZ    ->Scale(0.0);
  histo_ZZ    ->Scale(0.0);
  histo_Wjets ->Scale(0.0);
  histo_Wgamma->Scale(0.0);
  histo_VVV   ->Scale(0.0);
  histo_WH_htt_SM->Scale(0.0);
  histo_WH_hww_SM->Scale(0.0);
  histoWW->Scale(0.0);
  histoTT->Scale(0.0);

  TH1D* histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVUp	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAWH_httStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVDown = (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAWH_httStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVUp	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVDown = (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_WZ_CMS_MVAWZStatBounding_8TeVUp	        = (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAWZStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_WZ_CMS_MVAWZStatBounding_8TeVDown	        = (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAWZStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_ZZ_CMS_MVAZZStatBounding_8TeVUp	        = (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAZZStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_ZZ_CMS_MVAZZStatBounding_8TeVDown	        = (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAZZStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp     = (TH1D*) histo0->Clone(Form("histo_Wjets_CMS_MVAWjetsStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown   = (TH1D*) histo0->Clone(Form("histo_Wjets_CMS_MVAWjetsStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAWgammaStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown = (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAWgammaStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_VVV_CMS_MVAVVVStatBounding_8TeVUp         = (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAVVVStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_VVV_CMS_MVAVVVStatBounding_8TeVDown       = (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAVVVStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVUp   = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVDown = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeV_vh3l%dDown",nsel));
  TH1D* histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVUp   = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeV_vh3l%dUp",nsel));
  TH1D* histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVDown = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeV_vh3l%dDown",nsel));

  TH1D* histo_WH_htt_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_WH_htt_CMS_MVALepEffBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVALepEffBoundingDown"));
  TH1D* histo_WH_hww_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_WH_hww_CMS_MVALepEffBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVALepEffBoundingDown"));
  TH1D* histo_WZ_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_WZ_CMS_MVALepEffBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVALepEffBoundingDown"));
  TH1D* histo_ZZ_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_ZZ_CMS_MVALepEffBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVALepEffBoundingDown"));
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVALepEffBoundingDown"));
  TH1D* histo_VVV_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_VVV_CMS_MVALepEffBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVALepEffBoundingDown"));
  TH1D* histo_WH_htt_SM_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_WH_htt_SM_CMS_MVALepEffBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVALepEffBoundingDown"));
  TH1D* histo_WH_hww_SM_CMS_MVALepEffBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVALepEffBoundingUp"));  
  TH1D* histo_WH_hww_SM_CMS_MVALepEffBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVALepEffBoundingDown"));

  TH1D* histo_WH_htt_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVALepResBoundingUp"));  
  TH1D* histo_WH_htt_CMS_MVALepResBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVALepResBoundingDown"));
  TH1D* histo_WH_hww_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVALepResBoundingUp"));  
  TH1D* histo_WH_hww_CMS_MVALepResBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVALepResBoundingDown"));
  TH1D* histo_WZ_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVALepResBoundingUp"));  
  TH1D* histo_WZ_CMS_MVALepResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVALepResBoundingDown"));
  TH1D* histo_ZZ_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVALepResBoundingUp"));  
  TH1D* histo_ZZ_CMS_MVALepResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVALepResBoundingDown"));
  TH1D* histo_Wgamma_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVALepResBoundingUp"));  
  TH1D* histo_Wgamma_CMS_MVALepResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVALepResBoundingDown"));
  TH1D* histo_VVV_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVALepResBoundingUp"));  
  TH1D* histo_VVV_CMS_MVALepResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVALepResBoundingDown"));
  TH1D* histo_WH_htt_SM_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVALepResBoundingUp"));  
  TH1D* histo_WH_htt_SM_CMS_MVALepResBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVALepResBoundingDown"));
  TH1D* histo_WH_hww_SM_CMS_MVALepResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVALepResBoundingUp"));  
  TH1D* histo_WH_hww_SM_CMS_MVALepResBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVALepResBoundingDown"));

  TH1D* histo_WH_htt_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_WH_htt_CMS_MVAMETResBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAMETResBoundingDown"));
  TH1D* histo_WH_hww_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_WH_hww_CMS_MVAMETResBoundingDown       	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAMETResBoundingDown"));
  TH1D* histo_WZ_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_WZ_CMS_MVAMETResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAMETResBoundingDown"));
  TH1D* histo_ZZ_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_ZZ_CMS_MVAMETResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAMETResBoundingDown"));
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAMETResBoundingDown"));
  TH1D* histo_VVV_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_VVV_CMS_MVAMETResBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAMETResBoundingDown"));
  TH1D* histo_WH_htt_SM_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_WH_htt_SM_CMS_MVAMETResBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAMETResBoundingDown"));
  TH1D* histo_WH_hww_SM_CMS_MVAMETResBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAMETResBoundingUp"));  
  TH1D* histo_WH_hww_SM_CMS_MVAMETResBoundingDown       = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAMETResBoundingDown"));

  TH1D* histo_WH_htt_CMS_MVAJESBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAJESBoundingUp"));  
  TH1D* histo_WH_htt_CMS_MVAJESBoundingDown          	= (TH1D*) histo0->Clone(Form("histo_WH_htt_CMS_MVAJESBoundingDown"));
  TH1D* histo_WH_hww_CMS_MVAJESBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAJESBoundingUp"));  
  TH1D* histo_WH_hww_CMS_MVAJESBoundingDown          	= (TH1D*) histo0->Clone(Form("histo_WH_hww_CMS_MVAJESBoundingDown"));
  TH1D* histo_WZ_CMS_MVAJESBoundingUp	             	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAJESBoundingUp"));  
  TH1D* histo_WZ_CMS_MVAJESBoundingDown	             	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAJESBoundingDown"));
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp	             	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAJESBoundingUp"));  
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown	             	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAJESBoundingDown"));
  TH1D* histo_Wgamma_CMS_MVAJESBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAJESBoundingUp"));  
  TH1D* histo_Wgamma_CMS_MVAJESBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_Wgamma_CMS_MVAJESBoundingDown"));
  TH1D* histo_VVV_CMS_MVAJESBoundingUp	             	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAJESBoundingUp"));  
  TH1D* histo_VVV_CMS_MVAJESBoundingDown	     	= (TH1D*) histo0->Clone(Form("histo_VVV_CMS_MVAJESBoundingDown"));
  TH1D* histo_WH_htt_SM_CMS_MVAJESBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAJESBoundingUp"));  
  TH1D* histo_WH_htt_SM_CMS_MVAJESBoundingDown          = (TH1D*) histo0->Clone(Form("histo_WH_htt_SM_CMS_MVAJESBoundingDown"));
  TH1D* histo_WH_hww_SM_CMS_MVAJESBoundingUp	     	= (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAJESBoundingUp"));  
  TH1D* histo_WH_hww_SM_CMS_MVAJESBoundingDown          = (TH1D*) histo0->Clone(Form("histo_WH_hww_SM_CMS_MVAJESBoundingDown"));

  TH1D* histo_Wjets_CMS_MVAWBoundingUp	             	= (TH1D*) histo0->Clone(Form("histo_Wjets_CMS_MVAWBoundingUp"));
  TH1D* histo_Wjets_CMS_MVAWBoundingDown             	= (TH1D*) histo0->Clone(Form("histo_Wjets_CMS_MVAWBoundingDown"));

  TH1D* histo_WZ_CMS_MVAWZNLOBoundingUp              	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAWZNLOBoundingUp"));
  TH1D* histo_WZ_CMS_MVAWZNLOBoundingDown            	= (TH1D*) histo0->Clone(Form("histo_WZ_CMS_MVAWZNLOBoundingDown"));

  TH1D* histo_ZZ_CMS_MVAZZNLOBoundingUp              	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAZZNLOBoundingUp"));
  TH1D* histo_ZZ_CMS_MVAZZNLOBoundingDown            	= (TH1D*) histo0->Clone(Form("histo_ZZ_CMS_MVAZZNLOBoundingDown"));

  if(useTemplates == true){
    histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVUp   ->Rebin(rebinMVAHist);
    histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVDown ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVUp   ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVDown ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAWZStatBounding_8TeVUp	    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAWZStatBounding_8TeVDown	    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAZZStatBounding_8TeVUp	    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAZZStatBounding_8TeVDown	    ->Rebin(rebinMVAHist);
    histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp     ->Rebin(rebinMVAHist);
    histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown   ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp   ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAVVVStatBounding_8TeVUp	    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAVVVStatBounding_8TeVDown	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVUp   ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVDown ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVUp   ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVDown ->Rebin(rebinMVAHist);
    
    histo_WH_htt_CMS_MVALepEffBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_CMS_MVALepEffBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVALepEffBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVALepEffBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVALepEffBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVALepEffBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVALepEffBoundingUp		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVALepEffBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVALepEffBoundingUp	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVALepEffBoundingDown	    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVALepEffBoundingUp		    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVALepEffBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVALepEffBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVALepEffBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVALepEffBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVALepEffBoundingDown	    ->Rebin(rebinMVAHist);

    histo_WH_htt_CMS_MVALepResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_CMS_MVALepResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVALepResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVALepResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVALepResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVALepResBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVALepResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVALepResBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVALepResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVALepResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVALepResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVALepResBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVALepResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVALepResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVALepResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVALepResBoundingDown	    ->Rebin(rebinMVAHist);
    
    histo_WH_htt_CMS_MVAMETResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_CMS_MVAMETResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAMETResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAMETResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAMETResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAMETResBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAMETResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAMETResBoundingDown  	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAMETResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAMETResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAMETResBoundingUp		    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAMETResBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAMETResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAMETResBoundingDown	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAMETResBoundingUp	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAMETResBoundingDown	    ->Rebin(rebinMVAHist);
    
    histo_WH_htt_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WH_htt_CMS_MVAJESBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WH_hww_CMS_MVAJESBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAJESBoundingDown		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAJESBoundingDown		    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_Wgamma_CMS_MVAJESBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAJESBoundingUp		    ->Rebin(rebinMVAHist);
    histo_VVV_CMS_MVAJESBoundingDown		    ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAJESBoundingUp            ->Rebin(rebinMVAHist);
    histo_WH_htt_SM_CMS_MVAJESBoundingDown 	    ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAJESBoundingUp            ->Rebin(rebinMVAHist);
    histo_WH_hww_SM_CMS_MVAJESBoundingDown 	    ->Rebin(rebinMVAHist);

    histo_Wjets_CMS_MVAWBoundingUp		    ->Rebin(rebinMVAHist);
    histo_Wjets_CMS_MVAWBoundingDown		    ->Rebin(rebinMVAHist);

    histo_WZ_CMS_MVAWZNLOBoundingUp		    ->Rebin(rebinMVAHist);
    histo_WZ_CMS_MVAWZNLOBoundingDown		    ->Rebin(rebinMVAHist);

    histo_ZZ_CMS_MVAZZNLOBoundingUp		    ->Rebin(rebinMVAHist);
    histo_ZZ_CMS_MVAZZNLOBoundingDown		    ->Rebin(rebinMVAHist);
  }

  double bgdDecayFake[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double bgdDecayFakeE[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  double bgdDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double weiDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double nSigCut[7]  = {0,0,0,0,0,0,0}; // Last channel is for h->tt
  double nSigCutE[7] = {0,0,0,0,0,0,0};
  unsigned int patternTopVeto         = SmurfTree::TopVeto;
  //unsigned int patternTopTag          = SmurfTree::TopTag;
  //unsigned int patternTopTagNotInJets = SmurfTree::TopTagNotInJets;
  Float_t bdtg = 0.0;
  Float_t bdtg_aux0 = 0.0;
  Float_t bdtg_aux1 = 0.0;
  Float_t bdtg_aux2 = 0.0;

  int nBgd=bgdEvent.tree_->GetEntries();
  if(channel==1101) {
    //bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l"       ,(int)120), &bdtg);
    //bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux0"  ,(int)120), &bdtg_aux0);
    //bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux1"  ,(int)120), &bdtg_aux1);
    //bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux2"  ,(int)120), &bdtg_aux2);
  }
  for (int i=0; i<nBgd; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);
    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::other &&
            TMath::Abs(bgdEvent.scale1fb_-0.0003856)/0.0003856<0.0001) fDecay = 41;
    else if(bgdEvent.processId_==121 ||
    	    bgdEvent.processId_==122 ||	                     
            bgdEvent.processId_==24  ||
            bgdEvent.processId_==26  ||
            bgdEvent.processId_==10001 ||
            bgdEvent.processId_==10010)                      fDecay = 42;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;assert(0);}
    if(fDecay == -1 || fDecay > 100) fDecay = 0;//44;
    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
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

    double usedMet = TMath::Min(pmetA,pmetB);
    Float_t mTWMin   = bgdEvent.mt1_;
    Float_t mTWMax   = bgdEvent.mt1_;
    if(mTWMin > bgdEvent.mt2_                        ) mTWMin = bgdEvent.mt2_;
    if(mTWMin > bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMin = bgdEvent.mt3_;
    if(mTWMax < bgdEvent.mt2_	 		     ) mTWMax = bgdEvent.mt2_;
    if(mTWMax < bgdEvent.mt3_ && bgdEvent.lid3_ != 0.) mTWMax = bgdEvent.mt3_; 

    double Mjj = (bgdEvent.jet1_+bgdEvent.jet2_).M();
    if(bgdEvent.njets_ >= 3){
      if(TMath::Abs((bgdEvent.jet1_+bgdEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (bgdEvent.jet1_+bgdEvent.jet3_).M();
      if(TMath::Abs((bgdEvent.jet2_+bgdEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (bgdEvent.jet2_+bgdEvent.jet3_).M();
    }
    if(bgdEvent.njets_ >= 4){
      if(TMath::Abs((bgdEvent.jet1_+bgdEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (bgdEvent.jet1_+bgdEvent.jet4_).M();
      if(TMath::Abs((bgdEvent.jet2_+bgdEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (bgdEvent.jet2_+bgdEvent.jet4_).M();
      if(TMath::Abs((bgdEvent.jet3_+bgdEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (bgdEvent.jet3_+bgdEvent.jet4_).M();
    }

    double massZMin = 999.0; double massMin = 999.0; double massminSFOS = 999.0; double massminSFSS = 999.0; double mtQQLN = -1.0; double mQQLN = -1.0;
    double dRMin = 999.0; double mTW3 = 0.0;double mTW = 0.0; double type3l = 0.0; double mTFromW = 0.0; double ptlw = 0.0; double which3rdl = -1; double new3LVar = 0.0;
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
      mTW   = trilepton_info(3,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      mTW3  = trilepton_info(10,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      type3l= trilepton_info(4,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      massminSFOS  = trilepton_info(5,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      massminSFSS  = trilepton_info(6,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      mTFromW     = trilepton_info(7,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      ptlw        = trilepton_info(8,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      which3rdl   = trilepton_info(9,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_);
      new3LVar   = trilepton_info(thePlot,bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lep3_,
                                  bgdEvent.lq1_ ,bgdEvent.lq2_ ,bgdEvent.lq3_,
		                  bgdEvent.lid1_,bgdEvent.lid2_,bgdEvent.lid3_,
				  bgdEvent.mt1_ ,bgdEvent.mt2_ ,bgdEvent.mt3_,bgdEvent.met_, bgdEvent.metPhi_);
      double pQQLN[5] = {bgdEvent.jet1_.Px() + bgdEvent.jet2_.Px() + newMet*cos( bgdEvent.metPhi_),
                         bgdEvent.jet1_.Py() + bgdEvent.jet2_.Py() + newMet*sin( bgdEvent.metPhi_),
			 bgdEvent.jet1_.Pz() + bgdEvent.jet2_.Pz(),
			 bgdEvent.jet1_.P()  + bgdEvent.jet2_.P() + newMet,bgdEvent.jet1_.Pt()  + bgdEvent.jet2_.Pt() + newMet};
      if     (which3rdl == 1) {pQQLN[0]+=bgdEvent.lep1_.Px();pQQLN[1]+=bgdEvent.lep1_.Py();pQQLN[2]+=bgdEvent.lep1_.Pz();pQQLN[3]+=bgdEvent.lep1_.P();pQQLN[4]+=bgdEvent.lep1_.Pt();}
      else if(which3rdl == 2) {pQQLN[0]+=bgdEvent.lep2_.Px();pQQLN[1]+=bgdEvent.lep2_.Py();pQQLN[2]+=bgdEvent.lep2_.Pz();pQQLN[3]+=bgdEvent.lep2_.P();pQQLN[4]+=bgdEvent.lep2_.Pt();}
      else if(which3rdl == 3) {pQQLN[0]+=bgdEvent.lep3_.Px();pQQLN[1]+=bgdEvent.lep3_.Py();pQQLN[2]+=bgdEvent.lep3_.Pz();pQQLN[3]+=bgdEvent.lep3_.P();pQQLN[4]+=bgdEvent.lep3_.Pt();}
      mQQLN  = pQQLN[3]*pQQLN[3]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1]-pQQLN[2]*pQQLN[2]; if(mQQLN  > 0) mQQLN  = sqrt(mQQLN);  else mQQLN   = 0.0;
      mtQQLN = pQQLN[4]*pQQLN[4]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1];                   if(mtQQLN > 0) mtQQLN = sqrt(mtQQLN); else mtQQLN  = 0.0;
    }    
    bdtg = ((dRMin/4.0)-0.5)*2.0; /*bdtg = ((massMin/100.0)-0.5)*2.0; bdtg = Unroll2VarTo1ForWH(dRMin,mTW3,0); bdtg = Unroll2VarTo1ForWH(massMin,mTW3,1);*/ if(bdtg >= 1.0) bdtg = 0.999; if(bdtg <= -1) bdtg = -0.999;

    double deltaPhiQQL = -1;
    if     (which3rdl == 1) deltaPhiQQL = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep1_.Phi());
    else if(which3rdl == 2) deltaPhiQQL = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep2_.Phi());
    else if(which3rdl == 3) deltaPhiQQL = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep3_.Phi());

    bool cleanMode = true;
    if(nsel!=0&&
      (type3l==1||type3l==8||
      (type3l==2&&bgdEvent.lq1_*bgdEvent.lq2_<0)||(type3l==3&&bgdEvent.lq1_*bgdEvent.lq3_<0)||(type3l==4&&bgdEvent.lq2_*bgdEvent.lq3_<0)||
      (type3l==5&&bgdEvent.lq2_*bgdEvent.lq3_<0)||(type3l==6&&bgdEvent.lq1_*bgdEvent.lq2_<0)||(type3l==7&&bgdEvent.lq1_*bgdEvent.lq3_<0))) cleanMode = false;
    bool passMode = cleanMode; if(nsel==2) passMode = !passMode;
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 27){ // WZ selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 40. &&
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.jet1_.Pt() <= 40. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 1 == 1
	){
	passCuts = true;
	if(bgdEvent.dstype_ == SmurfTree::wz) isSignalDecay = true;
      }
    } // WZ selection
    if(channel == 28){ // 3lnq selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         bgdEvent.jet1_.Pt() > 30. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 usedMet < 50 &&
	 bgdEvent.jet1Btag_ < 2.1 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // 3lnq selection
    if(channel == 1098){ // top region
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 30 &&
         (bgdEvent.jet1_.Pt() > 40. || bgdEvent.jet1Btag_ > 2.1) &&
         massZMin > 15 &&
         massMin > 12 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // top region
    if(channel == 1099){ // Zjets region
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
         bgdEvent.lep3_.Pt() < 15. &&
         usedMet < 25 &&
         bgdEvent.jet1_.Pt() < 40. &&
         //massZMin < 15 &&
	 abs(bgdEvent.dilep_.M()-91.1876) < 15. &&
         massMin > 12 &&
	 //(type3l == 1 || type3l == 3) &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // Zjets region
    if(channel == 1097 || channel == 1096){ // ttX selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      int nBtag0 = (bgdEvent.jet1Btag_ > 2.1) + (bgdEvent.jet2Btag_ > 2.1);
      int nBtag1 = (bgdEvent.jet1Btag_ > 2.1) + (bgdEvent.jet2Btag_ > 2.1) + (bgdEvent.jet3Btag_ > 2.1) + (bgdEvent.jet4Btag_ > 2.1);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
         bgdEvent.lep3_.Pt() > 20. &&
         usedMet > 20. &&
       ((channel == 1097 && bgdEvent.njets_ >= 2 && bgdEvent.njets_ <= 3 &&
         nBtag0 == 2 && massZMin > 10) ||
	(channel == 1096 && bgdEvent.njets_ >= 4 && massZMin < 15 &&
         nBtag1 >= 2)) &&
         massMin > 12 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // ttX selection
    if(channel == 1100){ // HW->3l selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
         dRMin < 2.0 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection
    if(channel == 1101){ // HW->3l selection-BDT
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         bgdEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection-BDT
    if(channel == 1300 || channel == 1301){ // ZG selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 10. &&
	 massMin > 10 &&
	 massZMin > 10 &&
	 TMath::Abs((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.lep3_).M()-91.1976) < 15.0 &&
	 (type3l==2||type3l==3||type3l==4||type3l==8) &&
         //(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 1 == 1
	){
	if(type3l==1||type3l==8||
	  (type3l==2&&bgdEvent.lq1_*bgdEvent.lq2_<0)||(type3l==3&&bgdEvent.lq1_*bgdEvent.lq3_<0)||(type3l==4&&bgdEvent.lq2_*bgdEvent.lq3_<0)||
	  (type3l==5&&bgdEvent.lq2_*bgdEvent.lq3_<0)||(type3l==6&&bgdEvent.lq1_*bgdEvent.lq2_<0)||(type3l==7&&bgdEvent.lq1_*bgdEvent.lq3_<0)){
	passCuts = true;
	if(bgdEvent.dstype_ == SmurfTree::wgamma) isSignalDecay = true;
	}
      }
    } // ZG selection
    if(channel == 1200){ // WG* selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
         bgdEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
         bgdEvent.lep3_.Pt() > 3. &&
	 mTWMin > 25. && mTFromW > 45. &&
	(type3l == 1 || type3l == 2 || type3l == 3 || type3l == 4) &&
	 massminSFOS > 0 && massminSFOS < 12 && TMath::Abs(massminSFOS-3.1) > 0.1 &&
	(bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.jet1Btag_ < 2.1 && bgdEvent.jet2Btag_ < 2.1 && bgdEvent.jet3Btag_ < 2.1) &&
	 bgdEvent.njets_ <= 3 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WG* selection
    if(channel == 1102){ // HZ->3lqq selection
      charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_ + bgdEvent.lq3_);
      if(
           bgdEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           bgdEvent.lep1_.Pt() > 20. &&
           bgdEvent.lep2_.Pt() > 10. &&
           bgdEvent.lep3_.Pt() > 10. &&
	   (bgdEvent.jet1_+bgdEvent.jet2_).M() > 30. && (bgdEvent.jet1_+bgdEvent.jet2_).M() < 200. &&
           //(bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           bgdEvent.njets_ >= 2 &&
	   massZMin < 15 &&
	   //bgdEvent.met_ > 25 &&
	   //mTW > 30 &&
	   massMin > 12 &&
	   //mtQQLN > mhAna-70 && mtQQLN < mhAna+30 &&
           deltaPhiQQL*180.0/TMath::Pi() < 150.0 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HZ->3lqq selection
  
    if(passCuts == true){
      double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	 							fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        						       TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	 							fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
        						       TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     	 							fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        						       TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      double trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
      bool isRealLepton = false;
      if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13) &&
         (TMath::Abs(bgdEvent.lep3McId_) == 11 || TMath::Abs(bgdEvent.lep3McId_) == 13)) isRealLepton = true;
      double theWeight = 0.0;
      double add       = 1.0;
      //if(fDecay == 41){ // SMHtt
      //  add = add * 0.969212;
      //}
      //else if(fDecay == 42){ // SMHww
      //  add = add * 1.027976;
      //}
      int nFake = 0;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake > 1){
	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	fDecay = 14;
	theWeight	       = -1.0 * add;
	bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += theWeight*theWeight;
	printf("WARNING, double fakes, you should know what you are doing. Number of fakes: %d\n",nFake);
      }
      else if(nFake == 1){
        if(bgdEvent.dstype_ == SmurfTree::data){
	  //if(channel == 1200) add = 3.35;
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	  if(applyFRCorr == true){
	    if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)||
	       ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection)||
	       ((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection))
	       add = add * 0.78;
	    else 
	    if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4)  == SmurfTree::Lep1LooseEleV4  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)||
	       ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4)  == SmurfTree::Lep2LooseEleV4  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection)||
	       ((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4)  == SmurfTree::Lep3LooseEleV4  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection))
	       add = add * 0.78;
	  }
	  theWeight              = add;
	  bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += theWeight*theWeight;
	}
	//else if(fDecay == 5 || fDecay == 13){  
	else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){  
	  if(channel == 1200) add = 0.0;
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
	  add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);
	  add = add*nPUScaleFactor2012(fhDPU,(float)bgdEvent.npu_);

	  if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,bgdEvent.met_);
	  fDecay                 = 3;
	  isSignalDecay          = false;
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add*trigEff;
          bgdDecay[(int)fDecay] += theWeight;
          weiDecay[(int)fDecay] += TMath::Abs(theWeight)*TMath::Abs(theWeight);
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven || bgdEvent.dstype_ == SmurfTree::qcd) {
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
      else if((bgdEvent.dstype_ != SmurfTree::data && isRealLepton == true) || bgdEvent.dstype_ == SmurfTree::wgamma){
      //else if(fDecay == 5 || fDecay == 13){
        if(bgdEvent.dstype_ == SmurfTree::wz) add = 1.04;
        if(bgdEvent.dstype_ == SmurfTree::wgamma && channel != 1300) add = 1.40;
        if(bgdEvent.dstype_ == SmurfTree::zz && channel == 1301) add = 2.2;
	add = add*nPUScaleFactor2012(fhDPU,(float)bgdEvent.npu_);
        if((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
	  add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	if((bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
	if((bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)
	  add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

  	//if(bgdEvent.dstype_ == SmurfTree::www) add = bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_;
        //if(bgdEvent.dstype_ == SmurfTree::www && channel == 1098) add = 0.0;

        if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,bgdEvent.met_);
	if(bgdEvent.dstype_ == SmurfTree::wgstar && channel == 1200) add = 0.0; // to measure cross-section
	// CAREFUL HERE, no data-driven corrections, just Higgs k-factors
	// add = 1.0;

        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add*trigEff;
        bgdDecay[(int)fDecay] += theWeight;
        weiDecay[(int)fDecay] += theWeight*theWeight;
      }
      else {
        theWeight = 0.0;
      }

      double myVar = newMet;
      if     (thePlot == 1) myVar = TMath::Min((double)bgdEvent.lep1_.Pt(),199.999);
      else if(thePlot == 2) myVar = TMath::Min((double)bgdEvent.lep2_.Pt(),199.999);
      else if(thePlot == 3) myVar = TMath::Min((double)bgdEvent.lep3_.Pt(),199.999);
      else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = bgdEvent.dilep_.M();
      else if(thePlot == 8) myVar = bdtg;
      else if(thePlot == 9) myVar = TMath::Min((double)bgdEvent.mt1_,199.999);
      else if(thePlot ==10) myVar = TMath::Min((double)bgdEvent.mt2_,199.999);
      else if(thePlot ==11) myVar = TMath::Min((double)bgdEvent.mt3_,199.999);
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = newMet/bgdEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = bgdEvent.njets_;
      else if(thePlot ==18) myVar = bgdEvent.nvtx_;
      else if(thePlot ==19) myVar = (bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.lep3_).M();
      else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(bgdEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(bgdEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(bgdEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(bgdEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(bgdEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(bgdEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = mQQLN;
      else if(thePlot ==31) myVar = mtQQLN;
      else if(thePlot ==32) myVar = massminSFOS;
      else if(thePlot ==33) myVar = massminSFSS;
      else if(thePlot ==34) myVar = (bgdEvent.jet1_+bgdEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = (bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.lep3_+bgdEvent.jet1_).M();
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = TMath::Min(dRMin,3.999);
      else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = Mjj;
      else if(thePlot ==46) myVar = ptlw;
      else if(thePlot ==47) myVar = mTFromW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = newMet*cos(bgdEvent.metPhi_);
      else if(thePlot ==50) myVar = newMet*sin(bgdEvent.metPhi_);
      else if(thePlot ==51) myVar = newTrackMet*cos(bgdEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = newTrackMet*sin(bgdEvent.trackMetPhi_);
      else if(thePlot ==53) {if(bgdEvent.jet1_.Pt()>15&&bgdEvent.jet2_.Pt()>15)myVar = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi();else myVar=0.0;}
      else if(thePlot ==54) myVar = mTW;
      else if(thePlot>=61&&thePlot<=66) myVar = new3LVar;

      int nBinC = 0;
      if     (fDecay ==  0 || fDecay ==  1 || fDecay ==  2 || fDecay == 3 || //wjets
              fDecay ==  5 || fDecay == 13 || //ttbar/wt
	      fDecay == 14 || fDecay == 29 || fDecay == 30 || // ww
	      fDecay ==  6 || fDecay ==  7 || fDecay ==  8 || fDecay == 9 || fDecay == 10 || fDecay == 31){ // zjets
        histo4->Fill(myVar,theWeight);
	histo_Wjets->Fill(myVar,theWeight);
        nBinC = 3;
      }
      else if(fDecay == 19){ // W/Zgamma
        if(isSignalDecay == true) histo6->Fill(myVar,theWeight);
	else                      histo4->Fill(myVar,theWeight);
	histo_Wgamma->Fill(myVar,theWeight);
        nBinC = 4;
      }
      else if(fDecay == 27 || fDecay == 20){ // WZ/Wgamma*
        histo3->Fill(myVar,theWeight);
	histo_WZ->Fill(myVar,theWeight);
        nBinC = 1;
      }
      else if(fDecay == 21){ // VVV
        histo3->Fill(myVar,theWeight);
	histo_VVV->Fill(myVar,theWeight);
        nBinC = 5;
      }
      else if(fDecay == 28){ // ZZ
        histo2->Fill(myVar,theWeight);
	histo_ZZ->Fill(myVar,theWeight);
        nBinC = 2;
      }
      else if(fDecay == 41){ // SMHtt
        histo_WH_htt_SM->Fill(myVar,theWeight);
        nBinC = 6;
      }
      else if(fDecay == 42){ // SMHww
        histo_WH_hww_SM->Fill(myVar,theWeight);
        nBinC = 7;
      }
      else if(theWeight != 0){
        histo4->Fill(myVar,theWeight);
        printf("NOOOOOOOOOOOOOOOOOOOO %d\n",fDecay);
      }

      if(useTemplates == true){
        bool passJetCut[3] = {bgdEvent.jet1_.Pt() <= 40., bgdEvent.jet1_.Pt()*1.05 <= 40., bgdEvent.jet1_.Pt()*0.95 <= 40.};
        double addLepEff   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
                             leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_)*
			     leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);
        double addLepEffUp = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,+1)*
                             leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,+1)*
			     leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_,+1);
        double addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
                               leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1)*
			       leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_,-1);
        if     (nBinC == 1){
	  histo_WZ_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_WZ_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_WZ_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WZ_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WZ_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_WZ_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_WZ_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
        else if(nBinC == 2){
	  histo_ZZ_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_ZZ_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_ZZ_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_ZZ_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_ZZ_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_ZZ_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
        else if(nBinC == 3){
	  double addFR = fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
		        								       (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFR = addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                								      (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          addFR = addFR*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                								      (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	  double addFRS = fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											                (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											                (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          addFRS = addFRS*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
											                (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
          histo_Wjets_CMS_MVAWBoundingUp       ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addFRS/addFR);
	}
        else if(nBinC == 4){
	  histo_Wgamma_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_Wgamma_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_Wgamma_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_Wgamma_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_Wgamma_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_Wgamma_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_Wgamma_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
        else if(nBinC == 5){
	  histo_VVV_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_VVV_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_VVV_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_VVV_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_VVV_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_VVV_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_VVV_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
        else if(nBinC == 6){
	  histo_WH_htt_SM_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_WH_htt_SM_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_WH_htt_SM_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_htt_SM_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_htt_SM_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_WH_htt_SM_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_WH_htt_SM_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
        else if(nBinC == 7){
	  histo_WH_hww_SM_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_WH_hww_SM_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_WH_hww_SM_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_hww_SM_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_hww_SM_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_WH_hww_SM_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_WH_hww_SM_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
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

  if(systInputFile != "" && useTemplates == true){
    int nSyst=systEvent.tree_->GetEntries();
    if(channel==1101) {
      //systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l"	    ,(int)120), &bdtg);
      //systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux0"  ,(int)120), &bdtg_aux0);
      //systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux1"  ,(int)120), &bdtg_aux1);
      //systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux2"  ,(int)120), &bdtg_aux2);
    }
    for (int i=0; i<nSyst; ++i) {

      if (i%100000 == 0 && verboseLevel > 0)
	printf("--- reading syst event %5d of %5d\n",i,nSyst);
      systEvent.tree_->GetEntry(i);
      int fDecay = 0;
      if     (systEvent.dstype_ == SmurfTree::wjets 	    ) fDecay = 3;
      else if(systEvent.dstype_ == SmurfTree::ttbar 	    ) fDecay = 5;
      else if(systEvent.dstype_ == SmurfTree::dyee  	    ) fDecay = 9;
      else if(systEvent.dstype_ == SmurfTree::dymm  	    ) fDecay = 9;
      else if(systEvent.dstype_ == SmurfTree::dytt  	    ) fDecay = 10;
      else if(systEvent.dstype_ == SmurfTree::tw    	    ) fDecay = 13;
      else if(systEvent.dstype_ == SmurfTree::qqww  	    ) fDecay = 29;
      else if(systEvent.dstype_ == SmurfTree::wz    	    ) fDecay = 27;
      else if(systEvent.dstype_ == SmurfTree::zz    	    ) fDecay = 28;
      else if(systEvent.dstype_ == SmurfTree::ggww  	    ) fDecay = 30;
      else if(systEvent.dstype_ == SmurfTree::wgamma	    ) fDecay = 19;
      else if(systEvent.dstype_ == SmurfTree::data  	    ) fDecay =  1;
      else if(systEvent.dstype_ == SmurfTree::ggzz  	    ) fDecay = 18;
      else if(systEvent.dstype_ == SmurfTree::dyttDataDriven) fDecay = 10;
      else if(systEvent.dstype_ == SmurfTree::wgstar        ) fDecay = 20;
      else if(systEvent.dstype_ == SmurfTree::www           ) fDecay = 21;
      else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;assert(0);}
      if(fDecay == -1 || fDecay > 100) fDecay = 0;//44;
      int charge = (int)(systEvent.lq1_ + systEvent.lq2_);
      double newMet      = systEvent.met_;
      double newTrackMet = systEvent.trackMet_;
      bool applyCor = false;
      if(applyCor == true && systEvent.dstype_ != SmurfTree::data){
	double metx=0.0;double mety=0.0;double trkmetx=0.0;double trkmety=0.0;
	if      (systEvent.njets_ == 0){
          metx	= newMet*cos(systEvent.metPhi_)+gRandom->Gaus(0.0,4.8);
          mety	= newMet*sin(systEvent.metPhi_)+gRandom->Gaus(0.0,4.8);
          trkmetx = newTrackMet*cos(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,1.4);
          trkmety = newTrackMet*sin(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,1.4);
	}
	else if(systEvent.njets_ == 1){
          metx	= newMet*cos(systEvent.metPhi_)+gRandom->Gaus(0.0,4.9);
          mety	= newMet*sin(systEvent.metPhi_)+gRandom->Gaus(0.0,4.9);
          trkmetx = newTrackMet*cos(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.4);
          trkmety = newTrackMet*sin(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.4);
	}
	else if(systEvent.njets_ >= 2){
          metx	= newMet*cos(systEvent.metPhi_)+gRandom->Gaus(0.0,5.0);
          mety	= newMet*sin(systEvent.metPhi_)+gRandom->Gaus(0.0,5.0);
          trkmetx = newTrackMet*cos(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.8);
          trkmety = newTrackMet*sin(systEvent.trackMetPhi_)+gRandom->Gaus(0.0,3.8);
	}
	newMet      = sqrt(metx*metx+mety*mety);
	newTrackMet = sqrt(trkmetx*trkmetx+trkmety*trkmety);
      }
      double deltaPhiA[3] = {TMath::Abs(systEvent.lep1_.Phi()-systEvent.metPhi_),TMath::Abs(systEvent.lep2_.Phi()-systEvent.metPhi_),0.0};
      while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
      while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
      deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
      double pmetA = newMet;
      if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

      double deltaPhiB[3] = {TMath::Abs(systEvent.lep1_.Phi()-systEvent.trackMetPhi_),TMath::Abs(systEvent.lep2_.Phi()-systEvent.trackMetPhi_),0.0};
      while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
      while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
      deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
      double pmetB = newTrackMet;
      if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

      double deltaPhiC[3] = {TMath::Abs(systEvent.jet1_.Phi()-systEvent.metPhi_),TMath::Abs(systEvent.jet2_.Phi()-systEvent.metPhi_),0.0};
      while(deltaPhiC[0]>TMath::Pi()) deltaPhiC[0] = TMath::Abs(deltaPhiC[0] - 2*TMath::Pi());
      while(deltaPhiC[1]>TMath::Pi()) deltaPhiC[1] = TMath::Abs(deltaPhiC[1] - 2*TMath::Pi());
      double pmetC = newMet;

      double deltaPhiD[3] = {TMath::Abs(systEvent.jet1_.Phi()-systEvent.trackMetPhi_),TMath::Abs(systEvent.jet2_.Phi()-systEvent.trackMetPhi_),0.0};
      while(deltaPhiD[0]>TMath::Pi()) deltaPhiD[0] = TMath::Abs(deltaPhiD[0] - 2*TMath::Pi());
      while(deltaPhiD[1]>TMath::Pi()) deltaPhiD[1] = TMath::Abs(deltaPhiD[1] - 2*TMath::Pi());
      double pmetD = newTrackMet;

      if     (systEvent.njets_ == 0){
	pmetC = pmetA;
	pmetD = pmetB;
      }
      else if(systEvent.njets_ == 1){
	deltaPhiC[2] = deltaPhiC[0];
	deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
	if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
	deltaPhiD[2] = deltaPhiD[0];
	deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
	if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
      }
      else if(systEvent.njets_ >= 2){
	deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiC[1]);
	deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
	if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
	deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiD[1]);
	deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
	if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
      }

      double usedMet = TMath::Min(pmetA,pmetB);
      Float_t mTWMin   = systEvent.mt1_;
      Float_t mTWMax   = systEvent.mt1_;
      if(mTWMin > systEvent.mt2_                        ) mTWMin = systEvent.mt2_;
      if(mTWMin > systEvent.mt3_ && systEvent.lid3_ != 0.) mTWMin = systEvent.mt3_;
      if(mTWMax < systEvent.mt2_	 		     ) mTWMax = systEvent.mt2_;
      if(mTWMax < systEvent.mt3_ && systEvent.lid3_ != 0.) mTWMax = systEvent.mt3_; 

      double Mjj = (systEvent.jet1_+systEvent.jet2_).M();
      if(systEvent.njets_ >= 3){
        if(TMath::Abs((systEvent.jet1_+systEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (systEvent.jet1_+systEvent.jet3_).M();
        if(TMath::Abs((systEvent.jet2_+systEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (systEvent.jet2_+systEvent.jet3_).M();
      }
      if(systEvent.njets_ >= 4){
        if(TMath::Abs((systEvent.jet1_+systEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (systEvent.jet1_+systEvent.jet4_).M();
        if(TMath::Abs((systEvent.jet2_+systEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (systEvent.jet2_+systEvent.jet4_).M();
        if(TMath::Abs((systEvent.jet3_+systEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (systEvent.jet3_+systEvent.jet4_).M();
      }

      double massZMin = 999.0; double massMin = 999.0; double mtQQLN = -1.0; double mQQLN = -1.0;
      double dRMin = 999.0; /*double mTW3 = 0.0;double mTW = 0.0;*/ double which3rdl = -1; double type3l = 0.0;
      if(systEvent.lid3_ != 0){
	usedMet = TMath::Min(newMet,newTrackMet);
	dRMin = trilepton_info(2,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                 systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                 systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				 systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
	/*mTW   = trilepton_info(3,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                    systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                    systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				    systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
	mTW3  = trilepton_info(10,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                    systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                    systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				    systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);*/
	massZMin = trilepton_info(0,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                    systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                    systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				    systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
	massMin  = trilepton_info(1,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                    systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                    systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				    systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
        which3rdl   = trilepton_info(9,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                    systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                    systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				    systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
        type3l= trilepton_info(4,systEvent.lep1_,systEvent.lep2_,systEvent.lep3_,
                                 systEvent.lq1_ ,systEvent.lq2_ ,systEvent.lq3_,
		                 systEvent.lid1_,systEvent.lid2_,systEvent.lid3_,
				 systEvent.mt1_ ,systEvent.mt2_ ,systEvent.mt3_);
	double pQQLN[5] = {systEvent.jet1_.Px() + systEvent.jet2_.Px() + newMet*cos( systEvent.metPhi_),
                           systEvent.jet1_.Py() + systEvent.jet2_.Py() + newMet*sin( systEvent.metPhi_),
			   systEvent.jet1_.Pz() + systEvent.jet2_.Pz(),
			   systEvent.jet1_.P()  + systEvent.jet2_.P() + newMet,systEvent.jet1_.Pt()  + systEvent.jet2_.Pt() + newMet};
	if     (which3rdl == 1) {pQQLN[0]+=systEvent.lep1_.Px();pQQLN[1]+=systEvent.lep1_.Py();pQQLN[2]+=systEvent.lep1_.Pz();pQQLN[3]+=systEvent.lep1_.P();pQQLN[4]+=systEvent.lep1_.Pt();}
	else if(which3rdl == 2) {pQQLN[0]+=systEvent.lep2_.Px();pQQLN[1]+=systEvent.lep2_.Py();pQQLN[2]+=systEvent.lep2_.Pz();pQQLN[3]+=systEvent.lep2_.P();pQQLN[4]+=systEvent.lep2_.Pt();}
	else if(which3rdl == 3) {pQQLN[0]+=systEvent.lep3_.Px();pQQLN[1]+=systEvent.lep3_.Py();pQQLN[2]+=systEvent.lep3_.Pz();pQQLN[3]+=systEvent.lep3_.P();pQQLN[4]+=systEvent.lep3_.Pt();}
	mQQLN  = pQQLN[3]*pQQLN[3]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1]-pQQLN[2]*pQQLN[2]; if(mQQLN  > 0) mQQLN  = sqrt(mQQLN);  else mQQLN   = 0.0;
	mtQQLN = pQQLN[4]*pQQLN[4]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1];                   if(mtQQLN > 0) mtQQLN = sqrt(mtQQLN); else mtQQLN  = 0.0;
      }
      bdtg = ((dRMin/4.0)-0.5)*2.0; /*bdtg = ((massMin/100.0)-0.5)*2.0; bdtg = Unroll2VarTo1ForWH(dRMin,mTW3,0); bdtg = Unroll2VarTo1ForWH(massMin,mTW3,1);*/ if(bdtg >= 1.0) bdtg = 0.999; if(bdtg <= -1) bdtg = -0.999;

      bool cleanMode = true;
      if(nsel!=0&&
        (type3l==1||type3l==8||
        (type3l==2&&systEvent.lq1_*systEvent.lq2_<0)||(type3l==3&&systEvent.lq1_*systEvent.lq3_<0)||(type3l==4&&systEvent.lq2_*systEvent.lq3_<0)||
        (type3l==5&&systEvent.lq2_*systEvent.lq3_<0)||(type3l==6&&systEvent.lq1_*systEvent.lq2_<0)||(type3l==7&&systEvent.lq1_*systEvent.lq3_<0))) cleanMode = false;
      bool passMode = cleanMode; if(nsel==2) passMode = !passMode;
      bool passCuts = false;
      charge = (int)(systEvent.lq1_ + systEvent.lq2_ + systEvent.lq3_);
      if(
         systEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         systEvent.lep1_.Pt() > 20. &&
         systEvent.lep2_.Pt() > 10. &&
         systEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (systEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         systEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
	 passMode == true &&
         1 == 1
        ){
        passCuts = true;
      }

      if(passCuts == true){
	double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
     	 							  fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
        							 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
	double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
     	 							  fabs(systEvent.lep3_.Eta()), systEvent.lep3_.Pt(), 
        							 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid3_));
	double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep3_.Eta()), systEvent.lep3_.Pt() , 
     	 							  fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
        							 TMath::Abs( systEvent.lid3_), TMath::Abs(systEvent.lid2_));
	double trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
	bool isRealLepton = false;
	if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
           (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13) &&
           (TMath::Abs(systEvent.lep3McId_) == 11 || TMath::Abs(systEvent.lep3McId_) == 13)) isRealLepton = true;
	double theWeight = 0.0;
	double add       = 1.0;
	int nFake = 0;
	if(((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
	if(((systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
	if(((systEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
	if(((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
	if(((systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
	if(((systEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
	if(nFake == 1){
          if(systEvent.dstype_ == SmurfTree::data){
	    theWeight              = 0.0;
	  }
	  else {  
	    add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
            add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
            add = add*fakeRate(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
            add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	    add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
	    add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);
	    add = add*nPUScaleFactor2012(fhDPU,(float)systEvent.npu_);
	    fDecay                 = 3;
	    theWeight              = +1.0 * systEvent.scale1fb_*lumi*add*trigEff;
	  }
	}
	else if((systEvent.dstype_ != SmurfTree::data && isRealLepton == true) || systEvent.dstype_ == SmurfTree::wgamma){
	//else if(systEvent.dstype_ != SmurfTree::data){
          if(systEvent.dstype_ == SmurfTree::wz) add = 1.04;
	  add = add*nPUScaleFactor2012(fhDPU,(float)systEvent.npu_);
          if((systEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
	    add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	  if((systEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
	    add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
	  if((systEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)
	    add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

          if(systEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(systEvent.type_,systEvent.met_);
	  // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
	  // add = 1.0;
	  theWeight              = systEvent.scale1fb_*lumi*add*trigEff;
	}
	else {
          theWeight = 0.0;
	}

	int nBinC = 0;
	if     (fDecay ==  0 || fDecay ==  1 || fDecay ==  2 || fDecay == 3 || fDecay == 19 || //wjets/gamma
        	fDecay ==  5 || fDecay == 13 || //ttbar/wt
		fDecay == 14 || fDecay == 29 || fDecay == 30 || // ww
		fDecay ==  6 || fDecay ==  7 || fDecay ==  8 || fDecay == 9 || fDecay == 10 || fDecay == 31){ // zjets
          nBinC = 3;
	}
	else if(fDecay == 27 || fDecay == 20){
          nBinC = 1;
	}
	else if(fDecay == 28){
          nBinC = 2;
	}
	else if(theWeight != 0){
          printf("NOOOOOOOOOOOOOOOOOOOO %d\n",fDecay);
	}
	if(useTemplates == true){
          if     (nBinC == 1) histo_WZ_CMS_MVAWZNLOBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          else if(nBinC == 2) histo_ZZ_CMS_MVAZZNLOBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
      }
    } // end systematics loop
  }
  
  if((channel >= 0 && channel <= 8000)){
    int nSig=sigEvent.tree_->GetEntries();
    if(channel==1101) {
      //sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l"	   ,(int)120), &bdtg);
      //sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux0"  ,(int)120), &bdtg_aux0);
      //sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux1"  ,(int)120), &bdtg_aux1);
      //sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux2"  ,(int)120), &bdtg_aux2);
    }
    for (int i=0; i<nSig; ++i) {

      if (i%100000 == 0 && verboseLevel > 0)
	printf("--- reading Signal event %5d of %5d\n",i,nSig);
      sigEvent.tree_->GetEntry(i);
      int charge = (int)(sigEvent.lq1_ + sigEvent.lq2_);
      bool lid = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
	         (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
	         (sigEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection;
      if(!lid) continue;

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

      double usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
      Float_t mTWMin   = sigEvent.mt1_;
      Float_t mTWMax   = sigEvent.mt1_;
      if(mTWMin > sigEvent.mt2_			     ) mTWMin = sigEvent.mt2_;
      if(mTWMin > sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMin = sigEvent.mt3_;
      if(mTWMax < sigEvent.mt2_			     ) mTWMax = sigEvent.mt2_;
      if(mTWMax < sigEvent.mt3_ && sigEvent.lid3_ != 0.) mTWMax = sigEvent.mt3_;

      double Mjj = (sigEvent.jet1_+sigEvent.jet2_).M();
      if(sigEvent.njets_ >= 3){
        if(TMath::Abs((sigEvent.jet1_+sigEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (sigEvent.jet1_+sigEvent.jet3_).M();
        if(TMath::Abs((sigEvent.jet2_+sigEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (sigEvent.jet2_+sigEvent.jet3_).M();
      }
      if(sigEvent.njets_ >= 4){
        if(TMath::Abs((sigEvent.jet1_+sigEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (sigEvent.jet1_+sigEvent.jet4_).M();
        if(TMath::Abs((sigEvent.jet2_+sigEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (sigEvent.jet2_+sigEvent.jet4_).M();
        if(TMath::Abs((sigEvent.jet3_+sigEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (sigEvent.jet3_+sigEvent.jet4_).M();
      }

      double massZMin = 999.0; double massMin = 999.0; double massminSFOS = 999.0; double massminSFSS = 999.0; double mtQQLN = -1.0; double mQQLN = -1.0;
      double dRMin = 999.0; double mTW3 = 0.0;double mTW = 0.0; double type3l = 0.0; double mTFromW = 0.0; double ptlw = 0.0; double which3rdl = -1; double new3LVar = 0.0;
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
	mTW   = trilepton_info(3,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	mTW3  = trilepton_info(10,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	type3l= trilepton_info(4,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	massminSFOS  = trilepton_info(5,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	massminSFSS  = trilepton_info(6,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	mTFromW     = trilepton_info(7,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	ptlw        = trilepton_info(8,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
	which3rdl   = trilepton_info(9,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                    sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                    sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				    sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_);
       new3LVar   = trilepton_info(thePlot,sigEvent.lep1_,sigEvent.lep2_,sigEvent.lep3_,
                                  sigEvent.lq1_ ,sigEvent.lq2_ ,sigEvent.lq3_,
		                  sigEvent.lid1_,sigEvent.lid2_,sigEvent.lid3_,
				  sigEvent.mt1_ ,sigEvent.mt2_ ,sigEvent.mt3_,sigEvent.met_, sigEvent.metPhi_);
	double pQQLN[5] = {sigEvent.jet1_.Px() + sigEvent.jet2_.Px() + newMet*cos( sigEvent.metPhi_),
                           sigEvent.jet1_.Py() + sigEvent.jet2_.Py() + newMet*sin( sigEvent.metPhi_),
			   sigEvent.jet1_.Pz() + sigEvent.jet2_.Pz(),
			   sigEvent.jet1_.P()  + sigEvent.jet2_.P() + newMet,sigEvent.jet1_.Pt()  + sigEvent.jet2_.Pt() + newMet};
	if     (which3rdl == 1) {pQQLN[0]+=sigEvent.lep1_.Px();pQQLN[1]+=sigEvent.lep1_.Py();pQQLN[2]+=sigEvent.lep1_.Pz();pQQLN[3]+=sigEvent.lep1_.P();pQQLN[4]+=sigEvent.lep1_.Pt();}
	else if(which3rdl == 2) {pQQLN[0]+=sigEvent.lep2_.Px();pQQLN[1]+=sigEvent.lep2_.Py();pQQLN[2]+=sigEvent.lep2_.Pz();pQQLN[3]+=sigEvent.lep2_.P();pQQLN[4]+=sigEvent.lep2_.Pt();}
	else if(which3rdl == 3) {pQQLN[0]+=sigEvent.lep3_.Px();pQQLN[1]+=sigEvent.lep3_.Py();pQQLN[2]+=sigEvent.lep3_.Pz();pQQLN[3]+=sigEvent.lep3_.P();pQQLN[4]+=sigEvent.lep3_.Pt();}
	mQQLN  = pQQLN[3]*pQQLN[3]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1]-pQQLN[2]*pQQLN[2]; if(mQQLN  > 0) mQQLN  = sqrt(mQQLN);  else mQQLN   = 0.0;
	mtQQLN = pQQLN[4]*pQQLN[4]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1];                   if(mtQQLN > 0) mtQQLN = sqrt(mtQQLN); else mtQQLN  = 0.0;
      }
      bdtg = ((dRMin/4.0)-0.5)*2.0; /*bdtg = ((massMin/100.0)-0.5)*2.0; bdtg = Unroll2VarTo1ForWH(dRMin,mTW3,0); bdtg = Unroll2VarTo1ForWH(massMin,mTW3,1);*/ if(bdtg >= 1.0) bdtg = 0.999; if(bdtg <= -1) bdtg = -0.999;

      double deltaPhiQQL = -1;
      if     (which3rdl == 1) deltaPhiQQL = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep1_.Phi());
      else if(which3rdl == 2) deltaPhiQQL = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep2_.Phi());
      else if(which3rdl == 3) deltaPhiQQL = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep3_.Phi());
      bool cleanMode = true;
      if(nsel!=0&&
        (type3l==1||type3l==8||
        (type3l==2&&sigEvent.lq1_*sigEvent.lq2_<0)||(type3l==3&&sigEvent.lq1_*sigEvent.lq3_<0)||(type3l==4&&sigEvent.lq2_*sigEvent.lq3_<0)||
        (type3l==5&&sigEvent.lq2_*sigEvent.lq3_<0)||(type3l==6&&sigEvent.lq1_*sigEvent.lq2_<0)||(type3l==7&&sigEvent.lq1_*sigEvent.lq3_<0))) cleanMode = false;
      bool passMode = cleanMode; if(nsel==2) passMode = !passMode;
      bool passCuts = false;
      bool isSignalDecay = false;
      if(channel == 27){ // WZ selection
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
           sigEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           sigEvent.lep1_.Pt() > 20. &&
           sigEvent.lep2_.Pt() > 10. &&
           sigEvent.lep3_.Pt() > 10. &&
           usedMet > 40. &&
           (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           sigEvent.jet1_.Pt() <= 40. &&
	   massZMin < 15 &&
	   massMin > 12 &&
	   1 == 1
	  ){
	  passCuts = true;
	}
      } // WZ selection
      if(channel == 28){ // 3lnq selection
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
           sigEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           sigEvent.lep1_.Pt() > 20. &&
           sigEvent.lep2_.Pt() > 10. &&
           sigEvent.lep3_.Pt() > 10. &&
           sigEvent.jet1_.Pt() > 30. &&
	   massZMin < 15 &&
	   massMin > 12 &&
  	   usedMet < 50 &&
	   sigEvent.jet1Btag_ < 2.1 &&
	   1 == 1
	  ){
	  passCuts = true; isSignalDecay = true;
	}
      } // 3lnq selection
      if(channel == 1100){ // HW->3l selection
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
           sigEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           sigEvent.lep1_.Pt() > 20. &&
           sigEvent.lep2_.Pt() > 10. &&
           sigEvent.lep3_.Pt() > 10. &&
           usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
           (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           sigEvent.jet1_.Pt() <= 40. &&
	   massZMin > 25 &&
	   massMin > 12 && massMin < 100 &&
	   dRMin < 2.0 &&
	   passMode == true &&
	   1 == 1
	  ){
	  passCuts = true;
          isSignalDecay = true;
	}
      } // HW->3l selection
      if(channel == 1101){ // HW->3l selection-BDT
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
           sigEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           sigEvent.lep1_.Pt() > 20. &&
           sigEvent.lep2_.Pt() > 10. &&
           sigEvent.lep3_.Pt() > 10. &&
           usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
           (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           sigEvent.jet1_.Pt() <= 40. &&
	   massZMin > 15 &&
	   massMin > 12 && massMin < 100 &&
	   passMode == true &&
	   1 == 1
	  ){
	  passCuts = true;
          isSignalDecay = true;
	}
      } // HW->3l selection-BDT
      if(channel == 1200){ // WG* selection
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
         sigEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 10. &&
         sigEvent.lep3_.Pt() > 3. &&
	 mTWMin > 25. && mTFromW > 45. &&
	(type3l == 1 || type3l == 2 || type3l == 3 || type3l == 4) &&
	 massminSFOS > 0 && massminSFOS < 12 && TMath::Abs(massminSFOS-3.1) > 0.1 &&
	(sigEvent.jetLowBtag_ < 2.1 && sigEvent.jet1Btag_ < 2.1 && sigEvent.jet2Btag_ < 2.1 && sigEvent.jet3Btag_ < 2.1) &&
	 sigEvent.njets_ <= 3 &&
	 1 == 1
	  ){
	  passCuts = true;
          isSignalDecay = true;
	}
      } // WG* selection
      if(channel == 1102){ // HZ->3lqq selection
	charge = (int)(sigEvent.lq1_ + sigEvent.lq2_ + sigEvent.lq3_);
	if(
           sigEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           sigEvent.lep1_.Pt() > 20. &&
           sigEvent.lep2_.Pt() > 10. &&
           sigEvent.lep3_.Pt() > 10. &&
	   (sigEvent.jet1_+sigEvent.jet2_).M() > 30. && (sigEvent.jet1_+sigEvent.jet2_).M() < 200. &&
           //(sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           sigEvent.njets_ >= 2 &&
	   massZMin < 15 &&
	   //sigEvent.met_ > 25 &&
	   //mTW > 30 &&
	   massMin > 12 &&
	   //mtQQLN > mhAna-70 && mtQQLN < mhAna+30 &&
           deltaPhiQQL*180.0/TMath::Pi() < 150.0 &&
	   1 == 1
	  ){
	  passCuts = true;
          isSignalDecay = true;
	}
      } // HZ->3lqq selection

      //if(passCuts == true && sigEvent.processId_==10001){
      if(passCuts == true){
	double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
     	 							  fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        							 TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
	double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
     	 							  fabs(sigEvent.lep3_.Eta()), sigEvent.lep3_.Pt(), 
        							 TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid3_));
	double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep3_.Eta()), sigEvent.lep3_.Pt() , 
     	 							  fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        							 TMath::Abs( sigEvent.lid3_), TMath::Abs(sigEvent.lid2_));
	double trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));

	double add = 1.;
	add = add*nPUScaleFactor2012(fhDPU,(float)sigEvent.npu_);
	if((sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
          add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
	if((sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
          add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
	if((sigEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)
          add = add*leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_);
	// CAREFUL HERE, no data-driven corrections, just Higgs k-factors
	// add = 1.0;
        if(TMath::Abs((sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_)-add*trigEff)/add/trigEff>0.0001) {
        printf("PROBLEM: %f %f = %f - %f %f %f %f = %f\n",add,trigEff,add*trigEff,sigEvent.sfWeightFR_,sigEvent.sfWeightPU_,sigEvent.sfWeightEff_,sigEvent.sfWeightTrig_,sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_);
        }
	if(isFermioPhobic == true) {
          add = add * enhancementFactor(mhAna,2); // FF BR(H->WW) enhancement factor
	  if(sigEvent.processId_==121 || sigEvent.processId_==122) add = 0.0;
	}

	if(isSM4 == true) add = add * enhancementFactor(mhAna,1); // SM4 BR(H->WW) enhancement factor

	if (sigEvent.processId_ == 10010) {
          //add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(sigEvent.higgsPt_));
	  if(isFermioPhobic == true) add = 0.0;
	  if(isSM4 == true) add = add * enhancementFactor(mhAna,0); // SM4 ggH enhancement factor
	}

	if(channel == 1200 && sigEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(sigEvent.type_,sigEvent.met_);

	double theWeight = sigEvent.scale1fb_*lumi*add*trigEff;

	if(useTemplates == true){
          bool passJetCut[3] = {sigEvent.jet1_.Pt() <= 40., sigEvent.jet1_.Pt()*1.05 <= 40., sigEvent.jet1_.Pt()*0.95 <= 40.};
          double addLepEff   = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_)*
                               leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_)*
			       leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_);
          double addLepEffUp = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,+1)*
                               leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,+1)*
			       leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_,+1);
          double addLepEffDown = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,-1)*
                        	 leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,-1)*
				 leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_,-1);
	  histo_WH_hww_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	  histo_WH_hww_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	  histo_WH_hww_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_hww_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	  histo_WH_hww_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[1] == true) histo_WH_hww_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
          if(passJetCut[2] == true) histo_WH_hww_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
	}
	int nSigBin = -1;
	// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
	if     (sigEvent.processId_==121 ||
    		sigEvent.processId_==122)   nSigBin = 1;
	else if(sigEvent.processId_==24)    nSigBin = 2;
	else if(sigEvent.processId_==26)    nSigBin = 3;
	else if(sigEvent.processId_==10001) nSigBin = 4;
	else if(sigEvent.processId_==10010) nSigBin = 5;
	else  {if(channel != 1200) assert(0);}
	nSigCut[0]  = nSigCut[0]  + theWeight;
	nSigCutE[0] = nSigCutE[0] + theWeight*theWeight;
	nSigCut[nSigBin]  = nSigCut[nSigBin]  + theWeight;
	nSigCutE[nSigBin] = nSigCutE[nSigBin] + theWeight*theWeight;

	double myVar = newMet;
	if     (thePlot == 1) myVar = TMath::Min((double)sigEvent.lep1_.Pt(),199.999);
	else if(thePlot == 2) myVar = TMath::Min((double)sigEvent.lep2_.Pt(),199.999);
	else if(thePlot == 3) myVar = TMath::Min((double)sigEvent.lep3_.Pt(),199.999);
	else if(thePlot == 4) myVar = sigEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = sigEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = sigEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = sigEvent.dilep_.M();
	else if(thePlot == 8) myVar = bdtg;
	else if(thePlot == 9) myVar = TMath::Min((double)sigEvent.mt1_,199.999);
	else if(thePlot ==10) myVar = TMath::Min((double)sigEvent.mt2_,199.999);
	else if(thePlot ==11) myVar = TMath::Min((double)sigEvent.mt3_,199.999);
	else if(thePlot ==12) myVar = usedMet;
	else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = newMet/sigEvent.dilep_.Pt()/2.0;
	else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = sigEvent.njets_;
	else if(thePlot ==18) myVar = sigEvent.nvtx_;
	else if(thePlot ==19) myVar = (sigEvent.lep1_+sigEvent.lep2_+sigEvent.lep3_).M();
	else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = deltaPhiQQL*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = fabs(sigEvent.dilep_.Eta());
	else if(thePlot ==24) myVar = fabs(sigEvent.lep1_.Eta());
	else if(thePlot ==25) myVar = fabs(sigEvent.lep2_.Eta());
	else if(thePlot ==26) myVar = fabs(sigEvent.lep3_.Eta());
	else if(thePlot ==27) myVar = fabs(sigEvent.jet1_.Eta());
	else if(thePlot ==28) myVar = fabs(sigEvent.jet2_.Eta());
	else if(thePlot ==29) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = mQQLN;
	else if(thePlot ==31) myVar = mtQQLN;
	else if(thePlot ==32) myVar = massminSFOS;
	else if(thePlot ==33) myVar = massminSFSS;
	else if(thePlot ==34) myVar = (sigEvent.jet1_+sigEvent.jet2_).M();
	else if(thePlot ==35) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
	else if(thePlot ==36) myVar = (sigEvent.lep1_+sigEvent.lep2_+sigEvent.lep3_+sigEvent.jet1_).M();
	else if(thePlot ==37) myVar = massZMin;
	else if(thePlot ==38) myVar = massMin;
	else if(thePlot ==39) myVar = TMath::Min(dRMin,3.999);
	else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.metPhi_)*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,sigEvent.metPhi_)*180.0/TMath::Pi();
	else if(thePlot ==42) myVar = mTWMax;
	else if(thePlot ==43) myVar = mTWMin;
	else if(thePlot ==44) myVar = mTW3;
        else if(thePlot ==45) myVar = Mjj;
	else if(thePlot ==46) myVar = ptlw;
	else if(thePlot ==47) myVar = mTFromW;
	else if(thePlot ==48) myVar = type3l;
	else if(thePlot ==49) myVar = newMet*cos(sigEvent.metPhi_);
	else if(thePlot ==50) myVar = newMet*sin(sigEvent.metPhi_);
	else if(thePlot ==51) myVar = newTrackMet*cos(sigEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = newTrackMet*sin(sigEvent.trackMetPhi_);
	else if(thePlot ==53) {if(sigEvent.jet1_.Pt()>15&&sigEvent.jet2_.Pt()>15)myVar = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180.0/TMath::Pi();else myVar=0.0;}
	else if(thePlot ==54) myVar = mTW;
        else if(thePlot>=61&&thePlot<=66) myVar = new3LVar;
	histos->Fill(myVar,theWeight);
	histoWW->Fill(myVar,theWeight);

	myVar = myVar / xmaxPlot;

	if     (myVar <= 0) myVar = 0.001;
	else if(myVar >= 1) myVar = 0.999;
	for(int n1=0; n1<nBin; n1++){
          if(myVar > 1.0*n1/nBin){
            if(isSignalDecay == true ) S0[n1] = S0[n1] + theWeight;
          }
          if(myVar < 1.0*n1/nBin){
            if(isSignalDecay == true ) S1[n1] = S1[n1] + theWeight;
          }
	}
	if(isSignalDecay == true ) hDSigOpt->Fill(myVar,theWeight);
      }
    } // Loop over signal
  }


  if(signalInputFile2 != "" && isFermioPhobic == false && channel != 1102){
    int nSig=ttEvent.tree_->GetEntries();
    if(channel==1101) {
      //ttEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l"	  ,(int)120), &bdtg);
      //ttEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux0"  ,(int)120), &bdtg_aux0);
      //ttEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux1"  ,(int)120), &bdtg_aux1);
      //ttEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux2"  ,(int)120), &bdtg_aux2);
    }
    for (int i=0; i<nSig; ++i) {

      if (i%100000 == 0 && verboseLevel > 0)
	printf("--- reading tt-Signal event %5d of %5d\n",i,nSig);
    ttEvent.tree_->GetEntry(i);
    bool lid = (ttEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection &&
	       (ttEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
	       (ttEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection;
    if(!lid) continue;
    int charge = (int)(ttEvent.lq1_ + ttEvent.lq2_);

    double newMet      = ttEvent.met_;
    double newTrackMet = ttEvent.trackMet_;
    double deltaPhiA[3] = {TMath::Abs(ttEvent.lep1_.Phi()-ttEvent.metPhi_),TMath::Abs(ttEvent.lep2_.Phi()-ttEvent.metPhi_),0.0};
    while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
    while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
    deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
    double pmetA = newMet;
    if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

    double deltaPhiB[3] = {TMath::Abs(ttEvent.lep1_.Phi()-ttEvent.trackMetPhi_),TMath::Abs(ttEvent.lep2_.Phi()-ttEvent.trackMetPhi_),0.0};
    while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
    while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
    deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
    double pmetB = newTrackMet;
    if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

    double deltaPhiC[3] = {TMath::Abs(ttEvent.jet1_.Phi()-ttEvent.metPhi_),TMath::Abs(ttEvent.jet2_.Phi()-ttEvent.metPhi_),0.0};
    while(deltaPhiC[0]>TMath::Pi()) deltaPhiC[0] = TMath::Abs(deltaPhiC[0] - 2*TMath::Pi());
    while(deltaPhiC[1]>TMath::Pi()) deltaPhiC[1] = TMath::Abs(deltaPhiC[1] - 2*TMath::Pi());
    double pmetC = newMet;

    double deltaPhiD[3] = {TMath::Abs(ttEvent.jet1_.Phi()-ttEvent.trackMetPhi_),TMath::Abs(ttEvent.jet2_.Phi()-ttEvent.trackMetPhi_),0.0};
    while(deltaPhiD[0]>TMath::Pi()) deltaPhiD[0] = TMath::Abs(deltaPhiD[0] - 2*TMath::Pi());
    while(deltaPhiD[1]>TMath::Pi()) deltaPhiD[1] = TMath::Abs(deltaPhiD[1] - 2*TMath::Pi());
    double pmetD = newTrackMet;

    if     (ttEvent.njets_ == 0){
      pmetC = pmetA;
      pmetD = pmetB;
    }
    else if(ttEvent.njets_ == 1){
      deltaPhiC[2] = deltaPhiC[0];
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = deltaPhiD[0];
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }
    else if(ttEvent.njets_ >= 2){
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiC[1]);
      deltaPhiC[2] = TMath::Min(deltaPhiC[0],deltaPhiA[2]);
      if(deltaPhiC[2]<TMath::Pi()/2) pmetC = pmetC * sin(deltaPhiC[2]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiD[1]);
      deltaPhiD[2] = TMath::Min(deltaPhiD[0],deltaPhiA[2]);
      if(deltaPhiD[2]<TMath::Pi()/2) pmetD = pmetD * sin(deltaPhiD[2]);
    }

    double usedMet = TMath::Min(ttEvent.pmet_,ttEvent.pTrackMet_);
    Float_t mTWMin   = ttEvent.mt1_;
    Float_t mTWMax   = ttEvent.mt1_;
    if(mTWMin > ttEvent.mt2_			     ) mTWMin = ttEvent.mt2_;
    if(mTWMin > ttEvent.mt3_ && ttEvent.lid3_ != 0.) mTWMin = ttEvent.mt3_;
    if(mTWMax < ttEvent.mt2_			     ) mTWMax = ttEvent.mt2_;
    if(mTWMax < ttEvent.mt3_ && ttEvent.lid3_ != 0.) mTWMax = ttEvent.mt3_;

    double Mjj = (ttEvent.jet1_+ttEvent.jet2_).M();
    if(ttEvent.njets_ >= 3){
      if(TMath::Abs((ttEvent.jet1_+ttEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (ttEvent.jet1_+ttEvent.jet3_).M();
      if(TMath::Abs((ttEvent.jet2_+ttEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (ttEvent.jet2_+ttEvent.jet3_).M();
    }
    if(ttEvent.njets_ >= 4){
      if(TMath::Abs((ttEvent.jet1_+ttEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (ttEvent.jet1_+ttEvent.jet4_).M();
      if(TMath::Abs((ttEvent.jet2_+ttEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (ttEvent.jet2_+ttEvent.jet4_).M();
      if(TMath::Abs((ttEvent.jet3_+ttEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (ttEvent.jet3_+ttEvent.jet4_).M();
    }

    double massZMin = 999.0; double massMin = 999.0; double massminSFOS = 999.0; double massminSFSS = 999.0; double mtQQLN = -1.0; double mQQLN = -1.0;
    double dRMin = 999.0; double mTW3 = 0.0;double mTW = 0.0; double type3l = 0.0; double mTFromW = 0.0; double ptlw = 0.0; double which3rdl = -1; double new3LVar = 0.0;
    if(ttEvent.lid3_ != 0){
      usedMet = TMath::Min(newMet,newTrackMet);
      massZMin = trilepton_info(0,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      massMin  = trilepton_info(1,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      dRMin = trilepton_info(2,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      mTW   = trilepton_info(3,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      mTW3  = trilepton_info(10,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      type3l= trilepton_info(4,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      massminSFOS  = trilepton_info(5,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      massminSFSS  = trilepton_info(6,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      mTFromW     = trilepton_info(7,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      ptlw        = trilepton_info(8,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      which3rdl   = trilepton_info(9,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_);
      new3LVar   = trilepton_info(thePlot,ttEvent.lep1_,ttEvent.lep2_,ttEvent.lep3_,
                                  ttEvent.lq1_ ,ttEvent.lq2_ ,ttEvent.lq3_,
		                  ttEvent.lid1_,ttEvent.lid2_,ttEvent.lid3_,
				  ttEvent.mt1_ ,ttEvent.mt2_ ,ttEvent.mt3_,ttEvent.met_, ttEvent.metPhi_);
      double pQQLN[5] = {ttEvent.jet1_.Px() + ttEvent.jet2_.Px() + newMet*cos( ttEvent.metPhi_),
                         ttEvent.jet1_.Py() + ttEvent.jet2_.Py() + newMet*sin( ttEvent.metPhi_),
			 ttEvent.jet1_.Pz() + ttEvent.jet2_.Pz(),
			 ttEvent.jet1_.P()  + ttEvent.jet2_.P() + newMet,ttEvent.jet1_.Pt()  + ttEvent.jet2_.Pt() + newMet};
      if     (which3rdl == 1) {pQQLN[0]+=ttEvent.lep1_.Px();pQQLN[1]+=ttEvent.lep1_.Py();pQQLN[2]+=ttEvent.lep1_.Pz();pQQLN[3]+=ttEvent.lep1_.P();pQQLN[4]+=ttEvent.lep1_.Pt();}
      else if(which3rdl == 2) {pQQLN[0]+=ttEvent.lep2_.Px();pQQLN[1]+=ttEvent.lep2_.Py();pQQLN[2]+=ttEvent.lep2_.Pz();pQQLN[3]+=ttEvent.lep2_.P();pQQLN[4]+=ttEvent.lep2_.Pt();}
      else if(which3rdl == 3) {pQQLN[0]+=ttEvent.lep3_.Px();pQQLN[1]+=ttEvent.lep3_.Py();pQQLN[2]+=ttEvent.lep3_.Pz();pQQLN[3]+=ttEvent.lep3_.P();pQQLN[4]+=ttEvent.lep3_.Pt();}
      mQQLN  = pQQLN[3]*pQQLN[3]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1]-pQQLN[2]*pQQLN[2]; if(mQQLN  > 0) mQQLN  = sqrt(mQQLN);  else mQQLN   = 0.0;
      mtQQLN = pQQLN[4]*pQQLN[4]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1];                   if(mtQQLN > 0) mtQQLN = sqrt(mtQQLN); else mtQQLN  = 0.0;
    }
    bdtg = ((dRMin/4.0)-0.5)*2.0; /*bdtg = ((massMin/100.0)-0.5)*2.0; bdtg = Unroll2VarTo1ForWH(dRMin,mTW3,0); bdtg = Unroll2VarTo1ForWH(massMin,mTW3,1);*/ if(bdtg >= 1.0) bdtg = 0.999; if(bdtg <= -1) bdtg = -0.999;

    double deltaPhiQQL = -1;
    if     (which3rdl == 1) deltaPhiQQL = DeltaPhi((ttEvent.jet1_+ttEvent.jet2_).Phi(),ttEvent.lep1_.Phi());
    else if(which3rdl == 2) deltaPhiQQL = DeltaPhi((ttEvent.jet1_+ttEvent.jet2_).Phi(),ttEvent.lep2_.Phi());
    else if(which3rdl == 3) deltaPhiQQL = DeltaPhi((ttEvent.jet1_+ttEvent.jet2_).Phi(),ttEvent.lep3_.Phi());
    bool cleanMode = true;
    if(nsel!=0&&
      (type3l==1||type3l==8||
      (type3l==2&&ttEvent.lq1_*ttEvent.lq2_<0)||(type3l==3&&ttEvent.lq1_*ttEvent.lq3_<0)||(type3l==4&&ttEvent.lq2_*ttEvent.lq3_<0)||
      (type3l==5&&ttEvent.lq2_*ttEvent.lq3_<0)||(type3l==6&&ttEvent.lq1_*ttEvent.lq2_<0)||(type3l==7&&ttEvent.lq1_*ttEvent.lq3_<0))) cleanMode = false;
    bool passMode = cleanMode; if(nsel==2) passMode = !passMode;
    bool passCuts = false;
    bool isSignalDecay = false;
    if(channel == 27){ // WZ selection
      charge = (int)(ttEvent.lq1_ + ttEvent.lq2_ + ttEvent.lq3_);
      if(
         ttEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         ttEvent.lep1_.Pt() > 20. &&
         ttEvent.lep2_.Pt() > 10. &&
         ttEvent.lep3_.Pt() > 10. &&
         usedMet > 40. &&
         (ttEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         ttEvent.jet1_.Pt() <= 40. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WZ selection
    if(channel == 28){ // 3lnq selection
      charge = (int)(ttEvent.lq1_ + ttEvent.lq2_ + ttEvent.lq3_);
      if(
         ttEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         ttEvent.lep1_.Pt() > 20. &&
         ttEvent.lep2_.Pt() > 10. &&
         ttEvent.lep3_.Pt() > 10. &&
         ttEvent.jet1_.Pt() > 30. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 usedMet < 50 &&
         ttEvent.jet1Btag_ < 2.1 &&
	 1 == 1
	){
	passCuts = true; isSignalDecay = true;
      }
    } // 3lnq selection
    if(channel == 1100){ // HW->3l selection
      charge = (int)(ttEvent.lq1_ + ttEvent.lq2_ + ttEvent.lq3_);
      if(
         ttEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         ttEvent.lep1_.Pt() > 20. &&
         ttEvent.lep2_.Pt() > 10. &&
         ttEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (ttEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         ttEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
         dRMin < 2.0 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // HW->3l selection
    if(channel == 1101){ // HW->3l selection-BDT
      charge = (int)(ttEvent.lq1_ + ttEvent.lq2_ + ttEvent.lq3_);
      if(
         ttEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         ttEvent.lep1_.Pt() > 20. &&
         ttEvent.lep2_.Pt() > 10. &&
         ttEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (ttEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         ttEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
        isSignalDecay = true;
      }
    } // HW->3l selection-BDT
    //if(passCuts == true && ttEvent.processId_==10001){
    if(passCuts == true){
      double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(ttEvent.lep1_.Eta()), ttEvent.lep1_.Pt() , 
     	 							fabs(ttEvent.lep2_.Eta()), ttEvent.lep2_.Pt(), 
        						       TMath::Abs( ttEvent.lid1_), TMath::Abs(ttEvent.lid2_));
      double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(ttEvent.lep1_.Eta()), ttEvent.lep1_.Pt() , 
     	 							fabs(ttEvent.lep3_.Eta()), ttEvent.lep3_.Pt(), 
        						       TMath::Abs( ttEvent.lid1_), TMath::Abs(ttEvent.lid3_));
      double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(ttEvent.lep3_.Eta()), ttEvent.lep3_.Pt() , 
     	 							fabs(ttEvent.lep2_.Eta()), ttEvent.lep2_.Pt(), 
        						       TMath::Abs( ttEvent.lid3_), TMath::Abs(ttEvent.lid2_));
      double trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
      double add = 1.;
      /*
      add = add*nPUScaleFactor2012(fhDPU,(float)ttEvent.npu_);
      if((ttEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
        add = add*leptonEfficiency(ttEvent.lep1_.Pt(), ttEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid1_);
      if((ttEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)
        add = add*leptonEfficiency(ttEvent.lep2_.Pt(), ttEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid2_);
      if((ttEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)
        add = add*leptonEfficiency(ttEvent.lep3_.Pt(), ttEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid3_);
      */
      add = ttEvent.sfWeightPU_*ttEvent.sfWeightEff_;

      if(isSM4 == true) add = add * enhancementFactor(mhAna,3); // SM4 BR(H->tautau) enhancement factor

      double theWeight = ttEvent.scale1fb_*lumi*add*trigEff;

      if(useTemplates == true){
        bool passJetCut[3] = {ttEvent.jet1_.Pt() <= 40., ttEvent.jet1_.Pt()*1.05 <= 40., ttEvent.jet1_.Pt()*0.95 <= 40.};
        double addLepEff   = leptonEfficiency(ttEvent.lep1_.Pt(), ttEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid1_)*
                             leptonEfficiency(ttEvent.lep2_.Pt(), ttEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid2_)*
			     leptonEfficiency(ttEvent.lep3_.Pt(), ttEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid3_);
        double addLepEffUp = leptonEfficiency(ttEvent.lep1_.Pt(), ttEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid1_,+1)*
                             leptonEfficiency(ttEvent.lep2_.Pt(), ttEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid2_,+1)*
			     leptonEfficiency(ttEvent.lep3_.Pt(), ttEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid3_,+1);
        double addLepEffDown = leptonEfficiency(ttEvent.lep1_.Pt(), ttEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid1_,-1)*
                               leptonEfficiency(ttEvent.lep2_.Pt(), ttEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid2_,-1)*
			       leptonEfficiency(ttEvent.lep3_.Pt(), ttEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, ttEvent.lid3_,-1);
	histo_WH_htt_CMS_MVALepEffBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffUp  /addLepEff);
	histo_WH_htt_CMS_MVALepEffBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight*addLepEffDown/addLepEff);
	histo_WH_htt_CMS_MVALepResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux0,maxHis-0.001),minHis+0.001), theWeight);
	histo_WH_htt_CMS_MVALepResBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg_aux1,maxHis-0.001),minHis+0.001), theWeight);
	histo_WH_htt_CMS_MVAMETResBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg_aux2,maxHis-0.001),minHis+0.001), theWeight);
        if(passJetCut[1] == true) histo_WH_htt_CMS_MVAJESBoundingUp  ->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
        if(passJetCut[2] == true) histo_WH_htt_CMS_MVAJESBoundingDown->Fill(TMath::Max(TMath::Min((double)bdtg,maxHis-0.001),minHis+0.001), theWeight);
      }
      int nSigBin = 6;
      nSigCut[nSigBin]  = nSigCut[nSigBin]  + theWeight;
      nSigCutE[nSigBin] = nSigCutE[nSigBin] + theWeight*theWeight;

      double myVar = newMet;
      if     (thePlot == 1) myVar = TMath::Min((double)ttEvent.lep1_.Pt(),199.999);
      else if(thePlot == 2) myVar = TMath::Min((double)ttEvent.lep2_.Pt(),199.999);
      else if(thePlot == 3) myVar = TMath::Min((double)ttEvent.lep3_.Pt(),199.999);
      else if(thePlot == 4) myVar = ttEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = ttEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = ttEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = ttEvent.dilep_.M();
      else if(thePlot == 8) myVar = bdtg;
      else if(thePlot == 9) myVar = TMath::Min((double)ttEvent.mt1_,199.999);
      else if(thePlot ==10) myVar = TMath::Min((double)ttEvent.mt2_,199.999);
      else if(thePlot ==11) myVar = TMath::Min((double)ttEvent.mt3_,199.999);
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = ttEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(ttEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = newMet/ttEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = ttEvent.lep2_.Pt()/ttEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = ttEvent.njets_;
      else if(thePlot ==18) myVar = ttEvent.nvtx_;
      else if(thePlot ==19) myVar = (ttEvent.lep1_+ttEvent.lep2_+ttEvent.lep3_).M();
      else if(thePlot ==20) myVar = ttEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(ttEvent.dPhiLep1MET_,ttEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(ttEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(ttEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(ttEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(ttEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(ttEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(ttEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(ttEvent.jet1_.Eta()),fabs(ttEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = mQQLN;
      else if(thePlot ==31) myVar = mtQQLN;
      else if(thePlot ==32) myVar = massminSFOS;
      else if(thePlot ==33) myVar = massminSFSS;
      else if(thePlot ==34) myVar = (ttEvent.jet1_+ttEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(ttEvent.jet1_.Eta()-ttEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = (ttEvent.lep1_+ttEvent.lep2_+ttEvent.lep3_+ttEvent.jet1_).M();
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = TMath::Min(dRMin,3.999);
      else if(thePlot ==40) myVar = DeltaPhi(ttEvent.jet1_.Phi() ,ttEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(ttEvent.trackMetPhi_,ttEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = Mjj;
      else if(thePlot ==46) myVar = ptlw;
      else if(thePlot ==47) myVar = mTFromW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = newMet*cos(ttEvent.metPhi_);
      else if(thePlot ==50) myVar = newMet*sin(ttEvent.metPhi_);
      else if(thePlot ==51) myVar = newTrackMet*cos(ttEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = newTrackMet*sin(ttEvent.trackMetPhi_);
      else if(thePlot ==53) {if(ttEvent.jet1_.Pt()>15&&ttEvent.jet2_.Pt()>15)myVar = DeltaPhi((ttEvent.jet1_+ttEvent.jet2_).Phi(),ttEvent.dilep_.Phi())*180.0/TMath::Pi();else myVar=0.0;}
      else if(thePlot ==54) myVar = mTW;
      else if(thePlot>=61&&thePlot<=66) myVar = new3LVar;
      histos->Fill(myVar,theWeight);
      histoTT->Fill(myVar,theWeight);

      myVar = myVar / xmaxPlot;

      if     (myVar <= 0) myVar = 0.001;
      else if(myVar >= 1) myVar = 0.999;
      for(int n1=0; n1<nBin; n1++){
        if(myVar > 1.0*n1/nBin){
          if(isSignalDecay == true ) S0[n1] = S0[n1] + theWeight;
        }
        if(myVar < 1.0*n1/nBin){
          if(isSignalDecay == true ) S1[n1] = S1[n1] + theWeight;
        }
      }
      if(isSignalDecay == true ) hDSigOpt->Fill(myVar,theWeight);
    }
    } // Loop over h->tt signal
  }

  int nData=dataEvent.tree_->GetEntries();

  if(channel==1101) {
    //dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l"	    ,(int)120), &bdtg);
    //dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux0"  ,(int)120), &bdtg_aux0);
    //dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux1"  ,(int)120), &bdtg_aux1);
    //dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_wh3l_aux2"  ,(int)120), &bdtg_aux2);
  }
  int nSelectedData = 0;
  for (int i=0; i<nData; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);

    int FullLid  = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection;
	FullLid += (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
	FullLid += (dataEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection;
    int Photonid  = (dataEvent.cuts_ & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2;
	Photonid += (dataEvent.cuts_ & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2;
	Photonid += (dataEvent.cuts_ & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2;
    
    if(channel != 1300 && FullLid != 3) continue;
    if(channel == 1300 && (FullLid != 2 || Photonid != 1)) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;
    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

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

    double usedMet = TMath::Min(pmetA,pmetB);
    Float_t mTWMin   = dataEvent.mt1_;
    Float_t mTWMax   = dataEvent.mt1_;
    if(mTWMin > dataEvent.mt2_  		       ) mTWMin = dataEvent.mt2_;
    if(mTWMin > dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMin = dataEvent.mt3_;
    if(mTWMax < dataEvent.mt2_  		       ) mTWMax = dataEvent.mt2_;
    if(mTWMax < dataEvent.mt3_ && dataEvent.lid3_ != 0.) mTWMax = dataEvent.mt3_;

    double Mjj = (dataEvent.jet1_+dataEvent.jet2_).M();
    if(dataEvent.njets_ >= 3){
      if(TMath::Abs((dataEvent.jet1_+dataEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (dataEvent.jet1_+dataEvent.jet3_).M();
      if(TMath::Abs((dataEvent.jet2_+dataEvent.jet3_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (dataEvent.jet2_+dataEvent.jet3_).M();
    }
    if(dataEvent.njets_ >= 4){
      if(TMath::Abs((dataEvent.jet1_+dataEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (dataEvent.jet1_+dataEvent.jet4_).M();
      if(TMath::Abs((dataEvent.jet2_+dataEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (dataEvent.jet2_+dataEvent.jet4_).M();
      if(TMath::Abs((dataEvent.jet3_+dataEvent.jet4_).M()-85) < TMath::Abs(Mjj-85)) Mjj = (dataEvent.jet3_+dataEvent.jet4_).M();
    }   

    double massZMin = 999.0; double massMin = 999.0; double massminSFOS = 999.0; double massminSFSS = 999.0; double mtQQLN = -1.0; double mQQLN = -1.0;
    double dRMin = 999.0; double mTW3 = 0.0;double mTW = 0.0; double type3l = 0.0; double mTFromW = 0.0; double ptlw = 0.0; double which3rdl = -1; double new3LVar = 0.0;
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
      mTW   = trilepton_info(3,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      mTW3  = trilepton_info(10,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      type3l= trilepton_info(4,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      massminSFOS  = trilepton_info(5,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      massminSFSS  = trilepton_info(6,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      mTFromW     = trilepton_info(7,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      ptlw        = trilepton_info(8,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      which3rdl   = trilepton_info(9,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_);
      new3LVar   = trilepton_info(thePlot,dataEvent.lep1_,dataEvent.lep2_,dataEvent.lep3_,
                                  dataEvent.lq1_ ,dataEvent.lq2_ ,dataEvent.lq3_,
		                  dataEvent.lid1_,dataEvent.lid2_,dataEvent.lid3_,
				  dataEvent.mt1_ ,dataEvent.mt2_ ,dataEvent.mt3_,dataEvent.met_, dataEvent.metPhi_);
      double pQQLN[5] = {dataEvent.jet1_.Px() + dataEvent.jet2_.Px() + dataEvent.met_*cos( dataEvent.metPhi_),
                         dataEvent.jet1_.Py() + dataEvent.jet2_.Py() + dataEvent.met_*sin( dataEvent.metPhi_),
			 dataEvent.jet1_.Pz() + dataEvent.jet2_.Pz(),
			 dataEvent.jet1_.P()  + dataEvent.jet2_.P() + dataEvent.met_,dataEvent.jet1_.Pt()  + dataEvent.jet2_.Pt() + dataEvent.met_};
      if     (which3rdl == 1) {pQQLN[0]+=dataEvent.lep1_.Px();pQQLN[1]+=dataEvent.lep1_.Py();pQQLN[2]+=dataEvent.lep1_.Pz();pQQLN[3]+=dataEvent.lep1_.P();pQQLN[4]+=dataEvent.lep1_.Pt();}
      else if(which3rdl == 2) {pQQLN[0]+=dataEvent.lep2_.Px();pQQLN[1]+=dataEvent.lep2_.Py();pQQLN[2]+=dataEvent.lep2_.Pz();pQQLN[3]+=dataEvent.lep2_.P();pQQLN[4]+=dataEvent.lep2_.Pt();}
      else if(which3rdl == 3) {pQQLN[0]+=dataEvent.lep3_.Px();pQQLN[1]+=dataEvent.lep3_.Py();pQQLN[2]+=dataEvent.lep3_.Pz();pQQLN[3]+=dataEvent.lep3_.P();pQQLN[4]+=dataEvent.lep3_.Pt();}
      mQQLN  = pQQLN[3]*pQQLN[3]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1]-pQQLN[2]*pQQLN[2]; if(mQQLN  > 0) mQQLN  = sqrt(mQQLN);  else mQQLN   = 0.0;
      mtQQLN = pQQLN[4]*pQQLN[4]-pQQLN[0]*pQQLN[0]-pQQLN[1]*pQQLN[1];                   if(mtQQLN > 0) mtQQLN = sqrt(mtQQLN); else mtQQLN  = 0.0;
    }
    bdtg = ((dRMin/4.0)-0.5)*2.0; /*bdtg = ((massMin/100.0)-0.5)*2.0; bdtg = Unroll2VarTo1ForWH(dRMin,mTW3,0); bdtg = Unroll2VarTo1ForWH(massMin,mTW3,1);*/ if(bdtg >= 1.0) bdtg = 0.999; if(bdtg <= -1) bdtg = -0.999;

    double deltaPhiQQL = -1;
    if     (which3rdl == 1) deltaPhiQQL = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep1_.Phi());
    else if(which3rdl == 2) deltaPhiQQL = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep2_.Phi());
    else if(which3rdl == 3) deltaPhiQQL = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep3_.Phi());
    bool cleanMode = true;
    if(nsel!=0&&
      (type3l==1||type3l==8||
      (type3l==2&&dataEvent.lq1_*dataEvent.lq2_<0)||(type3l==3&&dataEvent.lq1_*dataEvent.lq3_<0)||(type3l==4&&dataEvent.lq2_*dataEvent.lq3_<0)||
      (type3l==5&&dataEvent.lq2_*dataEvent.lq3_<0)||(type3l==6&&dataEvent.lq1_*dataEvent.lq2_<0)||(type3l==7&&dataEvent.lq1_*dataEvent.lq3_<0))) cleanMode = false;
    bool passMode = cleanMode; if(nsel==2) passMode = !passMode;
    bool passCuts = false;
    if(channel == 27){ // WZ selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 40. &&
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.jet1_.Pt() <= 40. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WZ selection
    if(channel == 28){ // 3lnq selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         dataEvent.jet1_.Pt() > 30. &&
	 massZMin < 15 &&
	 massMin > 12 &&
	 usedMet < 50 &&
	 dataEvent.jet1Btag_ < 2.1 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // 3lnq selection
    if(channel == 1098){ // top region
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 30 &&
         (dataEvent.jet1_.Pt() > 40. || dataEvent.jet1Btag_ > 2.1) &&
         massZMin > 15 &&
         massMin > 12 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // top region
    if(channel == 1099){ // Zjets region
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20. &&
         dataEvent.lep3_.Pt() < 15. &&
         usedMet < 25 &&
         dataEvent.jet1_.Pt() < 40. &&
         //massZMin < 15 &&
	 abs(dataEvent.dilep_.M()-91.1876) < 15. &&
         massMin > 12 &&
	 //(type3l == 1 || type3l == 3) &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
	//if(type3l ==1 && usedMet > 20   && usedMet < 24.0)type3l = 2;
	//if(type3l ==1 && usedMet > 19   && usedMet < 20.5)type3l = 3;
      }
    } // Zjets region
    if(channel == 1097 || channel == 1096){ // ttX selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      int nBtag0 = (dataEvent.jet1Btag_ > 2.1) + (dataEvent.jet2Btag_ > 2.1);
      int nBtag1 = (dataEvent.jet1Btag_ > 2.1) + (dataEvent.jet2Btag_ > 2.1) + (dataEvent.jet3Btag_ > 2.1) + (dataEvent.jet4Btag_ > 2.1);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20. &&
         dataEvent.lep3_.Pt() > 20. &&
         usedMet > 20. &&
       ((channel == 1097 && dataEvent.njets_ >= 2 && dataEvent.njets_ <= 3 &&
         nBtag0 == 2 && massZMin > 10) ||
	(channel == 1096 && dataEvent.njets_ >= 4 && massZMin < 15 &&
         nBtag1 >= 2)) &&
         massMin > 12 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // ttX selection
    if(channel == 1100){ // HW->3l selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
         dRMin < 2.0 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection
    if(channel == 1101){ // HW->3l selection-BDT
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
         usedMet > 30. && (cleanMode == true || usedMet > 40.) &&
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
         dataEvent.jet1_.Pt() <= 40. &&
         massZMin > 25 &&
         massMin > 12 && massMin < 100 &&
	 passMode == true &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HW->3l selection-BDT
    if(channel == 1300 || channel == 1301){ // ZG selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 10. &&
	 massMin > 10 &&
	 massZMin > 10 &&
	 TMath::Abs((dataEvent.lep1_+dataEvent.lep2_+dataEvent.lep3_).M()-91.1976) < 15.0 &&
	 (type3l==2||type3l==3||type3l==4||type3l==8) &&
         //(dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 1 == 1
	){
	if(type3l==1||type3l==8||
	  (type3l==2&&dataEvent.lq1_*dataEvent.lq2_<0)||(type3l==3&&dataEvent.lq1_*dataEvent.lq3_<0)||(type3l==4&&dataEvent.lq2_*dataEvent.lq3_<0)||
	  (type3l==5&&dataEvent.lq2_*dataEvent.lq3_<0)||(type3l==6&&dataEvent.lq1_*dataEvent.lq2_<0)||(type3l==7&&dataEvent.lq1_*dataEvent.lq3_<0)){
	passCuts = true;
	}
      }
    } // ZG selection
    if(channel == 1200){ // WG* selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
         dataEvent.lid3_ != 0 &&
         abs(charge) == 1 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 10. &&
         dataEvent.lep3_.Pt() > 3. &&
	 mTWMin > 25. && mTFromW > 45. &&
	(type3l == 1 || type3l == 2 || type3l == 3 || type3l == 4) &&
	 massminSFOS > 0 && massminSFOS < 12 && TMath::Abs(massminSFOS-3.1) > 0.1 &&
	(dataEvent.jetLowBtag_ < 2.1 && dataEvent.jet1Btag_ < 2.1 && dataEvent.jet2Btag_ < 2.1 && dataEvent.jet3Btag_ < 2.1) &&
	 dataEvent.njets_ <= 3 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // WG* selection
    if(channel == 1102){ // HZ->3lqq selection
      charge = (int)(dataEvent.lq1_ + dataEvent.lq2_ + dataEvent.lq3_);
      if(
           dataEvent.lid3_ != 0 &&
           abs(charge) == 1 &&
           dataEvent.lep1_.Pt() > 20. &&
           dataEvent.lep2_.Pt() > 10. &&
           dataEvent.lep3_.Pt() > 10. &&
	   (dataEvent.jet1_+dataEvent.jet2_).M() > 30. && (dataEvent.jet1_+dataEvent.jet2_).M() < 200. &&
           //(dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
           dataEvent.njets_ >= 2 &&
	   massZMin < 15 &&
	   //dataEvent.met_ > 25 &&
	   //mTW > 30 &&
	   massMin > 12 &&
	   //mtQQLN > mhAna-70 && mtQQLN < mhAna+30 &&
           deltaPhiQQL*180.0/TMath::Pi() < 150.0 &&
	 1 == 1
	){
	passCuts = true;
      }
    } // HZ->3lqq selection
    
    if(passCuts == true){
      double myVar = dataEvent.met_;
      if     (thePlot == 1) myVar = TMath::Min((double)dataEvent.lep1_.Pt(),199.999);
      else if(thePlot == 2) myVar = TMath::Min((double)dataEvent.lep2_.Pt(),199.999);
      else if(thePlot == 3) myVar = TMath::Min((double)dataEvent.lep3_.Pt(),199.999);
      else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = dataEvent.dilep_.M();
      else if(thePlot == 8) myVar = bdtg;
      else if(thePlot == 9) myVar = TMath::Min((double)dataEvent.mt1_,199.999);
      else if(thePlot ==10) myVar = TMath::Min((double)dataEvent.mt2_,199.999);
      else if(thePlot ==11) myVar = TMath::Min((double)dataEvent.mt3_,199.999);
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = dataEvent.met_/dataEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = dataEvent.njets_;
      else if(thePlot ==18) myVar = dataEvent.nvtx_;
      else if(thePlot ==19) myVar = (dataEvent.lep1_+dataEvent.lep2_+dataEvent.lep3_).M();
      else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(dataEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(dataEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(dataEvent.lep2_.Eta());
      else if(thePlot ==26) myVar = fabs(dataEvent.lep3_.Eta());
      else if(thePlot ==27) myVar = fabs(dataEvent.jet1_.Eta());
      else if(thePlot ==28) myVar = fabs(dataEvent.jet2_.Eta());
      else if(thePlot ==29) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = mQQLN;
      else if(thePlot ==31) myVar = mtQQLN;
      else if(thePlot ==32) myVar = massminSFOS;
      else if(thePlot ==33) myVar = massminSFSS;
      else if(thePlot ==34) myVar = (dataEvent.jet1_+dataEvent.jet2_).M();
      else if(thePlot ==35) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
      else if(thePlot ==36) myVar = (dataEvent.lep1_+dataEvent.lep2_+dataEvent.lep3_+dataEvent.jet1_).M();
      else if(thePlot ==37) myVar = massZMin;
      else if(thePlot ==38) myVar = massMin;
      else if(thePlot ==39) myVar = TMath::Min(dRMin,3.999);
      else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==42) myVar = mTWMax;
      else if(thePlot ==43) myVar = mTWMin;
      else if(thePlot ==44) myVar = mTW3;
      else if(thePlot ==45) myVar = Mjj;
      else if(thePlot ==46) myVar = ptlw;
      else if(thePlot ==47) myVar = mTFromW;
      else if(thePlot ==48) myVar = type3l;
      else if(thePlot ==49) myVar = dataEvent.met_*cos(dataEvent.metPhi_);
      else if(thePlot ==50) myVar = dataEvent.met_*sin(dataEvent.metPhi_);
      else if(thePlot ==51) myVar = dataEvent.trackMet_*cos(dataEvent.trackMetPhi_);
      else if(thePlot ==52) myVar = dataEvent.trackMet_*sin(dataEvent.trackMetPhi_);
      else if(thePlot ==53) {if(dataEvent.jet1_.Pt()>15&&dataEvent.jet2_.Pt()>15)myVar = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi();else myVar=0.0;}
      else if(thePlot ==54) myVar = mTW;
      else if(thePlot>=61&&thePlot<=66) myVar = new3LVar;

      nSelectedData = nSelectedData + 1;
      histo5->Fill(myVar,1.0);
    }
  } // End loop data
  printf("data: %d\n",nSelectedData);

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
      histo6->Write();
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

  if(channel != 1200){
    bgdDecay[27] = bgdDecay[27] + bgdDecay[20];
    weiDecay[27] = weiDecay[27] + weiDecay[20];
  }
  //bgdDecay[20] = 0.0;
  //weiDecay[20] = 0.0;
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
    for(int i=0; i<7; i++) if(nSigCut[i] > 0) nSigCutE[i] = sqrt(nSigCutE[i]);
    for(int i=0; i<7; i++) printf("nSig(%1d) = %6.3f +/- %6.3f\n",i,nSigCut[i],nSigCutE[i]);
    bgdCombined[0]  = bgdCombined[0] + bgdDecay[29] + bgdDecay[30] + bgdCombined[1] + bgdCombined[2];
    bgdCombinedE[0] = bgdCombinedE[0]*bgdCombinedE[0] + weiDecay[29] + weiDecay[30] + 
                      bgdCombinedE[1]*bgdCombinedE[1] + bgdCombinedE[2]*bgdCombinedE[1];
    bgdCombinedE[0] = sqrt(bgdCombinedE[0]); 

    double totalB   = bgdDecay[27]+bgdDecay[28]+bgdCombined[0]+bgdCombined[3]+bgdDecay[21];
    double totalB_E = sqrt(weiDecay[27]+weiDecay[28]+bgdCombinedE[0]*bgdCombinedE[0]+bgdCombinedE[3]*bgdCombinedE[3]+weiDecay[21]);
    printf("%7.2f \\pm %5.2f & %7.2f \\pm %5.2f & %4d & %7.2f \\pm %5.2f & %7.2f \\pm %5.2f & %7.2f \\pm %5.2f & %7.2f \\pm %5.2f & %7.2f \\pm %5.2f & %7.2f \\pm %5.2f\\\\\n",
           nSigCut[6],nSigCutE[6],nSigCut[0],nSigCutE[0],nSelectedData,totalB,totalB_E,bgdDecay[27],sqrt(weiDecay[27]),bgdDecay[28],sqrt(weiDecay[28]),
	   bgdCombined[0],bgdCombinedE[0],bgdCombined[3],bgdCombinedE[3],bgdDecay[21],sqrt(weiDecay[21]));

    for(int i=0; i<45; i++) if(bgdDecay[i] < 0) bgdDecay[i] = 0.0;
    if(bgdDecay[27] > 0) weiDecay[27] = sqrt(weiDecay[27])/bgdDecay[27];
    if(bgdDecay[28] > 0) weiDecay[28] = sqrt(weiDecay[28])/bgdDecay[28];
    if(bgdCombined[0] > 0) bgdCombinedE[0] =  bgdCombinedE[0] / bgdCombined[0];
    if(bgdCombined[3] > 0) bgdCombinedE[3] =  bgdCombinedE[3] / bgdCombined[3];
    if(bgdDecay[21] > 0) weiDecay[21] = sqrt(weiDecay[21])/bgdDecay[21];
    if(bgdDecay[41] > 0) weiDecay[41] = sqrt(weiDecay[41])/bgdDecay[41];
    if(bgdDecay[42] > 0) weiDecay[42] = sqrt(weiDecay[42])/bgdDecay[42];
    
    char outputLimitsCut[200];
    double signalError[2] = {1.0, 1.0};
    if(nSigCut[0] > 0) signalError[0] = nSigCutE[0]/nSigCut[0]+1.0;
    if(nSigCut[6] > 0) signalError[1] = nSigCutE[6]/nSigCut[6]+1.0;
    
    double QCDscale_VH = 1.07;
    if(isFermioPhobic == true) QCDscale_VH += 0.05;
    if(channel == 1102) QCDscale_VH  -= 0.05;

    double PDF_VH = 1.04;
    if(isFermioPhobic == true) PDF_VH += 0.00;
    double CMS_met = 1.02; double CMS_jes = 1.02; double CMS_jes_ZZ = 1.02; double CMS_jes_WZ = 1.00;
    if(channel == 1102) {CMS_jes_ZZ = 1.10; CMS_jes_WZ = 1.10;}

    if(channel == 1100 || channel == 1102){
      if     (channel == 1100) sprintf(outputLimitsCut,"histo_limits_vh3l%d_mh%3.0f_cut_8TeV.txt",nsel,mhAna);
      else if(channel == 1102) sprintf(outputLimitsCut,"histo_limits_zh3l2q_mh%3.0f_cut_8TeV.txt",mhAna);
      else assert(0);
      ofstream newcardCut;
      newcardCut.open(outputLimitsCut);
      newcardCut << Form("imax 1 number of channels\n");
      newcardCut << Form("jmax * number of background\n");
      newcardCut << Form("kmax * number of nuisance parameters\n");
      newcardCut << Form("Observation %d\n",nSelectedData);
      newcardCut << Form("bin 1 1 1 1 1 1 1 1 1\n");
      newcardCut << Form("process WH_htt WH_hww WZ ZZ Wjets Wgamma VVV WH_htt_SM WH_hww_SM\n");
      newcardCut << Form("process -1 0 1 2 3 4 5 6 7\n");
      newcardCut << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f %6.3f %6.3f\n",nSigCut[6],nSigCut[0],bgdDecay[27],bgdDecay[28],bgdCombined[0],bgdCombined[3],bgdDecay[21],bgdDecay[41],bgdDecay[42]);
      newcardCut << Form("lumi_8TeV                  lnN 1.026 1.026	  -   1.026   -	  1.026   1.026 1.026 1.026\n");		     
      newcardCut << Form("CMS_eff_m                  lnN 1.040 1.040	  -   1.040   -	  1.040   1.040 1.040 1.040\n");		     
      newcardCut << Form("CMS_eff_e                  lnN 1.040 1.040	  -   1.040   -	  1.040   1.040 1.040 1.040\n");		     
      newcardCut << Form("CMS_p_scale_m              lnN 1.015 1.015	  -   1.015   -	  1.015   1.015 1.015 1.015\n");		     
      newcardCut << Form("CMS_p_scale_e              lnN 1.020 1.020	  -   1.020   -	  1.020   1.020 1.020 1.020\n");		     
      newcardCut << Form("CMS_res_met                lnN %5.3f %5.3f      -   %5.3f   -	  %5.3f   %5.3f %5.3f %5.3f\n",CMS_met,CMS_met,CMS_met,CMS_met,CMS_met,CMS_met,CMS_met);
      newcardCut << Form("CMS_res_j                  lnN %5.3f %5.3f    %5.3f %5.3f   -   %5.3f   %5.3f %5.3f %5.3f\n",CMS_jes,CMS_jes,CMS_jes_WZ,CMS_jes_ZZ,CMS_jes,CMS_jes,CMS_jes,CMS_jes);
      newcardCut << Form("CMS_vh3l_pu                lnN 1.010 1.010      -   1.010   -   1.010   1.010 1.010 1.010\n");		
      newcardCut << Form("FakeRate                   lnN   -	  -	  -     -   1.400   -       -     -	 - \n");
      newcardCut << Form("UEPS                       lnN 1.030 1.030	  -     -     -	    -       -   1.030 1.030\n");
      newcardCut << Form("pdf_qqbar                  lnN %5.3f %5.3f    1.015 1.040   -	  1.300   1.500 %5.3f %5.3f\n",PDF_VH,PDF_VH,PDF_VH,PDF_VH);
      newcardCut << Form("QCDscale_VH                lnN %5.3f %5.3f	  -     -     -	    -       -   %5.3f %5.3f\n",QCDscale_VH,QCDscale_VH,QCDscale_VH,QCDscale_VH);	  
      newcardCut << Form("CMS_vh3l_WZ_8TeV           lnN   -	 -      1.050   -     -     -       -     -	-  \n");	  
      newcardCut << Form("CMS_vh3l%d_stat_WH_htt_8TeV  lnN %5.3f   -	 -      -     -     -       -     -     -  \n",nsel,signalError[1]);	
      newcardCut << Form("CMS_vh3l%d_stat_WH_hww_8TeV  lnN   -   %5.3f	 -      -     -     -       -     -     -  \n",nsel,signalError[0]);	
      newcardCut << Form("CMS_vh3l%d_stat_WZ_8TeV      lnN   -	 -     %5.3f    -     -     -       -     -     -  \n",nsel,weiDecay[27]+1.0);	
      newcardCut << Form("CMS_vh3l%d_stat_ZZ_8TeV      lnN   -	 -	 -    %5.3f   -     -       -     -     -  \n",nsel,weiDecay[28]+1.0);
      newcardCut << Form("CMS_vh3l%d_stat_Wjets_8TeV   lnN   -	 -	 -      -   %5.3f   -       -     -     -  \n",nsel,bgdCombinedE[0]+1.0);
      newcardCut << Form("CMS_vh3l%d_stat_Wgamma_8TeV  lnN   -	 -	 -      -     -   %5.3f     -     -     -  \n",nsel,bgdCombinedE[3]+1.0);
      newcardCut << Form("CMS_vh3l%d_stat_VVV_8TeV     lnN   -	 -	 -      -     -     -     %5.3f   -     -  \n",nsel,weiDecay[21]+1.0);
      newcardCut << Form("CMS_vh3l%d_stat_WH_htt_SM_8TeV lnN   -   -	 -      -     -     -       -   %5.3f   -  \n",nsel,weiDecay[41]+1.0);   
      newcardCut << Form("CMS_vh3l%d_stat_WH_hww_SM_8TeV lnN   -   -	 -      -     -     -       -     -   %5.3f\n",nsel,weiDecay[42]+1.0);   
    }

    if(channel == 1101){
      //----------------------------------------------------------------------------
      // Nominal Shapes
      //----------------------------------------------------------------------------
      TH1D *histo_Data    = (TH1D*) histo5 ->Clone("histo_Data");
      TH1D *histo_WH_htt  = (TH1D*) histoTT->Clone("histo_WH_htt");
      TH1D *histo_WH_hww  = (TH1D*) histoWW->Clone("histo_WH_hww");
      histo_Data  ->Rebin(rebinMVAHist);
      histo_WH_htt->Rebin(rebinMVAHist);
      histo_WH_hww->Rebin(rebinMVAHist);
      histo_WZ    ->Rebin(rebinMVAHist);
      histo_ZZ    ->Rebin(rebinMVAHist);
      histo_Wjets ->Rebin(rebinMVAHist);
      histo_Wgamma->Rebin(rebinMVAHist);
      histo_VVV   ->Rebin(rebinMVAHist);
      histo_WH_htt_SM->Rebin(rebinMVAHist);
      histo_WH_hww_SM->Rebin(rebinMVAHist);
      double old_wjets_norm = histo_Wjets->GetSumOfWeights();
      for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
        if(histo_Wjets->GetBinContent(i) < 0) histo_Wjets->SetBinContent(i,0.000001);
      }
      // We need to renormalize
      if(old_wjets_norm > 0) histo_Wjets->Scale(old_wjets_norm/histo_Wjets ->GetSumOfWeights());
      //----------------------------------------------------------------------------
      // Systematics Shapes
      //----------------------------------------------------------------------------
      if(histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVUp->GetNbinsX() != histo_WH_htt->GetNbinsX()) {printf("PROBLEMMMMM\n");return;}
      for(int i=1; i<=histo_WH_htt->GetNbinsX(); i++){
        double factorUp = +1.0; double factorDown = -1.0;
	if(useAlternativeStatTemplates == true){
	  if     (i<1.0*histo_WH_htt->GetNbinsX()/3.0) {factorUp = -1.0; factorDown = +1.0;}
	  else if(i<2.0*histo_WH_htt->GetNbinsX()/3.0) {factorUp = +0.0; factorDown = +0.0;}
	  else                                       {factorUp = +1.0; factorDown = -1.0;}
	}
    	histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_WH_htt   ->GetBinContent(i)+factorUp  *histo_WH_htt  ->GetBinError(i),0.000001));
    	histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVDown      ->SetBinContent(i,TMath::Max(histo_WH_htt   ->GetBinContent(i)+factorDown*histo_WH_htt  ->GetBinError(i),0.000001));
    	histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_WH_hww   ->GetBinContent(i)+factorUp  *histo_WH_hww  ->GetBinError(i),0.000001));
    	histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVDown      ->SetBinContent(i,TMath::Max(histo_WH_hww   ->GetBinContent(i)+factorDown*histo_WH_hww  ->GetBinError(i),0.000001));
    	histo_WZ_CMS_MVAWZStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorUp  *histo_WZ    ->GetBinError(i),0.000001));
    	histo_WZ_CMS_MVAWZStatBounding_8TeVDown           ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorDown*histo_WZ    ->GetBinError(i),0.000001));
    	histo_ZZ_CMS_MVAZZStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorUp  *histo_ZZ    ->GetBinError(i),0.000001));
    	histo_ZZ_CMS_MVAZZStatBounding_8TeVDown           ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorDown*histo_ZZ    ->GetBinError(i),0.000001));
    	histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_Wjets  ->GetBinContent(i)+factorUp  *histo_Wjets ->GetBinError(i),0.000001));
    	histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown     ->SetBinContent(i,TMath::Max(histo_Wjets  ->GetBinContent(i)+factorDown*histo_Wjets ->GetBinError(i),0.000001));
    	histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_Wgamma  ->GetBinContent(i)+factorUp  *histo_Wgamma ->GetBinError(i),0.000001));
    	histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown     ->SetBinContent(i,TMath::Max(histo_Wgamma  ->GetBinContent(i)+factorDown*histo_Wgamma ->GetBinError(i),0.000001));
    	histo_VVV_CMS_MVAVVVStatBounding_8TeVUp	     ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorUp  *histo_VVV ->GetBinError(i),0.000001));
    	histo_VVV_CMS_MVAVVVStatBounding_8TeVDown     ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorDown*histo_VVV ->GetBinError(i),0.000001));
    	histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVUp	   ->SetBinContent(i,TMath::Max(histo_WH_htt_SM   ->GetBinContent(i)+factorUp  *histo_WH_htt_SM  ->GetBinError(i),0.000001));
    	histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVDown      ->SetBinContent(i,TMath::Max(histo_WH_htt_SM   ->GetBinContent(i)+factorDown*histo_WH_htt_SM  ->GetBinError(i),0.000001));
    	histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVUp	   ->SetBinContent(i,TMath::Max(histo_WH_hww_SM   ->GetBinContent(i)+factorUp  *histo_WH_hww_SM  ->GetBinError(i),0.000001));
    	histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVDown      ->SetBinContent(i,TMath::Max(histo_WH_hww_SM   ->GetBinContent(i)+factorDown*histo_WH_hww_SM  ->GetBinError(i),0.000001));
      }
      for(int i=1; i<=histo_WH_htt->GetNbinsX(); i++){	 
        double mean = histo_WH_htt			->GetBinContent(i);
        double up   = histo_WH_htt_CMS_MVAMETResBoundingUp->GetBinContent(i);
        double diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WH_htt_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WH_htt_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_WH_hww			 ->GetBinContent(i);
        up   = histo_WH_hww_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WH_hww_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WH_hww_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_WZ 		       ->GetBinContent(i);
        up   = histo_WZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_ZZ 		       ->GetBinContent(i);
        up   = histo_ZZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ZZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ZZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
	
        mean = histo_Wgamma 		       ->GetBinContent(i);
        up   = histo_Wgamma_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_VVV 		       ->GetBinContent(i);
        up   = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_WH_htt_SM			 ->GetBinContent(i);
        up   = histo_WH_htt_SM_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WH_htt_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WH_htt_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_WH_hww_SM			 ->GetBinContent(i);
        up   = histo_WH_hww_SM_CMS_MVAMETResBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WH_hww_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WH_hww_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

        mean = histo_Wjets 		     ->GetBinContent(i);
        up   = histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
	
        mean = histo_WZ 		      ->GetBinContent(i);
        up   = histo_WZ_CMS_MVAWZNLOBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_WZ_CMS_MVAWZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_WZ_CMS_MVAWZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
	
        mean = histo_ZZ 		      ->GetBinContent(i);
        up   = histo_ZZ_CMS_MVAZZNLOBoundingUp->GetBinContent(i);
        diff = TMath::Abs(mean-up);
        if     (mean-up >0) histo_ZZ_CMS_MVAZZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
        else		    histo_ZZ_CMS_MVAZZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
      }
      histo_WZ_CMS_MVALepEffBoundingUp     ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVALepEffBoundingUp   ->GetSumOfWeights());
      histo_WZ_CMS_MVALepEffBoundingDown   ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVALepEffBoundingDown ->GetSumOfWeights());
      histo_WZ_CMS_MVALepResBoundingUp     ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVALepResBoundingUp   ->GetSumOfWeights());
      histo_WZ_CMS_MVALepResBoundingDown   ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVALepResBoundingDown ->GetSumOfWeights());
      histo_WZ_CMS_MVAMETResBoundingUp     ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAMETResBoundingUp   ->GetSumOfWeights());
      histo_WZ_CMS_MVAMETResBoundingDown   ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAMETResBoundingDown ->GetSumOfWeights());
      histo_WZ_CMS_MVAJESBoundingUp	   ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAJESBoundingUp	  ->GetSumOfWeights());
      histo_WZ_CMS_MVAJESBoundingDown	   ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAJESBoundingDown    ->GetSumOfWeights());
      histo_Wjets_CMS_MVAWBoundingUp	   ->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingUp  ->GetSumOfWeights());
      histo_Wjets_CMS_MVAWBoundingDown     ->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingDown->GetSumOfWeights());
      histo_WZ_CMS_MVAWZNLOBoundingUp      ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAWZNLOBoundingUp  ->GetSumOfWeights());
      histo_WZ_CMS_MVAWZNLOBoundingDown    ->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_MVAWZNLOBoundingDown->GetSumOfWeights());
      histo_ZZ_CMS_MVAZZNLOBoundingUp	   ->Scale(histo_ZZ->GetSumOfWeights()/histo_ZZ_CMS_MVAZZNLOBoundingUp  ->GetSumOfWeights());
      histo_ZZ_CMS_MVAZZNLOBoundingDown    ->Scale(histo_ZZ->GetSumOfWeights()/histo_ZZ_CMS_MVAZZNLOBoundingDown->GetSumOfWeights());

      histo_Wgamma_CMS_MVALepEffBoundingUp  ->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVALepEffBoundingUp   ->GetSumOfWeights());
      histo_Wgamma_CMS_MVALepEffBoundingDown->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVALepEffBoundingDown ->GetSumOfWeights());
      histo_Wgamma_CMS_MVALepResBoundingUp  ->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVALepResBoundingUp   ->GetSumOfWeights());
      histo_Wgamma_CMS_MVALepResBoundingDown->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVALepResBoundingDown ->GetSumOfWeights());
      histo_Wgamma_CMS_MVAMETResBoundingUp  ->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVAMETResBoundingUp   ->GetSumOfWeights());
      histo_Wgamma_CMS_MVAMETResBoundingDown->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVAMETResBoundingDown ->GetSumOfWeights());
      histo_Wgamma_CMS_MVAJESBoundingUp	    ->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVAJESBoundingUp	   ->GetSumOfWeights());
      histo_Wgamma_CMS_MVAJESBoundingDown   ->Scale(histo_Wgamma->GetSumOfWeights()/histo_Wgamma_CMS_MVAJESBoundingDown    ->GetSumOfWeights());

      histo_VVV_CMS_MVALepEffBoundingUp     ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVALepEffBoundingUp   ->GetSumOfWeights());
      histo_VVV_CMS_MVALepEffBoundingDown   ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVALepEffBoundingDown ->GetSumOfWeights());
      histo_VVV_CMS_MVALepResBoundingUp     ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVALepResBoundingUp   ->GetSumOfWeights());
      histo_VVV_CMS_MVALepResBoundingDown   ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVALepResBoundingDown ->GetSumOfWeights());
      histo_VVV_CMS_MVAMETResBoundingUp     ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVAMETResBoundingUp   ->GetSumOfWeights());
      histo_VVV_CMS_MVAMETResBoundingDown   ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVAMETResBoundingDown ->GetSumOfWeights());
      histo_VVV_CMS_MVAJESBoundingUp	    ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVAJESBoundingUp	     ->GetSumOfWeights());
      histo_VVV_CMS_MVAJESBoundingDown	    ->Scale(histo_VVV->GetSumOfWeights()/histo_VVV_CMS_MVAJESBoundingDown    ->GetSumOfWeights());

      char outputLimits[200];
      sprintf(outputLimits,"vh3l%d_%3d.input_8TeV.root",nsel,(int)mhAna);     
      TFile* outFileLimits = new TFile(outputLimits,"recreate");
      outFileLimits->cd();
      histo_Data  ->Write();
      histo_WH_htt->Write();
      histo_WH_hww->Write();
      histo_WZ    ->Write();
      histo_ZZ    ->Write();
      histo_Wjets ->Write();
      histo_Wgamma->Write();
      histo_VVV   ->Write();
      histo_WH_htt_SM->Write();
      histo_WH_hww_SM->Write();
      cout << histo_Data ->GetSumOfWeights()   << " ";
      cout << histo_WH_htt ->GetSumOfWeights() << " ";
      cout << histo_WH_hww ->GetSumOfWeights() << " ";
      cout << histo_WZ   ->GetSumOfWeights()   << " ";
      cout << histo_ZZ   ->GetSumOfWeights()   << " ";
      cout << histo_Wjets->GetSumOfWeights()   << " ";
      cout << histo_Wgamma->GetSumOfWeights()  << " ";
      cout << histo_VVV->GetSumOfWeights()     << " ";
      cout << histo_WH_htt_SM ->GetSumOfWeights() << " ";
      cout << histo_WH_hww_SM ->GetSumOfWeights() << endl;
      histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVUp  ->Write();
      histo_WH_htt_CMS_MVAWH_httStatBounding_8TeVDown->Write();
      histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVUp  ->Write();
      histo_WH_hww_CMS_MVAWH_hwwStatBounding_8TeVDown->Write();
      histo_WZ_CMS_MVAWZStatBounding_8TeVUp  	     ->Write();
      histo_WZ_CMS_MVAWZStatBounding_8TeVDown	     ->Write();
      histo_ZZ_CMS_MVAZZStatBounding_8TeVUp  	     ->Write();
      histo_ZZ_CMS_MVAZZStatBounding_8TeVDown	     ->Write();
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVUp    ->Write();
      histo_Wjets_CMS_MVAWjetsStatBounding_8TeVDown  ->Write();
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVUp  ->Write();
      histo_Wgamma_CMS_MVAWgammaStatBounding_8TeVDown->Write();
      histo_VVV_CMS_MVAVVVStatBounding_8TeVUp        ->Write();
      histo_VVV_CMS_MVAVVVStatBounding_8TeVDown      ->Write();
      histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVUp  ->Write();
      histo_WH_htt_SM_CMS_MVAWH_htt_SMStatBounding_8TeVDown->Write();
      histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVUp  ->Write();
      histo_WH_hww_SM_CMS_MVAWH_hww_SMStatBounding_8TeVDown->Write();

      histo_WH_htt_CMS_MVALepEffBoundingUp	 ->Write();
      histo_WH_htt_CMS_MVALepEffBoundingDown	 ->Write();
      histo_WH_hww_CMS_MVALepEffBoundingUp	 ->Write();
      histo_WH_hww_CMS_MVALepEffBoundingDown	 ->Write();
      histo_WZ_CMS_MVALepEffBoundingUp  	 ->Write();
      histo_WZ_CMS_MVALepEffBoundingDown	 ->Write();
      histo_ZZ_CMS_MVALepEffBoundingUp  	 ->Write();
      histo_ZZ_CMS_MVALepEffBoundingDown	 ->Write();
      histo_Wgamma_CMS_MVALepEffBoundingUp  	 ->Write();
      histo_Wgamma_CMS_MVALepEffBoundingDown	 ->Write();
      histo_VVV_CMS_MVALepEffBoundingUp  	 ->Write();
      histo_VVV_CMS_MVALepEffBoundingDown	 ->Write();
      histo_WH_htt_SM_CMS_MVALepEffBoundingUp	 ->Write();
      histo_WH_htt_SM_CMS_MVALepEffBoundingDown	 ->Write();
      histo_WH_hww_SM_CMS_MVALepEffBoundingUp	 ->Write();
      histo_WH_hww_SM_CMS_MVALepEffBoundingDown	 ->Write();

      histo_WH_htt_CMS_MVALepResBoundingUp	 ->Write();
      histo_WH_htt_CMS_MVALepResBoundingDown	 ->Write();
      histo_WH_hww_CMS_MVALepResBoundingUp	 ->Write();
      histo_WH_hww_CMS_MVALepResBoundingDown	 ->Write();
      histo_WZ_CMS_MVALepResBoundingUp  	 ->Write();
      histo_WZ_CMS_MVALepResBoundingDown	 ->Write();
      histo_ZZ_CMS_MVALepResBoundingUp  	 ->Write();
      histo_ZZ_CMS_MVALepResBoundingDown	 ->Write();
      histo_Wgamma_CMS_MVALepResBoundingUp  	 ->Write();
      histo_Wgamma_CMS_MVALepResBoundingDown	 ->Write();
      histo_VVV_CMS_MVALepResBoundingUp  	 ->Write();
      histo_VVV_CMS_MVALepResBoundingDown	 ->Write();
      histo_WH_htt_SM_CMS_MVALepResBoundingUp	 ->Write();
      histo_WH_htt_SM_CMS_MVALepResBoundingDown	 ->Write();
      histo_WH_hww_SM_CMS_MVALepResBoundingUp	 ->Write();
      histo_WH_hww_SM_CMS_MVALepResBoundingDown	 ->Write();

      histo_WH_htt_CMS_MVAMETResBoundingUp	 ->Write();
      histo_WH_htt_CMS_MVAMETResBoundingDown	 ->Write();
      histo_WH_hww_CMS_MVAMETResBoundingUp	 ->Write();
      histo_WH_hww_CMS_MVAMETResBoundingDown	 ->Write();
      histo_WZ_CMS_MVAMETResBoundingUp  	 ->Write();
      histo_WZ_CMS_MVAMETResBoundingDown	 ->Write();
      histo_ZZ_CMS_MVAMETResBoundingUp  	 ->Write();
      histo_ZZ_CMS_MVAMETResBoundingDown	 ->Write();
      histo_Wgamma_CMS_MVAMETResBoundingUp  	 ->Write();
      histo_Wgamma_CMS_MVAMETResBoundingDown	 ->Write();
      histo_VVV_CMS_MVAMETResBoundingUp  	 ->Write();
      histo_VVV_CMS_MVAMETResBoundingDown	 ->Write();
      histo_WH_htt_SM_CMS_MVAMETResBoundingUp	 ->Write();
      histo_WH_htt_SM_CMS_MVAMETResBoundingDown	 ->Write();
      histo_WH_hww_SM_CMS_MVAMETResBoundingUp	 ->Write();
      histo_WH_hww_SM_CMS_MVAMETResBoundingDown	 ->Write();

      histo_WH_htt_CMS_MVAJESBoundingUp		 ->Write();
      histo_WH_htt_CMS_MVAJESBoundingDown 	 ->Write();
      histo_WH_hww_CMS_MVAJESBoundingUp		 ->Write();
      histo_WH_hww_CMS_MVAJESBoundingDown 	 ->Write();
      histo_WZ_CMS_MVAJESBoundingUp		 ->Write();
      histo_WZ_CMS_MVAJESBoundingDown		 ->Write();
      histo_ZZ_CMS_MVAJESBoundingUp		 ->Write();
      histo_ZZ_CMS_MVAJESBoundingDown		 ->Write();
      histo_Wgamma_CMS_MVAJESBoundingUp  	 ->Write();
      histo_Wgamma_CMS_MVAJESBoundingDown	 ->Write();
      histo_VVV_CMS_MVAJESBoundingUp  	         ->Write();
      histo_VVV_CMS_MVAJESBoundingDown	         ->Write();
      histo_WH_htt_SM_CMS_MVAJESBoundingUp       ->Write();
      histo_WH_htt_SM_CMS_MVAJESBoundingDown 	 ->Write();
      histo_WH_hww_SM_CMS_MVAJESBoundingUp       ->Write();
      histo_WH_hww_SM_CMS_MVAJESBoundingDown 	 ->Write();

      histo_Wjets_CMS_MVAWBoundingUp    	 ->Write();
      histo_Wjets_CMS_MVAWBoundingDown  	 ->Write();
      histo_WZ_CMS_MVAWZNLOBoundingUp    	 ->Write();
      histo_WZ_CMS_MVAWZNLOBoundingDown  	 ->Write();
      histo_ZZ_CMS_MVAZZNLOBoundingUp  	         ->Write();
      histo_ZZ_CMS_MVAZZNLOBoundingDown  	 ->Write();
      outFileLimits->Close();

      char theWH_httString[20];
      if(histo_WH_htt->GetSumOfWeights() > 0) sprintf(theWH_httString,"1.000");
      else                                    sprintf(theWH_httString,"  -  ");
      char theWH_WgammaString[20];
      if(histo_Wgamma->GetSumOfWeights() > 0) sprintf(theWH_WgammaString,"1.000");
      else                                    sprintf(theWH_WgammaString,"  -  ");
      char theWH_htt_SMString[20];
      if(histo_WH_htt_SM->GetSumOfWeights() > 0) sprintf(theWH_htt_SMString,"1.000");
      else                                       sprintf(theWH_htt_SMString,"  -  ");
      char theWH_hww_SMString[20];
      if(histo_WH_hww_SM->GetSumOfWeights() > 0) sprintf(theWH_hww_SMString,"1.000");
      else                                       sprintf(theWH_hww_SMString,"  -  ");
      char outputLimitsShape[200];
      sprintf(outputLimitsShape,"histo_limits_vh3l%d_mh%3.0f_shape_8TeV.txt",nsel,mhAna);     
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",nSelectedData);
      if(useTemplates == true)
        newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
      else
        newcardShape << Form("shapes *   *   %s  histo_$PROCESS\n",outputLimits);
      newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
      newcardShape << Form("bin 1 1 1 1 1 1 1 1 1\n");
      newcardShape << Form("process WH_htt WH_hww WZ ZZ Wjets Wgamma VVV WH_htt_SM WH_hww_SM\n");
      newcardShape << Form("process -1 0 1 2 3 4 5 6 7\n");
      newcardShape << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigCut[6],nSigCut[0],bgdDecay[27],bgdDecay[28],bgdCombined[0],bgdCombined[3],bgdDecay[21],bgdDecay[41],bgdDecay[42]);
      newcardShape << Form("lumi_8TeV				   lnN   1.026 1.026   -   1.026   -   1.026   1.026 1.026 1.026\n");		       
      newcardShape << Form("CMS_MVALepEffBounding		 shape     %s  1.000 1.000 1.000   -     %s    1.000   %s    %s \n",theWH_httString,theWH_WgammaString,theWH_htt_SMString,theWH_hww_SMString);		  
      newcardShape << Form("CMS_p_scale_m			   lnN   1.015 1.015   -   1.015   -   1.015   1.015 1.015 1.015\n");		
      newcardShape << Form("CMS_p_scale_e			   lnN   1.020 1.020   -   1.020   -   1.020   1.020 1.020 1.020\n");		
      newcardShape << Form("CMS_res_met 			   lnN   %5.3f %5.3f   -   %5.3f   -   %5.3f   %5.3f %5.3f %5.3f\n",CMS_met,CMS_met,CMS_met,CMS_met,CMS_met,CMS_met,CMS_met);
      newcardShape << Form("CMS_MVAJESBounding  		 shape     %s  1.000 1.000 1.000   -     %s    1.000   %s    %s \n",theWH_httString,theWH_WgammaString,theWH_htt_SMString,theWH_hww_SMString);			
      newcardShape << Form("CMS_vh3l_pu 			   lnN   1.010 1.010   -   1.010   -   1.000   1.000 1.010 1.010\n");		    
      newcardShape << Form("CMS_fake				   lnN     -	 -     -     -   1.400   -	 -      -     - \n");
      newcardShape << Form("CMS_MVAWBounding			 shape     -	 -     -     -   1.000   -	 -      -     - \n");
      newcardShape << Form("CMS_MVAWZNLOBounding		 shape     -	 -   1.000   -     -	 -	 -      -     - \n");
      newcardShape << Form("CMS_MVAZZNLOBounding		 shape     -	 -     -   1.000   -	 -	 -      -     - \n");
      newcardShape << Form("UEPS				   lnN   1.030 1.030   -     -     -	 -	 -   1.030 1.030\n");
      newcardShape << Form("pdf_qqbar				   lnN   %5.3f %5.3f 1.015 1.040   -	1.300  1.500 %5.3f %5.3f\n",PDF_VH,PDF_VH,PDF_VH,PDF_VH);
      newcardShape << Form("QCDscale_VH 			   lnN   %5.3f %5.3f   -     -     -	 -	 -   %5.3f %5.3f\n",QCDscale_VH,QCDscale_VH,QCDscale_VH,QCDscale_VH);			  
      newcardShape << Form("CMS_vh3l_WZ_8TeV			   lnN     -	 -   1.090   -     -	 -	 -     -     - \n");
      if(histo_WH_htt->GetSumOfWeights() > 0)	 
      newcardShape << Form("CMS_MVAWH_httStatBounding_8TeV_vh3l%d     shape   1.000   -     -     -     -	 -	 -     -     - \n",nsel);
      newcardShape << Form("CMS_MVAWH_hwwStatBounding_8TeV_vh3l%d     shape     -   1.000   -     -     -	 -	 -     -     - \n",nsel);
      newcardShape << Form("CMS_MVAWZStatBounding_8TeV_vh3l%d	      shape     -	 -    1.000  -     -	 -	 -     -     - \n",nsel);
      newcardShape << Form("CMS_MVAZZStatBounding_8TeV_vh3l%d	      shape     -	 -	-   1.000  -	 -	 -     -     - \n",nsel);
      newcardShape << Form("CMS_MVAWjetsStatBounding_8TeV_vh3l%d      shape     -	 -	-    -   1.000   -	 -     -     - \n",nsel);
      if(histo_Wgamma->GetSumOfWeights() > 0)
      newcardShape << Form("CMS_MVAWgammaStatBounding_8TeV_vh3l%d     shape     -	 -	-    -     -    1.000    -     -     - \n",nsel);
      newcardShape << Form("CMS_MVAVVVStatBounding_8TeV_vh3l%d        shape     -	 -	-    -     -     -     1.000   -     - \n",nsel);
      if(histo_WH_htt_SM->GetSumOfWeights() > 0)	 
      newcardShape << Form("CMS_MVAWH_htt_SMStatBounding_8TeV_vh3l%d  shape  -      -     -     -     -	 -	 -   1.000   - \n",nsel);
      if(histo_WH_hww_SM->GetSumOfWeights() > 0)	 
      newcardShape << Form("CMS_MVAWH_hww_SMStatBounding_8TeV_vh3l%d  shape  -      -     -     -     -	 -	 -     -    1.000\n",nsel);
    }
  }
}
