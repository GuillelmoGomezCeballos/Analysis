#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Ana/nt_scripts/trilepton.h"
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
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/HiggsQCDScaleSystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/PSUESystematics_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/PDFgHHSystematics_8TeV.h"

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;

// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCuts2l_53x
(
 int     mH  	 = 300,
 int thePlot = 7,
 TString bgdInputFile    = "ntuples_53x/backgroundA_skim8.root",
 TString signalInputFile = "ntuples_53x/hww125.root",
 TString dataInputFile   = "ntuples_53x/data_skim8.root",
 TString systInputFile   = "",
 bool fillInfoNote = false,
 int lDecay = 4,
 int period = 3
 )
{

  double lumi = 1.0;
 
  bool fCheckProblem = true;
  double FRCorrMu = 1.0; double FRCorrEl = 1.0; double WSSF = 1.0;
  if(mH == 301) {FRCorrMu = 0.65; FRCorrEl = 1.00; WSSF = 1.25;}

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  if(systInputFile != ""){
    systEvent.LoadTree(systInputFile,-1);
    systEvent.InitTree(0);
  }

  char finalStateName[2];
  sprintf(finalStateName,"ll");
  if	 (lDecay == 0) sprintf(finalStateName,"mm");
  else if(lDecay == 1) sprintf(finalStateName,"me");
  else if(lDecay == 2) sprintf(finalStateName,"em");
  else if(lDecay == 3) sprintf(finalStateName,"ee");
  else if(lDecay == 5) sprintf(finalStateName,"sf");
  else if(lDecay == 6) sprintf(finalStateName,"of");

  TString ECMsb  = "";
  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double lumiE = 1.099;
  if	 (period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;minRun =      0;maxRun = 999999;ECMsb="8TeV";lumiE = 1.026;
  }
  else if(period == 4){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
    lumi     = 4.924;minRun =	 0;maxRun = 999999;ECMsb="7TeV"; lumiE = 1.022;
    UseDyttDataDriven = false;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }
  
  //----------------------------------------------------------------------------
  // radio photon to electron
  //----------------------------------------------------------------------------
  TFile *fRatioPhotonElectron = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_photon_electron.root");
  TH1D *fhDRatioPhotonElectron = (TH1D*)(fRatioPhotonElectron->Get("hDRatioPhotonElectron"));
  assert(fhDRatioPhotonElectron);
  fhDRatioPhotonElectron->SetDirectory(0);
  fRatioPhotonElectron->Close();
  delete fRatioPhotonElectron;

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
 
  LeptonScaleLookup trigLookup(effPath.Data());

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int channel = mH;
  if     (channel > 1000 && channel < 2000) mH = channel-1000;
  else if(channel == 300                  ) mH = channel;
  else if(channel == 301                  ) mH = channel;
  else assert(0);

  Float_t bdtg = 0.0;
  Float_t bdtg_aux0 = 0.0;
  Float_t bdtg_aux1 = 0.0;
  Float_t bdtg_aux2 = 0.0;
  Float_t mll_lepup,mt_lepup,mll_lepdown,mt_lepdown,mll_metup,mt_metup;  
  int nBin    = 100;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 30; xminPlot = -0.5; xmaxPlot = 29.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 20 && thePlot <= 22) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 26 && thePlot <= 26) {nBinPlot = 100; xminPlot =-2.5; xmaxPlot = 2.5;}
  else if(thePlot >= 23 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 30 && thePlot <= 30) {nBinPlot = 200; xminPlot = -5.0; xmaxPlot = 15.0;}
  else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 300; xminPlot = 0.0; xmaxPlot = 600.0;}
  else if(thePlot >= 33 && thePlot <= 33) {nBinPlot = 90; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  800.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  2000.0;}
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  5.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 49) {nBinPlot = 200; xminPlot = 0.; xmaxPlot = 10.;}
  else if(thePlot >= 50 && thePlot <= 52) {nBinPlot = 200; xminPlot = 0.; xmaxPlot = 200.;}
  else if(thePlot >= 53 && thePlot <= 55) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 56 && thePlot <= 56) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 57 && thePlot <= 57) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}
  else if(thePlot >= 58 && thePlot <= 59) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 60 && thePlot <= 60) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 61 && thePlot <= 61) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot =  2000.0;}
  else if(thePlot >= 62 && thePlot <= 62) {nBinPlot = 144; xminPlot =  0.0; xmaxPlot =  12.0;}
  else if(thePlot >= 63 && thePlot <= 63) {nBinPlot = 120; xminPlot =  0.0; xmaxPlot =  12.0;}
  nBin = nBinPlot;

  double nSigCut[6]  = {0,0,0,0,0,0};
  double nSigECut[6] = {0,0,0,0,0,0};
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
  TH1D *histo_Data   = (TH1D*) histos->Clone("histo_Data");
  TH1D *histo_ttH    = (TH1D*) histos->Clone("histo_ttH");
  TH1D *histo_ZH     = (TH1D*) histos->Clone("histo_ZH");
  TH1D *histo_WH     = (TH1D*) histos->Clone("histo_WH");
  TH1D *histo_WZ     = (TH1D*) histos->Clone("histo_WZ");
  TH1D *histo_WS     = (TH1D*) histos->Clone("histo_WS");
  TH1D *histo_WjetsE = (TH1D*) histos->Clone("histo_WjetsE");
  TH1D *histo_WjetsM = (TH1D*) histos->Clone("histo_WjetsM");
  TH1D *histo_Wgamma = (TH1D*) histos->Clone("histo_Wgamma");
  TH1D *histo_Wg3l   = (TH1D*) histos->Clone("histo_Wg3l");
  TH1D *histo_VVV    = (TH1D*) histos->Clone("histo_VVV");
  TH1D *histo_ttH_SM = (TH1D*) histos->Clone("histo_ttH_SM");
  TH1D *histo_ZH_SM  = (TH1D*) histos->Clone("histo_ZH_SM");
  TH1D *histo_WH_SM  = (TH1D*) histos->Clone("histo_WH_SM");

  TH1D* histo_ttH_CMS_MVAttHStatBoundingUp   = new TH1D( Form("histo_ttH_CMS_vhss%s_MVAttHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ttH_CMS_vhss%s_MVAttHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAttHStatBoundingUp  ->Sumw2();
  TH1D* histo_ttH_CMS_MVAttHStatBoundingDown = new TH1D( Form("histo_ttH_CMS_vhss%s_MVAttHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ttH_CMS_vhss%s_MVAttHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAttHStatBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_MVAZHStatBoundingUp   = new TH1D( Form("histo_ZH_CMS_vhss%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZH_CMS_vhss%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_ZH_CMS_MVAZHStatBoundingDown = new TH1D( Form("histo_ZH_CMS_vhss%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZH_CMS_vhss%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAZHStatBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_MVAWHStatBoundingUp   = new TH1D( Form("histo_WH_CMS_vhss%s_MVAWHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WH_CMS_vhss%s_MVAWHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAWHStatBoundingUp  ->Sumw2();
  TH1D* histo_WH_CMS_MVAWHStatBoundingDown = new TH1D( Form("histo_WH_CMS_vhss%s_MVAWHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WH_CMS_vhss%s_MVAWHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAWHStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_vhss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_vhss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingUp   = new TH1D( Form("histo_WS_CMS_vhss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WS_CMS_vhss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAWSStatBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingDown = new TH1D( Form("histo_WS_CMS_vhss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WS_CMS_vhss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAWSStatBoundingDown->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingUp   = new TH1D( Form("histo_WjetsE_CMS_vhss%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_vhss%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingDown = new TH1D( Form("histo_WjetsE_CMS_vhss%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_vhss%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingUp   = new TH1D( Form("histo_WjetsM_CMS_vhss%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_vhss%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingDown = new TH1D( Form("histo_WjetsM_CMS_vhss%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_vhss%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingUp   = new TH1D( Form("histo_Wgamma_CMS_vhss%s_MVAWgammaStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wgamma_CMS_vhss%s_MVAWgammaStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAWgammaStatBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingDown = new TH1D( Form("histo_Wgamma_CMS_vhss%s_MVAWgammaStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wgamma_CMS_vhss%s_MVAWgammaStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAWgammaStatBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingUp   = new TH1D( Form("histo_Wg3l_CMS_vhss%s_MVAWg3lStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wg3l_CMS_vhss%s_MVAWg3lStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAWg3lStatBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingDown = new TH1D( Form("histo_Wg3l_CMS_vhss%s_MVAWg3lStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wg3l_CMS_vhss%s_MVAWg3lStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAWg3lStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp   = new TH1D( Form("histo_VVV_CMS_vhss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_vhss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown = new TH1D( Form("histo_VVV_CMS_vhss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_vhss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAttH_SMStatBoundingUp   = new TH1D( Form("histo_ttH_SM_CMS_vhss%s_MVAttH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ttH_SM_CMS_vhss%s_MVAttH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAttH_SMStatBoundingUp  ->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAttH_SMStatBoundingDown = new TH1D( Form("histo_ttH_SM_CMS_vhss%s_MVAttH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ttH_SM_CMS_vhss%s_MVAttH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAttH_SMStatBoundingDown->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAZH_SMStatBoundingUp   = new TH1D( Form("histo_ZH_SM_CMS_vhss%s_MVAZH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZH_SM_CMS_vhss%s_MVAZH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAZH_SMStatBoundingUp  ->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAZH_SMStatBoundingDown = new TH1D( Form("histo_ZH_SM_CMS_vhss%s_MVAZH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZH_SM_CMS_vhss%s_MVAZH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAZH_SMStatBoundingDown->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAWH_SMStatBoundingUp   = new TH1D( Form("histo_WH_SM_CMS_vhss%s_MVAWH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WH_SM_CMS_vhss%s_MVAWH_SMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAWH_SMStatBoundingUp  ->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAWH_SMStatBoundingDown = new TH1D( Form("histo_WH_SM_CMS_vhss%s_MVAWH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WH_SM_CMS_vhss%s_MVAWH_SMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAWH_SMStatBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_ttH_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_ttH_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ttH_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_ttH_CMS_vhss_MVALepEffBoundingDown"), Form("histo_ttH_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_ZH_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_ZH_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ZH_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_ZH_CMS_vhss_MVALepEffBoundingDown"), Form("histo_ZH_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_WH_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_WH_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WH_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_WH_CMS_vhss_MVALepEffBoundingDown"), Form("histo_WH_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_WZ_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss_MVALepEffBoundingDown"), Form("histo_WZ_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_WS_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_WS_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_WS_CMS_vhss_MVALepEffBoundingDown"), Form("histo_WS_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_Wgamma_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_Wgamma_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_Wgamma_CMS_vhss_MVALepEffBoundingDown"), Form("histo_Wgamma_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_Wg3l_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_Wg3l_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_Wg3l_CMS_vhss_MVALepEffBoundingDown"), Form("histo_Wg3l_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_VVV_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_VVV_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_VVV_CMS_vhss_MVALepEffBoundingDown"), Form("histo_VVV_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_ttH_SM_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVALepEffBoundingDown"), Form("histo_ttH_SM_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_ZH_SM_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVALepEffBoundingDown"), Form("histo_ZH_SM_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WH_SM_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_WH_SM_CMS_vhss_MVALepEffBoundingUp")  , Form("histo_WH_SM_CMS_vhss_MVALepEffBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WH_SM_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_WH_SM_CMS_vhss_MVALepEffBoundingDown"), Form("histo_WH_SM_CMS_vhss_MVALepEffBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVALepEffBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_ttH_CMS_vhss_MVALepResBoundingUp")  , Form("histo_ttH_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ttH_CMS_MVALepResBoundingDown = new TH1D( Form("histo_ttH_CMS_vhss_MVALepResBoundingDown"), Form("histo_ttH_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_ZH_CMS_vhss_MVALepResBoundingUp")  , Form("histo_ZH_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_CMS_MVALepResBoundingDown = new TH1D( Form("histo_ZH_CMS_vhss_MVALepResBoundingDown"), Form("histo_ZH_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_WH_CMS_vhss_MVALepResBoundingUp")  , Form("histo_WH_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WH_CMS_MVALepResBoundingDown = new TH1D( Form("histo_WH_CMS_vhss_MVALepResBoundingDown"), Form("histo_WH_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss_MVALepResBoundingUp")  , Form("histo_WZ_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss_MVALepResBoundingDown"), Form("histo_WZ_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_WS_CMS_vhss_MVALepResBoundingUp")  , Form("histo_WS_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepResBoundingDown = new TH1D( Form("histo_WS_CMS_vhss_MVALepResBoundingDown"), Form("histo_WS_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_Wgamma_CMS_vhss_MVALepResBoundingUp")  , Form("histo_Wgamma_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepResBoundingDown = new TH1D( Form("histo_Wgamma_CMS_vhss_MVALepResBoundingDown"), Form("histo_Wgamma_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_Wg3l_CMS_vhss_MVALepResBoundingUp")  , Form("histo_Wg3l_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepResBoundingDown = new TH1D( Form("histo_Wg3l_CMS_vhss_MVALepResBoundingDown"), Form("histo_Wg3l_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_VVV_CMS_vhss_MVALepResBoundingUp")  , Form("histo_VVV_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingDown = new TH1D( Form("histo_VVV_CMS_vhss_MVALepResBoundingDown"), Form("histo_VVV_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVALepResBoundingUp")  , Form("histo_ttH_SM_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVALepResBoundingDown = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVALepResBoundingDown"), Form("histo_ttH_SM_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVALepResBoundingUp")  , Form("histo_ZH_SM_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVALepResBoundingDown = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVALepResBoundingDown"), Form("histo_ZH_SM_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WH_SM_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_WH_SM_CMS_vhss_MVALepResBoundingUp")  , Form("histo_WH_SM_CMS_vhss_MVALepResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WH_SM_CMS_MVALepResBoundingDown = new TH1D( Form("histo_WH_SM_CMS_vhss_MVALepResBoundingDown"), Form("histo_WH_SM_CMS_vhss_MVALepResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVALepResBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_ttH_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_ttH_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ttH_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_ttH_CMS_vhss_MVAMETResBoundingDown"), Form("histo_ttH_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_ZH_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_ZH_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_ZH_CMS_vhss_MVAMETResBoundingDown"), Form("histo_ZH_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_WH_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_WH_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WH_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_WH_CMS_vhss_MVAMETResBoundingDown"), Form("histo_WH_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_WZ_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss_MVAMETResBoundingDown"), Form("histo_WZ_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_WS_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_WS_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_WS_CMS_vhss_MVAMETResBoundingDown"), Form("histo_WS_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_Wgamma_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_Wgamma_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_Wgamma_CMS_vhss_MVAMETResBoundingDown"), Form("histo_Wgamma_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_Wg3l_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_Wg3l_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_Wg3l_CMS_vhss_MVAMETResBoundingDown"), Form("histo_Wg3l_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_VVV_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_VVV_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_VVV_CMS_vhss_MVAMETResBoundingDown"), Form("histo_VVV_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_ttH_SM_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVAMETResBoundingDown"), Form("histo_ttH_SM_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_ZH_SM_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVAMETResBoundingDown"), Form("histo_ZH_SM_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_WH_SM_CMS_vhss_MVAMETResBoundingUp")  , Form("histo_WH_SM_CMS_vhss_MVAMETResBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_WH_SM_CMS_vhss_MVAMETResBoundingDown"), Form("histo_WH_SM_CMS_vhss_MVAMETResBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAMETResBoundingDown->Sumw2();

  TH1D* histo_ttH_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_ttH_CMS_vhss_MVAJESBoundingUp")  , Form("histo_ttH_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ttH_CMS_MVAJESBoundingDown = new TH1D( Form("histo_ttH_CMS_vhss_MVAJESBoundingDown"), Form("histo_ttH_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZH_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_ZH_CMS_vhss_MVAJESBoundingUp")  , Form("histo_ZH_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZH_CMS_MVAJESBoundingDown = new TH1D( Form("histo_ZH_CMS_vhss_MVAJESBoundingDown"), Form("histo_ZH_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WH_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_WH_CMS_vhss_MVAJESBoundingUp")  , Form("histo_WH_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WH_CMS_MVAJESBoundingDown = new TH1D( Form("histo_WH_CMS_vhss_MVAJESBoundingDown"), Form("histo_WH_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss_MVAJESBoundingUp")  , Form("histo_WZ_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss_MVAJESBoundingDown"), Form("histo_WZ_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_WS_CMS_vhss_MVAJESBoundingUp")  , Form("histo_WS_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingDown = new TH1D( Form("histo_WS_CMS_vhss_MVAJESBoundingDown"), Form("histo_WS_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_Wgamma_CMS_vhss_MVAJESBoundingUp")  , Form("histo_Wgamma_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAJESBoundingDown = new TH1D( Form("histo_Wgamma_CMS_vhss_MVAJESBoundingDown"), Form("histo_Wgamma_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wgamma_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_Wg3l_CMS_vhss_MVAJESBoundingUp")  , Form("histo_Wg3l_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAJESBoundingDown = new TH1D( Form("histo_Wg3l_CMS_vhss_MVAJESBoundingDown"), Form("histo_Wg3l_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_Wg3l_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_VVV_CMS_vhss_MVAJESBoundingUp")  , Form("histo_VVV_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown = new TH1D( Form("histo_VVV_CMS_vhss_MVAJESBoundingDown"), Form("histo_VVV_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVAJESBoundingUp")  , Form("histo_ttH_SM_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ttH_SM_CMS_MVAJESBoundingDown = new TH1D( Form("histo_ttH_SM_CMS_vhss_MVAJESBoundingDown"), Form("histo_ttH_SM_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ttH_SM_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVAJESBoundingUp")  , Form("histo_ZH_SM_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZH_SM_CMS_MVAJESBoundingDown = new TH1D( Form("histo_ZH_SM_CMS_vhss_MVAJESBoundingDown"), Form("histo_ZH_SM_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_ZH_SM_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAJESBoundingUp   = new TH1D( Form("histo_WH_SM_CMS_vhss_MVAJESBoundingUp")  , Form("histo_WH_SM_CMS_vhss_MVAJESBoundingUp")  , nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WH_SM_CMS_MVAJESBoundingDown = new TH1D( Form("histo_WH_SM_CMS_vhss_MVAJESBoundingDown"), Form("histo_WH_SM_CMS_vhss_MVAJESBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WH_SM_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_WjetsE_CMS_MVAWEBoundingUp   = new TH1D( Form("histo_WjetsE_CMS_vhss_MVAWEBoundingUp"),   Form("histo_WjetsE_CMS_vhss_MVAWEBoundingUp"),   nBinPlot,xminPlot,xmaxPlot); histo_WjetsE_CMS_MVAWEBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWEBoundingDown = new TH1D( Form("histo_WjetsE_CMS_vhss_MVAWEBoundingDown"), Form("histo_WjetsE_CMS_vhss_MVAWEBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WjetsE_CMS_MVAWEBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWMBoundingUp   = new TH1D( Form("histo_WjetsM_CMS_vhss_MVAWMBoundingUp"),   Form("histo_WjetsM_CMS_vhss_MVAWMBoundingUp"),   nBinPlot,xminPlot,xmaxPlot); histo_WjetsM_CMS_MVAWMBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWMBoundingDown = new TH1D( Form("histo_WjetsM_CMS_vhss_MVAWMBoundingDown"), Form("histo_WjetsM_CMS_vhss_MVAWMBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WjetsM_CMS_MVAWMBoundingDown->Sumw2();

  TH1D* histo_WZ_CMS_WZNLOBoundingUp   = new TH1D( Form("histo_WZ_CMS_vhss_WZNLOBoundingUp"),   Form("histo_WZ_CMS_vhss_WZNLOBoundingUp"),   nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_WZNLOBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_WZNLOBoundingDown = new TH1D( Form("histo_WZ_CMS_vhss_WZNLOBoundingDown"), Form("histo_WZ_CMS_vhss_WZNLOBoundingDown"), nBinPlot,xminPlot,xmaxPlot); histo_WZ_CMS_WZNLOBoundingDown->Sumw2();

  double bgdDecayFake[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double bgdDecayFakeE[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  double bgdDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double weiDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  unsigned int patternTopVeto         = SmurfTree::TopVeto;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);
    if(channel != 300 && channel != 301 && channel != 1125){
      bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww"       ,(int)125), &bdtg);
      bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux0"  ,(int)125), &bdtg_aux0);
      bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux1"  ,(int)125), &bdtg_aux1);
      bgdEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux2"  ,(int)125), &bdtg_aux2);
      bgdEvent.tree_->SetBranchAddress( "mll_lepup"      , &mll_lepup      );
      bgdEvent.tree_->SetBranchAddress( "mt_lepup"       , &mt_lepup       );
      bgdEvent.tree_->SetBranchAddress( "mll_lepdown"    , &mll_lepdown    );
      bgdEvent.tree_->SetBranchAddress( "mt_lepdown"     , &mt_lepdown     );
      bgdEvent.tree_->SetBranchAddress( "mll_metup"      , &mll_metup      );
      bgdEvent.tree_->SetBranchAddress( "mt_metup"       , &mt_metup       );
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
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.processId_==121 ||
            bgdEvent.processId_==122)   fDecay = 41;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 43;
    else if(bgdEvent.processId_==10001) fDecay = 44;
    else if(bgdEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
    if(bgdEvent.lep1MotherMcId_ == 23 && bgdEvent.lep2MotherMcId_ == 23) {
      //fDecay = 31;
    }

    double deltaPhiQQL[3] = {DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep1_.Phi()),DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool   passMET = usedMet > 30.;
    double ptww[3] = {bgdEvent.met_*cos(bgdEvent.metPhi_)+bgdEvent.lep1_.Px()+bgdEvent.lep2_.Px(),
                      bgdEvent.met_*sin(bgdEvent.metPhi_)+bgdEvent.lep1_.Py()+bgdEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);
    bool dPhiDiLepJetCut = true;
    if(bgdEvent.njets_ <= 1) dPhiDiLepJetCut = bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
    					       bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    else		     dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
        				       bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    if(bgdEvent.njets_ >= 1) dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    if((channel > 1000 && channel < 2000) || channel == 300) {passMET = usedMet > 30.; if(bgdEvent.type_ == SmurfTree::ee) passMET = usedMet > 50.;}
    double MJJL[2] = {0.0,0.0};
    if(DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep1_.Phi()) > DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.lep2_.Phi())) {
      MJJL[0] =  (bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.lep1_).M();
      MJJL[1] =  (bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.lep2_).M();
    } else {
      MJJL[0] =  (bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.lep2_).M();
      MJJL[1] =  (bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.lep1_).M();
    }

    int centrality = 0;
    if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
       ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    bool passCuts = false;
    bool isSignalDecay = false;
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;
    double HTMiss = (bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_+bgdEvent.jet3_+bgdEvent.jet4_).Pt();
    //double LD = bgdEvent.met_ * 0.00397 + HTMiss * 0.00265;

    if     (channel > 1000 && channel < 2000){ // HW->2l selection
      if(
          bgdEvent.dilep_.M() > 12 && bgdEvent.dilep_.M() < 200 &&
         (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 bgdEvent.njets_ <= 3. && bgdEvent.njets_ >= 2. &&
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
	if(isRealLepton == false &&
	   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
	    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
	    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // HW->2l selection
    else if(channel == 300){ // WW->2l selection
      if(
         bgdEvent.dilep_.M() > 200 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 bgdEvent.njets_ >= 2 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
         passMET == true &&
	 (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 &&
	 TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 &&
	 centrality == 1 &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
         //(dPhiDiLepJetCut == true || bgdEvent.type_ != SmurfTree::ee) && 
         (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	
	passCuts = true;
	if(isRealLepton == false &&
	   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
	    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
	    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // WW->2l selection
    else if(channel == 301){ // ttH selection
      int nBtag[3] = {0, 0, 0};
      if(bgdEvent.jet1Btag_     > 2.100) nBtag[0]++; if(bgdEvent.jet2Btag_     > 2.100) nBtag[0]++; if(bgdEvent.jet3Btag_     > 2.100) nBtag[0]++; if(bgdEvent.jet4Btag_     > 2.100) nBtag[0]++;
      if(bgdEvent.jet1ProbBtag_ > 0.244) nBtag[1]++; if(bgdEvent.jet2ProbBtag_ > 0.244) nBtag[1]++; if(bgdEvent.jet3ProbBtag_ > 0.244) nBtag[1]++; if(bgdEvent.jet4ProbBtag_ > 0.244) nBtag[1]++;
      if(bgdEvent.jet1ProbBtag_ > 0.679) nBtag[2]++; if(bgdEvent.jet2ProbBtag_ > 0.679) nBtag[2]++; if(bgdEvent.jet3ProbBtag_ > 0.679) nBtag[2]++; if(bgdEvent.jet4ProbBtag_ > 0.679) nBtag[2]++;
      if(
         bgdEvent.dilep_.M() > 200 &&
        (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20. &&
	 //bgdEvent.jet4_.Pt() > 25. && LD > 0.2 && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 bgdEvent.njets_ >= 1 && bgdEvent.met_ > 30. && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 bgdEvent.lep1_.Pt()+bgdEvent.lep2_.Pt()+bgdEvent.met_ > 100.0 &&
	 fabs(bgdEvent.jet1_.Eta()) < 2.5 && fabs(bgdEvent.jet2_.Eta()) < 2.5 &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 10. || bgdEvent.type_ != SmurfTree::ee) && 
	 (bgdEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	
	passCuts = true;
	if(isRealLepton == false &&
	   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
	    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
	    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www  || bgdEvent.dstype_ == SmurfTree::qqww2j)) passCuts = false;
      }
    } // ttH selection

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
      if(nFake < 0) assert(0);
 
      if(nFake > 1){
	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	fDecay = 22;
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
          if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);

	  if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)||
	     ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection))
	     add = add * FRCorrMu;
	  else 
	  if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4)  == SmurfTree::Lep1LooseEleV4  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)||
	     ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4)  == SmurfTree::Lep2LooseEleV4  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection))
	     add = add * FRCorrEl;

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
	  if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  fDecay                 = 3;
	  isSignalDecay          = false;
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
        
        if(fCheckProblem == true && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	  printf("PROBLEMC(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",bgdEvent.event_,add1,add2,trigEff,add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);

        int correctQ = 0;
        if(bgdEvent.lep1McId_ ==  11 && bgdEvent.lq1_ < 0) correctQ++;
        if(bgdEvent.lep1McId_ == -11 && bgdEvent.lq1_ > 0) correctQ++;
        if(bgdEvent.lep1McId_ ==  13 && bgdEvent.lq1_ < 0) correctQ++;
        if(bgdEvent.lep1McId_ == -13 && bgdEvent.lq1_ > 0) correctQ++;
        if(bgdEvent.lep2McId_ ==  11 && bgdEvent.lq2_ < 0) correctQ++;
        if(bgdEvent.lep2McId_ == -11 && bgdEvent.lq2_ > 0) correctQ++;
        if(bgdEvent.lep2McId_ ==  13 && bgdEvent.lq2_ < 0) correctQ++;
        if(bgdEvent.lep2McId_ == -13 && bgdEvent.lq2_ > 0) correctQ++;
	if(correctQ != 2 && channel == 301) add = add*WSSF;

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
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = bgdEvent.met_/bgdEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = bgdEvent.njets_;
      else if(thePlot ==18) myVar = bgdEvent.nvtx_;
      else if(thePlot ==19) myVar = HTMiss;
      else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(bgdEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(bgdEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(bgdEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(bgdEvent.lid1_)==13) myVar = bgdEvent.lep1_.Eta();else if(fabs(bgdEvent.lid2_)==13) myVar = bgdEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = TMath::Max(bgdEvent.jet1Btag_,bgdEvent.jet2Btag_);
      else if(thePlot ==37) myVar = (bgdEvent.jet1_+bgdEvent.jet2_).M();
      else if(thePlot ==38) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
      else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,bgdEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==44) myVar = bgdEvent.jet1_.Pt()+ bgdEvent.jet2_.Pt()+bgdEvent.jet3_.Pt();
      else if(thePlot ==48) myVar = bgdEvent.type_;
      else if(thePlot ==49) myVar = (bgdEvent.jet1_.Pt()*bgdEvent.jet2_.Pt()*bgdEvent.lep1_.Pt()*bgdEvent.lep2_.Pt())/10000000.0;
      else if(thePlot ==50) myVar = TMath::Min(ptww[2],199.999);
      else if(thePlot ==51) myVar = TMath::Min(TMath::Min(MJJL[0],MJJL[1]),199.999);
      else if(thePlot ==52) myVar = TMath::Min(MJJL[1],199.999);
      else if(thePlot ==53) myVar = DeltaPhi(bgdEvent.jet3_.Phi(),bgdEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = bgdEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==57) myVar = bgdEvent.dR_;
      else if(thePlot ==60) myVar = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.jet2_.Phi(),bgdEvent.jet2_.Eta());
      else if(thePlot ==61) myVar = (bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M();
      else if(thePlot ==62) {
        myVar     = Unroll2VarTo1ForWH2l(bgdEvent.dilep_.M(), bgdEvent.mt_)+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_;
	bdtg      = Unroll2VarTo1ForWH2l(bgdEvent.dilep_.M(), bgdEvent.mt_);
	bdtg_aux0 = Unroll2VarTo1ForWH2l(mll_lepup  , mt_lepup  );
	bdtg_aux1 = Unroll2VarTo1ForWH2l(mll_lepdown, mt_lepdown);
	bdtg_aux2 = Unroll2VarTo1ForWH2l(mll_metup  , mt_metup  );
      }
      else if(thePlot ==63) myVar = (bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_;

      double addLepEff     = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
        		     leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
      double addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
        		     leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
      double addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
        		     leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
      double NjetSyst[2] = {0., 0.};
      if(bgdEvent.jet1_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(bgdEvent.jet2_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(bgdEvent.jet3_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(bgdEvent.jet4_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(bgdEvent.jet1_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(bgdEvent.jet2_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(bgdEvent.jet3_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(bgdEvent.jet4_.Pt()*0.95 > 30) NjetSyst[1]++;

      if     (fDecay == 27){
	histo_WZ                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WZ_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_WZ_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_WZ_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WZ_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WZ_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WZ_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_WZ_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 14 || fDecay == 29 || fDecay == 30 || fDecay == 6 || fDecay ==  7 || fDecay ==  8 || fDecay ==  9 ||
              fDecay == 10 || fDecay == 31 || fDecay ==  4 || fDecay == 5 || fDecay == 11 || fDecay == 12 || fDecay == 13 || fDecay == 28){
	histo_WS                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WS_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_WS_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_WS_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WS_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WS_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WS_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_WS_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 0 || fDecay == 1 || fDecay == 2 || fDecay == 3){
        double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        histo_WjetsE                    ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_,theWeight);
        histo_WjetsE_CMS_MVAWEBoundingUp->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_,theWeight*addFRS/addFR);
      }
      else if(fDecay == 23){
        double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        histo_WjetsM                    ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_,theWeight);
        histo_WjetsM_CMS_MVAWMBoundingUp->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_,theWeight*addFRS/addFR);
      }
      else if(fDecay == 19){
        //histo_Wgamma                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        //histo_Wgamma_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        //histo_Wgamma_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        //histo_Wgamma_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        //histo_Wgamma_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        //histo_Wgamma_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        //histo_Wgamma_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        //histo_Wgamma_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 20){
        histo_Wg3l                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_Wg3l_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_Wg3l_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_Wg3l_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_Wg3l_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_Wg3l_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_Wg3l_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_Wg3l_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 21){
	histo_VVV                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_VVV_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_VVV_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_VVV_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_VVV_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_VVV_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_VVV_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_VVV_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 41){
        histo_ttH_SM                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ttH_SM_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_ttH_SM_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_ttH_SM_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ttH_SM_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ttH_SM_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ttH_SM_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_ttH_SM_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 42){
        histo_ZH_SM                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ZH_SM_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_ZH_SM_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_ZH_SM_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ZH_SM_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ZH_SM_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_ZH_SM_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_ZH_SM_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else if(fDecay == 43){
        histo_WH_SM                          ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WH_SM_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_WH_SM_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_WH_SM_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WH_SM_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WH_SM_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(bgdEvent.njets_)+3.0*(double)bgdEvent.type_, theWeight);
        histo_WH_SM_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)bgdEvent.type_, theWeight);
        histo_WH_SM_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)bgdEvent.type_, theWeight);
      }
      else {
        printf("NOOOOOOOOOOOOOOOOOOOO (%d)\n", fDecay); assert(0);
      }

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
      else if(fDecay == 41 || fDecay == 42 || fDecay == 43){
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

  if(systInputFile != ""){
  int nSyst=systEvent.tree_->GetEntries();
  for (int i=0; i<nSyst; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nSyst);
    systEvent.tree_->GetEntry(i);
    if(channel != 300 && channel != 301 && channel != 1125){
      systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww"       ,(int)125), &bdtg);
      systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux0"  ,(int)125), &bdtg_aux0);
      systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux1"  ,(int)125), &bdtg_aux1);
      systEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux2"  ,(int)125), &bdtg_aux2);
      systEvent.tree_->SetBranchAddress( "mll_lepup"      , &mll_lepup      );
      systEvent.tree_->SetBranchAddress( "mt_lepup"	    , &mt_lepup       );
      systEvent.tree_->SetBranchAddress( "mll_lepdown"    , &mll_lepdown    );
      systEvent.tree_->SetBranchAddress( "mt_lepdown"     , &mt_lepdown     );
      systEvent.tree_->SetBranchAddress( "mll_metup"      , &mll_metup      );
      systEvent.tree_->SetBranchAddress( "mt_metup"	    , &mt_metup       );
    }
    if(!(((systEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  systEvent.dstype_ != SmurfTree::data)) continue;
    if(systEvent.dstype_ == SmurfTree::data &&
      (systEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ <  minRun) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (systEvent.dstype_ == SmurfTree::wjets 	   ) fDecay = 3;
    else if(systEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay = 5;
    else if(systEvent.dstype_ == SmurfTree::dyee  	   ) fDecay = 9;
    else if(systEvent.dstype_ == SmurfTree::dymm  	   ) fDecay = 9;
    else if(systEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(systEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(systEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(systEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(systEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(systEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(systEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(systEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(systEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(systEvent.processId_==121 ||
            systEvent.processId_==122)   fDecay = 41;
    else if(systEvent.processId_==24)    fDecay = 42;
    else if(systEvent.processId_==26)    fDecay = 43;
    else if(systEvent.processId_==10001) fDecay = 44;
    else if(systEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;}

    int charge = (int)(systEvent.lq1_ + systEvent.lq2_);
    if((fDecay == 27 || fDecay == 28) && channel != 128 && 
       (channel <= 1000 || channel >= 2000)) {
      if(systEvent.lep1MotherMcId_ == 23 && systEvent.lep2MotherMcId_ == 23) {
        fDecay = 31;
      }
    }

    double deltaPhiQQL[3] = {DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.lep1_.Phi()),DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double usedMet = TMath::Min(systEvent.pmet_,systEvent.pTrackMet_);
    bool   passMET = usedMet > 30.;
    bool dPhiDiLepJetCut = true;
    if(systEvent.njets_ <= 1) dPhiDiLepJetCut = systEvent.jet1_.Pt() <= 15. || systEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
    					       systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;
    else		     dPhiDiLepJetCut = DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
        				       systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;
    if(systEvent.njets_ >= 1) dPhiDiLepJetCut = DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;
    if((channel > 1000 && channel < 2000) || channel == 300) {passMET = usedMet > 30.; if(systEvent.type_ == SmurfTree::ee) passMET = usedMet > 50.;}

    int centrality = 0;
    if(((systEvent.jet1_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep1_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep1_.Eta() < 0)) &&
       ((systEvent.jet1_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep2_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep2_.Eta() < 0))) centrality = 1; 
    bool passCuts = false;
    bool isRealLepton = false;
    if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
       (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) isRealLepton = true;
    //double HTMiss = (systEvent.lep1_+systEvent.lep2_+systEvent.jet1_+systEvent.jet2_+systEvent.jet3_+systEvent.jet4_).Pt();
    //double LD = systEvent.met_ * 0.00397 + HTMiss * 0.00265;

    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
          systEvent.dilep_.M() > 12 && systEvent.dilep_.M() < 200 &&
         (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 systEvent.njets_ <= 3. && systEvent.njets_ >= 2. &&
         systEvent.lep1_.Pt() > 20. &&
         systEvent.lep2_.Pt() > 10. &&
         (systEvent.jet1_.Pt() <= 30 || TMath::Abs(systEvent.jet1_.Eta()) < 2.5) && (systEvent.jet2_.Pt() <= 30 || TMath::Abs(systEvent.jet2_.Eta()) < 2.5) && (systEvent.jet3_.Pt() <= 30 || TMath::Abs(systEvent.jet3_.Eta()) < 2.5) &&
         passMET == true &&
         systEvent.mt_ > 0. && systEvent.mt_ < 200. && 
         (fabs(systEvent.dilep_.M()-91.1876) > 10. || systEvent.type_ != SmurfTree::ee) && 
         (dPhiDiLepJetCut == true || systEvent.type_ != SmurfTree::ee) && 
         (systEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (systEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	
	passCuts = true;
	if(isRealLepton == false &&
	   (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
	    systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
	    systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // HW->2l selection
    else if(channel == 300 && channel != 301){ // WW->2l selection
      if(
         systEvent.dilep_.M() > 200 &&
        (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 systEvent.njets_ >= 2 &&
         systEvent.lep1_.Pt() > 20. &&
         systEvent.lep2_.Pt() > 20. &&
         passMET == true &&
	 (systEvent.jet1_+systEvent.jet2_).M() > 500 &&
	 TMath::Abs(systEvent.jet1_.Eta()-systEvent.jet2_.Eta()) > 3.5 &&
	 centrality == 1 &&
         (fabs(systEvent.dilep_.M()-91.1876) > 10. || systEvent.type_ != SmurfTree::ee) && 
         //(dPhiDiLepJetCut == true || systEvent.type_ != SmurfTree::ee) && 
         (systEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (systEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	if(isRealLepton == false &&
	   (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
	    systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
	    systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // WW->2l selection
    else if(channel == 301){ // ttH selection
      int nBtag[3] = {0, 0, 0};
      if(systEvent.jet1Btag_     > 2.100) nBtag[0]++; if(systEvent.jet2Btag_     > 2.100) nBtag[0]++; if(systEvent.jet3Btag_     > 2.100) nBtag[0]++; if(systEvent.jet4Btag_     > 2.100) nBtag[0]++;
      if(systEvent.jet1ProbBtag_ > 0.244) nBtag[1]++; if(systEvent.jet2ProbBtag_ > 0.244) nBtag[1]++; if(systEvent.jet3ProbBtag_ > 0.244) nBtag[1]++; if(systEvent.jet4ProbBtag_ > 0.244) nBtag[1]++;
      if(systEvent.jet1ProbBtag_ > 0.679) nBtag[2]++; if(systEvent.jet2ProbBtag_ > 0.679) nBtag[2]++; if(systEvent.jet3ProbBtag_ > 0.679) nBtag[2]++; if(systEvent.jet4ProbBtag_ > 0.679) nBtag[2]++;
      if(
         systEvent.dilep_.M() > 200 &&
        (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
         systEvent.lep1_.Pt() > 20. &&
         systEvent.lep2_.Pt() > 20. &&
	 //systEvent.jet4_.Pt() > 25. && LD > 0.2 && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 systEvent.njets_ >= 1 && systEvent.met_ > 30. && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 systEvent.lep1_.Pt()+systEvent.lep2_.Pt()+systEvent.met_ > 100.0 &&
	 fabs(systEvent.jet1_.Eta()) < 2.5 && fabs(systEvent.jet2_.Eta()) < 2.5 &&
         (fabs(systEvent.dilep_.M()-91.1876) > 10. || systEvent.type_ != SmurfTree::ee) && 
	 (systEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	passCuts = true;
	if(isRealLepton == false &&
	   (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
	    systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
	    systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) passCuts = false;
      }
    } // ttH selection

    if(passCuts == true){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake < 0) assert(0);
 
      if(nFake > 1){
	add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	fDecay = 22;
	theWeight	       = -1.0*add;
      }
      else if(nFake == 1){
        if(systEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          if(fCheckProblem == true && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	}
	else if(isRealLepton == true || systEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,systEvent.npu_);
          add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	  add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt(), 
								   fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						   TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
          add = add*trigEff;
	  if(fCheckProblem == true && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  fDecay                 = 3;
	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * systEvent.scale1fb_*lumi*add;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(systEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
	        					        fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
							        TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        double sf_eff = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_)*
        	        leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
        theWeight = ZttScaleFactor(period,systEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(systEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,systEvent.npu_);
        double add2 = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	add2 = add2*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
								 fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        add = add1*add2*trigEff;

        if(fCheckProblem == true && add != 0 && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	  printf("PROBLEMC(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",systEvent.event_,add1,add2,trigEff,add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);

	if(systEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(systEvent.type_,systEvent.met_);

        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (systEvent.dstype_ == SmurfTree::dymm || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dytt) &&
          (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me)) add = 0.0;

        //----------------------------------------------------------------------------      
        // Apply weighting factor to wgamma (gamma->electron ratio)
        //----------------------------------------------------------------------------
        if(systEvent.dstype_ == SmurfTree::wgamma) {
          if(!(TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep1_.Eta()));
          if(!(TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep2_.Eta()));      
        }

	theWeight              = systEvent.scale1fb_*lumi*add;
      }

      double addLepEff     = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_, 0)*
        		     leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_, 0);
      double addLepEffUp   = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_, 1)*
        		     leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_, 1);
      double addLepEffDown = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_,-1)*
        		     leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_,-1);
      double NjetSyst[2] = {0., 0.};
      if(systEvent.jet1_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(systEvent.jet2_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(systEvent.jet3_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(systEvent.jet4_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(systEvent.jet1_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(systEvent.jet2_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(systEvent.jet3_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(systEvent.jet4_.Pt()*0.95 > 30) NjetSyst[1]++;

      if(thePlot ==62) {
	bdtg      = Unroll2VarTo1ForWH2l(systEvent.dilep_.M(), systEvent.mt_);
	bdtg_aux0 = Unroll2VarTo1ForWH2l(mll_lepup  , mt_lepup  );
	bdtg_aux1 = Unroll2VarTo1ForWH2l(mll_lepdown, mt_lepdown);
	bdtg_aux2 = Unroll2VarTo1ForWH2l(mll_metup  , mt_metup  );
      }

      if     (fDecay == 27){
	histo_WZ_CMS_WZNLOBoundingUp->Fill((bdtg+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_,theWeight);
      }
      else if(fDecay == 19){
        histo_Wgamma                          ->Fill((bdtg+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight);
        histo_Wgamma_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight*addLepEffUp  /addLepEff);
        histo_Wgamma_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight*addLepEffDown/addLepEff);
        histo_Wgamma_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight);
        histo_Wgamma_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight);
        histo_Wgamma_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(systEvent.njets_)+3.0*(double)systEvent.type_, theWeight);
        histo_Wgamma_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)systEvent.type_, theWeight);
        histo_Wgamma_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)systEvent.type_, theWeight);
      }

    }
  } // end syst loop
  }

  int nSig=sigEvent.tree_->GetEntries();
  for (int i=0; i<nSig; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
       printf("--- reading Signal event %5d of %5d\n",i,nSig);
    sigEvent.tree_->GetEntry(i);
    if(channel != 300 && channel != 301 && channel != 1125){
      sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww"     ,(int)125), &bdtg);
      sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux0"  ,(int)125), &bdtg_aux0);
      sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux1"  ,(int)125), &bdtg_aux1);
      sigEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux2"  ,(int)125), &bdtg_aux2);
      sigEvent.tree_->SetBranchAddress( "mll_lepup"      , &mll_lepup      );
      sigEvent.tree_->SetBranchAddress( "mt_lepup"       , &mt_lepup       );
      sigEvent.tree_->SetBranchAddress( "mll_lepdown"    , &mll_lepdown    );
      sigEvent.tree_->SetBranchAddress( "mt_lepdown"     , &mt_lepdown     );
      sigEvent.tree_->SetBranchAddress( "mll_metup"      , &mll_metup      );
      sigEvent.tree_->SetBranchAddress( "mt_metup"       , &mt_metup       );
    }
    bool lId = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    int charge = (int)(sigEvent.lq1_ + sigEvent.lq2_);

    double deltaPhiQQL[3] = {DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep1_.Phi()),DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
    bool   passMET = usedMet > 30.;
    double ptww[3] = {sigEvent.met_*cos(sigEvent.metPhi_)+sigEvent.lep1_.Px()+sigEvent.lep2_.Px(),
                      sigEvent.met_*sin(sigEvent.metPhi_)+sigEvent.lep1_.Py()+sigEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);
    bool dPhiDiLepJetCut = true;
    if(sigEvent.njets_ <= 1) dPhiDiLepJetCut = sigEvent.jet1_.Pt() <= 15. || sigEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
   					       sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
    else 		    dPhiDiLepJetCut = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
   					       sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
    if(sigEvent.njets_ >= 1) dPhiDiLepJetCut = DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
   							 sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me;
    if((channel > 1000 && channel < 2000) || channel == 300) {passMET = usedMet > 30.; if(sigEvent.type_ == SmurfTree::ee) passMET = usedMet > 50.;}
    double MJJL[2] = {0.0,0.0};
    if(DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep1_.Phi()) > DeltaPhi((sigEvent.jet1_+sigEvent.jet2_).Phi(),sigEvent.lep2_.Phi())) {
      MJJL[0] =  (sigEvent.jet1_+sigEvent.jet2_+sigEvent.lep1_).M();
      MJJL[1] =  (sigEvent.jet1_+sigEvent.jet2_+sigEvent.lep2_).M();
    } else {
      MJJL[0] =  (sigEvent.jet1_+sigEvent.jet2_+sigEvent.lep2_).M();
      MJJL[1] =  (sigEvent.jet1_+sigEvent.jet2_+sigEvent.lep1_).M();
    }

    int centrality = 0;
    if(((sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() < 0)) &&
       ((sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() < 0))) centrality = 1; 
    bool passCuts = false;
    bool isSignalDecay = false;
    bool isRealLepton = false;
    if((TMath::Abs(sigEvent.lep1McId_) == 11 || TMath::Abs(sigEvent.lep1McId_) == 13) &&
       (TMath::Abs(sigEvent.lep2McId_) == 11 || TMath::Abs(sigEvent.lep2McId_) == 13)) isRealLepton = true;
    double HTMiss = (sigEvent.lep1_+sigEvent.lep2_+sigEvent.jet1_+sigEvent.jet2_+sigEvent.jet3_+sigEvent.jet4_).Pt();
    //double LD = sigEvent.met_ * 0.00397 + HTMiss * 0.00265;

    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
          sigEvent.dilep_.M() > 12 && sigEvent.dilep_.M() < 200 &&
         (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 sigEvent.njets_ <= 3. && sigEvent.njets_ >= 2. &&
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
        if(isRealLepton == false) {passCuts = false; isSignalDecay = false;}
        if(sigEvent.processId_ == 10001 || sigEvent.processId_ == 10010) {passCuts = false; isSignalDecay = false;}
      }
    } // HW->2l selection
    else if(channel == 300){ // WW->2l selection
      if(
         sigEvent.dilep_.M() > 200 &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 sigEvent.njets_ >= 2 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 20. &&
         passMET == true &&
	 (sigEvent.jet1_+sigEvent.jet2_).M() > 500 &&
	 TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta()) > 3.5 &&
	 centrality == 1 &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 10. || sigEvent.type_ != SmurfTree::ee) && 
         //(dPhiDiLepJetCut == true || sigEvent.type_ != SmurfTree::ee) && 
         (sigEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
        passCuts = true;
        isSignalDecay = true;
        if(isRealLepton == false) {passCuts = false; isSignalDecay = false;}
        //if(sigEvent.processId_ == 10001 || sigEvent.processId_ == 10010) {passCuts = false; isSignalDecay = false;}
      }
    } // WW->2l selection
    else if(channel == 301){ // ttH selection
      int nBtag[3] = {0, 0, 0};
      if(sigEvent.jet1Btag_     > 2.100) nBtag[0]++; if(sigEvent.jet2Btag_     > 2.100) nBtag[0]++; if(sigEvent.jet3Btag_     > 2.100) nBtag[0]++; if(sigEvent.jet4Btag_     > 2.100) nBtag[0]++;
      if(sigEvent.jet1ProbBtag_ > 0.244) nBtag[1]++; if(sigEvent.jet2ProbBtag_ > 0.244) nBtag[1]++; if(sigEvent.jet3ProbBtag_ > 0.244) nBtag[1]++; if(sigEvent.jet4ProbBtag_ > 0.244) nBtag[1]++;
      if(sigEvent.jet1ProbBtag_ > 0.679) nBtag[2]++; if(sigEvent.jet2ProbBtag_ > 0.679) nBtag[2]++; if(sigEvent.jet3ProbBtag_ > 0.679) nBtag[2]++; if(sigEvent.jet4ProbBtag_ > 0.679) nBtag[2]++;
      if(
         sigEvent.dilep_.M() > 200 &&
        (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 20. &&
	 //sigEvent.jet4_.Pt() > 25. && LD > 0.2 && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 sigEvent.njets_ >= 1 && sigEvent.met_ > 30. && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 sigEvent.lep1_.Pt()+sigEvent.lep2_.Pt()+sigEvent.met_ > 100.0 &&
	 fabs(sigEvent.jet1_.Eta()) < 2.5 && fabs(sigEvent.jet2_.Eta()) < 2.5 &&
         (fabs(sigEvent.dilep_.M()-91.1876) > 10. || sigEvent.type_ != SmurfTree::ee) && 
	 (sigEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
        passCuts = true;
        isSignalDecay = true;
        if(isRealLepton == false) {passCuts = false; isSignalDecay = false;}
      }
    } // ttH selection

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

      double myWeight = sigEvent.scale1fb_*lumi*add;

      if(fCheckProblem == true && TMath::Abs((sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_)-add)/add>0.0001) {
	printf("PROBLEM: %f - %f %f %f %f %f = %f\n",add,sigEvent.sfWeightFR_,sigEvent.sfWeightPU_,sigEvent.sfWeightEff_,sigEvent.sfWeightTrig_,sigEvent.sfWeightHPt_,sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_*sigEvent.sfWeightHPt_);
	printf("  		%f %f %f %f\n",1.0,addPU,leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_)*
  	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_), trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  																			     fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
  																			     TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_)));
      }
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
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = sigEvent.met_/sigEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = sigEvent.njets_;
      else if(thePlot ==18) myVar = sigEvent.nvtx_;
      else if(thePlot ==19) myVar = HTMiss;
      else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(sigEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(sigEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(sigEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(sigEvent.lid1_)==13) myVar = sigEvent.lep1_.Eta();else if(fabs(sigEvent.lid2_)==13) myVar = sigEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = TMath::Max(sigEvent.jet1Btag_,sigEvent.jet2Btag_);
      else if(thePlot ==37) myVar = (sigEvent.jet1_+sigEvent.jet2_).M();
      else if(thePlot ==38) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
      else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,sigEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==44) myVar = sigEvent.jet1_.Pt()+ sigEvent.jet2_.Pt()+sigEvent.jet3_.Pt();
      else if(thePlot ==48) myVar = sigEvent.type_;
      else if(thePlot ==49) myVar = (sigEvent.jet1_.Pt()*sigEvent.jet2_.Pt()*sigEvent.lep1_.Pt()*sigEvent.lep2_.Pt())/10000000.0;
      else if(thePlot ==50) myVar = TMath::Min(ptww[2],199.999);
      else if(thePlot ==51) myVar = TMath::Min(TMath::Min(MJJL[0],MJJL[1]),199.999);
      else if(thePlot ==52) myVar = TMath::Min(MJJL[1],199.999);
      else if(thePlot ==53) myVar = DeltaPhi(sigEvent.jet3_.Phi(),sigEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = sigEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==57) myVar = sigEvent.dR_;
      else if(thePlot ==60) myVar = DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.jet2_.Phi(),sigEvent.jet2_.Eta());
      else if(thePlot ==61) myVar = (sigEvent.lep1_+sigEvent.lep2_+sigEvent.jet1_+sigEvent.jet2_).M();
      else if(thePlot ==62) {
        myVar     = Unroll2VarTo1ForWH2l(sigEvent.dilep_.M(), sigEvent.mt_)+(sigEvent.njets_)+3.0*(double)sigEvent.type_;
	bdtg      = Unroll2VarTo1ForWH2l(sigEvent.dilep_.M(), sigEvent.mt_);
	bdtg_aux0 = Unroll2VarTo1ForWH2l(mll_lepup  , mt_lepup  );
	bdtg_aux1 = Unroll2VarTo1ForWH2l(mll_lepdown, mt_lepdown);
	bdtg_aux2 = Unroll2VarTo1ForWH2l(mll_metup  , mt_metup  );
      }
      else if(thePlot ==63) myVar = (bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_;

      int nSigBin = -1;
      if     (sigEvent.processId_==121 ||
  	      sigEvent.processId_==122)	nSigBin = 1;
      else if(sigEvent.processId_==24)	nSigBin = 2;
      //else if(sigEvent.processId_==26)	nSigBin = 3;
      //else assert(0);
      else  nSigBin = 3;
      nSigCut[0]  = nSigCut[0]  + myWeight;
      nSigECut[0] = nSigECut[0] + myWeight*myWeight;
      nSigCut[nSigBin]  = nSigCut[nSigBin]  + myWeight;
      nSigECut[nSigBin] = nSigECut[nSigBin] + myWeight*myWeight;

      double addLepEff     = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 0)*
        		     leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 0);
      double addLepEffUp   = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 1)*
        		     leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 1);
      double addLepEffDown = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,-1)*
        		     leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,-1);
      double NjetSyst[2] = {0., 0.};
      if(sigEvent.jet1_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(sigEvent.jet2_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(sigEvent.jet3_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(sigEvent.jet4_.Pt()*1.05 > 30) NjetSyst[0]++;
      if(sigEvent.jet1_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(sigEvent.jet2_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(sigEvent.jet3_.Pt()*0.95 > 30) NjetSyst[1]++;
      if(sigEvent.jet4_.Pt()*0.95 > 30) NjetSyst[1]++;

      histos->Fill(myVar,myWeight);
      if     (sigEvent.processId_==121 || sigEvent.processId_==122) {
        histo_ttH                          ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ttH_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffUp  /addLepEff);
        histo_ttH_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffDown/addLepEff);
        histo_ttH_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ttH_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ttH_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ttH_CMS_MVAJESBoundingUp	   ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)sigEvent.type_, myWeight);
        histo_ttH_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)sigEvent.type_, myWeight);
      }
      else if(sigEvent.processId_==24) {
        histo_ZH                          ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ZH_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffUp  /addLepEff);
        histo_ZH_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffDown/addLepEff);
        histo_ZH_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ZH_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ZH_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_ZH_CMS_MVAJESBoundingUp     ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)sigEvent.type_, myWeight);
        histo_ZH_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)sigEvent.type_, myWeight);
      }
      else if(sigEvent.processId_==26) {
        histo_WH                          ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_WH_CMS_MVALepEffBoundingUp  ->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffUp  /addLepEff);
        histo_WH_CMS_MVALepEffBoundingDown->Fill((bdtg+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight*addLepEffDown/addLepEff);
        histo_WH_CMS_MVALepResBoundingUp  ->Fill((bdtg_aux0+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_WH_CMS_MVALepResBoundingDown->Fill((bdtg_aux1+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_WH_CMS_MVAMETResBoundingUp  ->Fill((bdtg_aux2+1.)/2.+(sigEvent.njets_)+3.0*(double)sigEvent.type_, myWeight);
        histo_WH_CMS_MVAJESBoundingUp	  ->Fill((bdtg+1.)/2.+NjetSyst[0]+3.0*(double)sigEvent.type_, myWeight);
        histo_WH_CMS_MVAJESBoundingDown   ->Fill((bdtg+1.)/2.+NjetSyst[1]+3.0*(double)sigEvent.type_, myWeight);
      } //else {assert(0);};

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
  
  int nData=dataEvent.tree_->GetEntries();
  double nSelectedData = 0;
  for (int i=0; i<nData; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nData);
    dataEvent.tree_->GetEntry(i);
    if(channel != 300 && channel != 301 && channel != 1125){
      dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww"      ,(int)125), &bdtg);
      dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux0"  ,(int)125), &bdtg_aux0);
      dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux1"  ,(int)125), &bdtg_aux1);
      dataEvent.tree_->SetBranchAddress(Form("bdtg_hww%i_999jet_ww_aux2"  ,(int)125), &bdtg_aux2);
      dataEvent.tree_->SetBranchAddress( "mll_lepup"      , &mll_lepup      );
      dataEvent.tree_->SetBranchAddress( "mt_lepup"       , &mt_lepup       );
      dataEvent.tree_->SetBranchAddress( "mll_lepdown"    , &mll_lepdown    );
      dataEvent.tree_->SetBranchAddress( "mt_lepdown"     , &mt_lepdown     );
      dataEvent.tree_->SetBranchAddress( "mll_metup"      , &mll_metup      );
      dataEvent.tree_->SetBranchAddress( "mt_metup"       , &mt_metup       );
    }
    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;
    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    double deltaPhiQQL[3] = {DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep1_.Phi()),DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep2_.Phi()),0.0};
    deltaPhiQQL[2] = TMath::Min(deltaPhiQQL[0],deltaPhiQQL[1]);
    double usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool   passMET = usedMet > 30.;
    double ptww[3] = {dataEvent.met_*cos(dataEvent.metPhi_)+dataEvent.lep1_.Px()+dataEvent.lep2_.Px(),
                      dataEvent.met_*sin(dataEvent.metPhi_)+dataEvent.lep1_.Py()+dataEvent.lep2_.Py(),0.0};
    ptww[2] = sqrt(ptww[0]*ptww[0]+ptww[1]*ptww[1]);
    bool dPhiDiLepJetCut = true;
    if(dataEvent.njets_ <= 1) dPhiDiLepJetCut = dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
    					       dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    else		     dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
        				       dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    if(dataEvent.njets_ >= 1) dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    if((channel > 1000 && channel < 2000) || channel == 300) {passMET = usedMet > 30.; if(dataEvent.type_ == SmurfTree::ee) passMET = usedMet > 50.;}
    double MJJL[2] = {0.0,0.0};
    if(DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep1_.Phi()) > DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.lep2_.Phi())) {
      MJJL[0] =  (dataEvent.jet1_+dataEvent.jet2_+dataEvent.lep1_).M();
      MJJL[1] =  (dataEvent.jet1_+dataEvent.jet2_+dataEvent.lep2_).M();
    } else {
      MJJL[0] =  (dataEvent.jet1_+dataEvent.jet2_+dataEvent.lep2_).M();
      MJJL[1] =  (dataEvent.jet1_+dataEvent.jet2_+dataEvent.lep1_).M();
    }

    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    double HTMiss = (dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_+dataEvent.jet3_+dataEvent.jet4_).Pt();
    //double LD = dataEvent.met_ * 0.00397 + HTMiss * 0.00265;

    bool passCuts = false;
    if(channel > 1000 && channel < 2000){ // HW->2l selection
      if(
          dataEvent.dilep_.M() > 12 && dataEvent.dilep_.M() < 200 &&
         (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 dataEvent.njets_ <= 3. && dataEvent.njets_ >= 2. &&
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
    else if(channel == 300){ // WW->2l selection
      if(
         dataEvent.dilep_.M() > 200 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
	 dataEvent.njets_ >= 2 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20. &&
         passMET == true &&
	 (dataEvent.jet1_+dataEvent.jet2_).M() > 500 &&
	 TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 &&
	 centrality == 1 &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
         //(dPhiDiLepJetCut == true || dataEvent.type_ != SmurfTree::ee) && 
         (dataEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
        passCuts = true;
      }
    } // WW->2l selection
    else if(channel == 301){ // ttH selection
      int nBtag[3] = {0, 0, 0};
      if(dataEvent.jet1Btag_     > 2.100) nBtag[0]++; if(dataEvent.jet2Btag_     > 2.100) nBtag[0]++; if(dataEvent.jet3Btag_     > 2.100) nBtag[0]++; if(dataEvent.jet4Btag_     > 2.100) nBtag[0]++;
      if(dataEvent.jet1ProbBtag_ > 0.244) nBtag[1]++; if(dataEvent.jet2ProbBtag_ > 0.244) nBtag[1]++; if(dataEvent.jet3ProbBtag_ > 0.244) nBtag[1]++; if(dataEvent.jet4ProbBtag_ > 0.244) nBtag[1]++;
      if(dataEvent.jet1ProbBtag_ > 0.679) nBtag[2]++; if(dataEvent.jet2ProbBtag_ > 0.679) nBtag[2]++; if(dataEvent.jet3ProbBtag_ > 0.679) nBtag[2]++; if(dataEvent.jet4ProbBtag_ > 0.679) nBtag[2]++;
      if(
         dataEvent.dilep_.M() > 200 &&
        (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto &&
         charge != 0 &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20. &&
	 //dataEvent.jet4_.Pt() > 25. && LD > 0.2 && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 dataEvent.njets_ >= 1 && dataEvent.met_ > 30. && (nBtag[1] >= 2 && nBtag[2] >= 1) &&
	 dataEvent.lep1_.Pt()+dataEvent.lep2_.Pt()+dataEvent.met_ > 100.0 &&
	 fabs(dataEvent.jet1_.Eta()) < 2.5 && fabs(dataEvent.jet2_.Eta()) < 2.5 &&
         (fabs(dataEvent.dilep_.M()-91.1876) > 10. || dataEvent.type_ != SmurfTree::ee) && 
	 (dataEvent.type_ == lDecay || lDecay == 4 || (lDecay == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) || (lDecay == 6 && (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me))) &&
	 1 == 1
	){
	if(dataEvent.event_!=543549984&&dataEvent.event_!=1140495707&&dataEvent.event_!=290847449&&dataEvent.event_!=76088896&&dataEvent.event_!=691685534)
	passCuts = true;
      }
    } // ttH selection

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
      else if(thePlot ==12) myVar = usedMet;
      else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
      else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
      else if(thePlot ==15) myVar = dataEvent.met_/dataEvent.dilep_.Pt()/2.0;
      else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
      else if(thePlot ==17) myVar = dataEvent.njets_;
      else if(thePlot ==18) myVar = dataEvent.nvtx_;
      else if(thePlot ==19) myVar = HTMiss;
      else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = deltaPhiQQL[2]*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = fabs(dataEvent.dilep_.Eta());
      else if(thePlot ==24) myVar = fabs(dataEvent.lep1_.Eta());
      else if(thePlot ==25) myVar = fabs(dataEvent.lep2_.Eta());
      else if(thePlot ==26) {if(fabs(dataEvent.lid1_)==13) myVar = dataEvent.lep1_.Eta();else if(fabs(dataEvent.lid2_)==13) myVar = dataEvent.lep2_.Eta(); else myVar=9;}
      else if(thePlot ==27) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==28) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==29) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
      else if(thePlot ==30) myVar = TMath::Max(dataEvent.jet1Btag_,dataEvent.jet2Btag_);
      else if(thePlot ==37) myVar = (dataEvent.jet1_+dataEvent.jet2_).M();
      else if(thePlot ==38) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
      else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.jet2_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==44) myVar = dataEvent.jet1_.Pt()+ dataEvent.jet2_.Pt()+dataEvent.jet3_.Pt();
      else if(thePlot ==48) myVar = dataEvent.type_;
      else if(thePlot ==49) myVar = (dataEvent.jet1_.Pt()*dataEvent.jet2_.Pt()*dataEvent.lep1_.Pt()*dataEvent.lep2_.Pt())/10000000.0;
      else if(thePlot ==50) myVar = TMath::Min(ptww[2],199.999);
      else if(thePlot ==51) myVar = TMath::Min(TMath::Min(MJJL[0],MJJL[1]),199.999);
      else if(thePlot ==52) myVar = TMath::Min(MJJL[1],199.999);
      else if(thePlot ==53) myVar = DeltaPhi(dataEvent.jet3_.Phi(),dataEvent.jet4_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==55) myVar = dataEvent.dPhiDiLepMET_*180.0/TMath::Pi();
      else if(thePlot ==57) myVar = dataEvent.dR_;
      else if(thePlot ==60) myVar = DeltaR(dataEvent.jet1_.Phi(),dataEvent.jet1_.Eta(),dataEvent.jet2_.Phi(),dataEvent.jet2_.Eta());
      else if(thePlot ==61) myVar = (dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_).M();
      else if(thePlot ==62) {
        myVar     = Unroll2VarTo1ForWH2l(dataEvent.dilep_.M(), dataEvent.mt_)+(dataEvent.njets_)+3.0*(double)dataEvent.type_;
      }
      else if(thePlot ==63) myVar = (bdtg+1.)/2.+(dataEvent.njets_)+3.0*(double)dataEvent.type_;

      double wei = 1.0;
      nSelectedData = nSelectedData + wei;
      histo5->Fill(myVar,1.0*wei);
      histo_Data->Fill(myVar,1.0*wei);
    }
  } // End loop data

  printf("data: %f\n",nSelectedData);

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
    double totalB   = bgdDecay[21]+bgdDecay[27]+bgdDecay[28]+bgdDecay[29]+bgdDecay[30]+bgdCombined[0]+bgdCombined[1]+bgdCombined[2]+bgdCombined[3];
    double totalB_E = sqrt(weiDecay[21]+weiDecay[27]+weiDecay[28]+weiDecay[29]+weiDecay[27]+weiDecay[30]+bgdCombinedE[0]*bgdCombinedE[0]+bgdCombinedE[1]*bgdCombinedE[1]+bgdCombinedE[2]*bgdCombinedE[2]+bgdCombinedE[3]*bgdCombinedE[3]);
    printf(" Higgs           & data & SumBkg           & qqWW+ggWW        & Top              & Wjets            & VV               & Zjets            & Vgamma(*)        \\\\\n"); 
    printf("%6.2f \\pm %5.2f & %3d  & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f & %6.2f \\pm %5.2f \\\\\n",
           nSigCut[0],nSigECut[0],(int)nSelectedData,totalB,totalB_E,bgdDecay[29]+bgdDecay[30],sqrt(weiDecay[29]+weiDecay[30]),bgdCombined[2],bgdCombinedE[2],bgdCombined[0],bgdCombinedE[0],
	   bgdDecay[27]+bgdDecay[28]+bgdDecay[21],sqrt(weiDecay[27]+weiDecay[28]+weiDecay[21]),bgdCombined[1],bgdCombinedE[1],bgdDecay[19]+bgdDecay[20],sqrt(weiDecay[19]+weiDecay[20]));

    bgdDecay[27] = bgdDecay[27];
    weiDecay[27] = sqrt(weiDecay[27]);
    if(bgdDecay[27] > 0) weiDecay[27] = weiDecay[27]/bgdDecay[27];
    if(bgdDecay[29] > 0) weiDecay[29] = weiDecay[29]/bgdDecay[29];
    if(bgdDecay[30] > 0) weiDecay[30] = weiDecay[30]/bgdDecay[30];
    if(bgdDecay[21] > 0) weiDecay[21] = weiDecay[21]/bgdDecay[21];
    if(bgdCombined[0] > 0) bgdCombinedE[0] =  bgdCombinedE[0] / bgdCombined[0];
    if(bgdCombined[1] > 0) bgdCombinedE[1] =  bgdCombinedE[1] / bgdCombined[1];
    if(bgdCombined[2] > 0) bgdCombinedE[2] =  bgdCombinedE[2] / bgdCombined[2];
    if(bgdCombined[3] > 0) bgdCombinedE[3] =  bgdCombinedE[3] / bgdCombined[3];
    histo_Wgamma_CMS_MVALepEffBoundingUp  ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVALepEffBoundingDown->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVALepResBoundingUp  ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVALepResBoundingDown->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVAMETResBoundingUp  ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVAJESBoundingUp	  ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    histo_Wgamma_CMS_MVAJESBoundingDown   ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());	
    histo_Wgamma                          ->Scale(bgdDecay[19]/histo_Wgamma->GetSumOfWeights());
    for(int i=1; i<=histo_ttH->GetNbinsX(); i++){
      double factorUp = +1.0; double factorDown = -1.0;
      histo_ttH_CMS_MVAttHStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_ttH   ->GetBinContent(i)+factorUp  *histo_ttH   ->GetBinError(i),0.000001));
      histo_ttH_CMS_MVAttHStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_ttH   ->GetBinContent(i)+factorDown*histo_ttH   ->GetBinError(i),0.000001));
      histo_ZH_CMS_MVAZHStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_ZH    ->GetBinContent(i)+factorUp  *histo_ZH    ->GetBinError(i),0.000001));
      histo_ZH_CMS_MVAZHStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_ZH    ->GetBinContent(i)+factorDown*histo_ZH    ->GetBinError(i),0.000001));
      histo_WH_CMS_MVAWHStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_WH    ->GetBinContent(i)+factorUp  *histo_WH    ->GetBinError(i),0.000001));
      histo_WH_CMS_MVAWHStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WH    ->GetBinContent(i)+factorDown*histo_WH    ->GetBinError(i),0.000001));
      histo_WZ_CMS_MVAWZStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_WZ    ->GetBinContent(i)+factorUp  *histo_WZ    ->GetBinError(i),0.000001));
      histo_WZ_CMS_MVAWZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WZ    ->GetBinContent(i)+factorDown*histo_WZ    ->GetBinError(i),0.000001));
      histo_WS_CMS_MVAWSStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_WS    ->GetBinContent(i)+factorUp  *histo_WS    ->GetBinError(i),0.000001));
      histo_WS_CMS_MVAWSStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WS    ->GetBinContent(i)+factorDown*histo_WS    ->GetBinError(i),0.000001));
      histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_WjetsE->GetBinContent(i)+factorUp  *histo_WjetsE->GetBinError(i),0.000001));
      histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->SetBinContent(i,TMath::Max(histo_WjetsE->GetBinContent(i)+factorDown*histo_WjetsE->GetBinError(i),0.000001));
      histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_WjetsM->GetBinContent(i)+factorUp  *histo_WjetsM->GetBinError(i),0.000001));
      histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->SetBinContent(i,TMath::Max(histo_WjetsM->GetBinContent(i)+factorDown*histo_WjetsM->GetBinError(i),0.000001));
      histo_Wgamma_CMS_MVAWgammaStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_Wgamma->GetBinContent(i)+factorUp  *histo_Wgamma->GetBinError(i),0.000001));
      histo_Wgamma_CMS_MVAWgammaStatBoundingDown->SetBinContent(i,TMath::Max(histo_Wgamma->GetBinContent(i)+factorDown*histo_Wgamma->GetBinError(i),0.000001));
      histo_Wg3l_CMS_MVAWg3lStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_Wg3l  ->GetBinContent(i)+factorUp  *histo_Wg3l  ->GetBinError(i),0.000001));
      histo_Wg3l_CMS_MVAWg3lStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_Wg3l  ->GetBinContent(i)+factorDown*histo_Wg3l  ->GetBinError(i),0.000001));
      histo_VVV_CMS_MVAVVVStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_VVV   ->GetBinContent(i)+factorUp  *histo_VVV   ->GetBinError(i),0.000001));
      histo_VVV_CMS_MVAVVVStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_VVV   ->GetBinContent(i)+factorDown*histo_VVV   ->GetBinError(i),0.000001));
      histo_ttH_SM_CMS_MVAttH_SMStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_ttH_SM   ->GetBinContent(i)+factorUp  *histo_ttH_SM   ->GetBinError(i),0.000001));
      histo_ttH_SM_CMS_MVAttH_SMStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_ttH_SM   ->GetBinContent(i)+factorDown*histo_ttH_SM   ->GetBinError(i),0.000001));
      histo_ZH_SM_CMS_MVAZH_SMStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_ZH_SM   ->GetBinContent(i)+factorUp  *histo_ZH_SM   ->GetBinError(i),0.000001));
      histo_ZH_SM_CMS_MVAZH_SMStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_ZH_SM   ->GetBinContent(i)+factorDown*histo_ZH_SM   ->GetBinError(i),0.000001));
      histo_WH_SM_CMS_MVAWH_SMStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_WH_SM   ->GetBinContent(i)+factorUp  *histo_WH_SM   ->GetBinError(i),0.000001));
      histo_WH_SM_CMS_MVAWH_SMStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_WH_SM   ->GetBinContent(i)+factorDown*histo_WH_SM   ->GetBinError(i),0.000001));
    }

    for(int i=1; i<=histo_WjetsE->GetNbinsX(); i++){
      if(histo_WjetsE->GetBinContent(i) <= 0) {histo_WjetsE->SetBinContent(i,0.000001);histo_WjetsE->SetBinError(i,0.000001);}
      if(histo_WjetsM->GetBinContent(i) <= 0) {histo_WjetsM->SetBinContent(i,0.000001);histo_WjetsM->SetBinError(i,0.000001);}
    }

    double mean,up,diff;

    if(histo_WjetsE_CMS_MVAWEBoundingUp->GetNbinsX() != histo_WjetsE->GetNbinsX()) {printf("Different binning in W!\n"); return;}
    histo_WjetsE_CMS_MVAWEBoundingUp  ->Scale(histo_WjetsE->GetSumOfWeights()/histo_WjetsE_CMS_MVAWEBoundingUp  ->GetSumOfWeights());
    for(int i=1; i<=histo_WjetsE_CMS_MVAWEBoundingUp->GetNbinsX(); i++){
      mean = histo_WjetsE		    ->GetBinContent(i);
      up   = histo_WjetsE_CMS_MVAWEBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WjetsE_CMS_MVAWEBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WjetsE_CMS_MVAWEBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    }
    histo_WjetsE_CMS_MVAWEBoundingDown->Scale(histo_WjetsE->GetSumOfWeights()/histo_WjetsE_CMS_MVAWEBoundingDown->GetSumOfWeights());

    if(histo_WjetsM_CMS_MVAWMBoundingUp->GetNbinsX() != histo_WjetsM->GetNbinsX()) {printf("Different binning in W!\n"); return;}
    histo_WjetsM_CMS_MVAWMBoundingUp  ->Scale(histo_WjetsM->GetSumOfWeights()/histo_WjetsM_CMS_MVAWMBoundingUp  ->GetSumOfWeights());
    for(int i=1; i<=histo_WjetsM_CMS_MVAWMBoundingUp->GetNbinsX(); i++){
      mean = histo_WjetsM		    ->GetBinContent(i);
      up   = histo_WjetsM_CMS_MVAWMBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WjetsM_CMS_MVAWMBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WjetsM_CMS_MVAWMBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    }
    histo_WjetsM_CMS_MVAWMBoundingDown->Scale(histo_WjetsM->GetSumOfWeights()/histo_WjetsM_CMS_MVAWMBoundingDown->GetSumOfWeights());

    if(histo_WZ_CMS_WZNLOBoundingUp->GetSumOfWeights()>0) printf("nuisance WZ/Wg: %f %f\n",histo_WZ_CMS_WZNLOBoundingUp->GetSumOfWeights(),histo_Wgamma->GetSumOfWeights());
    histo_WZ_CMS_WZNLOBoundingUp->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingUp->GetSumOfWeights());
    for(int i=1; i<=histo_ttH_CMS_MVAMETResBoundingUp->GetNbinsX(); i++){
      mean = histo_ttH  		          ->GetBinContent(i);
      up   = histo_ttH_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_ttH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_ttH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_ZH  		                 ->GetBinContent(i);
      up   = histo_ZH_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_ZH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_ZH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_WH  		                 ->GetBinContent(i);
      up   = histo_WH_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_WZ  		                 ->GetBinContent(i);
      up   = histo_WZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_WS  		                 ->GetBinContent(i);
      up   = histo_WS_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WS_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WS_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_Wgamma  		             ->GetBinContent(i);
      up   = histo_Wgamma_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_WH  		                 ->GetBinContent(i);
      up   = histo_WH_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WH_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_Wg3l  		           ->GetBinContent(i);
      up   = histo_Wg3l_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_Wg3l_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_Wg3l_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_VVV  		                 ->GetBinContent(i);
      up   = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_ttH_SM  		             ->GetBinContent(i);
      up   = histo_ttH_SM_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_ttH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_ttH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    
      mean = histo_ZH_SM  		            ->GetBinContent(i);
      up   = histo_ZH_SM_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_ZH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_ZH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

      mean = histo_WH_SM  		        ->GetBinContent(i);
      up   = histo_WH_SM_CMS_MVAMETResBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WH_SM_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));    
	
      mean = histo_WZ	     		 ->GetBinContent(i);
      up   = histo_WZ_CMS_WZNLOBoundingUp->GetBinContent(i);
      diff = TMath::Abs(mean-up);
      if     (mean-up >0) histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
      else		  histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
    }
    histo_WZ_CMS_WZNLOBoundingDown->Scale(histo_WZ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingDown->GetSumOfWeights());

    //----------------------------------------------------------------------------
    // Produce output cards for cut-based analysis
    //----------------------------------------------------------------------------
    char outputLimits[200];
    sprintf(outputLimits,"vhss%s_%3d.input_8TeV.root",finalStateName,mH);
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    histo_Data  ->Write();
    histo_ttH	->Write();
    histo_ZH	->Write();
    histo_WH	->Write();
    histo_WZ	->Write();
    histo_WS	->Write();
    histo_WjetsE->Write();
    histo_WjetsM->Write();
    histo_Wgamma->Write();
    histo_Wg3l  ->Write();
    histo_VVV	->Write();
    histo_ttH_SM->Write();
    histo_ZH_SM	->Write();
    histo_WH_SM	->Write();

    cout << histo_Data  ->GetSumOfWeights() << " ";
    cout << histo_ttH	->GetSumOfWeights() << " ";
    cout << histo_ZH	->GetSumOfWeights() << " ";
    cout << histo_WH	->GetSumOfWeights() << " ";
    cout << histo_WZ	->GetSumOfWeights() << " ";
    cout << histo_WS	->GetSumOfWeights() << " ";
    cout << histo_WjetsE->GetSumOfWeights() << " ";
    cout << histo_WjetsM->GetSumOfWeights() << " ";
    cout << histo_Wgamma->GetSumOfWeights() << " ";
    cout << histo_Wg3l  ->GetSumOfWeights() << " ";
    cout << histo_VVV	->GetSumOfWeights() << " ";
    cout << histo_ttH_SM->GetSumOfWeights() << " ";
    cout << histo_ZH_SM	->GetSumOfWeights() << " ";
    cout << histo_WH_SM	->GetSumOfWeights() << " ";
    cout << endl;

    histo_ttH_CMS_MVAttHStatBoundingUp->Write();
    histo_ttH_CMS_MVAttHStatBoundingDown->Write();
    histo_ZH_CMS_MVAZHStatBoundingUp->Write();
    histo_ZH_CMS_MVAZHStatBoundingDown->Write();
    histo_WH_CMS_MVAWHStatBoundingUp->Write();
    histo_WH_CMS_MVAWHStatBoundingDown->Write();
    histo_WZ_CMS_MVAWZStatBoundingUp->Write();
    histo_WZ_CMS_MVAWZStatBoundingDown->Write();
    histo_WS_CMS_MVAWSStatBoundingUp->Write();
    histo_WS_CMS_MVAWSStatBoundingDown->Write();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingUp->Write();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Write();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingUp->Write();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Write();
    histo_Wgamma_CMS_MVAWgammaStatBoundingUp->Write();
    histo_Wgamma_CMS_MVAWgammaStatBoundingDown->Write();
    histo_Wg3l_CMS_MVAWg3lStatBoundingUp->Write();
    histo_Wg3l_CMS_MVAWg3lStatBoundingDown->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown->Write();
    histo_ttH_SM_CMS_MVAttH_SMStatBoundingUp->Write();
    histo_ttH_SM_CMS_MVAttH_SMStatBoundingDown->Write();
    histo_ZH_SM_CMS_MVAZH_SMStatBoundingUp->Write();
    histo_ZH_SM_CMS_MVAZH_SMStatBoundingDown->Write();
    histo_WH_SM_CMS_MVAWH_SMStatBoundingUp->Write();
    histo_WH_SM_CMS_MVAWH_SMStatBoundingDown->Write();
    histo_ttH_CMS_MVALepEffBoundingUp->Write();
    histo_ttH_CMS_MVALepEffBoundingDown->Write();
    histo_ZH_CMS_MVALepEffBoundingUp->Write();
    histo_ZH_CMS_MVALepEffBoundingDown->Write();
    histo_WH_CMS_MVALepEffBoundingUp->Write();
    histo_WH_CMS_MVALepEffBoundingDown->Write();
    histo_WZ_CMS_MVALepEffBoundingUp->Write();
    histo_WZ_CMS_MVALepEffBoundingDown->Write();
    histo_WS_CMS_MVALepEffBoundingUp->Write();
    histo_WS_CMS_MVALepEffBoundingDown->Write();
    histo_Wgamma_CMS_MVALepEffBoundingUp->Write();
    histo_Wgamma_CMS_MVALepEffBoundingDown->Write();
    histo_Wg3l_CMS_MVALepEffBoundingUp->Write();
    histo_Wg3l_CMS_MVALepEffBoundingDown->Write();
    histo_VVV_CMS_MVALepEffBoundingUp->Write();
    histo_VVV_CMS_MVALepEffBoundingDown->Write();
    histo_ttH_SM_CMS_MVALepEffBoundingUp->Write();
    histo_ttH_SM_CMS_MVALepEffBoundingDown->Write();
    histo_ZH_SM_CMS_MVALepEffBoundingUp->Write();
    histo_ZH_SM_CMS_MVALepEffBoundingDown->Write();
    histo_WH_SM_CMS_MVALepEffBoundingUp->Write();
    histo_WH_SM_CMS_MVALepEffBoundingDown->Write();
    histo_ttH_CMS_MVALepResBoundingUp->Write();
    histo_ttH_CMS_MVALepResBoundingDown->Write();
    histo_ZH_CMS_MVALepResBoundingUp->Write();
    histo_ZH_CMS_MVALepResBoundingDown->Write();
    histo_WH_CMS_MVALepResBoundingUp->Write();
    histo_WH_CMS_MVALepResBoundingDown->Write();
    histo_WZ_CMS_MVALepResBoundingUp->Write();
    histo_WZ_CMS_MVALepResBoundingDown->Write();
    histo_WS_CMS_MVALepResBoundingUp->Write();
    histo_WS_CMS_MVALepResBoundingDown->Write();
    histo_Wgamma_CMS_MVALepResBoundingUp->Write();
    histo_Wgamma_CMS_MVALepResBoundingDown->Write();
    histo_Wg3l_CMS_MVALepResBoundingUp->Write();
    histo_Wg3l_CMS_MVALepResBoundingDown->Write();
    histo_VVV_CMS_MVALepResBoundingUp->Write();
    histo_VVV_CMS_MVALepResBoundingDown->Write();
    histo_ttH_SM_CMS_MVALepResBoundingUp->Write();
    histo_ttH_SM_CMS_MVALepResBoundingDown->Write();
    histo_ZH_SM_CMS_MVALepResBoundingUp->Write();
    histo_ZH_SM_CMS_MVALepResBoundingDown->Write();
    histo_WH_SM_CMS_MVALepResBoundingUp->Write();
    histo_WH_SM_CMS_MVALepResBoundingDown->Write();
    histo_ttH_CMS_MVAMETResBoundingUp->Write();
    histo_ttH_CMS_MVAMETResBoundingDown->Write();
    histo_ZH_CMS_MVAMETResBoundingUp->Write();
    histo_ZH_CMS_MVAMETResBoundingDown->Write();
    histo_WH_CMS_MVAMETResBoundingUp->Write();
    histo_WH_CMS_MVAMETResBoundingDown->Write();
    histo_WZ_CMS_MVAMETResBoundingUp->Write();
    histo_WZ_CMS_MVAMETResBoundingDown->Write();
    histo_WS_CMS_MVAMETResBoundingUp->Write();
    histo_WS_CMS_MVAMETResBoundingDown->Write();
    histo_Wgamma_CMS_MVAMETResBoundingUp->Write();
    histo_Wgamma_CMS_MVAMETResBoundingDown->Write();
    histo_Wg3l_CMS_MVAMETResBoundingUp->Write();
    histo_Wg3l_CMS_MVAMETResBoundingDown->Write();
    histo_VVV_CMS_MVAMETResBoundingUp->Write();
    histo_VVV_CMS_MVAMETResBoundingDown->Write();
    histo_ttH_SM_CMS_MVAMETResBoundingUp->Write();
    histo_ttH_SM_CMS_MVAMETResBoundingDown->Write();
    histo_ZH_SM_CMS_MVAMETResBoundingUp->Write();
    histo_ZH_SM_CMS_MVAMETResBoundingDown->Write();
    histo_WH_SM_CMS_MVAMETResBoundingUp->Write();
    histo_WH_SM_CMS_MVAMETResBoundingDown->Write();
    histo_ttH_CMS_MVAJESBoundingUp->Write();
    histo_ttH_CMS_MVAJESBoundingDown->Write();
    histo_ZH_CMS_MVAJESBoundingUp->Write();
    histo_ZH_CMS_MVAJESBoundingDown->Write();
    histo_WH_CMS_MVAJESBoundingUp->Write();
    histo_WH_CMS_MVAJESBoundingDown->Write();
    histo_WZ_CMS_MVAJESBoundingUp->Write();
    histo_WZ_CMS_MVAJESBoundingDown->Write();
    histo_WS_CMS_MVAJESBoundingUp->Write();
    histo_WS_CMS_MVAJESBoundingDown->Write();
    histo_Wgamma_CMS_MVAJESBoundingUp->Write();
    histo_Wgamma_CMS_MVAJESBoundingDown->Write();
    histo_Wg3l_CMS_MVAJESBoundingUp->Write();
    histo_Wg3l_CMS_MVAJESBoundingDown->Write();
    histo_VVV_CMS_MVAJESBoundingUp->Write();
    histo_VVV_CMS_MVAJESBoundingDown->Write();
    histo_ttH_SM_CMS_MVAJESBoundingUp->Write();
    histo_ttH_SM_CMS_MVAJESBoundingDown->Write();
    histo_ZH_SM_CMS_MVAJESBoundingUp->Write();
    histo_ZH_SM_CMS_MVAJESBoundingDown->Write();
    histo_WH_SM_CMS_MVAJESBoundingUp->Write();
    histo_WH_SM_CMS_MVAJESBoundingDown->Write();
    histo_WjetsE_CMS_MVAWEBoundingUp->Write();
    histo_WjetsE_CMS_MVAWEBoundingDown->Write();
    histo_WjetsM_CMS_MVAWMBoundingUp->Write();
    histo_WjetsM_CMS_MVAWMBoundingDown->Write();
    histo_WZ_CMS_WZNLOBoundingUp->Write();
    histo_WZ_CMS_WZNLOBoundingDown->Write();

    char thettH_SMString[20];
    if(histo_ttH_SM->GetSumOfWeights() > 0) sprintf(thettH_SMString,"1.0");
    else				    sprintf(thettH_SMString," - ");
    char theZH_SMString[20];
    if(histo_ZH_SM->GetSumOfWeights() > 0) sprintf(theZH_SMString,"1.0");
    else				   sprintf(theZH_SMString," - ");
    char theWH_SMString[20];
    if(histo_WH_SM->GetSumOfWeights() > 0) sprintf(theWH_SMString,"1.0");
    else				   sprintf(theWH_SMString," - ");

    char outputLimitsCut[200];
    sprintf(outputLimitsCut,"histo_limits_whss%2s_mh%3d_cut.txt",finalStateName,mH);
    ofstream newcardShape;
    newcardShape.open(outputLimitsCut);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)nSelectedData);
    newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
    newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
    newcardShape << Form("bin vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s vhss%2s\n",finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName,finalStateName);
    newcardShape << Form("process ttH ZH WH WZ WS WjetsE WjetsM Wgamma Wg3l VVV ttH_SM ZH_SM WH_SM\n");
    newcardShape << Form("process -2 -1 0 1 2 3 4 5 6 7 8 9 10\n");
    newcardShape << Form("rate  %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",nSigCut[1],nSigCut[2],nSigCut[3],bgdDecay[27],bgdDecay[29]+bgdDecay[30]+bgdDecay[5]+bgdDecay[13]+bgdDecay[28]+bgdDecay[9]+bgdDecay[10],bgdDecay[0]+bgdDecay[1]+bgdDecay[2]+bgdDecay[3],bgdDecay[23],bgdDecay[19],bgdDecay[20],bgdDecay[21],bgdDecay[41],bgdDecay[42],bgdDecay[43]);
    newcardShape << Form("lumi_%4s                   lnN %5.3f %5.3f %5.3f %5.3f %5.3f   -     -   %5.3f %5.3f  %5.3f  %5.3f  %5.3f %5.3f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);			   
    newcardShape << Form("QCDscale_VH                lnN  -    1.050 1.050   -     -	 -     -     -     -      -    1.050  1.050 1.050\n");		     
    newcardShape << Form("QCDscale_ttH               lnN 1.070   -     -     -     -     -     -     -     -      -    1.070   -     -   \n");		     
    newcardShape << Form("UEPS                       lnN 1.030 1.030 1.030   -     -     -     -     -     -      -    1.030  1.030 1.030\n");
    newcardShape << Form("pdf_qqbar                  lnN   -   1.040 1.040 1.040 1.040   -     -   1.040 1.040 1.040    -     1.040 1.040\n");
    newcardShape << Form("pdf_gg                     lnN 1.080   -     -     -     -     -     -     -     -      -    1.080   -     -   \n");		     
    newcardShape << Form("CMS_vhss_WZ                lnN   -     -     -   1.100   -     -     -     -     -      -     -      -     -   \n");	
    newcardShape << Form("CMS_vhss_VV                lnN   -     -     -     -   1.200   -     -     -     -      -     -      -     -   \n");	
    newcardShape << Form("CMS_vhss_WjetsE            lnN   -     -     -     -	 -     1.360   -     -     -      -     -      -     -   \n");	
    newcardShape << Form("CMS_vhss_WjetsM            lnN   -     -     -     -   -       -   1.360   -     -      -     -      -     -   \n");	   
    newcardShape << Form("CMS_vhss_Wgamma            lnN   -     -     -     -	 -       -     -    1.300  -      -     -      -     -   \n");	
    newcardShape << Form("CMS_vhss_Wg3l              lnN   -     -     -     -	 -       -     -     -    1.400   -     -      -     -   \n");	
    newcardShape << Form("CMS_vhss_VVV               lnN   -     -     -     -	 -       -     -     -     -     1.500  -      -     -   \n");	
    newcardShape << Form("CMS_vhss_MVALepEffBounding	     shape    1.0   1.0   1.0	1.0    1.0	 -     -     1.0   1.0	1.0    %s     %s    %s\n",thettH_SMString,theZH_SMString,theWH_SMString);
    newcardShape << Form("CMS_vhss_MVALepResBounding	     shape    1.0   1.0   1.0	1.0    1.0	 -     -     1.0   1.0	1.0    %s     %s    %s\n",thettH_SMString,theZH_SMString,theWH_SMString);
    newcardShape << Form("CMS_vhss_MVAMETResBounding	     shape    1.0   1.0   1.0	1.0    1.0	 -     -     1.0   1.0	1.0    %s     %s    %s\n",thettH_SMString,theZH_SMString,theWH_SMString);
    newcardShape << Form("CMS_vhss_MVAJESBounding	     shape    1.0   1.0   1.0	1.0    1.0 	 -     -     1.0   1.0	1.0    %s     %s    %s\n",thettH_SMString,theZH_SMString,theWH_SMString);				
    newcardShape << Form("CMS_vhss_WZNLOBounding             shape    -      -     -	1.0	 -	 -     -     -     -	  -	-      -     -\n");
    newcardShape << Form("CMS_vhss_MVAWEBounding             shape    -      -     -	 -	 -	1.0    -     -     -	  -	-      -     -\n");
    newcardShape << Form("CMS_vhss_MVAWMBounding             shape    -      -     -	 -	 -	 -    1.0    -     -	  -	-      -     -\n");
    newcardShape << Form("CMS_vhss%s_MVAttHStatBounding_%s shape      1.0    -     -	 -	 -	 -     -     -     -	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAZHStatBounding_%s shape        -    1.0    -	 -	 -	 -     -     -     -	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWHStatBounding_%s shape        -     -    1.0    -	 -	 -     -     -     -	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWZStatBounding_%s shape        -     -     -	 1.0	 -	 -     -     -     -	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWSStatBounding_%s shape        -     -     -	 -	1.0	 -     -     -     -	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWjetsEStatBounding_%s shape    -     -     -	 -	 -      1.0    -     -	   -	  -     -      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWjetsMStatBounding_%s shape    -     -     -	 -	 -	 -    1.0    -	   -	  -     -      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWgammaStatBounding_%s shape    -     -     -	 -	 -	 -     -    1.0    -	  -     -      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAWg3lStatBounding_%s shape      -     -     -     -	 -	 -     -     -    1.0	  -	-      -     -\n",finalStateName,ECMsb.Data());
    newcardShape << Form("CMS_vhss%s_MVAVVVStatBounding_%s shape       -     -     -     -	 -	 -     -     -     -	 1.0	-      -     -\n",finalStateName,ECMsb.Data());
    if(histo_ttH_SM->GetSumOfWeights() > 0) 
    newcardShape << Form("CMS_vhss%s_MVAttH_SMStatBounding_%s shape    -     -     -     -	 -	 -     -     -     -	  -	%s     -     -\n",finalStateName,ECMsb.Data(),thettH_SMString);
    if(histo_ZH_SM->GetSumOfWeights() > 0) 
    newcardShape << Form("CMS_vhss%s_MVAZH_SMStatBounding_%s shape     -     -     -     -	 -	 -     -     -     -	  -	-      %s    -\n",finalStateName,ECMsb.Data(),theZH_SMString);
    if(histo_WH_SM->GetSumOfWeights() > 0) 
    newcardShape << Form("CMS_vhss%s_MVAWH_SMStatBounding_%s shape      -     -     -     -	 -       -     -     -     -      -     -      -    %s\n",finalStateName,ECMsb.Data(),theWH_SMString);
    newcardShape.close();
  }
  return;
}
