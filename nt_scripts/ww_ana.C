#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/makeSystematicEffects.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/LeptonEfficiencyZH.h"
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
#include "TLorentzVector.h"

Double_t DYBkgScaleFactor(Int_t jetBin);
Double_t DYBkgScaleFactorKappa(Int_t jetBin);
Double_t TopBkgScaleFactor(Int_t jetBin);
Double_t TopBkgScaleFactorKappa(Int_t jetBin);
double weightJetPt(Int_t nsel, Int_t jetpt0, Int_t jetpt1);

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 1;
const unsigned int nSelTypesSyst = 7;
const bool showSignalOnly = false;
const bool useDYMVA = true;
const bool useWeightEWKCorr = false;
const bool useWeightNNLOCorr = true;

enum selType {WWSEL};
TString selTypeName[nSelTypes*2] = {"WWSEL-SS",
                                    "WWSEL-OS"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-SS", "JESDOWN-SS", "LEPP-SS", "LEPM-SS", "MET-SS", "EFFP-SS", "EFFM-SS",
                                            "JESUP-OS", "JESDOWN-OS", "LEPP-OS", "LEPM-OS", "MET-OS", "EFFP-OS", "EFFM-OS"};

void ww_ana
(
 int thePlot = 9,
 int lSel = 4,
 unsigned int nJetsType = 0,
 int whichGen = 0, // 0 (powheg), 1 (madgraph), 2 (mcatnlo)
 TString bgdInputFile    = "ntuples_53x/backgroundD_skim6.root",
 TString signalInputFile = "ntuples_53x/hww125.root",
 TString dataInputFile   = "ntuples_53x/data_skim6.root",
 TString systInputFile   = "ntuples_53x/hww_syst_skim6.root",
 int period = 3
 )
{
  double events_0Jet[4] = {0.0, 0.0, 0.0, 0.0}; // qqWW, Higgs, Top, Top
  double lumi = 1.0;
  double ptJetMin = 30.0; double ptLLMin = 30.0; double metMin = 20.0; double ptLMin = 20.0;
  double useFullStatTemplates = true;

  bool fCheckProblem = false;

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

  TString ECMsb  = "";
  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double lumiE = 1.099; int year = 1;
  if	 (period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;minRun =      0;maxRun = 999999;ECMsb="8TeV";lumiE = 1.026; year = 2012;
  }
  else if(period == 4){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
    lumi     = 4.924;minRun =	 0;maxRun = 999999;ECMsb="7TeV"; lumiE = 1.022; year = 2011;
    UseDyttDataDriven = false;fCheckProblem = false;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //----------------------------------------------------------------------------
  // NNLO weights
  //----------------------------------------------------------------------------
  TFile *fRatioNNLOFile = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_output_nnlo.root");
  TH1D *fhDRatioNNLO[5];
  if     (whichGen == 0){
    printf("*** using powheg generator ***\n");
    fhDRatioNNLO[0] = (TH1D*)(fRatioNNLOFile->Get("ratio_powheg_central"));
    fhDRatioNNLO[1] = (TH1D*)(fRatioNNLOFile->Get("ratio_powheg_Qup"));
    fhDRatioNNLO[2] = (TH1D*)(fRatioNNLOFile->Get("ratio_powheg_Qdown"));
    fhDRatioNNLO[3] = (TH1D*)(fRatioNNLOFile->Get("ratio_powheg_Rup"));
    fhDRatioNNLO[4] = (TH1D*)(fRatioNNLOFile->Get("ratio_powheg_Rdown"));
  }
  else if(whichGen == 1){
    printf("*** using madgraph generator ***\n");
    fhDRatioNNLO[0] = (TH1D*)(fRatioNNLOFile->Get("ratio_madgraph_central"));
    fhDRatioNNLO[1] = (TH1D*)(fRatioNNLOFile->Get("ratio_madgraph_Qup"));
    fhDRatioNNLO[2] = (TH1D*)(fRatioNNLOFile->Get("ratio_madgraph_Qdown"));
    fhDRatioNNLO[3] = (TH1D*)(fRatioNNLOFile->Get("ratio_madgraph_Rup"));
    fhDRatioNNLO[4] = (TH1D*)(fRatioNNLOFile->Get("ratio_madgraph_Rdown"));
  }
  else if(whichGen == 2){
    printf("*** using mcatnlo generator ***\n");
    fhDRatioNNLO[0] = (TH1D*)(fRatioNNLOFile->Get("ratio_mcatnlo_central"));
    fhDRatioNNLO[1] = (TH1D*)(fRatioNNLOFile->Get("ratio_mcatnlo_Qup"));
    fhDRatioNNLO[2] = (TH1D*)(fRatioNNLOFile->Get("ratio_mcatnlo_Qdown"));
    fhDRatioNNLO[3] = (TH1D*)(fRatioNNLOFile->Get("ratio_mcatnlo_Rup"));
    fhDRatioNNLO[4] = (TH1D*)(fRatioNNLOFile->Get("ratio_mcatnlo_Rdown"));
  }
  else {assert(0);};
  for(int l=0; l<5; l++) assert(fhDRatioNNLO[l]);
  for(int l=0; l<5; l++) fhDRatioNNLO[l]->SetDirectory(0);
  fRatioNNLOFile->Close();
  delete fRatioNNLOFile;

  //----------------------------------------------------------------------------
  // radio photon to electron
  //------------------------;----------------------------------------------------
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

  // useful if using ZZ lepton selection
  LeptonEfficiencyZH theLeptonEfficiencyZH(year);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int nBinMVA = 1;
  Float_t xbins[nBinMVA+1] = {0,100};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Data      = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Higgs     = (TH1D*) histoMVA->Clone("histo_Higgs");
  TH1D *histo_qqWW      = (TH1D*) histoMVA->Clone("histo_qqWW");
  TH1D *histo_ggWW      = (TH1D*) histoMVA->Clone("histo_ggWW");
  TH1D *histo_VV        = (TH1D*) histoMVA->Clone("histo_VV");
  TH1D *histo_VVV       = (TH1D*) histoMVA->Clone("histo_VVV"); 
  TH1D *histo_Top       = (TH1D*) histoMVA->Clone("histo_Top");
  TH1D *histo_Zjets     = (TH1D*) histoMVA->Clone("histo_Zjets");
  TH1D *histo_Ztt       = (TH1D*) histoMVA->Clone("histo_Ztt");
  TH1D *histo_Wg3l      = (TH1D*) histoMVA->Clone("histo_Wg3l");
  TH1D *histo_Wgamma    = (TH1D*) histoMVA->Clone("histo_Wgamma");
  TH1D *histo_WjetsE    = (TH1D*) histoMVA->Clone("histo_WjetsE");
  TH1D *histo_WjetsM    = (TH1D*) histoMVA->Clone("histo_WjetsM");

  char finalStateName[2],effName[10],momName[10];sprintf(effName,"CMS_eff_l");sprintf(momName,"CMS_p_scale_l");
  if     (lSel == 0) {sprintf(finalStateName,"mm");}
  else if(lSel == 1) {sprintf(finalStateName,"me");}
  else if(lSel == 2) {sprintf(finalStateName,"em");}
  else if(lSel == 3) {sprintf(finalStateName,"mm");}
  else if(lSel == 4) {sprintf(finalStateName,"ll");}
  else if(lSel == 5) {sprintf(finalStateName,"sf");}
  else if(lSel == 6) {sprintf(finalStateName,"of");}
  else {printf("Wrong lSel: %d\n",lSel); assert(0);}

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot =100;  xminPlot =   0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 12 && thePlot <= 12) {nBinPlot =200;  xminPlot =   0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 13 && thePlot <= 13) {nBinPlot =100;  xminPlot =   0.0; xmaxPlot = 200.0;}
  else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =100; xminPlot =    0.0; xmaxPlot = 400.0;} // mll
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 20; xminPlot = 0.0; xmaxPlot = 2000.0;} // mlljj
  else if(thePlot >= 20 && thePlot <= 24) {nBinPlot = 18; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 25 && thePlot <= 25) {nBinPlot = 40; xminPlot = 0.0; xmaxPlot = 2.0;}
  else if(thePlot >= 26 && thePlot <= 26) {nBinPlot = 200;  xminPlot =  -1.0; xmaxPlot = 1.0;}
  else if(thePlot >= 27 && thePlot <= 28) {nBinPlot = 100;  xminPlot =   0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 29 && thePlot <= 29) {nBinPlot = 100;  xminPlot =   0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 30 && thePlot <= 30) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 300; xminPlot = 0.0; xmaxPlot = 600.0;}
  else if(thePlot >= 33 && thePlot <= 33) {nBinPlot = 90; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  800.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 8; xminPlot = 0.0; xmaxPlot =  2000;} // mjj
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 9; xminPlot = 0.0; xmaxPlot =  8.75;} // detajjs
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  5.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  else if(thePlot >= 53 && thePlot <= 55) {nBinPlot = 36; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 56 && thePlot <= 56) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 57 && thePlot <= 57) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot = 5.0;}

  TH1D* histos;
  histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
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

  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp   = new TH1D( Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown = new TH1D( Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingUp           = new TH1D( Form("histo_qqWW_CMS_wwana%s__MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_wwana%s__MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingDown         = new TH1D( Form("histo_qqWW_CMS_wwana%s_MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_wwana%s__MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingUp           = new TH1D( Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingDown         = new TH1D( Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingUp                 = new TH1D( Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingDown               = new TH1D( Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp               = new TH1D( Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown             = new TH1D( Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingUp                 = new TH1D( Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingDown               = new TH1D( Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp                 = new TH1D( Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown               = new TH1D( Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVAZttStatBoundingUp                 = new TH1D( Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Ztt_CMS_MVAZttStatBoundingUp  ->Sumw2();
  TH1D* histo_Ztt_CMS_MVAZttStatBoundingDown               = new TH1D( Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Ztt_CMS_MVAZttStatBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingUp                 = new TH1D( Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wg3l_CMS_MVAWg3lStatBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingDown               = new TH1D( Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wg3l_CMS_MVAWg3lStatBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingUp                 = new TH1D( Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wgamma_CMS_MVAWgammaStatBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingDown               = new TH1D( Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wgamma_CMS_MVAWgammaStatBoundingDown->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingUp           = new TH1D( Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingDown         = new TH1D( Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingUp           = new TH1D( Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingDown         = new TH1D( Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Sumw2();

  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp  [nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinUp          [nBinMVA];
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinDown        [nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinUp          [nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinDown        [nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinUp                [nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinDown              [nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp              [nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown            [nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinUp                [nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinDown              [nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinUp                [nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinDown              [nBinMVA];
  TH1D* histo_Ztt_CMS_MVAZttStatBoundingBinUp                [nBinMVA];
  TH1D* histo_Ztt_CMS_MVAZttStatBoundingBinDown              [nBinMVA];
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp                [nBinMVA];
  TH1D* histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown              [nBinMVA];
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp                [nBinMVA];
  TH1D* histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown              [nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp          [nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown        [nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp          [nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown        [nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%s_Bin%dUp"	     ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%s_Bin%dUp"	    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]	   ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb] = new TH1D(Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wwana%s_MVAHiggsStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]	  ->Sumw2();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_qqWW_CMS_wwana%s_MVAqqWWStatBounding_%s_Bin%dUp"	    ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_wwana%s_MVAqqWWStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]        ->Sumw2();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]	    = new TH1D(Form("histo_qqWW_CMS_wwana%s_MVAqqWWStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_wwana%s_MVAqqWWStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]	 ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]	  ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]         = new TH1D(Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_wwana%s_MVAggWWStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]	  ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinUp[nb]	            = new TH1D(Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinUp[nb]	  ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinDown[nb]	            = new TH1D(Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_wwana%s_MVAVVStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinDown[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	            = new TH1D(Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	            = new TH1D(Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wwana%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinUp[nb]	            = new TH1D(Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinUp[nb]	  ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinDown[nb]	            = new TH1D(Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_wwana%s_MVATopStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinDown[nb]	  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]	            = new TH1D(Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]	  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]	            = new TH1D(Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_wwana%s_MVAZjetsStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]	  ->Sumw2();
    histo_Ztt_CMS_MVAZttStatBoundingBinUp[nb]	            = new TH1D(Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Ztt_CMS_MVAZttStatBoundingBinUp[nb]	  ->Sumw2();
    histo_Ztt_CMS_MVAZttStatBoundingBinDown[nb]	            = new TH1D(Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Ztt_CMS_wwana%s_MVAZttStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Ztt_CMS_MVAZttStatBoundingBinDown[nb]	  ->Sumw2();
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp[nb]	            = new TH1D(Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp[nb]	  ->Sumw2();
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown[nb]	            = new TH1D(Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Wg3l_CMS_wwana%s_MVAWg3lStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown[nb]	  ->Sumw2();
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp[nb]	            = new TH1D(Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp[nb]	  ->Sumw2();
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown[nb]	            = new TH1D(Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Wgamma_CMS_wwana%s_MVAWgammaStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown[nb]	  ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]           = new TH1D(Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]  ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb]         = new TH1D(Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_wwana%s_MVAWjetsEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb]->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]           = new TH1D(Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]  ->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb]         = new TH1D(Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_wwana%s_MVAWjetsMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb]->Sumw2();
  }

  TH1D* histo_Higgs_CMS_MVALepEffBoundingUp    = new TH1D( Form("histo_Higgs_%sUp",effName)  , Form("histo_Higgs_%sUp",effName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffBoundingDown  = new TH1D( Form("histo_Higgs_%sDown",effName), Form("histo_Higgs_%sDown",effName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffBoundingUp        = new TH1D( Form("histo_qqWW_%sUp",effName)  , Form("histo_qqWW_%sUp",effName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffBoundingDown      = new TH1D( Form("histo_qqWW_%sDown",effName), Form("histo_qqWW_%sDown",effName), nBinMVA, xbins); histo_qqWW_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffBoundingUp        = new TH1D( Form("histo_ggWW_%sUp",effName)  , Form("histo_ggWW_%sUp",effName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffBoundingDown      = new TH1D( Form("histo_ggWW_%sDown",effName), Form("histo_ggWW_%sDown",effName), nBinMVA, xbins); histo_ggWW_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_VV_%sUp",effName)  , Form("histo_VV_%sUp",effName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_VV_%sDown",effName), Form("histo_VV_%sDown",effName), nBinMVA, xbins); histo_VV_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingUp          = new TH1D( Form("histo_VVV_%sUp",effName)  , Form("histo_VVV_%sUp",effName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingDown        = new TH1D( Form("histo_VVV_%sDown",effName), Form("histo_VVV_%sDown",effName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_Top_%sUp",effName)  , Form("histo_Top_%sUp",effName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_Top_%sDown",effName), Form("histo_Top_%sDown",effName), nBinMVA, xbins); histo_Top_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_Zjets_%sUp",effName)  , Form("histo_Zjets_%sUp",effName)  , nBinMVA, xbins); histo_Zjets_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_Zjets_%sDown",effName), Form("histo_Zjets_%sDown",effName), nBinMVA, xbins); histo_Zjets_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_Ztt_%sUp",effName)  , Form("histo_Ztt_%sUp",effName)  , nBinMVA, xbins); histo_Ztt_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Ztt_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_Ztt_%sDown",effName), Form("histo_Ztt_%sDown",effName), nBinMVA, xbins); histo_Ztt_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_Wgamma_%sUp",effName)  , Form("histo_Wgamma_%sUp",effName)  , nBinMVA, xbins); histo_Wgamma_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_Wgamma_%sDown",effName), Form("histo_Wgamma_%sDown",effName), nBinMVA, xbins); histo_Wgamma_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_Wg3l_%sUp",effName)  , Form("histo_Wg3l_%sUp",effName)  , nBinMVA, xbins); histo_Wg3l_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_Wg3l_%sDown",effName), Form("histo_Wg3l_%sDown",effName), nBinMVA, xbins); histo_Wg3l_CMS_MVALepEffBoundingDown->Sumw2();

  TH1D* histo_Higgs_CMS_MVALepResBoundingUp    = new TH1D( Form("histo_Higgs_%sUp",momName)  , Form("histo_Higgs_%sUp",momName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepResBoundingDown  = new TH1D( Form("histo_Higgs_%sDown",momName), Form("histo_Higgs_%sDown",momName), nBinMVA, xbins); histo_Higgs_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepResBoundingUp        = new TH1D( Form("histo_qqWW_%sUp",momName)  , Form("histo_qqWW_%sUp",momName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepResBoundingDown      = new TH1D( Form("histo_qqWW_%sDown",momName), Form("histo_qqWW_%sDown",momName), nBinMVA, xbins); histo_qqWW_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepResBoundingUp        = new TH1D( Form("histo_ggWW_%sUp",momName)  , Form("histo_ggWW_%sUp",momName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepResBoundingDown      = new TH1D( Form("histo_ggWW_%sDown",momName), Form("histo_ggWW_%sDown",momName), nBinMVA, xbins); histo_ggWW_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_VV_%sUp",momName)  , Form("histo_VV_%sUp",momName)  , nBinMVA, xbins); histo_VV_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_VV_%sDown",momName), Form("histo_VV_%sDown",momName), nBinMVA, xbins); histo_VV_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingUp          = new TH1D( Form("histo_VVV_%sUp",momName)  , Form("histo_VVV_%sUp",momName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingDown        = new TH1D( Form("histo_VVV_%sDown",momName), Form("histo_VVV_%sDown",momName), nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_Top_%sUp",momName)  , Form("histo_Top_%sUp",momName)  , nBinMVA, xbins); histo_Top_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_Top_%sDown",momName), Form("histo_Top_%sDown",momName), nBinMVA, xbins); histo_Top_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_Zjets_%sUp",momName)  , Form("histo_Zjets_%sUp",momName)  , nBinMVA, xbins); histo_Zjets_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_Zjets_%sDown",momName), Form("histo_Zjets_%sDown",momName), nBinMVA, xbins); histo_Zjets_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_Ztt_%sUp",momName)  , Form("histo_Ztt_%sUp",momName)  , nBinMVA, xbins); histo_Ztt_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Ztt_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_Ztt_%sDown",momName), Form("histo_Ztt_%sDown",momName), nBinMVA, xbins); histo_Ztt_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_Wgamma_%sUp",momName)  , Form("histo_Wgamma_%sUp",momName)  , nBinMVA, xbins); histo_Wgamma_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_Wgamma_%sDown",momName), Form("histo_Wgamma_%sDown",momName), nBinMVA, xbins); histo_Wgamma_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_Wg3l_%sUp",momName)  , Form("histo_Wg3l_%sUp",momName)  , nBinMVA, xbins); histo_Wg3l_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_Wg3l_%sDown",momName), Form("histo_Wg3l_%sDown",momName), nBinMVA, xbins); histo_Wg3l_CMS_MVALepResBoundingDown->Sumw2();

  TH1D* histo_Higgs_CMS_MVAMETResBoundingUp    = new TH1D( Form("histo_Higgs_%sUp","CMS_scale_met")  , Form("histo_Higgs_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETResBoundingDown  = new TH1D( Form("histo_Higgs_%sDown","CMS_scale_met"), Form("histo_Higgs_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVAMETResBoundingUp        = new TH1D( Form("histo_qqWW_%sUp","CMS_scale_met")  , Form("histo_qqWW_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_qqWW_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAMETResBoundingDown      = new TH1D( Form("histo_qqWW_%sDown","CMS_scale_met"), Form("histo_qqWW_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_qqWW_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETResBoundingUp        = new TH1D( Form("histo_ggWW_%sUp","CMS_scale_met")  , Form("histo_ggWW_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_ggWW_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETResBoundingDown      = new TH1D( Form("histo_ggWW_%sDown","CMS_scale_met"), Form("histo_ggWW_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_ggWW_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_VV_%sUp","CMS_scale_met")  , Form("histo_VV_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_VV_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_VV_%sDown","CMS_scale_met"), Form("histo_VV_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_VV_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingUp          = new TH1D( Form("histo_VVV_%sUp","CMS_scale_met")  , Form("histo_VVV_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingDown        = new TH1D( Form("histo_VVV_%sDown","CMS_scale_met"), Form("histo_VVV_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_Top_%sUp","CMS_scale_met")  , Form("histo_Top_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Top_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_Top_%sDown","CMS_scale_met"), Form("histo_Top_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Top_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_Zjets_%sUp","CMS_scale_met")  , Form("histo_Zjets_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Zjets_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_Zjets_%sDown","CMS_scale_met"), Form("histo_Zjets_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Zjets_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_Ztt_%sUp","CMS_scale_met")  , Form("histo_Ztt_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Ztt_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Ztt_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_Ztt_%sDown","CMS_scale_met"), Form("histo_Ztt_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Ztt_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_Wgamma_%sUp","CMS_scale_met")  , Form("histo_Wgamma_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Wgamma_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_Wgamma_%sDown","CMS_scale_met"), Form("histo_Wgamma_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Wgamma_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_Wg3l_%sUp","CMS_scale_met")  , Form("histo_Wg3l_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_Wg3l_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_Wg3l_%sDown","CMS_scale_met"), Form("histo_Wg3l_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_Wg3l_CMS_MVAMETResBoundingDown->Sumw2();

  TH1D* histo_Higgs_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_Higgs_%sUp","CMS_scale_j")  , Form("histo_Higgs_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_Higgs_%sDown","CMS_scale_j"), Form("histo_Higgs_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_qqWW_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_qqWW_%sUp","CMS_scale_j")  , Form("histo_qqWW_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_qqWW_%sDown","CMS_scale_j"), Form("histo_qqWW_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_ggWW_%sUp","CMS_scale_j")  , Form("histo_ggWW_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_ggWW_%sDown","CMS_scale_j"), Form("histo_ggWW_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_VV_%sUp","CMS_scale_j")  , Form("histo_VV_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_VV_%sDown","CMS_scale_j"), Form("histo_VV_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp          = new TH1D( Form("histo_VVV_%sUp","CMS_scale_j")  , Form("histo_VVV_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown        = new TH1D( Form("histo_VVV_%sDown","CMS_scale_j"), Form("histo_VVV_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Top_%sUp","CMS_scale_j")  , Form("histo_Top_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Top_%sDown","CMS_scale_j"), Form("histo_Top_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Zjets_%sUp","CMS_scale_j")  , Form("histo_Zjets_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Zjets_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Zjets_%sDown","CMS_scale_j"), Form("histo_Zjets_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Zjets_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Ztt_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Ztt_%sUp","CMS_scale_j")  , Form("histo_Ztt_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Ztt_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Ztt_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Ztt_%sDown","CMS_scale_j"), Form("histo_Ztt_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Ztt_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Wgamma_%sUp","CMS_scale_j")  , Form("histo_Wgamma_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Wgamma_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Wgamma_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Wgamma_%sDown","CMS_scale_j"), Form("histo_Wgamma_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Wgamma_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Wg3l_%sUp","CMS_scale_j")  , Form("histo_Wg3l_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_Wg3l_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Wg3l_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Wg3l_%sDown","CMS_scale_j"), Form("histo_Wg3l_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_Wg3l_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_VV_CMS_VVNLOBoundingUp        = new TH1D( Form("histo_VV_CMS_wwana_VVNLOBoundingUp"),   Form("histo_VV_CMS_wwana_VVNLOBoundingUp"),   nBinMVA, xbins); histo_VV_CMS_VVNLOBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_VVNLOBoundingDown      = new TH1D( Form("histo_VV_CMS_wwana_VVNLOBoundingDown"), Form("histo_VV_CMS_wwana_VVNLOBoundingDown"), nBinMVA, xbins); histo_VV_CMS_VVNLOBoundingDown->Sumw2();

  TH1D* histo_WjetsE_CMS_MVAWEBoundingUp    = new TH1D( Form("histo_WjetsE_CMS_wwana_MVAWEBoundingUp"),   Form("histo_WjetsE_CMS_wwana_MVAWEBoundingUp"),   nBinMVA, xbins); histo_WjetsE_CMS_MVAWEBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWEBoundingDown  = new TH1D( Form("histo_WjetsE_CMS_wwana_MVAWEBoundingDown"), Form("histo_WjetsE_CMS_wwana_MVAWEBoundingDown"), nBinMVA, xbins); histo_WjetsE_CMS_MVAWEBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWMBoundingUp    = new TH1D( Form("histo_WjetsM_CMS_wwana_MVAWMBoundingUp"),   Form("histo_WjetsM_CMS_wwana_MVAWMBoundingUp"),   nBinMVA, xbins); histo_WjetsM_CMS_MVAWMBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWMBoundingDown  = new TH1D( Form("histo_WjetsM_CMS_wwana_MVAWMBoundingDown"), Form("histo_WjetsM_CMS_wwana_MVAWMBoundingDown"), nBinMVA, xbins); histo_WjetsM_CMS_MVAWMBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWWBoundingUp = new TH1D( Form("histo_qqWW_CMS_wwana_MVAWWBoundingUp"), Form("histo_qqWW_CMS_wwana_MVAWWBoundingUp"), nBinMVA, xbins); histo_qqWW_CMS_MVAWWBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWWBoundingDown = new TH1D( Form("histo_qqWW_CMS_wwana_MVAWWBoundingDown"), Form("histo_qqWW_CMS_wwana_MVAWWBoundingDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAWWBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWWNLOBoundingUp = new TH1D( Form("histo_qqWW_CMS_wwana_MVAWWNLOBoundingUp"), Form("histo_qqWW_CMS_wwana_MVAWWNLOBoundingUp"), nBinMVA, xbins); histo_qqWW_CMS_MVAWWNLOBoundingUp->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWWNLOBoundingDown = new TH1D( Form("histo_qqWW_CMS_wwana_MVAWWNLOBoundingDown"), Form("histo_qqWW_CMS_wwana_MVAWWNLOBoundingDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAWWNLOBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_WWNLOQUp   = new TH1D( Form("histo_qqWW_CMS_wwana_WWNLOQUp"),   Form("histo_qqWW_wwana_CMS_WWNLOQUp"),   nBinMVA, xbins); histo_qqWW_CMS_WWNLOQUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_WWNLOQDown = new TH1D( Form("histo_qqWW_CMS_wwana_WWNLOQDown"), Form("histo_qqWW_wwana_CMS_WWNLOQDown"), nBinMVA, xbins); histo_qqWW_CMS_WWNLOQDown->Sumw2();
  TH1D* histo_qqWW_CMS_WWNLORUp   = new TH1D( Form("histo_qqWW_CMS_wwana_WWNLORUp"),   Form("histo_qqWW_wwana_CMS_WWNLORUp"),   nBinMVA, xbins); histo_qqWW_CMS_WWNLORUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_WWNLORDown = new TH1D( Form("histo_qqWW_CMS_wwana_WWNLORDown"), Form("histo_qqWW_wwana_CMS_WWNLORDown"), nBinMVA, xbins); histo_qqWW_CMS_WWNLORDown->Sumw2();

  double nSelectedData[nSelTypes*2];
  double nSigCut[nSelTypes*2],nSigECut[nSelTypes*2];
  double bgdDecay[nSelTypes*2][45],weiDecay[nSelTypes*2][45];
  double nSigCutSyst[nSelTypesSyst*2],nSigECutSyst[nSelTypesSyst*2];
  double bgdDecaySyst[nSelTypesSyst*2][45],weiDecaySyst[nSelTypesSyst*2][45];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    nSelectedData[i] = 0.0; nSigCut[i] = 0.0; nSigECut[i] = 0.0;
    for(int j=0; j<45; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    nSigCutSyst[i] = 0.0; nSigECutSyst[i] = 0.0;
    for(int j=0; j<45; j++) {
      bgdDecaySyst[i][j] = 0.0; weiDecaySyst[i][j] = 0.0; 
    }
  }

  unsigned int patternTopVeto = SmurfTree::TopVeto;
  double genLevelNorm[4] = {0.,0.,0.,0.};

  float dymvaMET=0.;float dymvaJESU=0;float dymvaJESD=0;
  bgdEvent.tree_->SetBranchAddress("dymvaMET", &dymvaMET );
  bgdEvent.tree_->SetBranchAddress("dymvaJESU",&dymvaJESU);
  bgdEvent.tree_->SetBranchAddress("dymvaJESD",&dymvaJESD);
  int nBgd=bgdEvent.tree_->GetEntries();
  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

    // generator level selection
    bool minGenCuts = !(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
                        ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
			((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection && bgdEvent.lid3_ != 0));
    bool genLevelSel = false;
    if(minGenCuts == true) {
      genLevelNorm[0]++;
      int nGenJets = 0;
      if(bgdEvent.genjet1_.Pt() > 30 && TMath::Abs(bgdEvent.genjet1_.Eta()) < 5.0) nGenJets++;
      if(bgdEvent.genjet2_.Pt() > 30 && TMath::Abs(bgdEvent.genjet2_.Eta()) < 5.0) nGenJets++;
      if(bgdEvent.genjet2_.Pt() > 30 && TMath::Abs(bgdEvent.genjet3_.Eta()) < 5.0) nGenJets++;

      if(nGenJets == 0) {
        genLevelNorm[1]++;
	genLevelSel = true;
      }
    }

    if(bgdEvent.lep1_.Pt() < 1.0) continue;

    if(!(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  bgdEvent.dstype_ != SmurfTree::data)) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    LorentzVector genjet_jer1;
    double deltaRJJGen = 999.9;
    if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet1_.Phi(),bgdEvent.genjet1_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet1_.Phi(),bgdEvent.genjet1_.Eta());
      genjet_jer1 = bgdEvent.genjet1_;
    }
    if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet2_.Phi(),bgdEvent.genjet2_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet2_.Phi(),bgdEvent.genjet2_.Eta());
      genjet_jer1 = bgdEvent.genjet2_;
    }
    if(DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet3_.Phi(),bgdEvent.genjet3_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(bgdEvent.jet1_.Phi(),bgdEvent.jet1_.Eta(),bgdEvent.genjet3_.Phi(),bgdEvent.genjet3_.Eta());
      genjet_jer1 = bgdEvent.genjet3_;
    }

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::wjets  	   ) fDecay =  3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.processId_==121 ||
            bgdEvent.processId_==122)   fDecay = 41;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 43;
    else if(bgdEvent.processId_==10001) fDecay = 44;
    else if(bgdEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false},
                                   {false}};
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = bgdEvent.met_; double theMETPHI = bgdEvent.metPhi_; 
    
    int lType = 0;
    if     (bgdEvent.lq1_ * bgdEvent.lq2_ < 0) lType = 1;

    double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool   passMET = usedMet > metMin;
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

    double ptLLCut = ptLLMin; if(bgdEvent.type_ == SmurfTree::ee || bgdEvent.type_ == SmurfTree::mm) ptLLCut = 45;
    bool passNjets         = bgdEvent.njets_ == nJetsType; if(nJetsType == 3) passNjets = bgdEvent.njets_ <= 1;
    bool preselCuts        = bgdEvent.lep1_.Pt() > ptLMin && bgdEvent.lep2_.Pt() > ptLMin && bgdEvent.dilep_.Pt() > ptLLCut;
    bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = bgdEvent.dilep_.M() > 12 && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);

    // 0      1      2       3     4   5      6        7           8  9            10            11     12  13    14
    // lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtZ,mlljj,mjj;
    double outputVarLepP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[6] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13]};
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
    		 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    unsigned int NjetSyst[2] = {0., 0.};
    if(bgdEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;

    bool   passMETSyst[3] = {TMath::Min(outputVarJESP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin,
                             TMath::Min(outputVarJESM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin,
			     TMath::Min(outputVarMET[4] /bgdEvent.met_*bgdEvent.pmet_,outputVarMET[6] /bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin};
    if(bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee){
      if(useDYMVA == false){
        if     (NjetSyst[0]     <= 1) passMETSyst[0] = passMETSyst[0] && TMath::Min(outputVarJESP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[0] = passMETSyst[0] &&	    outputVarJESP[4]											   > 45.0;

        if     (NjetSyst[1]     <= 1) passMETSyst[1] = passMETSyst[1] && TMath::Min(outputVarJESM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[1] = passMETSyst[1] &&	    outputVarJESM[4]											   > 45.0;

        if     (bgdEvent.njets_ <= 1) passMETSyst[2] = passMETSyst[2] && TMath::Min( outputVarMET[4]/bgdEvent.met_*bgdEvent.pmet_, outputVarMET[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[3] = passMETSyst[3] &&	     outputVarMET[4]											   > 45.0;


      } else {
        if     (NjetSyst[0]     == 0) passMETSyst[0] = passMETSyst[0] && dymvaJESU >  0.88;
        else if(NjetSyst[0]     == 1) passMETSyst[0] = passMETSyst[0] && dymvaJESU >  0.84;
        else                          passMETSyst[0] = passMETSyst[0] &&            outputVarJESP[4]                                                                                       > 45.0;

        if     (NjetSyst[0]     == 0) passMETSyst[1] = passMETSyst[1] && dymvaJESD >  0.88;
        else if(NjetSyst[0]     == 1) passMETSyst[1] = passMETSyst[1] && dymvaJESD >  0.84;
        else                          passMETSyst[1] = passMETSyst[1] &&            outputVarJESM[4]                                                                                       > 45.0;

        if     (bgdEvent.njets_ == 0) passMETSyst[2] = passMETSyst[2] && dymvaMET  >  0.88;
        else if(bgdEvent.njets_ == 1) passMETSyst[2] = passMETSyst[2] && dymvaMET  >  0.84;
        else                          passMETSyst[2] = passMETSyst[2] &&            outputVarMET[4]                                                                                        > 45.0;

      }
    }

    bool passLSel = false;
    if     (lSel == 0 && bgdEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && bgdEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && bgdEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && bgdEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (bgdEvent.type_ == SmurfTree::me || bgdEvent.type_ == SmurfTree::em)) passLSel = true;

    if(nJetsType == 3 && NjetSyst[0] <= 1) NjetSyst[0] = 3;
    if(nJetsType == 3 && NjetSyst[1] <= 1) NjetSyst[1] = 3;

    if(passLSel && NjetSyst[0] == nJetsType     && passMETSyst[0]                                                                                                                        && dPhiDiLepJetCut && outputVar[0]	> ptLMin && outputVar[1]     > ptLMin && outputVar[3]     > ptLLCut && passBtagVeto && pass3rLVeto && outputVar[2]	> 12.0 && (TMath::Abs(outputVar[2]    -91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) passSystCuts[lType][JESUP] = true;
    if(passLSel && NjetSyst[1] == nJetsType     && passMETSyst[1]                                                                                                                        && dPhiDiLepJetCut && outputVar[0]	> ptLMin && outputVar[1]     > ptLMin && outputVar[3]     > ptLLCut && passBtagVeto && pass3rLVeto && outputVar[2]	> 12.0 && (TMath::Abs(outputVar[2]    -91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) passSystCuts[lType][JESDOWN] = true;
    if(passLSel && passNjets                    && TMath::Min(outputVarLepP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && passMET && dPhiDiLepJetCut && outputVarLepP[0] > ptLMin && outputVarLepP[1] > ptLMin && outputVarLepP[3] > ptLLCut && passBtagVeto && pass3rLVeto && outputVarLepP[2] > 12.0 && (TMath::Abs(outputVarLepP[2]-91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) passSystCuts[lType][LEPP] = true;
    if(passLSel && passNjets                    && TMath::Min(outputVarLepM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && passMET && dPhiDiLepJetCut && outputVarLepM[0] > ptLMin && outputVarLepM[1] > ptLMin && outputVarLepM[3] > ptLLCut && passBtagVeto && pass3rLVeto && outputVarLepM[2] > 12.0 && (TMath::Abs(outputVarLepM[2]-91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) passSystCuts[lType][LEPM] = true;
    if(passLSel && passNjets                    && passMETSyst[2]															 && dPhiDiLepJetCut && outputVarMET[0]  > ptLMin && outputVarMET[1]  > ptLMin && outputVarMET[3]  > ptLLCut && passBtagVeto && pass3rLVeto && outputVarMET[2]  > 12.0 && (TMath::Abs(outputVarMET[2] -91.1876) > 15 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) passSystCuts[lType][MET] = true;

    if(passLSel) {
       
       if(passNjets && preselCuts && passBtagVeto && pass3rLVeto && passMass && dPhiDiLepJetCut && passMET) passCuts[lType][WWSEL] = true;

    }

    if(passLSel){
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
        add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	fDecay = 22;
	theWeight	       = -1.0*add;
      }
      else if(nFake == 1){
        if(bgdEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
          if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	}
	else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,bgdEvent.npu_);
          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt(), 
								   fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						   TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
            double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	 							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        						             TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	    double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	  							      fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	    double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     	  							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	    trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
         }
	  
	  add = add*trigEff;
	  if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  fDecay                 = 1;

	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					        fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							        TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
        	 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);

        theWeight = ZttScaleFactor(period,bgdEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add2 = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
								 fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						 TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
           double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
                						    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	   double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      								    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	   double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      								    TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	   trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
        }
        add = add1*add2*trigEff;

        if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	 printf("PROBLEMCB(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",bgdEvent.event_,add1,add2,trigEff,add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);

  	if(fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) {
  	  if(bgdEvent.njets_ == 0) add=add*DYBkgScaleFactor(0);
  	  if(bgdEvent.njets_ == 1) add=add*DYBkgScaleFactor(1);
  	  if(bgdEvent.njets_ >= 2) add=add*DYBkgScaleFactor(2);
  	}
  	if(fDecay == 5 || fDecay == 13) {
  	  if	 (bgdEvent.njets_ == 0) add=add*TopBkgScaleFactor(0);
  	  else if(bgdEvent.njets_ == 1) add=add*TopBkgScaleFactor(1); 
  	  else if(bgdEvent.njets_ >= 2) add=add*TopBkgScaleFactor(2);
  	}

	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,theMET);
        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add;
      }
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::qqww  )  theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,3);
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::ggww  )  theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,3);
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::qqww2j)  theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,3);
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::qqwwPWG) theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,3);

      if(useWeightNNLOCorr == true && (bgdEvent.dstype_ == SmurfTree::qqww||bgdEvent.dstype_ == SmurfTree::qqwwPWG)) theWeight = theWeight * weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,0);
      if(useWeightNNLOCorr == true && (bgdEvent.dstype_ == SmurfTree::qqww||bgdEvent.dstype_ == SmurfTree::qqwwPWG||
                                       bgdEvent.dstype_ == SmurfTree::qqww2j||bgdEvent.dstype_ == SmurfTree::ggww)) theWeight = theWeight * weightNLOToNNLOCorr(period);

      //if(bgdEvent.dstype_ == SmurfTree::qqww  )  theWeight = theWeight * weightJetPt(0,bgdEvent.jet3McId_,bgdEvent.jet4McId_);
      //if(bgdEvent.dstype_ == SmurfTree::qqww2j)  theWeight = theWeight * weightJetPt(0,bgdEvent.jet3McId_,bgdEvent.jet4McId_);
      //if(bgdEvent.dstype_ == SmurfTree::qqwwPWG) theWeight = theWeight * weightJetPt(0,bgdEvent.jet3McId_,bgdEvent.jet4McId_);

      if(minGenCuts == true && passCuts[1][WWSEL] && genLevelSel == false) genLevelNorm[2]++;
      if(minGenCuts == true && passCuts[1][WWSEL] && genLevelSel == true)  genLevelNorm[3]++;
 
       if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = theMET;
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
	else if(thePlot ==15) myVar = fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = bgdEvent.njets_;
	else if(thePlot ==18) myVar = bgdEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(bgdEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = DeltaPhi(bgdEvent.dilep_.Phi(),bgdEvent.trackMetPhi_)*180.0/TMath::Pi();
	else if(thePlot ==25) myVar = bgdEvent.lep1_.Pt()*bgdEvent.lep2_.Pt()/bgdEvent.jet1_.Pt()/bgdEvent.jet2_.Pt();
	else if(thePlot ==26) myVar = bgdEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = (bgdEvent.dilep_+bgdEvent.jet1_).Pt();
	else if(thePlot ==30) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,bgdEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((bgdEvent.jet1_+bgdEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = bgdEvent.jet1_.Pt()+ bgdEvent.jet2_.Pt()+bgdEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = bgdEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = bgdEvent.trackMet_*cos(bgdEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = bgdEvent.trackMet_*sin(bgdEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(bgdEvent.jet3_.Phi(),bgdEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = bgdEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = TMath::Min(deltaRJJGen,4.999);
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
      	  printf("NOOOOOOOOOOOOOOOOOOOO %d\n",fDecay);
      	}
      } // end making plots
      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            bgdDecay[i+j*nSelTypes][(int)fDecay] += theWeight;
            weiDecay[i+j*nSelTypes][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            bgdDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight;
            weiDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WWSEL]) {
        bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
        weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
        weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }

      if     (fDecay == 21){
        if(passCuts[1][WWSEL])  	     histo_VVV  			->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_VVV_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_VVV_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_VVV_CMS_MVAJESBoundingUp	->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_VVV_CMS_MVAJESBoundingDown	->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_VVV_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_VVV_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_VVV_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 29){
        if(bgdEvent.njets_ == 0 && passCuts[1][WWSEL]) events_0Jet[0] = events_0Jet[0] + theWeight;
        if(passCuts[1][WWSEL])  	     histo_qqWW			         ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_qqWW_CMS_MVAJESBoundingUp	 ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_qqWW_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_qqWW_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_qqWW_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_qqWW_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;

        if(useWeightNNLOCorr == true && (bgdEvent.dstype_ == SmurfTree::qqww||bgdEvent.dstype_ == SmurfTree::qqwwPWG)){
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLOQUp    	 ->Fill(MVAVar[0], theWeight*weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,1)/weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,0));
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLOQDown  	 ->Fill(MVAVar[0], theWeight*weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,2)/weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,0));
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLORUp    	 ->Fill(MVAVar[0], theWeight*weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,3)/weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,0));
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLORDown  	 ->Fill(MVAVar[0], theWeight*weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,4)/weightNNLOCorr(fhDRatioNNLO,bgdEvent.jet3McId_,0));
        } else {
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLOQUp    	 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLOQDown  	 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLORUp    	 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][WWSEL])  	     histo_qqWW_CMS_WWNLORDown  	 ->Fill(MVAVar[0], theWeight);
	}
      }
      else if(fDecay == 30){
        if(passCuts[1][WWSEL])  	     histo_ggWW		          ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_ggWW_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_ggWW_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_ggWW_CMS_MVAJESBoundingUp	  ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_ggWW_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_ggWW_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_ggWW_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_ggWW_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 27 || fDecay == 28){
        if(passCuts[1][WWSEL])  	     histo_VV			       ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_VV_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_VV_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_VV_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_VV_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_VV_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_VV_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_VV_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay ==  5 || fDecay == 13){
        if(bgdEvent.njets_ == 0 && passCuts[1][WWSEL]) events_0Jet[2] = events_0Jet[2] + theWeight;
        if(passCuts[1][WWSEL])  	     histo_Top		               ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_Top_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_Top_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Top_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Top_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Top_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Top_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Top_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);
      }
      else if(fDecay ==  9){
        if(bgdEvent.njets_ == 0 && passCuts[1][WWSEL]) events_0Jet[3] = events_0Jet[3] + theWeight;
        if(passCuts[1][WWSEL])  	     histo_Zjets		               ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_Zjets_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_Zjets_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Zjets_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Zjets_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Zjets_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Zjets_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Zjets_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);
      }
      else if(fDecay == 10){
        if(passCuts[1][WWSEL])  	     histo_Ztt		               ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_Ztt_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_Ztt_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Ztt_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Ztt_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Ztt_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Ztt_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Ztt_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);
      }
      else if(fDecay == 20){
        if(passCuts[1][WWSEL])  	     histo_Wg3l		               ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_Wg3l_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_Wg3l_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Wg3l_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Wg3l_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Wg3l_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Wg3l_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Wg3l_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);
      }
      else if(fDecay == 19){
        if(passCuts[1][WWSEL])  	     histo_Wgamma                          ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_Wgamma_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_Wgamma_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Wgamma_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Wgamma_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Wgamma_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Wgamma_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Wgamma_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);
      }
      else if(fDecay == 1){
        double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        if(passCuts[1][WWSEL]) 	     histo_WjetsE		        ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	     histo_WjetsE_CMS_MVAWEBoundingUp   ->Fill(MVAVar[0], theWeight*addFRS/addFR);
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
        if(passCuts[1][WWSEL]) 	     histo_WjetsM		        ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	     histo_WjetsM_CMS_MVAWMBoundingUp   ->Fill(MVAVar[0], theWeight*addFRS/addFR);
      }
      else assert(0);
    } // if passCuts
  } // end background loop

  if(systInputFile != ""){
  int nSyst=systEvent.tree_->GetEntries();
  for (int evt=0; evt<nSyst; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nSyst);
    systEvent.tree_->GetEntry(evt);

    if(systEvent.lep1_.Pt() < 1.0) continue;

    if(systEvent.dstype_ == SmurfTree::data &&
      (systEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ <  minRun) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (systEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(systEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
    else if(systEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(systEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::dyttDataDriven ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(systEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(systEvent.dstype_ == SmurfTree::wgstar         ) fDecay = 20;
    else if(systEvent.dstype_ == SmurfTree::www            ) fDecay = 21;
    else if(systEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(systEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(systEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 12;
    else if(systEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(systEvent.dstype_ == SmurfTree::other          ) fDecay = 40;
    else if(systEvent.processId_==121 ||
            systEvent.processId_==122)   fDecay = 41;
    else if(systEvent.processId_==24)    fDecay = 42;
    else if(systEvent.processId_==26)    fDecay = 43;
    else if(systEvent.processId_==10001) fDecay = 44;
    else if(systEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;}

    bool passCuts[3][nSelTypes] = {{false},
                                   {false}};
    bool isRealLepton = false;
    if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
       (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = systEvent.met_; double theMETPHI = systEvent.metPhi_; 
    
    int lType = 0;
    if     (systEvent.lq1_ * systEvent.lq2_ < 0) lType = 1;

    double usedMet = TMath::Min(systEvent.pmet_,systEvent.pTrackMet_);
    bool   passMET = usedMet > metMin;
    if(useDYMVA == false){
      if     (systEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
      else if(systEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
      else                           passMET = passMET && (systEvent.met_ > 45.0 || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
    } else {
      if     (systEvent.njets_ == 0) passMET = passMET && (systEvent.dymva_ >  0.88 || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
      else if(systEvent.njets_ == 1) passMET = passMET && (systEvent.dymva_ >  0.84 || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
      else                           passMET = passMET && (systEvent.met_   >  45.0 || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);
    }
    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(systEvent.njets_ <= 1) dPhiDiLepJetCut = systEvent.jet1_.Pt() <= 15. || systEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
      	                                         systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
    	                                         systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;
    }
    if(systEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((systEvent.jet1_+systEvent.jet2_).Phi(),systEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me;

    double ptLLCut = ptLLMin; if(systEvent.type_ == SmurfTree::ee || systEvent.type_ == SmurfTree::mm) ptLLCut = 45;
    bool passNjets         = systEvent.njets_ == nJetsType; if(nJetsType == 3) passNjets = systEvent.njets_ <= 1;
    bool preselCuts        = systEvent.lep1_.Pt() > ptLMin && systEvent.lep2_.Pt() > ptLMin && systEvent.dilep_.Pt() > ptLLCut;
    bool passBtagVeto      = (systEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = systEvent.dilep_.M() > 12 && (TMath::Abs(systEvent.dilep_.M()-91.1876) > 15 || systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me);

    bool passLSel = false;
    if     (lSel == 0 && systEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && systEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && systEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && systEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (systEvent.type_ == SmurfTree::me || systEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passLSel) {
       
       if(passNjets && preselCuts && passBtagVeto && pass3rLVeto && passMass && dPhiDiLepJetCut && passMET) passCuts[lType][WWSEL] = true;

    }

    if(passCuts[lType][WWSEL]){
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
          if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt(), 
								   fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						   TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
          add = add*trigEff;
	  if(fCheckProblem == true && (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMBSyst: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  fDecay                 = 1;

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
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_)*
        	 leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);

        theWeight = ZttScaleFactor(period,systEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(systEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,systEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	add2 = add2*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
        if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
								 fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        add = add1*add2*trigEff;

        if(fCheckProblem == true && add != 0 && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	printf("PROBLEMCSy(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",systEvent.event_,add1,add2,trigEff,add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);

	if(systEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(systEvent.type_,theMET);

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

      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][WWSEL]){
	if     (fDecay == 27 || fDecay == 28){
	  histo_VV_CMS_VVNLOBoundingUp->Fill(MVAVar[0], theWeight);
        }
        if     (fDecay == 30) histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Fill(MVAVar[0], theWeight);
	else if(fDecay == 12) histo_qqWW_CMS_MVAWWNLOBoundingDown->Fill(MVAVar[0], theWeight);
	else if(fDecay == 29) histo_qqWW_CMS_MVAWWBoundingUp     ->Fill(MVAVar[0], theWeight);
      }

    } // if passCuts
  } // end syst loop
  } // if want to use it at all
  if(showSignalOnly == false) {
    printf("nuisance WW: %0.1f/%0.1f/%0.1f/%0.1f\n",histo_qqWW                       ->GetSumOfWeights(),histo_qqWW_CMS_MVAWWBoundingUp     ->GetSumOfWeights(),
						    histo_qqWW_CMS_MVAWWNLOBoundingUp->GetSumOfWeights(),histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights());
  }

  sigEvent.tree_->SetBranchAddress("dymvaMET", &dymvaMET );
  sigEvent.tree_->SetBranchAddress("dymvaJESU",&dymvaJESU);
  sigEvent.tree_->SetBranchAddress("dymvaJESD",&dymvaJESD);
  int nSig=sigEvent.tree_->GetEntries();
  for (int evt=0; evt<nSig; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
       printf("--- reading Signal event %5d of %5d\n",evt,nSig);
    sigEvent.tree_->GetEntry(evt);

    if(sigEvent.lep1_.Pt() < 1.0) continue;

    bool lId = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
               (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    LorentzVector genjet_jer1;
    double deltaRJJGen = 999.9;
    if(DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet1_.Phi(),sigEvent.genjet1_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet1_.Phi(),sigEvent.genjet1_.Eta());
      genjet_jer1 = sigEvent.genjet1_;
    }
    if(DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet2_.Phi(),sigEvent.genjet2_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet2_.Phi(),sigEvent.genjet2_.Eta());
      genjet_jer1 = sigEvent.genjet2_;
    }
    if(DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet3_.Phi(),sigEvent.genjet3_.Eta()) < deltaRJJGen) {
      deltaRJJGen = DeltaR(sigEvent.jet1_.Phi(),sigEvent.jet1_.Eta(),sigEvent.genjet3_.Phi(),sigEvent.genjet3_.Eta());
      genjet_jer1 = sigEvent.genjet3_;
    }

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false},
                                   {false}};
    double theMET = sigEvent.met_; double theMETPHI = sigEvent.metPhi_; 
    
    int lType = 0;
    if     (sigEvent.lq1_ * sigEvent.lq2_ < 0) lType = 1;

    double usedMet = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_);
    bool   passMET = usedMet > metMin;
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

    double ptLLCut = ptLLMin; if(sigEvent.type_ == SmurfTree::ee || sigEvent.type_ == SmurfTree::mm) ptLLCut = 45;
    bool passNjets         = sigEvent.njets_ == nJetsType; if(nJetsType == 3) passNjets = sigEvent.njets_ <= 1;
    bool preselCuts        = sigEvent.lep1_.Pt() > ptLMin && sigEvent.lep2_.Pt() > ptLMin && sigEvent.dilep_.Pt() > ptLLCut;
    bool passBtagVeto      = (sigEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = sigEvent.dilep_.M() > 12 && (TMath::Abs(sigEvent.dilep_.M()-91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me);

    double outputVarLepP[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[6] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13]};
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 0)*
    		 leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 1)*
        	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,-1)*
        	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    unsigned int NjetSyst[2] = {0., 0.};
    if(sigEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
 
    bool   passMETSyst[3] = {TMath::Min(outputVarJESP[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESP[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin,
                             TMath::Min(outputVarJESM[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESM[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin,
			     TMath::Min(outputVarMET[4] /sigEvent.met_*sigEvent.pmet_,outputVarMET[6] /sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin};
    if(sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee){
      if(useDYMVA == false){
        if     (NjetSyst[0]     <= 1) passMETSyst[0] = passMETSyst[0] && TMath::Min(outputVarJESP[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESP[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[0] = passMETSyst[0] &&	    outputVarJESP[4]											   > 45.0;

        if     (NjetSyst[1]     <= 1) passMETSyst[1] = passMETSyst[1] && TMath::Min(outputVarJESM[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESM[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[1] = passMETSyst[1] &&	    outputVarJESM[4]											   > 45.0;

        if     (sigEvent.njets_ <= 1) passMETSyst[2] = passMETSyst[2] && TMath::Min( outputVarMET[4]/sigEvent.met_*sigEvent.pmet_, outputVarMET[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > 45.0;
        else                          passMETSyst[3] = passMETSyst[3] &&	     outputVarMET[4]											   > 45.0;


      } else {
        if     (NjetSyst[0]     == 0) passMETSyst[0] = passMETSyst[0] && dymvaJESU >  0.88;
        else if(NjetSyst[0]     == 1) passMETSyst[0] = passMETSyst[0] && dymvaJESU >  0.84;
        else                          passMETSyst[0] = passMETSyst[0] &&            outputVarJESP[4]                                                                                       > 45.0;

        if     (NjetSyst[0]     == 0) passMETSyst[1] = passMETSyst[1] && dymvaJESD >  0.88;
        else if(NjetSyst[0]     == 1) passMETSyst[1] = passMETSyst[1] && dymvaJESD >  0.84;
        else                          passMETSyst[1] = passMETSyst[1] &&            outputVarJESM[4]                                                                                       > 45.0;

        if     (sigEvent.njets_ == 0) passMETSyst[2] = passMETSyst[2] && dymvaMET  >  0.88;
        else if(sigEvent.njets_ == 1) passMETSyst[2] = passMETSyst[2] && dymvaMET  >  0.84;
        else                          passMETSyst[2] = passMETSyst[2] &&            outputVarMET[4]                                                                                        > 45.0;

      }
    }

    bool passLSel = false;
    if     (lSel == 0 && sigEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && sigEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && sigEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && sigEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (sigEvent.type_ == SmurfTree::me || sigEvent.type_ == SmurfTree::em)) passLSel = true;

    if(nJetsType == 3 && NjetSyst[0] <= 1) NjetSyst[0] = 3;
    if(nJetsType == 3 && NjetSyst[1] <= 1) NjetSyst[1] = 3;

    if(passLSel && NjetSyst[0] == nJetsType     && passMETSyst[0]                                                                                                                        && dPhiDiLepJetCut && outputVar[0]	> ptLMin && outputVar[1]     > ptLMin && outputVar[3]     > ptLLCut && passBtagVeto && pass3rLVeto && outputVar[2]	> 12.0 && (TMath::Abs(outputVar[2]    -91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) passSystCuts[lType][JESUP] = true;
    if(passLSel && NjetSyst[1] == nJetsType     && passMETSyst[1]                                                                                                                        && dPhiDiLepJetCut && outputVar[0]	> ptLMin && outputVar[1]     > ptLMin && outputVar[3]     > ptLLCut && passBtagVeto && pass3rLVeto && outputVar[2]	> 12.0 && (TMath::Abs(outputVar[2]    -91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) passSystCuts[lType][JESDOWN] = true;
    if(passLSel && passNjets                    && TMath::Min(outputVarLepP[4]/sigEvent.met_*sigEvent.pmet_,outputVarLepP[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && passMET && dPhiDiLepJetCut && outputVarLepP[0] > ptLMin && outputVarLepP[1] > ptLMin && outputVarLepP[3] > ptLLCut && passBtagVeto && pass3rLVeto && outputVarLepP[2] > 12.0 && (TMath::Abs(outputVarLepP[2]-91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) passSystCuts[lType][LEPP] = true;
    if(passLSel && passNjets                    && TMath::Min(outputVarLepM[4]/sigEvent.met_*sigEvent.pmet_,outputVarLepM[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && passMET && dPhiDiLepJetCut && outputVarLepM[0] > ptLMin && outputVarLepM[1] > ptLMin && outputVarLepM[3] > ptLLCut && passBtagVeto && pass3rLVeto && outputVarLepM[2] > 12.0 && (TMath::Abs(outputVarLepM[2]-91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) passSystCuts[lType][LEPM] = true;
    if(passLSel && passNjets                    && passMETSyst[2]                                                                                                                        && dPhiDiLepJetCut && outputVarMET[0]  > ptLMin && outputVarMET[1]  > ptLMin && outputVarMET[3]  > ptLLCut && passBtagVeto && pass3rLVeto && outputVarMET[2]  > 12.0 && (TMath::Abs(outputVarMET[2] -91.1876) > 15 || sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me)) passSystCuts[lType][MET] = true;

    if(passLSel) {
       
       if(passNjets && preselCuts && passBtagVeto && pass3rLVeto && passMass && dPhiDiLepJetCut && passMET) passCuts[lType][WWSEL] = true;

    }

    if(passLSel){
      double add = 1.;
      double addPU = 1.;
      addPU = nPUScaleFactor2012(fhDPU,sigEvent.npu_);
      add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
      add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
      if((sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
      add = add*leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_);

      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  							       fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        						       TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
      add = add*trigEff*addPU;

      double theWeight = sigEvent.scale1fb_*lumi*add;

      if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = theMET;
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
	else if(thePlot ==15) myVar = fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = sigEvent.njets_;
	else if(thePlot ==18) myVar = sigEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((sigEvent.lep1_+sigEvent.lep2_+sigEvent.jet1_+sigEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(sigEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(sigEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = DeltaPhi(sigEvent.dilep_.Phi(),sigEvent.trackMetPhi_)*180.0/TMath::Pi();
	else if(thePlot ==25) myVar = sigEvent.lep1_.Pt()*sigEvent.lep2_.Pt()/sigEvent.jet1_.Pt()/sigEvent.jet2_.Pt();
	else if(thePlot ==26) myVar = sigEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = (sigEvent.dilep_+sigEvent.jet1_).Pt();
	else if(thePlot ==30) myVar = DeltaPhi(sigEvent.dilep_.Phi() ,sigEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((sigEvent.jet1_+sigEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = sigEvent.jet1_.Pt()+ sigEvent.jet2_.Pt()+sigEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = sigEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = sigEvent.trackMet_*cos(sigEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = sigEvent.trackMet_*sin(sigEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(sigEvent.jet3_.Phi(),sigEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = sigEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = TMath::Min(deltaRJJGen,4.999);
      	histos->Fill(myVar,theWeight);
      } // end making plots

      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            nSigCut[i+j*nSelTypes]  += theWeight;
            nSigECut[i+j*nSelTypes] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            nSigCutSyst[i+j*nSelTypesSyst]  += theWeight;
            nSigECutSyst[i+j*nSelTypesSyst] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WWSEL]) {
        nSigCutSyst[EFFP+lType*nSelTypesSyst]  += theWeight          *addLepEffUp  /addLepEff;
        nSigECutSyst[EFFP+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        nSigCutSyst[EFFM+lType*nSelTypesSyst]  += theWeight          *addLepEffDown/addLepEff;
        nSigECutSyst[EFFM+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }
      if(passLSel){
        if(sigEvent.njets_ == 0 && passCuts[1][WWSEL]) events_0Jet[1] = events_0Jet[1] + theWeight;
	if(passCuts[1][WWSEL]) 	             histo_Higgs                           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	             histo_Higgs_CMS_MVALepEffBoundingUp   ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL]) 	             histo_Higgs_CMS_MVALepEffBoundingDown ->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_Higgs_CMS_MVAJESBoundingUp      ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_Higgs_CMS_MVAJESBoundingDown    ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_Higgs_CMS_MVALepResBoundingUp   ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_Higgs_CMS_MVALepResBoundingDown ->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_Higgs_CMS_MVAMETResBoundingUp   ->Fill(MVAVar[5], theWeight);;
      }
    } // if passCuts
  } // Loop over signal
  
  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    if(dataEvent.lep1_.Pt() < 1.0) continue;

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;

    bool passCuts[3][nSelTypes] = {{false},
                                   {false}};

    double theMET = dataEvent.met_; double theMETPHI = dataEvent.metPhi_; 

    int lType = 0;
    if     (dataEvent.lq1_ * dataEvent.lq2_ < 0) lType = 1;

    double usedMet = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_);
    bool   passMET = usedMet > metMin;
    if(useDYMVA == false){
      if     (dataEvent.njets_ == 0) passMET = passMET && (usedMet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (usedMet > 45. || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (dataEvent.met_ > 45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    } else {
      if     (dataEvent.njets_ == 0) passMET = passMET && (dataEvent.dymva_ >  0.88 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else if(dataEvent.njets_ == 1) passMET = passMET && (dataEvent.dymva_ >  0.84 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (dataEvent.met_   >  45.0 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);
    }
    bool dPhiDiLepJetCut = true;
    if(useDYMVA == false){
      if(dataEvent.njets_ <= 1) dPhiDiLepJetCut = dataEvent.jet1_.Pt() <= 15. || dataEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
      	                                         dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
    	                                         dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;
    }
    if(dataEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((dataEvent.jet1_+dataEvent.jet2_).Phi(),dataEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me;

    double ptLLCut = ptLLMin; if(dataEvent.type_ == SmurfTree::ee || dataEvent.type_ == SmurfTree::mm) ptLLCut = 45;
    bool passNjets         = dataEvent.njets_ == nJetsType; if(nJetsType == 3) passNjets = dataEvent.njets_ <= 1;
    bool preselCuts        = dataEvent.lep1_.Pt() > ptLMin && dataEvent.lep2_.Pt() > ptLMin && dataEvent.dilep_.Pt() > ptLLCut;
    bool passBtagVeto      = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = dataEvent.dilep_.M() > 12 && (TMath::Abs(dataEvent.dilep_.M()-91.1876) > 15 || dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me);

    bool passLSel = false;
    if     (lSel == 0 && dataEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && dataEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && dataEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && dataEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                     passLSel = true;
    else if(lSel == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (dataEvent.type_ == SmurfTree::me || dataEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passLSel) {
       
       if(passNjets && preselCuts && passBtagVeto && pass3rLVeto && passMass && dPhiDiLepJetCut && passMET) passCuts[lType][WWSEL] = true;

    }

    if(passLSel) {

      if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = theMET;
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
	else if(thePlot ==15) myVar = fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = dataEvent.njets_;
	else if(thePlot ==18) myVar = dataEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(dataEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = DeltaPhi(dataEvent.dilep_.Phi(),dataEvent.trackMetPhi_)*180.0/TMath::Pi();
	else if(thePlot ==25) myVar = dataEvent.lep1_.Pt()*dataEvent.lep2_.Pt()/dataEvent.jet1_.Pt()/dataEvent.jet2_.Pt();
	else if(thePlot ==26) myVar = dataEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = (dataEvent.dilep_+dataEvent.jet1_).Pt();
	else if(thePlot ==30) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,dataEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((dataEvent.jet1_+dataEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = dataEvent.jet1_.Pt()+ dataEvent.jet2_.Pt()+dataEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = dataEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = dataEvent.trackMet_*cos(dataEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = dataEvent.trackMet_*sin(dataEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(dataEvent.jet3_.Phi(),dataEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = dataEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = TMath::Min(0.0,4.999);
      	histo5->Fill(myVar,1.0);
      } // end making plots

      double outputVar[15];
      makeSystematicEffects(dataEvent.lid1_, dataEvent.lid2_, dataEvent.lep1_, dataEvent.lep2_, dataEvent.dilep_, 
                            dataEvent.mt_, theMET, theMETPHI, 
                            dataEvent.trackMet_, dataEvent.trackMetPhi_, 
			    dataEvent.njets_, dataEvent.jet1_, dataEvent.jet2_, 
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][WWSEL]){
	histo_Data->Fill(MVAVar[0], 1.0);
      }

      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            nSelectedData[i+j*nSelTypes]  += 1.0;
          }
        }
      }

    } // if passCuts
  } // End loop data

  printf("gen_eff: %f / %f = %f | rec_eff: %f / %f = %f\n",genLevelNorm[1],genLevelNorm[0],genLevelNorm[1]/genLevelNorm[0],
                                                           genLevelNorm[2],genLevelNorm[3],genLevelNorm[2]/genLevelNorm[3]);
  char output[200];
  sprintf(output,Form("histo_nice%s.root",ECMsb.Data()));	 
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
    double nOldH[6] = {histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights(),histos->GetSumOfWeights()};
    for(int i=1; i<=histo0->GetNbinsX(); i++){
      if(histo0->GetBinContent(i) < 0) {histo0->SetBinContent(i,0.000001);histo0->SetBinError(i,0.000001);}
      if(histo1->GetBinContent(i) < 0) {histo1->SetBinContent(i,0.000001);histo1->SetBinError(i,0.000001);}
      if(histo2->GetBinContent(i) < 0) {histo2->SetBinContent(i,0.000001);histo2->SetBinError(i,0.000001);}
      if(histo3->GetBinContent(i) < 0) {histo3->SetBinContent(i,0.000001);histo3->SetBinError(i,0.000001);}
      if(histo4->GetBinContent(i) < 0) {histo4->SetBinContent(i,0.000001);histo4->SetBinError(i,0.000001);}
      if(histos->GetBinContent(i) < 0) {histos->SetBinContent(i,0.000001);histos->SetBinError(i,0.000001);}
    }
    if(nOldH[0] > 0) histo0->Scale(nOldH[0]/histo0->GetSumOfWeights());
    if(nOldH[1] > 0) histo1->Scale(nOldH[1]/histo1->GetSumOfWeights());
    if(nOldH[2] > 0) histo2->Scale(nOldH[2]/histo2->GetSumOfWeights());
    if(nOldH[3] > 0) histo3->Scale(nOldH[3]/histo3->GetSumOfWeights());
    if(nOldH[4] > 0) histo4->Scale(nOldH[4]/histo4->GetSumOfWeights());
    if(nOldH[5] > 0) histos->Scale(nOldH[5]/histos->GetSumOfWeights());

    printf("histo -> s: %8.2f d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f\n",histos->GetSumOfWeights(),histo5->GetSumOfWeights(),
    histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights());

    histos->Write();
    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();

  outFilePlotsNote->Close();
  
  const unsigned int nBkg = 11;
  double nTot[nSelTypes*2]; double nETot[nSelTypes*2];
  double bgdCombined[nSelTypes*2][nBkg],bgdCombinedE[nSelTypes*2][nBkg];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) {bgdCombined[i][j] = 0.0; bgdCombinedE[i][j] = 0.0;}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("selection: %s\n",selTypeName[i].Data());
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("data(%2d): %f\n",i,nSelectedData[i]);
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("nSigCut(%2d): %11.3f +/- %8.3f\n",i,nSigCut[i],sqrt(nSigECut[i]));
    nTot[i] = 0.0; nETot[i] = 0.0;
    for(int j=0; j<45; j++){

      if(showSignalOnly == false || i%nSelTypes == WWSEL) if(bgdDecay[i][j] != 0) printf("bdg(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecay[i][j],sqrt(weiDecay[i][j]));

      nTot[i]  += bgdDecay[i][j];
      nETot[i] += weiDecay[i][j];

      if     (j == 29)            {bgdCombined[i][ 0] += bgdDecay[i][j]; bgdCombinedE[i][ 0] += weiDecay[i][j];}
      else if(j == 30)            {bgdCombined[i][ 1] += bgdDecay[i][j]; bgdCombinedE[i][ 1] += weiDecay[i][j];}
      else if(j == 27 || j== 28)  {bgdCombined[i][ 2] += bgdDecay[i][j]; bgdCombinedE[i][ 2] += weiDecay[i][j];}
      else if(j == 21)		  {bgdCombined[i][ 3] += bgdDecay[i][j]; bgdCombinedE[i][ 3] += weiDecay[i][j];}
      else if(j ==  5 || j == 13) {bgdCombined[i][ 4] += bgdDecay[i][j]; bgdCombinedE[i][ 4] += weiDecay[i][j];}
      else if(j ==  9)            {bgdCombined[i][ 5] += bgdDecay[i][j]; bgdCombinedE[i][ 5] += weiDecay[i][j];}
      else if(j == 10)            {bgdCombined[i][ 6] += bgdDecay[i][j]; bgdCombinedE[i][ 6] += weiDecay[i][j];}
      else if(j == 20)            {bgdCombined[i][ 7] += bgdDecay[i][j]; bgdCombinedE[i][ 7] += weiDecay[i][j];}
      else if(j == 19)            {bgdCombined[i][ 8] += bgdDecay[i][j]; bgdCombinedE[i][ 8] += weiDecay[i][j];}
      else if(j == 1)             {bgdCombined[i][ 9] += bgdDecay[i][j]; bgdCombinedE[i][ 9] += weiDecay[i][j];}
      else if(j == 23)            {bgdCombined[i][10] += bgdDecay[i][j]; bgdCombinedE[i][10] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]+bgdCombined[i][7]+bgdCombined[i][8]+bgdCombined[i][9]+bgdCombined[i][10]-nTot[i])/nTot[i] > 0.00001) 
                    {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]+bgdCombined[i][7]+bgdCombined[i][8]+bgdCombined[i][9]+bgdCombined[i][10]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(qqWW) = %11.3f +/- %8.3f\n",bgdCombined[i][ 0],sqrt(bgdCombinedE[i][ 0]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(ggWW) = %11.3f +/- %8.3f\n",bgdCombined[i][ 1],sqrt(bgdCombinedE[i][ 1]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xxVV) = %11.3f +/- %8.3f\n",bgdCombined[i][ 2],sqrt(bgdCombinedE[i][ 2]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xVVV) = %11.3f +/- %8.3f\n",bgdCombined[i][ 3],sqrt(bgdCombinedE[i][ 3]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(.Top) = %11.3f +/- %8.3f\n",bgdCombined[i][ 4],sqrt(bgdCombinedE[i][ 4]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(Zjet) = %11.3f +/- %8.3f\n",bgdCombined[i][ 5],sqrt(bgdCombinedE[i][ 5]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xZtt) = %11.3f +/- %8.3f\n",bgdCombined[i][ 6],sqrt(bgdCombinedE[i][ 6]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(Wg3l) = %11.3f +/- %8.3f\n",bgdCombined[i][ 7],sqrt(bgdCombinedE[i][ 7]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xxWg) = %11.3f +/- %8.3f\n",bgdCombined[i][ 8],sqrt(bgdCombinedE[i][ 8]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWjE) = %11.3f +/- %8.3f\n",bgdCombined[i][ 9],sqrt(bgdCombinedE[i][ 9]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWjM) = %11.3f +/- %8.3f\n",bgdCombined[i][10],sqrt(bgdCombinedE[i][10]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("*******************************\n");
  }

  if(showSignalOnly == false) printf("+++++++++++++++++++++++++++++++\n");
  double nTotSyst[nSelTypesSyst*2]; double nETotSyst[nSelTypesSyst*2];
  double bgdCombinedSyst[nSelTypesSyst*2][nBkg],bgdCombinedESyst[nSelTypesSyst*2][nBkg];
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    if(showSignalOnly == false) printf("selectionSyst: %s\n",selTypeNameSyst[i].Data());
    for(unsigned int j=0; j<nBkg; j++) {bgdCombinedSyst[i][j] = 0.0; bgdCombinedESyst[i][j] = 0.0;}
    if(showSignalOnly == false) printf("nSigCutSyst(%2d): %11.3f +/- %8.3f\n",i,nSigCutSyst[i],sqrt(nSigECutSyst[i]));
    nTotSyst[i] = 0.0; nETotSyst[i] = 0.0;
    for(int j=0; j<45; j++){
      if(showSignalOnly == false) if(bgdDecaySyst[i][j] != 0) printf("bdgSyst(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecaySyst[i][j],sqrt(weiDecaySyst[i][j]));
      nTotSyst[i]  += bgdDecaySyst[i][j];
      nETotSyst[i] += weiDecaySyst[i][j];

      if     (j == 29)            {bgdCombinedSyst[i][ 0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 0] += weiDecaySyst[i][j];}
      else if(j == 30)            {bgdCombinedSyst[i][ 1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 1] += weiDecaySyst[i][j];}
      else if(j == 27 || j== 28)  {bgdCombinedSyst[i][ 2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 2] += weiDecaySyst[i][j];}
      else if(j == 21)		  {bgdCombinedSyst[i][ 3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 3] += weiDecaySyst[i][j];}
      else if(j ==  5 || j == 13) {bgdCombinedSyst[i][ 4] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 4] += weiDecaySyst[i][j];}
      else if(j ==  9)            {bgdCombinedSyst[i][ 5] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 5] += weiDecaySyst[i][j];}
      else if(j == 10)            {bgdCombinedSyst[i][ 6] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 6] += weiDecaySyst[i][j];}
      else if(j == 20)            {bgdCombinedSyst[i][ 7] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 7] += weiDecaySyst[i][j];}
      else if(j == 19)            {bgdCombinedSyst[i][ 8] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 8] += weiDecaySyst[i][j];}
      else if(j == 1)             {bgdCombinedSyst[i][ 9] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][ 9] += weiDecaySyst[i][j];}
      else if(j == 23)            {bgdCombinedSyst[i][10] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][10] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]+bgdCombinedSyst[i][7]+bgdCombinedSyst[i][8]+bgdCombinedSyst[i][9]+bgdCombinedSyst[i][10]-nTotSyst[i])/nTotSyst[i] > 0.00001) 
                    {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]+bgdCombinedSyst[i][7]+bgdCombinedSyst[i][8]+bgdCombinedSyst[i][9]+bgdCombinedSyst[i][10]);assert(0);}
    if(showSignalOnly == false) printf("------\n");
    if(showSignalOnly == false) printf("bgdSyst(qqWW) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 0],sqrt(bgdCombinedESyst[i][ 0]));
    if(showSignalOnly == false) printf("bgdSyst(ggWW) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 1],sqrt(bgdCombinedESyst[i][ 1]));
    if(showSignalOnly == false) printf("bgdSyst(xxVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 2],sqrt(bgdCombinedESyst[i][ 2]));
    if(showSignalOnly == false) printf("bgdSyst(xVVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 3],sqrt(bgdCombinedESyst[i][ 3]));
    if(showSignalOnly == false) printf("bgdSyst(.Top) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 4],sqrt(bgdCombinedESyst[i][ 4]));
    if(showSignalOnly == false) printf("bgdSyst(Zjet) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 5],sqrt(bgdCombinedESyst[i][ 5]));
    if(showSignalOnly == false) printf("bgdSyst(xZtt) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 6],sqrt(bgdCombinedESyst[i][ 6]));
    if(showSignalOnly == false) printf("bgdSyst(Wg3l) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 7],sqrt(bgdCombinedESyst[i][ 7]));
    if(showSignalOnly == false) printf("bgdSyst(xxWg) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 8],sqrt(bgdCombinedESyst[i][ 8]));
    if(showSignalOnly == false) printf("bgdSyst(xWjE) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][ 9],sqrt(bgdCombinedESyst[i][ 9]));
    if(showSignalOnly == false) printf("bgdSyst(xWjM) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][10],sqrt(bgdCombinedESyst[i][10]));
    if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }
  double NFinal[nBkg+1]  = {bgdCombined[WWSEL+nSelTypes][0],bgdCombined[WWSEL+nSelTypes][1],bgdCombined[WWSEL+nSelTypes][2],
                            bgdCombined[WWSEL+nSelTypes][3],bgdCombined[WWSEL+nSelTypes][4],bgdCombined[WWSEL+nSelTypes][5],
                            bgdCombined[WWSEL+nSelTypes][6],bgdCombined[WWSEL+nSelTypes][7],bgdCombined[WWSEL+nSelTypes][8],
                            bgdCombined[WWSEL+nSelTypes][9],bgdCombined[WWSEL+nSelTypes][10],nSigCut[WWSEL+nSelTypes]};
  
  double systEffect[nSelTypesSyst][nBkg+1];
  for(unsigned int i=0 ; i<nSelTypesSyst; i++){
    for(unsigned int j=0 ; j<nBkg; j++){
      if(bgdCombinedE[WWSEL+nSelTypes][j] > 0){
        systEffect[i][j] = bgdCombinedSyst[i+nSelTypesSyst][j]/bgdCombined[WWSEL+nSelTypes][j];
        if(systEffect[i][j] < 1) systEffect[i][j] = 1.0/systEffect[i][j];
      } else {systEffect[i][j] = 1.0;}
    }
    systEffect[i][nBkg] = nSigCutSyst[i+nSelTypesSyst]/nSigCut[WWSEL+nSelTypes];    
    if(systEffect[i][nBkg] < 1) systEffect[i][nBkg] = 1.0/systEffect[i][nBkg];
  }
  if(showSignalOnly == false) printf("Syst(xHig) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][nBkg]-1,systEffect[JESDOWN][nBkg]-1,systEffect[LEPP][nBkg]-1,systEffect[LEPM][nBkg]-1,systEffect[MET][nBkg]-1,systEffect[EFFP][nBkg]-1,systEffect[EFFM][nBkg]-1);
  if(showSignalOnly == false) printf("Syst(qqWW) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 0]-1,systEffect[JESDOWN][ 0]-1,systEffect[LEPP][ 0]-1,systEffect[LEPM][ 0]-1,systEffect[MET][ 0]-1,systEffect[EFFP][ 0]-1,systEffect[EFFM][ 0]-1);
  if(showSignalOnly == false) printf("Syst(ggWW) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 1]-1,systEffect[JESDOWN][ 1]-1,systEffect[LEPP][ 1]-1,systEffect[LEPM][ 1]-1,systEffect[MET][ 1]-1,systEffect[EFFP][ 1]-1,systEffect[EFFM][ 1]-1);
  if(showSignalOnly == false) printf("Syst(xxVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 2]-1,systEffect[JESDOWN][ 2]-1,systEffect[LEPP][ 2]-1,systEffect[LEPM][ 2]-1,systEffect[MET][ 2]-1,systEffect[EFFP][ 2]-1,systEffect[EFFM][ 2]-1);
  if(showSignalOnly == false) printf("Syst(xVVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 3]-1,systEffect[JESDOWN][ 3]-1,systEffect[LEPP][ 3]-1,systEffect[LEPM][ 3]-1,systEffect[MET][ 3]-1,systEffect[EFFP][ 3]-1,systEffect[EFFM][ 3]-1);
  if(showSignalOnly == false) printf("Syst(.Top) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 4]-1,systEffect[JESDOWN][ 4]-1,systEffect[LEPP][ 4]-1,systEffect[LEPM][ 4]-1,systEffect[MET][ 4]-1,systEffect[EFFP][ 4]-1,systEffect[EFFM][ 4]-1);
  if(showSignalOnly == false) printf("Syst(Zjet) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 5]-1,systEffect[JESDOWN][ 5]-1,systEffect[LEPP][ 5]-1,systEffect[LEPM][ 5]-1,systEffect[MET][ 5]-1,systEffect[EFFP][ 5]-1,systEffect[EFFM][ 5]-1);
  if(showSignalOnly == false) printf("Syst(xZtt) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 6]-1,systEffect[JESDOWN][ 6]-1,systEffect[LEPP][ 6]-1,systEffect[LEPM][ 6]-1,systEffect[MET][ 6]-1,systEffect[EFFP][ 6]-1,systEffect[EFFM][ 6]-1);
  if(showSignalOnly == false) printf("Syst(Wg3l) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 7]-1,systEffect[JESDOWN][ 7]-1,systEffect[LEPP][ 7]-1,systEffect[LEPM][ 7]-1,systEffect[MET][ 7]-1,systEffect[EFFP][ 7]-1,systEffect[EFFM][ 7]-1);
  if(showSignalOnly == false) printf("Syst(xxWg) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 8]-1,systEffect[JESDOWN][ 8]-1,systEffect[LEPP][ 8]-1,systEffect[LEPM][ 8]-1,systEffect[MET][ 8]-1,systEffect[EFFP][ 8]-1,systEffect[EFFM][ 8]-1);
  if(showSignalOnly == false) printf("Syst(xWjE) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][ 9]-1,systEffect[JESDOWN][ 9]-1,systEffect[LEPP][ 9]-1,systEffect[LEPM][ 9]-1,systEffect[MET][ 9]-1,systEffect[EFFP][ 9]-1,systEffect[EFFM][ 9]-1);
  if(showSignalOnly == false) printf("Syst(xWjM) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][10]-1,systEffect[JESDOWN][10]-1,systEffect[LEPP][10]-1,systEffect[LEPM][10]-1,systEffect[MET][10]-1,systEffect[EFFP][10]-1,systEffect[EFFM][10]-1);

  double WjetsSyst = 1.36;
  double pdf_qqbar[5] = {1.013,1.040,1.040,1.040,1.040}; // qqWW,VV,VVV,Wg3l,Wgamma
  double pdf_gg[2] = {1.008,1.080}; // ggWW,Higgs
  double QCDscale_VV[2] = {1.03,1.05};

  double XS_QCDscale_ggH[3] = {1.160, 0.920, 1.000};
  //double XS_QCDscale_WW[3]  = {1.035, 0.987, 1.000};
  if     (nJetsType == 1){
    pdf_qqbar[0] = 1.033;
    pdf_gg[0] = 1.007;
    XS_QCDscale_ggH[0] = 1.000;
    XS_QCDscale_ggH[1] = 1.280;
    XS_QCDscale_ggH[2] = 0.970;
    //XS_QCDscale_WW[0] = 1.000; XS_QCDscale_WW[1] = 1.076; XS_QCDscale_WW[2] = 0.914;
  }
  else if(nJetsType == 2){
    pdf_qqbar[0] = 1.033;
    pdf_gg[0] = 1.007;
    XS_QCDscale_ggH[0] = 1.000;
    XS_QCDscale_ggH[1] = 1.000;
    XS_QCDscale_ggH[2] = 1.350;
    //XS_QCDscale_WW[0] = 1.000; XS_QCDscale_WW[1] = 1.000; XS_QCDscale_WW[2] = 1.420;
  }
  else if(nJetsType == 3){
    pdf_qqbar[0] = 1.019;
    pdf_gg[0] = 1.008;
    XS_QCDscale_ggH[0] = 1.000;
    XS_QCDscale_ggH[1] = 1.0+0.28*(1.0-events_0Jet[1]/histo_Higgs->GetSumOfWeights());
    XS_QCDscale_ggH[2] = 1.0-0.03*(1.0-events_0Jet[1]/histo_Higgs->GetSumOfWeights());
    //XS_QCDscale_WW[0] = 1.000; XS_QCDscale_WW[1] = 1.0+0.076*(1.0-events_0Jet[0]/histo_qqWW->GetSumOfWeights()); 
    //                           XS_QCDscale_WW[2] = 1.0-0.086*(1.0-events_0Jet[0]/histo_qqWW->GetSumOfWeights());
  }
  double topXS_E = 1.0;double ZXS_E = 1.0;
  if(nJetsType != 3) {
    topXS_E = TopBkgScaleFactorKappa(nJetsType);
    ZXS_E   = DYBkgScaleFactorKappa (nJetsType);
  } else {
    topXS_E   = 1.0+((TopBkgScaleFactorKappa(0)-1.0)*     events_0Jet[2]/histo_Top  ->GetSumOfWeights()+
                     (TopBkgScaleFactorKappa(1)-1.0)*(1.0-events_0Jet[2]/histo_Top  ->GetSumOfWeights()));
    if(histo_Zjets->GetSumOfWeights() > 0){
      ZXS_E   = 1.0+((DYBkgScaleFactorKappa (0)-1.0)*     events_0Jet[3]/histo_Zjets->GetSumOfWeights()+
                     (DYBkgScaleFactorKappa (1)-1.0)*(1.0-events_0Jet[3]/histo_Zjets->GetSumOfWeights()));
    }
  }
  
  double nOldNorm[7] = {histo_WjetsE->GetSumOfWeights(),histo_WjetsE_CMS_MVAWEBoundingUp->GetSumOfWeights(),
                        histo_WjetsM->GetSumOfWeights(),histo_WjetsM_CMS_MVAWMBoundingUp->GetSumOfWeights(),
                        histo_qqWW_CMS_MVAWWNLOBoundingUp->GetSumOfWeights(),histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights(),
		        histo_qqWW_CMS_MVAWWBoundingUp->GetSumOfWeights()};
  for(int i=1; i<=histo_WjetsE->GetNbinsX(); i++){
    if(histo_WjetsE->GetBinContent(i)                        < 0) {histo_WjetsE		     	      ->SetBinContent(i,0.000001);histo_WjetsE                       ->SetBinError(i,0.000001);}
    if(histo_WjetsE_CMS_MVAWEBoundingUp->GetBinContent(i)    < 0) {histo_WjetsE_CMS_MVAWEBoundingUp   ->SetBinContent(i,0.000001);histo_WjetsE_CMS_MVAWEBoundingUp   ->SetBinError(i,0.000001);}
    if(histo_WjetsM->GetBinContent(i)                        < 0) {histo_WjetsM		     	      ->SetBinContent(i,0.000001);histo_WjetsM                       ->SetBinError(i,0.000001);}
    if(histo_WjetsM_CMS_MVAWMBoundingUp->GetBinContent(i)    < 0) {histo_WjetsM_CMS_MVAWMBoundingUp   ->SetBinContent(i,0.000001);histo_WjetsM_CMS_MVAWMBoundingUp   ->SetBinError(i,0.000001);}
    if(histo_qqWW_CMS_MVAWWNLOBoundingUp->GetBinContent(i)   < 0) {histo_qqWW_CMS_MVAWWNLOBoundingUp  ->SetBinContent(i,0.000001);histo_qqWW_CMS_MVAWWNLOBoundingUp  ->SetBinError(i,0.000001);}
    if(histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(i) < 0) {histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,0.000001);histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinError(i,0.000001);}
    if(histo_qqWW_CMS_MVAWWBoundingUp->GetBinContent(i)      < 0) {histo_qqWW_CMS_MVAWWBoundingUp     ->SetBinContent(i,0.000001);histo_qqWW_CMS_MVAWWBoundingUp     ->SetBinError(i,0.000001);}
  }
  histo_WjetsE                       ->Scale(nOldNorm[0]/histo_WjetsE			    ->GetSumOfWeights());
  histo_WjetsE_CMS_MVAWEBoundingUp   ->Scale(nOldNorm[1]/histo_WjetsE_CMS_MVAWEBoundingUp   ->GetSumOfWeights());
  histo_WjetsM                       ->Scale(nOldNorm[2]/histo_WjetsM			    ->GetSumOfWeights());
  histo_WjetsM_CMS_MVAWMBoundingUp   ->Scale(nOldNorm[3]/histo_WjetsM_CMS_MVAWMBoundingUp   ->GetSumOfWeights());
  histo_qqWW_CMS_MVAWWNLOBoundingUp  ->Scale(nOldNorm[4]/histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetSumOfWeights());
  histo_qqWW_CMS_MVAWWNLOBoundingDown->Scale(nOldNorm[5]/histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights());
  histo_qqWW_CMS_MVAWWBoundingUp     ->Scale(nOldNorm[6]/histo_qqWW_CMS_MVAWWBoundingUp     ->GetSumOfWeights());
  
  for(int i=1; i<=histo_Higgs->GetNbinsX(); i++){
    double mean = histo_qqWW  ->GetBinContent(i);
    double up	= histo_qqWW_CMS_MVAWWBoundingUp->GetBinContent(i);
    double diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_qqWW_CMS_MVAWWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_qqWW_CMS_MVAWWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    double meanNLO    = histo_qqWW_CMS_MVAWWBoundingUp     ->GetBinContent(i);
    double theNLOUp   = histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetBinContent(i);
    //double theNLODown = histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(i);
    if(meanNLO  > 0 && histo_qqWW_CMS_MVAWWBoundingUp	->GetBinError(i)/meanNLO  < 0.5 &&
       theNLOUp > 0 && histo_qqWW_CMS_MVAWWNLOBoundingUp->GetBinError(i)/theNLOUp < 0.5 &&
    	meanNLO/histo_qqWW_CMS_MVAWWBoundingUp->GetSumOfWeights() > 0.000001) {
      histo_qqWW_CMS_MVAWWNLOBoundingUp->SetBinContent(i,histo_qqWW->GetBinContent(i)*histo_qqWW_CMS_MVAWWNLOBoundingUp->GetBinContent(i)/meanNLO);
    }
    else {
      histo_qqWW_CMS_MVAWWNLOBoundingUp->SetBinContent(i,histo_qqWW->GetBinContent(i));
    }
    
    mean = histo_qqWW  ->GetBinContent(i);
    up   = histo_qqWW_CMS_MVAWWNLOBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    //if(meanNLO    > 0 && histo_qqWW_CMS_MVAWWBoundingUp     ->GetBinError(i)/meanNLO	< 0.5 &&
    //   theNLODown > 0 && histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinError(i)/theNLODown < 0.5 &&
    //	meanNLO/histo_qqWW_CMS_MVAWWBoundingUp->GetSumOfWeights() > 0.000001) {
    //  histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,histo_qqWW->GetBinContent(i)*histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(i)/meanNLO);
    //}
    //else {
    //  histo_qqWW_CMS_MVAWWNLOBoundingDown->SetBinContent(i,histo_qqWW->GetBinContent(i));
    //}

    double factorUp = +1.0; double factorDown = -1.0;
    histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_Higgs->GetBinContent(i)+factorUp  *histo_Higgs	 ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown->SetBinContent(i,TMath::Max(histo_Higgs->GetBinContent(i)+factorDown*histo_Higgs	 ->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorUp  *histo_qqWW	 ->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorDown*histo_qqWW	 ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorUp  *histo_ggWW	 ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorDown*histo_ggWW	 ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_VV    	->GetBinContent(i)+factorUp  *histo_VV   	 ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_VV    	->GetBinContent(i)+factorDown*histo_VV   	 ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp        	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorUp  *histo_VVV  	 ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown      	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorDown*histo_VVV  	 ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_Top    	->GetBinContent(i)+factorUp  *histo_Top   	 ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_Top    	->GetBinContent(i)+factorDown*histo_Top   	 ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_Zjets    	->GetBinContent(i)+factorUp  *histo_Zjets   	 ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_Zjets    	->GetBinContent(i)+factorDown*histo_Zjets   	 ->GetBinError(i),0.000001));
    histo_Ztt_CMS_MVAZttStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_Ztt    	->GetBinContent(i)+factorUp  *histo_Ztt   	 ->GetBinError(i),0.000001));
    histo_Ztt_CMS_MVAZttStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_Ztt    	->GetBinContent(i)+factorDown*histo_Ztt   	 ->GetBinError(i),0.000001));
    histo_Wg3l_CMS_MVAWg3lStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_Wg3l    	->GetBinContent(i)+factorUp  *histo_Wg3l   	 ->GetBinError(i),0.000001));
    histo_Wg3l_CMS_MVAWg3lStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_Wg3l    	->GetBinContent(i)+factorDown*histo_Wg3l   	 ->GetBinError(i),0.000001));
    histo_Wgamma_CMS_MVAWgammaStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_Wgamma    	->GetBinContent(i)+factorUp  *histo_Wgamma   	 ->GetBinError(i),0.000001));
    histo_Wgamma_CMS_MVAWgammaStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_Wgamma    	->GetBinContent(i)+factorDown*histo_Wgamma   	 ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingUp    	    ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorUp  *histo_WjetsE        ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingDown  	    ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorDown*histo_WjetsE        ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingUp    	    ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorUp  *histo_WjetsM        ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingDown  	    ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorDown*histo_WjetsM        ->GetBinError(i),0.000001));

    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]  ->Add(histo_Higgs   ); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_Higgs   ->GetBinContent(i)+factorUp  *histo_Higgs->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]->Add(histo_Higgs   ); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_Higgs   ->GetBinContent(i)+factorDown*histo_Higgs->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[i-1]	    ->Add(histo_qqWW       ); histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_qqWW       ->GetBinContent(i)+factorUp  *histo_qqWW    ->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[i-1]	    ->Add(histo_qqWW       ); histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_qqWW       ->GetBinContent(i)+factorDown*histo_qqWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingBinUp[i-1]	    ->Add(histo_ggWW       ); histo_ggWW_CMS_MVAggWWStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_ggWW       ->GetBinContent(i)+factorUp  *histo_ggWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingBinDown[i-1]	    ->Add(histo_ggWW       ); histo_ggWW_CMS_MVAggWWStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_ggWW       ->GetBinContent(i)+factorDown*histo_ggWW    ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingBinUp[i-1]	            ->Add(histo_VV          ); histo_VV_CMS_MVAVVStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_VV          ->GetBinContent(i)+factorUp  *histo_VV       ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingBinDown[i-1]	            ->Add(histo_VV          ); histo_VV_CMS_MVAVVStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_VV          ->GetBinContent(i)+factorDown*histo_VV       ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	            ->Add(histo_VVV         ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]              ->SetBinContent(i,TMath::Max(histo_VVV         ->GetBinContent(i)+factorUp  *histo_VVV      ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]            ->Add(histo_VVV         ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]            ->SetBinContent(i,TMath::Max(histo_VVV         ->GetBinContent(i)+factorDown*histo_VVV      ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingBinUp[i-1]	            ->Add(histo_Top          ); histo_Top_CMS_MVATopStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_Top          ->GetBinContent(i)+factorUp  *histo_Top       ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingBinDown[i-1]	            ->Add(histo_Top          ); histo_Top_CMS_MVATopStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_Top          ->GetBinContent(i)+factorDown*histo_Top       ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[i-1]	            ->Add(histo_Zjets          ); histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_Zjets          ->GetBinContent(i)+factorUp  *histo_Zjets       ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[i-1]	            ->Add(histo_Zjets          ); histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_Zjets          ->GetBinContent(i)+factorDown*histo_Zjets       ->GetBinError(i),0.000001));
    histo_Ztt_CMS_MVAZttStatBoundingBinUp[i-1]	            ->Add(histo_Ztt          ); histo_Ztt_CMS_MVAZttStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_Ztt          ->GetBinContent(i)+factorUp  *histo_Ztt       ->GetBinError(i),0.000001));
    histo_Ztt_CMS_MVAZttStatBoundingBinDown[i-1]	            ->Add(histo_Ztt          ); histo_Ztt_CMS_MVAZttStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_Ztt          ->GetBinContent(i)+factorDown*histo_Ztt       ->GetBinError(i),0.000001));
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp[i-1]	            ->Add(histo_Wg3l          ); histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_Wg3l          ->GetBinContent(i)+factorUp  *histo_Wg3l       ->GetBinError(i),0.000001));
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown[i-1]	            ->Add(histo_Wg3l          ); histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_Wg3l          ->GetBinContent(i)+factorDown*histo_Wg3l       ->GetBinError(i),0.000001));
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp[i-1]	            ->Add(histo_Wgamma          ); histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_Wgamma          ->GetBinContent(i)+factorUp  *histo_Wgamma       ->GetBinError(i),0.000001));
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown[i-1]	            ->Add(histo_Wgamma          ); histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_Wgamma          ->GetBinContent(i)+factorDown*histo_Wgamma       ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[i-1]          ->Add(histo_WjetsE       ); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_WjetsE       ->GetBinContent(i)+factorUp  *histo_WjetsE    ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[i-1]        ->Add(histo_WjetsE       ); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_WjetsE       ->GetBinContent(i)+factorDown*histo_WjetsE    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[i-1]          ->Add(histo_WjetsM       ); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_WjetsM       ->GetBinContent(i)+factorUp  *histo_WjetsM    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[i-1]        ->Add(histo_WjetsM       ); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_WjetsM       ->GetBinContent(i)+factorDown*histo_WjetsM    ->GetBinError(i),0.000001));
  }
  double mean,up,diff;
  for(int i=1; i<=histo_Higgs->GetNbinsX(); i++){
    
    mean = histo_Higgs			   ->GetBinContent(i);
    up   = histo_Higgs_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Higgs_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Higgs_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_qqWW			   ->GetBinContent(i);
    up   = histo_qqWW_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_qqWW_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_qqWW_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_ggWW			   ->GetBinContent(i);
    up   = histo_ggWW_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ggWW_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ggWW_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_VV			   ->GetBinContent(i);
    up   = histo_VV_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VV		       ->GetBinContent(i);
    up   = histo_VV_CMS_VVNLOBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VV_CMS_VVNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VV_CMS_VVNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Top			   ->GetBinContent(i);
    up   = histo_Top_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Top_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Top_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_Zjets			   ->GetBinContent(i);
    up   = histo_Zjets_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Zjets_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Zjets_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));  

    mean = histo_Ztt			   ->GetBinContent(i);
    up   = histo_Ztt_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Ztt_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Ztt_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Wg3l			   ->GetBinContent(i);
    up   = histo_Wg3l_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wg3l_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wg3l_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Wgamma			   ->GetBinContent(i);
    up   = histo_Wgamma_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wgamma_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WjetsE 		         ->GetBinContent(i);
    up   = histo_WjetsE_CMS_MVAWEBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WjetsE_CMS_MVAWEBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WjetsE_CMS_MVAWEBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WjetsM 		         ->GetBinContent(i);
    up   = histo_WjetsM_CMS_MVAWMBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WjetsM_CMS_MVAWMBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WjetsM_CMS_MVAWMBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
  }

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"wwana%2s.input%d_%4s.root",finalStateName,nJetsType,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  histo_Higgs  ->Write();
  histo_qqWW   ->Write();
  histo_ggWW   ->Write();
  histo_VV     ->Write();
  histo_VVV    ->Write();
  histo_Top    ->Write();
  histo_Zjets  ->Write();
  histo_Ztt    ->Write();
  histo_Wg3l   ->Write();
  histo_Wgamma ->Write();
  histo_WjetsE ->Write();
  histo_WjetsM ->Write();

  cout << histo_Data   ->GetSumOfWeights() << " ";
  cout << histo_Higgs  ->GetSumOfWeights() << " ";
  cout << histo_qqWW   ->GetSumOfWeights() << " ";
  cout << histo_ggWW   ->GetSumOfWeights() << " ";
  cout << histo_VV     ->GetSumOfWeights() << " ";
  cout << histo_VVV    ->GetSumOfWeights() << " ";
  cout << histo_Top    ->GetSumOfWeights() << " ";
  cout << histo_Zjets  ->GetSumOfWeights() << " ";
  cout << histo_Ztt    ->GetSumOfWeights() << " ";
  cout << histo_Wg3l   ->GetSumOfWeights() << " ";
  cout << histo_Wgamma ->GetSumOfWeights() << " ";
  cout << histo_WjetsE ->GetSumOfWeights() << " ";
  cout << histo_WjetsM ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingUp 	 ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingDown->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingDown      ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW	->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingUp	     ->GetBinContent(i)/histo_qqWW   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW	->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingDown      ->GetBinContent(i)/histo_qqWW   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingUp  	->GetBinContent(i)/histo_ggWW   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingDown      ->GetBinContent(i)/histo_ggWW   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingUp	->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingDown	->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown   	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingUp	->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingDown	->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAZjetsStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAZjetsStatBoundingUp	->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAZjetsStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAZjetsStatBoundingDown	->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAZttStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAZttStatBoundingUp	->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAZttStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAZttStatBoundingDown	->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAWg3lStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAWg3lStatBoundingUp	->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAWg3lStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAWg3lStatBoundingDown	->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAWgammaStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAWgammaStatBoundingUp	->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAWgammaStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAWgammaStatBoundingDown	->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->GetBinContent(i)/histo_WjetsE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->GetBinContent(i)/histo_WjetsE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->GetBinContent(i)/histo_WjetsM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->GetBinContent(i)/histo_WjetsM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_Higgs_CMS_MVALepEffBoundingUp    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Higgs     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffBoundingUp	  ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffBoundingUp	     ->GetBinContent(i)/histo_ggWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffBoundingDown	->GetBinContent(i)/histo_ggWW     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingDown        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_Higgs_CMS_MVALepResBoundingUp    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepResBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepResBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Higgs     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepResBoundingUp	  ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepResBoundingUp	     ->GetBinContent(i)/histo_ggWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepResBoundingDown	->GetBinContent(i)/histo_ggWW     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingDown        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_Higgs_CMS_MVAMETResBoundingUp    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs	->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETResBoundingUp	  ->GetBinContent(i)/histo_Higgs	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETResBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs	->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Higgs	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETResBoundingUp	     ->GetBinContent(i)/histo_qqWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETResBoundingDown	->GetBinContent(i)/histo_qqWW     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETResBoundingUp	     ->GetBinContent(i)/histo_ggWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETResBoundingDown	->GetBinContent(i)/histo_ggWW     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingDown        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_Higgs_CMS_MVAJESBoundingUp       ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingUp	      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingDown     ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingDown       ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingUp	  ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW  ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingDown   ->GetBinContent(i)/histo_qqWW	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingUp     ->GetBinContent(i)/histo_ggWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW     ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingDown   ->GetBinContent(i)/histo_ggWW	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV	   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp             ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown           ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown	     ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_Top   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Zjets	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_Zjets   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Ztt_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Ztt	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Ztt_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_Ztt   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wg3l_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wg3l	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wg3l_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_Wg3l   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wgamma_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Wgamma	   ->GetBinContent(i)>0)printf("%5.1f ",histo_Wgamma_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_Wgamma   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties GEN\n");
  histo_VV_CMS_VVNLOBoundingUp            ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_VVNLOBoundingUp	     ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_VVNLOBoundingDown          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV   ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_VVNLOBoundingDown       ->GetBinContent(i)/histo_VV   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWEBoundingUp       ->GetBinContent(i)/histo_WjetsE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWEBoundingDown     ->GetBinContent(i)/histo_WjetsE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWMBoundingUp       ->GetBinContent(i)/histo_WjetsM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWMBoundingDown     ->GetBinContent(i)/histo_WjetsM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWWBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWWBoundingUp	 ->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWWBoundingDown        ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWWBoundingDown   ->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWWNLOBoundingUp       ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWWNLOBoundingDown     ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  histo_qqWW_CMS_WWNLOQUp	->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_WWNLOQUp  ->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_WWNLOQDown	->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_WWNLOQDown->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_WWNLORUp	->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_WWNLORUp  ->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_WWNLORDown	->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_WWNLORDown->GetBinContent(i)/histo_qqWW->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=0; nb<nBinMVA; nb++){
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]     ->Write();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]   ->Write();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]	      ->Write();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]     ->Write();
    histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]	      ->Write();
    histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]     ->Write();
    histo_VV_CMS_MVAVVStatBoundingBinUp[nb]	      ->Write();
    histo_VV_CMS_MVAVVStatBoundingBinDown[nb]	      ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	      ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	      ->Write();
    histo_Top_CMS_MVATopStatBoundingBinUp[nb]	      ->Write();
    histo_Top_CMS_MVATopStatBoundingBinDown[nb]	      ->Write();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]     ->Write();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]   ->Write();
    histo_Ztt_CMS_MVAZttStatBoundingBinUp[nb]	      ->Write();
    histo_Ztt_CMS_MVAZttStatBoundingBinDown[nb]	      ->Write();
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinUp[nb]	      ->Write();
    histo_Wg3l_CMS_MVAWg3lStatBoundingBinDown[nb]     ->Write();
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinUp[nb]   ->Write();
    histo_Wgamma_CMS_MVAWgammaStatBoundingBinDown[nb] ->Write();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]   ->Write();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb] ->Write();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]   ->Write();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb] ->Write();
  }

  if(showSignalOnly == false) {
    printf("nuisance VV | WjE | WjM | WWGen | WWNLO: %0.1f/%0.1f/%0.1f | %0.1f/%0.1f/%0.1f | %0.1f/%0.1f/%0.1f | %0.1f/%0.1f/%0.1f | %0.1f/%0.1f/%0.1f\n",
                                                              histo_VV    ->GetSumOfWeights(),histo_VV_CMS_VVNLOBoundingUp     ->GetSumOfWeights(),histo_VV_CMS_VVNLOBoundingDown     ->GetSumOfWeights(),
                                                              histo_WjetsE->GetSumOfWeights(),histo_WjetsE_CMS_MVAWEBoundingUp ->GetSumOfWeights(),histo_WjetsE_CMS_MVAWEBoundingDown ->GetSumOfWeights(),
                                                              histo_WjetsM->GetSumOfWeights(),histo_WjetsM_CMS_MVAWMBoundingUp ->GetSumOfWeights(),histo_WjetsM_CMS_MVAWMBoundingDown ->GetSumOfWeights(),
							      histo_qqWW  ->GetSumOfWeights(),histo_qqWW_CMS_MVAWWBoundingUp   ->GetSumOfWeights(),histo_qqWW_CMS_MVAWWBoundingDown   ->GetSumOfWeights(),
							      histo_qqWW  ->GetSumOfWeights(),histo_qqWW_CMS_MVAWWNLOBoundingUp->GetSumOfWeights(),histo_qqWW_CMS_MVAWWNLOBoundingDown->GetSumOfWeights());
  }
  histo_Data   ->Write();
  histo_qqWW   ->Write();
  histo_ggWW   ->Write();
  histo_Higgs  ->Write();
  histo_VV     ->Write();
  histo_VVV    ->Write();
  histo_Top    ->Write();
  histo_Zjets  ->Write();
  histo_Ztt    ->Write();
  histo_Wg3l   ->Write();
  histo_Wgamma ->Write();
  histo_WjetsE ->Write();
  histo_WjetsM ->Write();

  char outputLimitsShape[200];
  sprintf(outputLimitsShape,"histo_limits_wwana%2s_shape%d_%4s.txt",finalStateName,nJetsType,ECMsb.Data());
  ofstream newcardShape;
  newcardShape.open(outputLimitsShape);
  newcardShape << Form("imax 1 number of channels\n");
  newcardShape << Form("jmax * number of background\n");
  newcardShape << Form("kmax * number of nuisance parameters\n");
  newcardShape << Form("Observation %d\n",(int)nSelectedData[WWSEL+nSelTypes]);
  if(nBinMVA != 1){    
    newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
    newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
  }
  newcardShape << Form("bin wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s wwana%2s%4s\n",finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data());
  newcardShape << Form("process qqWW ggWW Higgs VV VVV Top Zjets Ztt Wg3l Wgamma WjetsE WjetsM\n");
  newcardShape << Form("process -1 0 1 2 3 4 5 6 7 8 9 10\n");
  newcardShape << Form("rate %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",NFinal[0],NFinal[1],NFinal[nBkg],NFinal[2],NFinal[3],NFinal[4],NFinal[5],NFinal[6],NFinal[7],NFinal[8],NFinal[9],NFinal[10]);
  newcardShape << Form("lumi_%4s                               lnN  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -     -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);		  
  if(nBinMVA != 1){
    newcardShape << Form("%s                                   shape  1.000 1.000 1.000 1.000 1.000   -     -	-   1.000 1.000   -	-  \n",effName);
    newcardShape << Form("%s                                   shape  1.000 1.000 1.000 1.000 1.000   -     -	-   1.000 1.000   -	-  \n",momName);
    newcardShape << Form("CMS_scale_met                        shape  1.000 1.000 1.000 1.000 1.000   -     -	-   1.000 1.000   -	-  \n");
    newcardShape << Form("CMS_scale_j                          shape  1.000 1.000 1.000 1.000 1.000   -     -	-   1.000 1.000   -	-  \n");		       
    newcardShape << Form("CMS_wwana_VVNLOBounding              shape    -     -   -   1.000   -     -     -     -     -     -     -     -  \n");
    newcardShape << Form("CMS_wwana_MVAWEBounding              shape    -     -	  -     -     -	    -	  -	-     -     -	1.000   -  \n");
    newcardShape << Form("CMS_wwana_MVAWMBounding              shape    -     -	  -     -     -	    -	  -	-     -     -	  -   1.000\n");
    //newcardShape << Form("CMS_wwana_MVAWWNLOBounding           shape  1.000   -	  -     -     -	    -     -     -     -     -	  -     -  \n");
    newcardShape << Form("CMS_wwana_WWNLOQ                     shape    -     -	  -     -     -	    -	  -	-     -     -	  -   1.000\n");
    newcardShape << Form("CMS_wwana_WWNLOR                     shape    -     -	  -     -     -	    -	  -	-     -     -	  -   1.000\n");
  } else {
    double systMVALepEff[7][3];
    systMVALepEff[0][0] = histo_qqWW_CMS_MVALepEffBoundingUp->GetBinContent(1)    /histo_qqWW  ->GetBinContent(1);
    systMVALepEff[0][1] = histo_qqWW_CMS_MVALepEffBoundingDown->GetBinContent(1)  /histo_qqWW  ->GetBinContent(1);
    systMVALepEff[1][0] = histo_ggWW_CMS_MVALepEffBoundingUp->GetBinContent(1)    /histo_ggWW  ->GetBinContent(1);
    systMVALepEff[1][1] = histo_ggWW_CMS_MVALepEffBoundingDown->GetBinContent(1)  /histo_ggWW  ->GetBinContent(1);
    systMVALepEff[2][0] = histo_Higgs_CMS_MVALepEffBoundingUp->GetBinContent(1)   /histo_Higgs ->GetBinContent(1);
    systMVALepEff[2][1] = histo_Higgs_CMS_MVALepEffBoundingDown->GetBinContent(1) /histo_Higgs ->GetBinContent(1);
    systMVALepEff[3][0] = histo_VV_CMS_MVALepEffBoundingUp->GetBinContent(1)      /histo_VV    ->GetBinContent(1);
    systMVALepEff[3][1] = histo_VV_CMS_MVALepEffBoundingDown->GetBinContent(1)    /histo_VV    ->GetBinContent(1);
    systMVALepEff[4][0] = histo_VVV_CMS_MVALepEffBoundingUp->GetBinContent(1)     /histo_VVV   ->GetBinContent(1);
    systMVALepEff[4][1] = histo_VVV_CMS_MVALepEffBoundingDown->GetBinContent(1)   /histo_VVV   ->GetBinContent(1);
    systMVALepEff[5][0] = histo_Wg3l_CMS_MVALepEffBoundingUp->GetBinContent(1)    /histo_Wg3l  ->GetBinContent(1);
    systMVALepEff[5][1] = histo_Wg3l_CMS_MVALepEffBoundingDown->GetBinContent(1)  /histo_Wg3l  ->GetBinContent(1);
    systMVALepEff[6][0] = histo_Wgamma_CMS_MVALepEffBoundingUp->GetBinContent(1)  /histo_Wgamma->GetBinContent(1);
    systMVALepEff[6][1] = histo_Wgamma_CMS_MVALepEffBoundingDown->GetBinContent(1)/histo_Wgamma->GetBinContent(1);
    for(int i=0; i<7; i++) for(int j=0; j<2; j++) if(systMVALepEff[i][j] < 1.0) systMVALepEff[i][j] = 1.0/systMVALepEff[i][j];
    for(int i=0; i<7; i++) systMVALepEff[i][2] = (systMVALepEff[i][0]+systMVALepEff[i][1])/2.0;

    double systMVALepRes[7][3];
    systMVALepRes[0][0] = histo_qqWW_CMS_MVALepResBoundingUp->GetBinContent(1)    /histo_qqWW  ->GetBinContent(1);
    systMVALepRes[0][1] = histo_qqWW_CMS_MVALepResBoundingDown->GetBinContent(1)  /histo_qqWW  ->GetBinContent(1);
    systMVALepRes[1][0] = histo_ggWW_CMS_MVALepResBoundingUp->GetBinContent(1)    /histo_ggWW  ->GetBinContent(1);
    systMVALepRes[1][1] = histo_ggWW_CMS_MVALepResBoundingDown->GetBinContent(1)  /histo_ggWW  ->GetBinContent(1);
    systMVALepRes[2][0] = histo_Higgs_CMS_MVALepResBoundingUp->GetBinContent(1)   /histo_Higgs ->GetBinContent(1);
    systMVALepRes[2][1] = histo_Higgs_CMS_MVALepResBoundingDown->GetBinContent(1) /histo_Higgs ->GetBinContent(1);
    systMVALepRes[3][0] = histo_VV_CMS_MVALepResBoundingUp->GetBinContent(1)      /histo_VV    ->GetBinContent(1);
    systMVALepRes[3][1] = histo_VV_CMS_MVALepResBoundingDown->GetBinContent(1)    /histo_VV    ->GetBinContent(1);
    systMVALepRes[4][0] = histo_VVV_CMS_MVALepResBoundingUp->GetBinContent(1)     /histo_VVV   ->GetBinContent(1);
    systMVALepRes[4][1] = histo_VVV_CMS_MVALepResBoundingDown->GetBinContent(1)   /histo_VVV   ->GetBinContent(1);
    systMVALepRes[5][0] = histo_Wg3l_CMS_MVALepResBoundingUp->GetBinContent(1)    /histo_Wg3l  ->GetBinContent(1);
    systMVALepRes[5][1] = histo_Wg3l_CMS_MVALepResBoundingDown->GetBinContent(1)  /histo_Wg3l  ->GetBinContent(1);
    systMVALepRes[6][0] = histo_Wgamma_CMS_MVALepResBoundingUp->GetBinContent(1)  /histo_Wgamma->GetBinContent(1);
    systMVALepRes[6][1] = histo_Wgamma_CMS_MVALepResBoundingDown->GetBinContent(1)/histo_Wgamma->GetBinContent(1);
    for(int i=0; i<7; i++) for(int j=0; j<2; j++) if(systMVALepRes[i][j] < 1.0) systMVALepRes[i][j] = 1.0/systMVALepRes[i][j];
    for(int i=0; i<7; i++) systMVALepRes[i][2] = (systMVALepRes[i][0]+systMVALepRes[i][1])/2.0;

    double systMVAMETRes[7][3];
    systMVAMETRes[0][0] = histo_qqWW_CMS_MVAMETResBoundingUp->GetBinContent(1)    /histo_qqWW  ->GetBinContent(1);
    systMVAMETRes[0][1] = histo_qqWW_CMS_MVAMETResBoundingDown->GetBinContent(1)  /histo_qqWW  ->GetBinContent(1);
    systMVAMETRes[1][0] = histo_ggWW_CMS_MVAMETResBoundingUp->GetBinContent(1)    /histo_ggWW  ->GetBinContent(1);
    systMVAMETRes[1][1] = histo_ggWW_CMS_MVAMETResBoundingDown->GetBinContent(1)  /histo_ggWW  ->GetBinContent(1);
    systMVAMETRes[2][0] = histo_Higgs_CMS_MVAMETResBoundingUp->GetBinContent(1)   /histo_Higgs ->GetBinContent(1);
    systMVAMETRes[2][1] = histo_Higgs_CMS_MVAMETResBoundingDown->GetBinContent(1) /histo_Higgs ->GetBinContent(1);
    systMVAMETRes[3][0] = histo_VV_CMS_MVAMETResBoundingUp->GetBinContent(1)      /histo_VV    ->GetBinContent(1);
    systMVAMETRes[3][1] = histo_VV_CMS_MVAMETResBoundingDown->GetBinContent(1)    /histo_VV    ->GetBinContent(1);
    systMVAMETRes[4][0] = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(1)     /histo_VVV   ->GetBinContent(1);
    systMVAMETRes[4][1] = histo_VVV_CMS_MVAMETResBoundingDown->GetBinContent(1)   /histo_VVV   ->GetBinContent(1);
    systMVAMETRes[5][0] = histo_Wg3l_CMS_MVAMETResBoundingUp->GetBinContent(1)    /histo_Wg3l  ->GetBinContent(1);
    systMVAMETRes[5][1] = histo_Wg3l_CMS_MVAMETResBoundingDown->GetBinContent(1)  /histo_Wg3l  ->GetBinContent(1);
    systMVAMETRes[6][0] = histo_Wgamma_CMS_MVAMETResBoundingUp->GetBinContent(1)  /histo_Wgamma->GetBinContent(1);
    systMVAMETRes[6][1] = histo_Wgamma_CMS_MVAMETResBoundingDown->GetBinContent(1)/histo_Wgamma->GetBinContent(1);
    for(int i=0; i<7; i++) for(int j=0; j<2; j++) if(systMVAMETRes[i][j] < 1.0) systMVAMETRes[i][j] = 1.0/systMVAMETRes[i][j];
    for(int i=0; i<7; i++) systMVAMETRes[i][2] = (systMVAMETRes[i][0]+systMVAMETRes[i][1])/2.0;

    double systMVAJES[7][3];
    systMVAJES[0][0] = histo_qqWW_CMS_MVAJESBoundingUp->GetBinContent(1)    /histo_qqWW  ->GetBinContent(1);
    systMVAJES[0][1] = histo_qqWW_CMS_MVAJESBoundingDown->GetBinContent(1)  /histo_qqWW  ->GetBinContent(1);
    systMVAJES[1][0] = histo_ggWW_CMS_MVAJESBoundingUp->GetBinContent(1)    /histo_ggWW  ->GetBinContent(1);
    systMVAJES[1][1] = histo_ggWW_CMS_MVAJESBoundingDown->GetBinContent(1)  /histo_ggWW  ->GetBinContent(1);
    systMVAJES[2][0] = histo_Higgs_CMS_MVAJESBoundingUp->GetBinContent(1)   /histo_Higgs ->GetBinContent(1);
    systMVAJES[2][1] = histo_Higgs_CMS_MVAJESBoundingDown->GetBinContent(1) /histo_Higgs ->GetBinContent(1);
    systMVAJES[3][0] = histo_VV_CMS_MVAJESBoundingUp->GetBinContent(1)      /histo_VV    ->GetBinContent(1);
    systMVAJES[3][1] = histo_VV_CMS_MVAJESBoundingDown->GetBinContent(1)    /histo_VV    ->GetBinContent(1);
    systMVAJES[4][0] = histo_VVV_CMS_MVAJESBoundingUp->GetBinContent(1)     /histo_VVV   ->GetBinContent(1);
    systMVAJES[4][1] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(1)   /histo_VVV   ->GetBinContent(1);
    systMVAJES[5][0] = histo_Wg3l_CMS_MVAJESBoundingUp->GetBinContent(1)    /histo_Wg3l  ->GetBinContent(1);
    systMVAJES[5][1] = histo_Wg3l_CMS_MVAJESBoundingDown->GetBinContent(1)  /histo_Wg3l  ->GetBinContent(1);
    systMVAJES[6][0] = histo_Wgamma_CMS_MVAJESBoundingUp->GetBinContent(1)  /histo_Wgamma->GetBinContent(1);
    systMVAJES[6][1] = histo_Wgamma_CMS_MVAJESBoundingDown->GetBinContent(1)/histo_Wgamma->GetBinContent(1);
    for(int i=0; i<7; i++) for(int j=0; j<2; j++) if(systMVAJES[i][j] < 1.0) systMVAJES[i][j] = 1.0/systMVAJES[i][j];
    for(int i=0; i<7; i++) systMVAJES[i][2] = (systMVAJES[i][0]+systMVAJES[i][1])/2.0;

    double systNLOBck[6][3];
    systNLOBck[0][0] = histo_VV_CMS_VVNLOBoundingUp       ->GetBinContent(1)/histo_VV    ->GetBinContent(1);
    systNLOBck[0][1] = histo_VV_CMS_VVNLOBoundingDown     ->GetBinContent(1)/histo_VV    ->GetBinContent(1);
    systNLOBck[1][0] = histo_WjetsE_CMS_MVAWEBoundingUp   ->GetBinContent(1)/histo_WjetsE->GetBinContent(1);
    systNLOBck[1][1] = histo_WjetsE_CMS_MVAWEBoundingDown ->GetBinContent(1)/histo_WjetsE->GetBinContent(1);
    systNLOBck[2][0] = histo_WjetsM_CMS_MVAWMBoundingUp   ->GetBinContent(1)/histo_WjetsM->GetBinContent(1);
    systNLOBck[2][1] = histo_WjetsM_CMS_MVAWMBoundingDown ->GetBinContent(1)/histo_WjetsM->GetBinContent(1);
    systNLOBck[3][0] = histo_qqWW_CMS_MVAWWNLOBoundingUp  ->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    systNLOBck[3][1] = histo_qqWW_CMS_MVAWWNLOBoundingDown->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    systNLOBck[4][0] = histo_qqWW_CMS_WWNLOQUp            ->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    systNLOBck[4][1] = histo_qqWW_CMS_WWNLOQDown          ->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    systNLOBck[5][0] = histo_qqWW_CMS_WWNLORUp            ->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    systNLOBck[5][1] = histo_qqWW_CMS_WWNLORDown          ->GetBinContent(1)/histo_qqWW  ->GetBinContent(1);
    for(int i=0; i<6; i++) for(int j=0; j<2; j++) if(systNLOBck[i][j] < 1.0) systNLOBck[i][j] = 1.0/systNLOBck[i][j];
    for(int i=0; i<6; i++) systNLOBck[i][2] = (systNLOBck[i][0]+systNLOBck[i][1])/2.0;

    newcardShape << Form("%s                                     lnN  %5.3f %5.3f %5.3f %5.3f %5.3f   -     -	-     %5.3f %5.3f  -	 -  \n",effName,systMVALepEff[0][2],systMVALepEff[1][2],systMVALepEff[2][2],systMVALepEff[3][2],systMVALepEff[4][2],systMVALepEff[5][2],systMVALepEff[6][2]);
    newcardShape << Form("%s                                     lnN  %5.3f %5.3f %5.3f %5.3f %5.3f   -     -	-     %5.3f %5.3f  -	 -  \n",momName,systMVALepRes[0][2],systMVALepRes[1][2],systMVALepRes[2][2],systMVALepRes[3][2],systMVALepRes[4][2],systMVALepRes[5][2],systMVALepRes[6][2]);
    newcardShape << Form("CMS_scale_met                          lnN  %5.3f %5.3f %5.3f %5.3f %5.3f   -     -	-     %5.3f %5.3f  -	 -  \n",	systMVAMETRes[0][2],systMVAMETRes[1][2],systMVAMETRes[2][2],systMVAMETRes[3][2],systMVAMETRes[4][2],systMVAMETRes[5][2],systMVAMETRes[6][2]);
    newcardShape << Form("CMS_scale_j                            lnN  %5.3f %5.3f %5.3f %5.3f %5.3f   -     -	-     %5.3f %5.3f  -	 -  \n",	systMVAJES[0][2],systMVAJES[1][2],systMVAJES[2][2],systMVAJES[3][2],systMVAJES[4][2],systMVAJES[5][2],systMVAJES[6][2]);
    newcardShape << Form("CMS_wwana_VVNLOBounding                lnN	-     -     -   %5.3f   -     -     -	-       -     -	   -	 -  \n",systNLOBck[0][2]);
    newcardShape << Form("CMS_wwana_MVAWEBounding                lnN	-     -     -	  -     -     -	    -	-       -     -  %5.3f	 -  \n",systNLOBck[1][2]);
    newcardShape << Form("CMS_wwana_MVAWMBounding                lnN	-     -     -	  -     -     -	    -	-       -     -    -   %5.3f\n",systNLOBck[2][2]);
    //newcardShape << Form("CMS_wwana_MVAWWNLOBounding             lnN  %5.3f   -     -	  -     -     -	    -	-       -     -    -	-   \n",systNLOBck[3][2]);
    newcardShape << Form("CMS_wwana_WWNLOQ                       lnN  %5.3f   -	    -     -     -     -     -	-	-     -    -	-   \n",systNLOBck[4][2]);
    newcardShape << Form("CMS_wwana_WWNLOR                       lnN  %5.3f   -	    -     -     -     -     -	-	-     -    -	-   \n",systNLOBck[5][2]);
  }
  //newcardShape << Form("QCDscale_WW			       lnN  %5.3f   -	  -	-     -     -	  -     -     -     -	  -     -  \n",XS_QCDscale_WW[0]);  
  //newcardShape << Form("QCDscale_WW1in  		       lnN  %5.3f   -	  -	-     -     -	  -     -     -     -	  -     -  \n",XS_QCDscale_WW[1]);  
  //newcardShape << Form("QCDscale_WW2in  		       lnN  %5.3f   -	  - 	-     -     -	  -     -     -     -	  -     -  \n",XS_QCDscale_WW[2]);  
  newcardShape << Form("pdf_qqbar                              lnN  %5.3f   -     -   %5.3f %5.3f   -     -     -   %5.3f %5.3f   -	-  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2],pdf_qqbar[3],pdf_qqbar[4]);
  newcardShape << Form("pdf_gg                                 lnN    -   %5.3f %5.3f   -     -     -     -     -     -     -     -	-  \n",pdf_gg[0],pdf_gg[1]);
  newcardShape << Form("QCDscale_ggH                           lnN    -     -   %5.3f   -     -     -	  -     -     -     -     -	-  \n",XS_QCDscale_ggH[0]);  
  newcardShape << Form("QCDscale_ggH1in                        lnN    -     -   %5.3f   -     -     -	  -     -     -     -     -	-  \n",XS_QCDscale_ggH[1]);  
  newcardShape << Form("QCDscale_ggH2in                        lnN    -     -   %5.3f   -     -     -	  -     -     -     -	  -	-  \n",XS_QCDscale_ggH[2]);  
  newcardShape << Form("QCDscale_VV		               lnN  %5.3f   -     -   %5.3f   -     -     -     -     -     -     -	-  \n",QCDscale_VV[0],QCDscale_VV[1]);  
  newcardShape << Form("CMS_UEPS                               lnN  %5.3f %5.3f   -     -     -     -     -     -     -     -     -	-  \n",1.035,1.035);  
  newcardShape << Form("QCDscale_ggVV		               lnN    -   1.300   -     -     -     -     -     -     -     -     -	-  \n");  
  newcardShape << Form("QCDscale_VVV		               lnN    -     -     -     -   1.500   -     -     -     -     -     -	-  \n");	     
  newcardShape << Form("CMS_ww_Wg3l		               lnN    -     -     -     -     -     -     -     -   1.400   -     -	-  \n");	     
  newcardShape << Form("QCDscale_Vgamma		               lnN    -     -     -     -     -     -     -     -     -   1.300   -     -  \n"); 	     
  newcardShape << Form("CMS_FakeRateE                          lnN    -     -     -     -     -	    -	  -	-     -     -	%5.3f   -  \n",WjetsSyst);  
  newcardShape << Form("CMS_FakeRateM                          lnN    -     -     -     -     -     -	  -	-     -     -	  -   %5.3f\n",WjetsSyst);  
  newcardShape << Form("CMS_ww_%1dj_top_%4s                    lnN    -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",nJetsType,ECMsb.Data(),topXS_E);      
  newcardShape << Form("CMS_ww_%1dj_Z_%4s                      lnN    -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",nJetsType,ECMsb.Data(),ZXS_E);		       
  newcardShape << Form("CMS_ww_Ztt                             lnN    -	    -     -     -     -	    -     -   %5.3f   -	    -     -     -  \n",1.10);
  if(useFullStatTemplates == false){
    if(NFinal[0]    > 0) newcardShape << Form("CMS_wwana%s_MVAqqWWStatBounding_%s   shape  1.000   -     -     -     -     -     -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[1]    > 0) newcardShape << Form("CMS_wwana%s_MVAggWWStatBounding_%s   shape    -   1.000   -     -     -     -	 -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[nBkg] > 0) newcardShape << Form("CMS_wwana%s_MVAHiggsStatBounding_%s  shape    -     -   1.000   -     -     -	 -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[2]    > 0) newcardShape << Form("CMS_wwana%s_MVAVVStatBounding_%s	    shape    -     -	 -   1.000   -     -	 -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[3]    > 0) newcardShape << Form("CMS_wwana%s_MVAVVVStatBounding_%s    shape    -     -	 -     -   1.000   -	 -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[4]    > 0) newcardShape << Form("CMS_wwana%s_MVATopStatBounding_%s    shape    -     -	 -     -     -   1.000   -     -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[5]    > 0) newcardShape << Form("CMS_wwana%s_MVAZjetsStatBounding_%s  shape    -     -     -     -     -     -   1.000   -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[6]    > 0) newcardShape << Form("CMS_wwana%s_MVAZttStatBounding_%s    shape    -     -     -     -     -     -	 -   1.000   -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[7]    > 0) newcardShape << Form("CMS_wwana%s_MVAWg3lStatBounding_%s   shape    -     -     -     -     -     -	 -     -   1.000   -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[8]    > 0) newcardShape << Form("CMS_wwana%s_MVAWgammaStatBounding_%s shape    -     -     -     -     -     -	 -     -     -   1.000   -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[9]    > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsEStatBounding_%s shape    -     -     -     -     -     -	 -     -     -     -   1.000   -  \n",finalStateName,ECMsb.Data());
    if(NFinal[10]   > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsMStatBounding_%s shape    -     -     -     -     -     -	 -     -     -     -     -   1.000\n",finalStateName,ECMsb.Data());
  } else {
    if(nBinMVA != 1){    
      for(int nb=1; nb<=nBinMVA; nb++){
	if(NFinal[0]    > 0 && histo_qqWW  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAqqWWStatBounding_%s_Bin%d   shape  1.000   -	-     -     -	  -	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[1]    > 0 && histo_ggWW  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAggWWStatBounding_%s_Bin%d   shape    -	1.000	-     -     -	  -	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[nBkg] > 0 && histo_Higgs ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAHiggsStatBounding_%s_Bin%d  shape    -	  -   1.000   -     -	  -	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[2]    > 0 && histo_VV	   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAVVStatBounding_%s_Bin%d     shape    -	  -	-   1.000   -	  -	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[3]    > 0 && histo_VVV   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAVVVStatBounding_%s_Bin%d    shape    -	  -	-     -   1.000   -	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[4]    > 0 && histo_Top   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVATopStatBounding_%s_Bin%d    shape    -	  -	-     -     -	1.000	-     -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[5]    > 0 && histo_Zjets ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAZjetsStatBounding_%s_Bin%d  shape    -	  -	-     -     -	  -   1.000   -     -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[6]    > 0 && histo_Ztt   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAZttStatBounding_%s_Bin%d    shape    -	  -	-     -     -	  -	-   1.000   -	  -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[7]    > 0 && histo_Wg3l  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAWg3lStatBounding_%s_Bin%d   shape    -	  -	-     -     -	  -	-     -   1.000   -	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[8]    > 0 && histo_Wgamma->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAWgammaStatBounding_%s_Bin%d shape    -	  -	-     -     -	  -	-     -     -	1.000	-     -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[9]    > 0 && histo_WjetsE->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsEStatBounding_%s_Bin%d shape    -	  -	-     -     -	  -	-     -     -	  -   1.000   -  \n",finalStateName,ECMsb.Data(),nb-1);
	if(NFinal[10]   > 0 && histo_WjetsM->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsMStatBounding_%s_Bin%d shape    -	  -	-     -     -	  -	-     -     -	  -	-   1.000\n",finalStateName,ECMsb.Data(),nb-1);
      }
    } else {
      if(NFinal[0]    > 0 && histo_qqWW  ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAqqWWStatBounding_%s   lnN  %5.3f    -    -      -	  -	-     -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_qqWW_CMS_MVAqqWWStatBoundingUp	->GetBinContent(1)/histo_qqWW  ->GetBinContent(1));
      if(NFinal[1]    > 0 && histo_ggWW  ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAggWWStatBounding_%s   lnN    -   %5.3f   -      -	  -	-     -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_ggWW_CMS_MVAggWWStatBoundingUp	->GetBinContent(1)/histo_ggWW  ->GetBinContent(1));
      if(NFinal[nBkg] > 0 && histo_Higgs ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAHiggsStatBounding_%s  lnN    -	-   %5.3f   -	  -	-     -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_Higgs_CMS_MVAHiggsStatBoundingUp	->GetBinContent(1)/histo_Higgs ->GetBinContent(1));
      if(NFinal[2]    > 0 && histo_VV    ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAVVStatBounding_%s	    lnN    -	-     -   %5.3f   -	-     -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_VV_CMS_MVAVVStatBoundingUp 	->GetBinContent(1)/histo_VV    ->GetBinContent(1));
      if(NFinal[3]    > 0 && histo_VVV   ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAVVVStatBounding_%s    lnN    -	-     -     -	%5.3f	-     -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_VVV_CMS_MVAVVVStatBoundingUp	->GetBinContent(1)/histo_VVV   ->GetBinContent(1));
      if(NFinal[4]    > 0 && histo_Top   ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVATopStatBounding_%s    lnN    -	-     -     -	  -   %5.3f   -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_Top_CMS_MVATopStatBoundingUp	->GetBinContent(1)/histo_Top   ->GetBinContent(1));
      if(NFinal[5]    > 0 && histo_Zjets ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAZjetsStatBounding_%s  lnN    -	-     -     -	  -	-   %5.3f   -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_Zjets_CMS_MVAZjetsStatBoundingUp	->GetBinContent(1)/histo_Zjets ->GetBinContent(1));
      if(NFinal[6]    > 0 && histo_Ztt   ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAZttStatBounding_%s    lnN    -	-     -     -	  -	-     -   %5.3f   -	-     -     -  \n",finalStateName,ECMsb.Data(),histo_Ztt_CMS_MVAZttStatBoundingUp	->GetBinContent(1)/histo_Ztt   ->GetBinContent(1));
      if(NFinal[7]    > 0 && histo_Wg3l  ->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAWg3lStatBounding_%s   lnN    -	-     -     -	  -	-     -     -	%5.3f	-     -     -  \n",finalStateName,ECMsb.Data(),histo_Wg3l_CMS_MVAWg3lStatBoundingUp	->GetBinContent(1)/histo_Wg3l  ->GetBinContent(1));
      if(NFinal[8]    > 0 && histo_Wgamma->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAWgammaStatBounding_%s lnN    -	-     -     -	  -	-     -     -	  -   %5.3f   -     -  \n",finalStateName,ECMsb.Data(),histo_Wgamma_CMS_MVAWgammaStatBoundingUp ->GetBinContent(1)/histo_Wgamma->GetBinContent(1));
      if(NFinal[9]    > 0 && histo_WjetsE->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsEStatBounding_%s lnN    -	-     -     -	  -	-     -     -	  -	-   %5.3f   -  \n",finalStateName,ECMsb.Data(),histo_WjetsE_CMS_MVAWjetsEStatBoundingUp ->GetBinContent(1)/histo_WjetsE->GetBinContent(1));
      if(NFinal[10]   > 0 && histo_WjetsM->GetBinContent(1) > 0) newcardShape << Form("CMS_wwana%s_MVAWjetsMStatBounding_%s lnN    -	-     -     -	  -	-     -     -	  -	-     -   %5.3f\n",finalStateName,ECMsb.Data(),histo_WjetsM_CMS_MVAWjetsMStatBoundingUp ->GetBinContent(1)/histo_WjetsM->GetBinContent(1));
    }
  }
  newcardShape.close();
  }

  return;
}

Double_t DYBkgScaleFactor(Int_t jetBin) {
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.05311, 3.84797, 2.18044  };
  return DYBkgScaleFactorWWPreselection[jetBin];
}

Double_t DYBkgScaleFactorKappa(Int_t jetBin) {
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31548, 1.31969, 1.3022  };
  return DYBkgScaleFactorWWPreselectionKappa[jetBin];
}

Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = {  1.11623, 1.0799, 1.15335  };
  return TopBkgScaleFactor[jetBin];
}

Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.12232, 1.0307, 1.02679   };
  return TopBkgScaleFactorKappa[jetBin];
}

double weightJetPt(Int_t nsel, Int_t jetpt0, Int_t jetpt1){
  const int nPoints = 40;

  double weightsLowPtMCNLO0[nPoints] = {0.65087, 0.70000, 0.76000, 0.84076, 0.91597, 0.99955, 1.07837, 1.08566, 1.10771, 1.12779, 1.14879, 1.10786, 1.11974, 1.08315, 1.09204, 1.07666, 1.07266, 1.06118, 1.04589, 1.02708, 1.02882, 1.02535, 1.01788, 1.01469, 1.03972, 1.03880, 1.03756, 1.04755, 1.03103, 1.01557, 1.00861, 1.00821, 1.00634, 1.03981, 1.03034, 1.02481, 0.98525, 1.00569, 0.97956, 0.97421};
  double funcHighPtMCNLO0[4]         = {+1.143727,-0.00359348,+1.29374e-05,-1.42165e-08};
                                      
  double weightsLowPtMCNLO1[nPoints] = {0.68588, 0.74000, 0.80000, 0.87553, 0.96178, 1.05265, 1.11364, 1.14806, 1.18738, 1.16817, 1.14539, 1.11272, 1.12017, 1.08063, 1.05242, 1.04867, 1.01538, 1.03444, 0.98613, 0.99706, 1.03438, 0.98340, 0.98930, 0.98947, 0.97965, 0.98921, 0.97735, 0.97862, 0.96612, 0.97913, 0.96389, 0.94993, 0.94549, 0.97022, 0.94213, 0.94853, 0.93399, 0.94286, 0.91894, 0.94716};
  double funcHighPtMCNLO1[4]         = {+1.162400,-0.00781737,+5.18438e-05,-1.01260e-07};
                                      
  double weightsLowPtMG0[nPoints]   = {1.16602, 1.15500, 0.145000, 1.13685, 1.14275, 1.14062, 1.13216, 1.11681, 1.10233, 1.11667, 1.12721, 1.08748, 1.10651, 1.08424, 1.09132, 1.07137, 1.08236, 1.06970, 1.06747, 1.06177, 1.07718, 1.06670, 1.08122, 1.08012, 1.05013, 1.09825, 1.06345, 1.08442, 1.06397, 1.04977, 1.04595, 1.06918, 1.06977, 1.06965, 1.05999, 1.03867, 1.03521, 1.04005, 1.02379, 1.04293};
  double funcHighPtMG0[4]           = {+1.37488,-0.00959298,+3.87912e-05,-5.27704e-08};
  
  double weightsLowPtMG1[nPoints]   = {1.11503, 1.11700, 1.120000, 1.12289, 1.11773, 1.11894, 1.11425, 1.09906, 1.10818, 1.09070, 1.08716, 1.06869, 1.08183, 1.04786, 1.02879, 1.02570, 1.01640, 0.99879, 0.98461, 0.96576, 0.97864, 0.96116, 0.94301, 0.95326, 0.93664, 0.91158, 0.90140, 0.89942, 0.89236, 0.86326, 0.87440, 0.86260, 0.83675, 0.85810, 0.85107, 0.83513, 0.82202, 0.81394, 0.81715, 0.77647};
  double funcHighPtMG1[4]           = {+1.06397,-0.00883160,+3.58076e-05,-5.13992e-08};
  
  double weightsLowPt0[nPoints];
  double funcHighPt0[4];

  double weightsLowPt1[nPoints];
  double funcHighPt1[4];
  
  if(nsel == 0){ // MC@NLO
    for(int i=0; i<nPoints; i++) weightsLowPt0[i] = weightsLowPtMCNLO0[i];
    for(int i=0; i<4; i++) funcHighPt0[i] = funcHighPtMCNLO0[i];
    for(int i=0; i<nPoints; i++) weightsLowPt1[i] = weightsLowPtMCNLO1[i];
    for(int i=0; i<4; i++) funcHighPt1[i] = funcHighPtMCNLO1[i];
  } else { // MG
    for(int i=0; i<nPoints; i++) weightsLowPt0[i] = weightsLowPtMG0[i];
    for(int i=0; i<4; i++) funcHighPt0[i] = funcHighPtMG0[i];
    for(int i=0; i<nPoints; i++) weightsLowPt1[i] = weightsLowPtMG1[i];
    for(int i=0; i<4; i++) funcHighPt1[i] = funcHighPtMG1[i];
  }

  double weight0 = 1.0;
  double weight1 = 1.0;
  
  if(jetpt0 >= 40) weight0 = (funcHighPt0[0] + funcHighPt0[1]*jetpt0 + funcHighPt0[2]*jetpt0*jetpt0 + funcHighPt0[3]*jetpt0*jetpt0*jetpt0);
  else             weight0 = weightsLowPt0[jetpt0];
  
  if(jetpt1 >= 40) weight1 = (funcHighPt1[0] + funcHighPt1[1]*jetpt1 + funcHighPt1[2]*jetpt1*jetpt1 + funcHighPt1[3]*jetpt1*jetpt1*jetpt1);
  else             weight1 = weightsLowPt1[jetpt1];

  return (weight0*weight1);
}
