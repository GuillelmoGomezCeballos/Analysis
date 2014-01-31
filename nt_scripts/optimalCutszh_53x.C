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

void metChange(double met, double metPhi, double nvtx, int year, bool isData, double metNew[2]);
double collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2, const int q1);

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 10;
const unsigned int nSelTypesSyst = 7;
const bool showSignalOnly = false;
const bool useGJForShapes = true;
const bool isZZLepSel     = true;
const bool properWjetsBkg = false; // taking W+jets background from EM background estimation otherwise

enum selType {ZSEL=0, WWSEL, BTAGSEL, WZSEL, SIGSEL, BCK1REG, BCK2REG, BCK12REG, WWLOOSESEL, ZLLSEL};
TString selTypeName[nSelTypes*2] = {"ZSEL-EM", "WWSEL-EM", "BTAGSEL-EM", "WZSEL-EM", "SIGSEL-EM", "BCK1REG-EM", "BCK2REG-EM", "BCK12REG-EM", "WWLOOSESEL-EM", "ZLLSEL-MM",
                                    "ZSEL-LL", "WWSEL-LL", "BTAGSEL-LL", "WZSEL-LL", "SIGSEL-LL", "BCK1REG-LL", "BCK2REG-LL", "BCK12REG-LL", "WWLOOSESEL-LL", "ZLLSEL-EE"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-EM", "JESDOWN-EM", "LEPP-EM", "LEPM-EM", "MET-EM", "EFFP-EM", "EFFM-EM",
                                            "JESUP-LL", "JESDOWN-LL", "LEPP-LL", "LEPM-LL", "MET-LL", "EFFP-LL", "EFFM-LL"};

// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCutszh_53x
(
 int     mH  	 = 1125,
 int thePlot = 7,
 TString bgdInputFile    = "ntuples_53x/backgroundA_skim10.root",
 TString signalInputFile = "ntuples_53x/zh125inv.root",
 TString dataInputFile   = "ntuples_53x/data_skim10.root",
 TString systInputFile   = "ntuples_53x/hww_syst_skim10.root",
 int period = 3,
 int lSel = 2, 
 int nJetsType = 0
 ,double var0 = 120, double var1 = 2.7, double var2 = 0.25, double var3 = 180
 )
{
  double lumi = 1.0;
  
  //                    MET,   DPhiZMET, PTFrac, dPhi
  double cutValue[4] = {120.0, 2.7,     0.25,   180};
  cutValue[0] = var0; cutValue[1] = var1; cutValue[2] = var2; cutValue[3] = var3;
  double ptJetMin = -9.0;
  double ptJet1st =  30.;
  double ptJet2nd =  30.;
  if     (nJetsType == 1) {ptJetMin = 30; ptJet1st = 999999999.;}
  else if(nJetsType >  1) assert(0);
  double metMin   = 60.;
  double useFullStatTemplates = true;
  bool useWeightEWKCorr       = true;

  bool fCheckProblem = true;
  if(isZZLepSel == true) fCheckProblem = false;

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

  // useful if using ZZ lepton selection
  LeptonEfficiencyZH theLeptonEfficiencyZH(year);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int channel = mH;
  if(channel > 1000 && channel < 2000) mH = channel-1000;
  else assert(0);

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 400.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 18;  xminPlot =   0.0; xmaxPlot = 900.0;}
  else if(thePlot >=  0 && thePlot <=  0) {nBinPlot = 40;  xminPlot =   0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 40;  xminPlot =   0.0; xmaxPlot = 400.0;}
  else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  60; xminPlot =  60.0; xmaxPlot = 120.0;}
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 20 && thePlot <= 23) {nBinPlot = 18; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 24 && thePlot <= 24) {nBinPlot = 50; xminPlot =-1.0; xmaxPlot = 1.0;}
  else if(thePlot >= 25 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 30 && thePlot <= 30) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 300; xminPlot = 0.0; xmaxPlot = 600.0;}
  else if(thePlot >= 33 && thePlot <= 33) {nBinPlot = 90; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  800.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  2.5;}
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  5.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  else if(thePlot >= 53 && thePlot <= 55) {nBinPlot = 36; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 56 && thePlot <= 56) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 57 && thePlot <= 57) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}

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

  const int nBinMVA = 11;
  Float_t xbins[nBinMVA+1] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_ZH_hinv= (TH1D*) histoMVA->Clone("histo_ZH_hinv");
  TH1D *histo_Zjets  = (TH1D*) histoMVA->Clone("histo_Zjets");
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EM     = (TH1D*) histoMVA->Clone("histo_EM");
  TH1D *histo_Wjets  = (TH1D*) histoMVA->Clone("histo_Wjets");

  TH1D *histo_GAMMA_JETS       = (TH1D*) histoMVA->Clone("histo_GAMMA_JETS");
  TH1D *histo_GAMMA_JETS_LOOSE = (TH1D*) histoMVA->Clone("histo_GAMMA_JETS_LOOSE");

  char finalStateName[2],effName[10],momName[10];
  if     (lSel == 2 && nJetsType == 0) {sprintf(finalStateName,"ll0j");sprintf(effName,"CMS_eff_l");sprintf(momName,"CMS_scale_l");}
  else if(lSel == 0 && nJetsType == 0) {sprintf(finalStateName,"mm0j");sprintf(effName,"CMS_eff_m");sprintf(momName,"CMS_scale_m");}
  else if(lSel == 1 && nJetsType == 0) {sprintf(finalStateName,"ee0j");sprintf(effName,"CMS_eff_e");sprintf(momName,"CMS_scale_e");}
  else if(lSel == 2 && nJetsType == 1) {sprintf(finalStateName,"ll1j");sprintf(effName,"CMS_eff_l");sprintf(momName,"CMS_scale_l");}
  else if(lSel == 0 && nJetsType == 1) {sprintf(finalStateName,"mm1j");sprintf(effName,"CMS_eff_m");sprintf(momName,"CMS_scale_m");}
  else if(lSel == 1 && nJetsType == 1) {sprintf(finalStateName,"ee1j");sprintf(effName,"CMS_eff_e");sprintf(momName,"CMS_scale_e");}
  else {printf("Wrong lSel: %d\n",lSel); assert(0);}

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingUp      = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingDown    = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingUp           = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingUp  ->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingDown         = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingDown->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingUp     = new TH1D( Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingDown   = new TH1D( Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinUp[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinDown[nBinMVA];
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]    = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]	  ->Sumw2();
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]  = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	  ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	  ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	  ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	  ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	    = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	  ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	    = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	  ->Sumw2();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb] = new TH1D(Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_Wjets_CMS_zllhinv%s_MVAWjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb]->Sumw2();
  }

  TH1D* histo_ZH_hinv_CMS_MVALepEffBoundingUp   = new TH1D( Form("histo_ZH_hinv_%sUp",effName)  , Form("histo_ZH_hinv_%sUp",effName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffBoundingDown = new TH1D( Form("histo_ZH_hinv_%sDown",effName), Form("histo_ZH_hinv_%sDown",effName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingUp   	= new TH1D( Form("histo_VVV_%sUp",effName)  , Form("histo_VVV_%sUp",effName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingDown 	= new TH1D( Form("histo_VVV_%sDown",effName), Form("histo_VVV_%sDown",effName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingUp    	= new TH1D( Form("histo_WZ_%sUp",effName)  , Form("histo_WZ_%sUp",effName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingDown  	= new TH1D( Form("histo_WZ_%sDown",effName), Form("histo_WZ_%sDown",effName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffBoundingUp    	= new TH1D( Form("histo_ZZ_%sUp",effName)  , Form("histo_ZZ_%sUp",effName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffBoundingDown  	= new TH1D( Form("histo_ZZ_%sDown",effName), Form("histo_ZZ_%sDown",effName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVALepResBoundingUp   = new TH1D( Form("histo_ZH_hinv_%sUp",momName)  , Form("histo_ZH_hinv_%sUp",momName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepResBoundingDown = new TH1D( Form("histo_ZH_hinv_%sDown",momName), Form("histo_ZH_hinv_%sDown",momName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingUp   	= new TH1D( Form("histo_VVV_%sUp",momName)  , Form("histo_VVV_%sUp",momName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingDown 	= new TH1D( Form("histo_VVV_%sDown",momName), Form("histo_VVV_%sDown",momName), nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingUp    	= new TH1D( Form("histo_WZ_%sUp",momName)  , Form("histo_WZ_%sUp",momName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingDown  	= new TH1D( Form("histo_WZ_%sDown",momName), Form("histo_WZ_%sDown",momName), nBinMVA, xbins); histo_WZ_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepResBoundingUp    	= new TH1D( Form("histo_ZZ_%sUp",momName)  , Form("histo_ZZ_%sUp",momName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepResBoundingDown  	= new TH1D( Form("histo_ZZ_%sDown",momName), Form("histo_ZZ_%sDown",momName), nBinMVA, xbins); histo_ZZ_CMS_MVALepResBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAMETResBoundingUp   = new TH1D( Form("histo_ZH_hinv_CMS_scale_metUp")  , Form("histo_ZH_hinv_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAMETResBoundingDown = new TH1D( Form("histo_ZH_hinv_CMS_scale_metDown"), Form("histo_ZH_hinv_CMS_scale_metDown"), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETResBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETResBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETResBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingUp      = new TH1D( Form("histo_ZH_hinv_CMS_scale_jUp")  , Form("histo_ZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingDown    = new TH1D( Form("histo_ZH_hinv_CMS_scale_jDown"), Form("histo_ZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp      	= new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown    	= new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_WZ_CMS_scale_jUp")  , Form("histo_WZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_WZ_CMS_scale_jDown"), Form("histo_WZ_CMS_scale_jDown"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_WZ_CMS_WZNLOBoundingUp        = new TH1D( Form("histo_WZ_CMS_zllhinv_WZNLOBoundingUp"),   Form("histo_WZ_CMS_zllhinv_WZNLOBoundingUp"),   nBinMVA, xbins); histo_WZ_CMS_WZNLOBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_WZNLOBoundingDown      = new TH1D( Form("histo_WZ_CMS_zllhinv_WZNLOBoundingDown"), Form("histo_WZ_CMS_zllhinv_WZNLOBoundingDown"), nBinMVA, xbins); histo_WZ_CMS_WZNLOBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_ZZNLOBoundingUp        = new TH1D( Form("histo_ZZ_CMS_zllhinv_ZZNLOBoundingUp"),   Form("histo_ZZ_CMS_zllhinv_ZZNLOBoundingUp"),   nBinMVA, xbins); histo_ZZ_CMS_ZZNLOBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_ZZNLOBoundingDown      = new TH1D( Form("histo_ZZ_CMS_zllhinv_ZZNLOBoundingDown"), Form("histo_ZZ_CMS_zllhinv_ZZNLOBoundingDown"), nBinMVA, xbins); histo_ZZ_CMS_ZZNLOBoundingDown->Sumw2();

  TH1D* histo_Wjets_CMS_MVAWBoundingUp      = new TH1D( Form("histo_Wjets_CMS_zllhinv_MVAWBoundingUp"),   Form("histo_Wjets_CMS_zllhinv_MVAWBoundingUp"),   nBinMVA, xbins); histo_Wjets_CMS_MVAWBoundingUp  ->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWBoundingDown    = new TH1D( Form("histo_Wjets_CMS_zllhinv_MVAWBoundingDown"), Form("histo_Wjets_CMS_zllhinv_MVAWBoundingDown"), nBinMVA, xbins); histo_Wjets_CMS_MVAWBoundingDown->Sumw2();

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

  //unsigned int patternTopVeto = SmurfTree::TopVeto;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

    if(!(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  bgdEvent.dstype_ != SmurfTree::data)) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
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

    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);

    if(bgdEvent.lep1MotherMcId_ == 23 && bgdEvent.lep2MotherMcId_ == 23) {
      if     (fDecay == 21) fDecay = 11;
      else if(fDecay == 27) fDecay = 17;
      else if(fDecay == 28) fDecay = 18;
    }

    bool passSystCuts[3][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[3][nSelTypes] = {{false, false, false, false, false, false, false, false, false, false},
                                   {false, false, false, false, false, false, false, false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;

    double metNew[2]; metChange(bgdEvent.met_,bgdEvent.metPhi_,bgdEvent.nvtx_,year,bgdEvent.dstype_ == SmurfTree::data,metNew);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 
    //bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
    //bool passBtagVeto      = bgdEvent.nSoftMuons_ == 0 && (bgdEvent.jet1_.Pt() <= 20 || bgdEvent.jet1Btag_     < 2.100) && (bgdEvent.jet2_.Pt() <= 20 || bgdEvent.jet2Btag_     < 2.100) && (bgdEvent.jet3_.Pt() <= 20 || bgdEvent.jet3Btag_     < 2.100) && (bgdEvent.jet4_.Pt() <= 20 || bgdEvent.jet4Btag_     < 2.100);
    bool passBtagVeto      = bgdEvent.nSoftMuons_ == 0 && (bgdEvent.jet1_.Pt() <= 20 || bgdEvent.jet1ProbBtag_ < 0.244) && (bgdEvent.jet2_.Pt() <= 20 || bgdEvent.jet2ProbBtag_ < 0.244) && (bgdEvent.jet3_.Pt() <= 20 || bgdEvent.jet3ProbBtag_ < 0.244) && (bgdEvent.jet4_.Pt() <= 20 || bgdEvent.jet4ProbBtag_ < 0.244);
    bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passZMass         = fabs(bgdEvent.dilep_.M()-91.1876) < 15.;
    bool passZMassLarge    = fabs(bgdEvent.dilep_.M()-91.1876) < 30.;
    bool passZMassSB       = (bgdEvent.dilep_.M() > 110.0 && bgdEvent.dilep_.M() < 200.0);
    bool passMET           = theMET > cutValue[0];
    bool passDPhiZMET      = DeltaPhi(bgdEvent.dilep_.Phi() ,theMETPHI) > cutValue[1];
    bool passPTFrac        = fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < cutValue[2];
    bool passDPhiLL        = bgdEvent.dPhi_*180.0/TMath::Pi() < cutValue[3];
    bool passPTLL          = bgdEvent.dilep_.Pt() > 60.;
    int lType = -1;
    if     (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) lType = 0;
    else if(lSel == 0 && bgdEvent.type_ == SmurfTree::mm) lType = 1;
    else if(lSel == 1 && bgdEvent.type_ == SmurfTree::ee) lType = 1;
    else if(lSel == 2 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) lType = 1;
    else lType = 2;
    
    bool passBCK1REG       = !passDPhiZMET &&  passPTFrac &&
 			      fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < 1.00;
    bool passBCK2REG	   =  passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < 1.00;
    bool passBCK12REG	   = !passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt() < 1.00;

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
    double MVAVar[6] = {outputVar[8],outputVar[8],outputVar[8],outputVarLepP[8],outputVarLepM[8],outputVarMET[8]};
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    if(isZZLepSel == false) {
      addLepEff	 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
                   leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
      if(addLepEff > 0) {
        addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
    			leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
        addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
    			leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
      } else {addLepEff = 1.0;}
    } else {
      std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lid1_);
      std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.lid2_);
      addLepEff	    = effZH1.first * effZH2.first;
      addLepEffUp   = (effZH1.first+effZH1.second)*(effZH2.first+effZH2.second);
      addLepEffDown = (effZH1.first-effZH1.second)*(effZH2.first-effZH2.second);
      if(addLepEffDown < 0) addLepEffDown = 0.;
    }
    double NjetSyst[2] = {0., 0.};
    if(bgdEvent.jet1_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;

    if(channel > 1000 && channel < 2000){ // ZH->2l+inv selection
      if(charge == 0 && NjetSyst[0] == nJetsType                                                                           && outputVar[4]     > metMin && outputVar[0]	  > 20.0 &&  outputVar[1]	  > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVar[2]-91.1876)	 < 15. && passDPhiLL && outputVar[4]	> cutValue[0] && outputVar[9]     > cutValue[1] && outputVar[11]     < cutValue[2]) passSystCuts[lType][JESUP] = true;
      if(charge == 0 && NjetSyst[1] == nJetsType                                                                           && outputVar[4]     > metMin && outputVar[0]	  > 20.0 &&  outputVar[1]	   > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVar[2]-91.1876)    < 15. && passDPhiLL && outputVar[4]	> cutValue[0] && outputVar[9]     > cutValue[1] && outputVar[11]	< cutValue[2]) passSystCuts[lType][JESDOWN] = true;
      if(charge == 0 && bgdEvent.jet1_.Pt() > ptJetMin && bgdEvent.jet1_.Pt() < ptJet1st && bgdEvent.jet2_.Pt() < ptJet2nd && outputVarLepP[4] > metMin && outputVarLepP[0] > 20.0 &&  outputVarLepP[1] > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarLepP[2]-91.1876) < 15. && passDPhiLL && outputVarLepP[4] > cutValue[0] && outputVarLepP[9] > cutValue[1] && outputVarLepP[11] < cutValue[2]) passSystCuts[lType][LEPP] = true;
      if(charge == 0 && bgdEvent.jet1_.Pt() > ptJetMin && bgdEvent.jet1_.Pt() < ptJet1st && bgdEvent.jet2_.Pt() < ptJet2nd && outputVarLepM[4] > metMin && outputVarLepM[0] > 20.0 &&  outputVarLepM[1] > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarLepM[2]-91.1876) < 15. && passDPhiLL && outputVarLepM[4] > cutValue[0] && outputVarLepM[9] > cutValue[1] && outputVarLepM[11] < cutValue[2]) passSystCuts[lType][LEPM] = true;
      if(charge == 0 && bgdEvent.jet1_.Pt() > ptJetMin && bgdEvent.jet1_.Pt() < ptJet1st && bgdEvent.jet2_.Pt() < ptJet2nd && outputVarMET[4]  > metMin && outputVarMET[0]  > 20.0 &&  outputVarMET[1]  > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarMET[2]-91.1876)  < 15. && passDPhiLL && outputVarMET[4]  > cutValue[0] && outputVarMET[9]  > cutValue[1] && outputVarMET[11]  < cutValue[2]) passSystCuts[lType][MET] = true;
      if(
         charge == 0 &&
	 bgdEvent.jet1_.Pt() > ptJetMin && bgdEvent.jet1_.Pt() < ptJet1st && bgdEvent.jet2_.Pt() < ptJet2nd &&
	 theMET > metMin &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 20.) {
	 
         passCuts[lType][ZSEL]   = true;	 
	 if(!passBtagVeto &&  pass3rLVeto && !passZMass && passZMassSB                                                                      ) passCuts[lType][WWLOOSESEL] = true;
	 if( passBtagVeto &&  pass3rLVeto && !passZMass && passZMassLarge && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][WWSEL]	  = true;
	 if(!passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][BTAGSEL]	  = true;
	 if( passBtagVeto && !pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac && bgdEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL]    = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGSEL]   = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK1REG			  ) passCuts[lType][BCK1REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK2REG			  ) passCuts[lType][BCK2REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK12REG			  ) passCuts[lType][BCK12REG] = true;
         if     (bgdEvent.type_ == SmurfTree::mm && pass3rLVeto && passPTLL && passZMassLarge) passCuts[0][ZLLSEL] = true;
         else if(bgdEvent.type_ == SmurfTree::ee && pass3rLVeto && passPTLL && passZMassLarge) passCuts[1][ZLLSEL] = true;

	if(isRealLepton == false &&
	   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
	    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
	    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) 
	  {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false; passCuts[0][ZLLSEL] = false; passCuts[1][ZLLSEL] = false;}
      }
    } // ZH->2l+inv selection

    if(1){
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
          if(isZZLepSel == false) {
            add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	    add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
            if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
            add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);
          } else {
            std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lid1_);
            std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.lid2_);
            std::pair<float,float> effZH3 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), bgdEvent.lid3_);
            add = add * effZH1.first * effZH2.first;
	    if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) add = add * effZH3.first;
	  }

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
	if(isZZLepSel == false) {
	  sf_eff = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
                   leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
	} else {
          std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lid1_);
          std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.lid2_);
          sf_eff = effZH1.first * effZH2.first;
	}
        theWeight = ZttScaleFactor(period,bgdEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);
        double add2 = 1.0;
        if(isZZLepSel == false) {
	  add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add2 = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add2 = add2*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);
        } else {
          std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lid1_);
          std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.lid2_);
          std::pair<float,float> effZH3 = theLeptonEfficiencyZH.getLeptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), bgdEvent.lid3_);
          add2 = effZH1.first * effZH2.first;
	  if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) add2 = add2 * effZH3.first;
	}
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

	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,theMET);
        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add;
      }

      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::wz) theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,1);
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::zz) theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,2);

      if(passCuts[1][SIGSEL] && theMET > metMin){ // begin making plots
	double myVar = theMET; // var0
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
	else if(thePlot ==12) myVar = bgdEvent.trackMet_;
	else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt(); // var2
	else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = bgdEvent.njets_;
	else if(thePlot ==18) myVar = bgdEvent.nvtx_;
	else if(thePlot ==19) myVar = bgdEvent.pTrackMet_;
	else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi(); // var1
	else if(thePlot ==23) myVar = DeltaPhi(bgdEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = collins_Soper(bgdEvent.lep1_,bgdEvent.lep2_,bgdEvent.lq1_);
	else if(thePlot ==25) myVar = fabs(bgdEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = bgdEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,bgdEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Min(TMath::Abs(bgdEvent.dilep_.Eta()),2.4999);
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
	else if(thePlot ==57) myVar = bgdEvent.dR_;
      	if     (fDecay == 9 || fDecay == 19){
      	  if(useGJForShapes == false) histo0->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 11){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 17){
      	  histo2->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 18){
      	  histo3->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 21 || fDecay == 27 || fDecay == 28 ||
                fDecay ==  1 || fDecay ==  3 || fDecay == 23 || 
                fDecay == 29 || fDecay == 30 || fDecay ==  5 ||
	        fDecay == 13 || fDecay == 20 || fDecay == 10){
      	  histo4->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 41 || fDecay == 42 || fDecay == 43){
      	}
      	else {
      	  printf("NOOOOOOOOOOOOOOOOOOOO\n");
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
      if(passCuts[lType][SIGSEL]) {
        bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
        weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
        weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;

        if     (fDecay == 9 || fDecay == 19){
	  if(passCuts[1][SIGSEL])              histo_Zjets                        ->Fill(MVAVar[0], theWeight);
        }
	else if(fDecay == 11){
	  if(passCuts[1][SIGSEL])	       histo_VVV			  ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][SIGSEL])              histo_VVV_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])              histo_VVV_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_VVV_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_VVV_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_VVV_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_VVV_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
          if(passSystCuts[1][MET]     == true) histo_VVV_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
        }
	else if(fDecay == 17){
	  if(passCuts[1][SIGSEL])	       histo_WZ 			 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][SIGSEL])	       histo_WZ_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])	       histo_WZ_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_WZ_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_WZ_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_WZ_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_WZ_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
          if(passSystCuts[1][MET]     == true) histo_WZ_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
        }
	else if(fDecay == 18){
	  if(passCuts[1][SIGSEL])	       histo_ZZ 			 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][SIGSEL])	       histo_ZZ_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])	       histo_ZZ_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_ZZ_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_ZZ_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_ZZ_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_ZZ_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
          if(passSystCuts[1][MET]     == true) histo_ZZ_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
        }
	else if(fDecay == 21 || fDecay == 27 || fDecay == 28 ||
                fDecay == 29 || fDecay == 30 || fDecay ==  5 ||
	        fDecay == 13 || fDecay == 20 || fDecay == 10 ||
		fDecay == 1  || fDecay == 23){
	  if(passCuts[1][SIGSEL])              histo_EM                          ->Fill(MVAVar[0], theWeight);
        }
	/*
	else if(fDecay == 1 || fDecay == 23){
          double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                       (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
                 addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                       (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                                       (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
                 addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                                       (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  if(passCuts[1][SIGSEL])	       histo_Wjets			 ->Fill(MVAVar[0], theWeight);
          if(passCuts[1][SIGSEL])              histo_Wjets_CMS_MVAWBoundingUp    ->Fill(MVAVar[0], theWeight*addFRS/addFR);
        }
	*/
      }
    } // if passCuts
  } // end background loop

  if(systInputFile != ""){
  int nSyst=systEvent.tree_->GetEntries();
  for (int evt=0; evt<nSyst; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nSyst);
    systEvent.tree_->GetEntry(evt);

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
    else if(systEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(systEvent.dstype_ == SmurfTree::other          ) fDecay = 40;
    else if(systEvent.processId_==121 ||
            systEvent.processId_==122)   fDecay = 41;
    else if(systEvent.processId_==24)    fDecay = 42;
    else if(systEvent.processId_==26)    fDecay = 43;
    else if(systEvent.processId_==10001) fDecay = 44;
    else if(systEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;}

    int charge = (int)(systEvent.lq1_ + systEvent.lq2_);

    if(systEvent.lep1MotherMcId_ == 23 && systEvent.lep2MotherMcId_ == 23) {
      if     (fDecay == 21) fDecay = 11;
      else if(fDecay == 27) fDecay = 17;
      else if(fDecay == 28) fDecay = 18;
    }

    bool passCuts[3][nSelTypes] = {{false, false, false, false, false, false, false, false, false, false},
                                   {false, false, false, false, false, false, false, false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
       (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) isRealLepton = true;

    double metNew[2]; metChange(systEvent.met_,systEvent.metPhi_,systEvent.nvtx_,year,systEvent.dstype_ == SmurfTree::data,metNew);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 
    //bool passBtagVeto      = (systEvent.cuts_ & patternTopVeto) == patternTopVeto;
    //bool passBtagVeto      = systEvent.nSoftMuons_ == 0 && (systEvent.jet1_.Pt() <= 20 || systEvent.jet1Btag_     < 2.100) && (systEvent.jet2_.Pt() <= 20 || systEvent.jet2Btag_     < 2.100) && (systEvent.jet3_.Pt() <= 20 || systEvent.jet3Btag_     < 2.100) && (systEvent.jet4_.Pt() <= 20 || systEvent.jet4Btag_     < 2.100);
    bool passBtagVeto      = systEvent.nSoftMuons_ == 0 && (systEvent.jet1_.Pt() <= 20 || systEvent.jet1ProbBtag_ < 0.244) && (systEvent.jet2_.Pt() <= 20 || systEvent.jet2ProbBtag_ < 0.244) && (systEvent.jet3_.Pt() <= 20 || systEvent.jet3ProbBtag_ < 0.244) && (systEvent.jet4_.Pt() <= 20 || systEvent.jet4ProbBtag_ < 0.244);
    bool pass3rLVeto       = (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passZMass         = fabs(systEvent.dilep_.M()-91.1876) < 15.;
    bool passMET           = theMET > cutValue[0];
    bool passDPhiZMET      = DeltaPhi(systEvent.dilep_.Phi() ,theMETPHI) > cutValue[1];
    bool passPTFrac        = fabs(theMET-systEvent.dilep_.Pt())/systEvent.dilep_.Pt() < cutValue[2];
    bool passDPhiLL        = systEvent.dPhi_*180.0/TMath::Pi() < cutValue[3];
    bool passPTLL          = systEvent.dilep_.Pt() > 60.;
    int lType = -1;
    if     (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me) lType = 0;
    else if(lSel == 0 && systEvent.type_ == SmurfTree::mm) lType = 1;
    else if(lSel == 1 && systEvent.type_ == SmurfTree::ee) lType = 1;
    else if(lSel == 2 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) lType = 1;
    else lType = 2;

    if(channel > 1000 && channel < 2000){ // ZH->2l+inv selection
      if(
         fDecay != 1 && // special treatment for gamma+jets
	 charge == 0 &&
	 systEvent.jet1_.Pt() > ptJetMin && systEvent.jet1_.Pt() < ptJet1st && systEvent.jet2_.Pt() < ptJet2nd &&
	 theMET > metMin &&
	 passPTLL &&
         systEvent.lep1_.Pt() > 20. &&
         systEvent.lep2_.Pt() > 20.) {
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass && passMET && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGSEL]  = true;

	if(isRealLepton == false &&
	   (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
	    systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
	    systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) 
	  {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false; passCuts[0][ZLLSEL] = false; passCuts[1][ZLLSEL] = false;}

      }
      else if(
         fDecay == 1 && // special treatment for gamma+jets
	 systEvent.jet1_.Pt() > ptJetMin && systEvent.jet1_.Pt() < ptJet1st && systEvent.jet2_.Pt() < ptJet2nd &&
	 theMET > metMin &&
	(systEvent.lid1_%16) != 15 && TMath::Abs(systEvent.dilep_.Eta()) < 2.3 &&
	 passPTLL) {
	 if( passMET && passDPhiZMET && passPTFrac) passCuts[0][SIGSEL]  = true;
	 passCuts[0][ZLLSEL] = true;
      }
      
    } // ZH->2l+inv selection

    if     (fDecay == 1 && passCuts[0][ZLLSEL]){
      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[8],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      double theWeight = systEvent.scale1fb_;
      histo_GAMMA_JETS_LOOSE->Fill(MVAVar[0], theWeight);
    }

    if     (fDecay == 1 && passCuts[0][SIGSEL]){
      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[8],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      double theWeight = systEvent.scale1fb_;
      histo_GAMMA_JETS->Fill(MVAVar[0], theWeight);
      if(1) {
	double myVar = theMET; // var0
	if     (thePlot == 1) myVar = systEvent.lep1_.Pt();
	else if(thePlot == 2) myVar = systEvent.lep2_.Pt();
	else if(thePlot == 3) myVar = systEvent.lep3_.Pt();
	else if(thePlot == 4) myVar = systEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = systEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = systEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = systEvent.dilep_.M();
	else if(thePlot == 8) myVar = systEvent.mt_;
	else if(thePlot == 9) myVar = systEvent.mt1_;
	else if(thePlot ==10) myVar = systEvent.mt2_;
	else if(thePlot ==12) myVar = systEvent.trackMet_;
	else if(thePlot ==13) myVar = systEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(systEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-systEvent.dilep_.Pt())/systEvent.dilep_.Pt(); // var2
	else if(thePlot ==16) myVar = systEvent.lep2_.Pt()/systEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = systEvent.njets_;
	else if(thePlot ==18) myVar = systEvent.nvtx_;
	else if(thePlot ==19) myVar = systEvent.pTrackMet_;
	else if(thePlot ==20) myVar = systEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(systEvent.dPhiLep1MET_,systEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(systEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi(); // var1
	else if(thePlot ==23) myVar = DeltaPhi(systEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = collins_Soper(systEvent.lep1_,systEvent.lep2_,systEvent.lq1_);
	else if(thePlot ==25) myVar = fabs(systEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = systEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(systEvent.jet1_.Eta()),fabs(systEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(systEvent.jet1_.Eta()),fabs(systEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(systEvent.jet1_.Eta()),fabs(systEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(systEvent.dilep_.Phi() ,systEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Min(TMath::Abs(systEvent.dilep_.Eta()),2.4999);
	else if(thePlot ==38) myVar = TMath::Abs(systEvent.jet1_.Eta()-systEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(systEvent.jet1_.Phi() ,systEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(systEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = systEvent.jet1_.Pt()+ systEvent.jet2_.Pt()+systEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = systEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = systEvent.trackMet_*cos(systEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = systEvent.trackMet_*sin(systEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(systEvent.jet3_.Phi(),systEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = systEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = systEvent.dR_;
      	if(useGJForShapes == true) histo0->Fill(myVar,theWeight*0.01);
      }
    }
    else if(passCuts[lType][SIGSEL]){
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
          if(isZZLepSel == false) {
            add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	    add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
            if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
            add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);
          } else {
            std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), systEvent.lid1_);
            std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), systEvent.lid2_);
            std::pair<float,float> effZH3 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), systEvent.lid3_);
            add = add * effZH1.first * effZH2.first;
	    if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) add = add * effZH3.first;
	  }
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
	if(isZZLepSel == false) {
	  sf_eff = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_)*
                   leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
	} else {
          std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), systEvent.lid1_);
          std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), systEvent.lid2_);
          sf_eff = effZH1.first * effZH2.first;
	}
        theWeight = ZttScaleFactor(period,systEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(systEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,systEvent.npu_);
        double add2 = 1.0;
        if(isZZLepSel == false) {
	  add2 = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	  add2 = add2*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
          if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add2 = add2*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);
        } else {
          std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), systEvent.lid1_);
          std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), systEvent.lid2_);
          std::pair<float,float> effZH3 = theLeptonEfficiencyZH.getLeptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), systEvent.lid3_);
          add2 = effZH1.first * effZH2.first;
	  if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) add2 = add2 * effZH3.first;
	}
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

      if(useWeightEWKCorr == true && systEvent.dstype_ == SmurfTree::wz) theWeight = theWeight * weightEWKCorr(systEvent.higgsPt_,1);
      if(useWeightEWKCorr == true && systEvent.dstype_ == SmurfTree::zz) theWeight = theWeight * weightEWKCorr(systEvent.higgsPt_,2);

      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[8],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][SIGSEL]){
	if     (fDecay == 17){
	  histo_WZ_CMS_WZNLOBoundingUp->Fill(MVAVar[0], theWeight);
        }
	else if(fDecay == 18){
	  histo_ZZ_CMS_ZZNLOBoundingUp->Fill(MVAVar[0], theWeight);
        }
      }

    } // if passCuts
  } // end syst loop
  } // if want to use it at all

  int nSig=sigEvent.tree_->GetEntries();
  for (int evt=0; evt<nSig; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
       printf("--- reading Signal event %5d of %5d\n",evt,nSig);
    sigEvent.tree_->GetEntry(evt);

    bool lId = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
               (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

//    bool passEvent =
//sigEvent.event_==10375||
//sigEvent.event_==11737;
//    if(!passEvent) continue;

    int charge = (int)(sigEvent.lq1_ + sigEvent.lq2_);

    bool passSystCuts[3][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[3][nSelTypes] = {{false, false, false, false, false, false, false, false, false, false},
                                   {false, false, false, false, false, false, false, false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(sigEvent.lep1McId_) == 11 || TMath::Abs(sigEvent.lep1McId_) == 13) &&
       (TMath::Abs(sigEvent.lep2McId_) == 11 || TMath::Abs(sigEvent.lep2McId_) == 13)) isRealLepton = true;

    double metNew[2]; metChange(sigEvent.met_,sigEvent.metPhi_,sigEvent.nvtx_,year,sigEvent.dstype_ == SmurfTree::data,metNew);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 
    //bool passBtagVeto      = (sigEvent.cuts_ & patternTopVeto) == patternTopVeto;
    //bool passBtagVeto      = sigEvent.nSoftMuons_ == 0 && (sigEvent.jet1_.Pt() <= 20 || sigEvent.jet1Btag_     < 2.100) && (sigEvent.jet2_.Pt() <= 20 || sigEvent.jet2Btag_     < 2.100) && (sigEvent.jet3_.Pt() <= 20 || sigEvent.jet3Btag_     < 2.100) && (sigEvent.jet4_.Pt() <= 20 || sigEvent.jet4Btag_     < 2.100);
    bool passBtagVeto      = sigEvent.nSoftMuons_ == 0 && (sigEvent.jet1_.Pt() <= 20 || sigEvent.jet1ProbBtag_ < 0.244) && (sigEvent.jet2_.Pt() <= 20 || sigEvent.jet2ProbBtag_ < 0.244) && (sigEvent.jet3_.Pt() <= 20 || sigEvent.jet3ProbBtag_ < 0.244) && (sigEvent.jet4_.Pt() <= 20 || sigEvent.jet4ProbBtag_ < 0.244);
    bool pass3rLVeto       = (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passZMass         = fabs(sigEvent.dilep_.M()-91.1876) < 15.;
    bool passZMassLarge    = fabs(sigEvent.dilep_.M()-91.1876) < 30.;
    bool passZMassSB       = (sigEvent.dilep_.M() > 110.0 && sigEvent.dilep_.M() < 200.0);
    bool passMET           = theMET > cutValue[0];
    bool passDPhiZMET      = DeltaPhi(sigEvent.dilep_.Phi() ,theMETPHI) > cutValue[1];
    bool passPTFrac        = fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt() < cutValue[2];
    bool passDPhiLL        = sigEvent.dPhi_*180.0/TMath::Pi() < cutValue[3];
    bool passPTLL          = sigEvent.dilep_.Pt() > 60.;
    int lType = -1;
    if     (sigEvent.type_ == SmurfTree::em || sigEvent.type_ == SmurfTree::me) lType = 0;
    else if(lSel == 0 && sigEvent.type_ == SmurfTree::mm) lType = 1;
    else if(lSel == 1 && sigEvent.type_ == SmurfTree::ee) lType = 1;
    else if(lSel == 2 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) lType = 1;
    else lType = 2;
    bool passBCK1REG       = !passDPhiZMET &&  passPTFrac &&
 			      fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt() < 1.00;
    bool passBCK2REG	   =  passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt() < 1.00;
    bool passBCK12REG	   = !passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt() < 1.00;

    //printf("PRESEL %d %6d %d %d %d %d %d %d %d %d %d %d | %6.2f %6.3f %6.2f %6.2f\n",sigEvent.type_,sigEvent.event_,sigEvent.njets_==nJetsType,sigEvent.lep1_.Pt()>20&&sigEvent.lep2_.Pt()>20,passBtagVeto,pass3rLVeto,passZMass,passMET,passPTLL,passDPhiLL,passDPhiZMET,passPTFrac,theMET,fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt(),sigEvent.jet1_.Pt(),sigEvent.jet2_.Pt());
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
    double MVAVar[6] = {outputVar[8],outputVar[8],outputVar[8],outputVarLepP[8],outputVarLepM[8],outputVarMET[8]};
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    if(isZZLepSel == false) {
      addLepEff	 = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 0)*
                   leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 0);
      if(addLepEff > 0) {
        addLepEffUp   = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 1)*
    			leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 1);
        addLepEffDown = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,-1)*
    			leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,-1);
      } else {addLepEff = 1.0;}
    } else {
      std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), sigEvent.lid1_);
      std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), sigEvent.lid2_);
      addLepEff	    = effZH1.first * effZH2.first;
      addLepEffUp   = (effZH1.first+effZH1.second)*(effZH2.first+effZH2.second);
      addLepEffDown = (effZH1.first-effZH1.second)*(effZH2.first-effZH2.second);
      if(addLepEffDown < 0) addLepEffDown = 0.;
    }
    double NjetSyst[2] = {0., 0.};
    if(sigEvent.jet1_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(sigEvent.jet2_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(sigEvent.jet3_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(sigEvent.jet4_.Pt()*1.05 > ptJet2nd) NjetSyst[0]++;
    if(sigEvent.jet1_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(sigEvent.jet2_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(sigEvent.jet3_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;
    if(sigEvent.jet4_.Pt()*0.95 > ptJet2nd) NjetSyst[1]++;

    if(channel > 1000 && channel < 2000){ // ZH->2l+inv selection
      if(charge == 0 && NjetSyst[0] == nJetsType									   && outputVar[4]     > metMin && passDPhiLL && outputVar[0]	   > 20.0 &&  outputVar[1]    > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVar[2]-91.1876)     < 15. && outputVar[4]	> cutValue[0] && outputVar[9]     > cutValue[1] && outputVar[11]     < cutValue[2]) passSystCuts[lType][JESUP] = true;
      if(charge == 0 && NjetSyst[1] == nJetsType									   && outputVar[4]     > metMin && passDPhiLL && outputVar[0]	   > 20.0 &&  outputVar[1]    > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVar[2]-91.1876)     < 15. && outputVar[4]	> cutValue[0] && outputVar[9]     > cutValue[1] && outputVar[11]     < cutValue[2]) passSystCuts[lType][JESDOWN] = true;
      if(charge == 0 && sigEvent.jet1_.Pt() > ptJetMin && sigEvent.jet1_.Pt() < ptJet1st && sigEvent.jet2_.Pt() < ptJet2nd && outputVarLepP[4] > metMin && passDPhiLL && outputVarLepP[0] > 20.0 &&  outputVarLepP[1] > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarLepP[2]-91.1876) < 15. && outputVarLepP[4] > cutValue[0] && outputVarLepP[9] > cutValue[1] && outputVarLepP[11] < cutValue[2]) passSystCuts[lType][LEPP] = true;
      if(charge == 0 && sigEvent.jet1_.Pt() > ptJetMin && sigEvent.jet1_.Pt() < ptJet1st && sigEvent.jet2_.Pt() < ptJet2nd && outputVarLepM[4] > metMin && passDPhiLL && outputVarLepM[0] > 20.0 &&  outputVarLepM[1] > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarLepM[2]-91.1876) < 15. && outputVarLepM[4] > cutValue[0] && outputVarLepM[9] > cutValue[1] && outputVarLepM[11] < cutValue[2]) passSystCuts[lType][LEPM] = true;
      if(charge == 0 && sigEvent.jet1_.Pt() > ptJetMin && sigEvent.jet1_.Pt() < ptJet1st && sigEvent.jet2_.Pt() < ptJet2nd && outputVarMET[4]  > metMin && passDPhiLL && outputVarMET[0]  > 20.0 &&  outputVarMET[1]  > 20.0 && passBtagVeto && pass3rLVeto && fabs(outputVarMET[2]-91.1876)  < 15. && outputVarMET[4]  > cutValue[0] && outputVarMET[9]  > cutValue[1] && outputVarMET[11]  < cutValue[2]) passSystCuts[lType][MET] = true;
      if(
         charge == 0 &&
	 sigEvent.jet1_.Pt() > ptJetMin && sigEvent.jet1_.Pt() < ptJet1st && sigEvent.jet2_.Pt() < ptJet2nd &&
	 theMET > metMin &&
         sigEvent.lep1_.Pt() > 20. &&
         sigEvent.lep2_.Pt() > 20.) {
	 
         passCuts[lType][ZSEL]   = true;	 
	 if(!passBtagVeto &&  pass3rLVeto && !passZMass && passZMassSB                                                                      ) passCuts[lType][WWLOOSESEL] = true;
	 if( passBtagVeto &&  pass3rLVeto && !passZMass && passZMassLarge && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][WWSEL]	  = true;
	 if(!passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][BTAGSEL]	  = true;
	 if( passBtagVeto && !pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac && sigEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL]    = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGSEL]   = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK1REG			  ) passCuts[lType][BCK1REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK2REG			  ) passCuts[lType][BCK2REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK12REG			  ) passCuts[lType][BCK12REG] = true;
         if     (sigEvent.type_ == SmurfTree::mm && pass3rLVeto && passPTLL && passZMassLarge) passCuts[0][ZLLSEL] = true;
         else if(sigEvent.type_ == SmurfTree::ee && pass3rLVeto && passPTLL && passZMassLarge) passCuts[1][ZLLSEL] = true;

	if(isRealLepton == false) {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false; passCuts[0][ZLLSEL] = false; passCuts[1][ZLLSEL] = false;}
      }
    } // ZH->2l+inv selection

    if(isRealLepton == true){
      double add = 1.;
      double addPU = 1.;
      addPU = nPUScaleFactor2012(fhDPU,sigEvent.npu_);
      if(isZZLepSel == false) {
        add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
        add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
        if((sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add = add*leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_);
      } else {
        std::pair<float,float> effZH1 = theLeptonEfficiencyZH.getLeptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), sigEvent.lid1_);
        std::pair<float,float> effZH2 = theLeptonEfficiencyZH.getLeptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), sigEvent.lid2_);
        std::pair<float,float> effZH3 = theLeptonEfficiencyZH.getLeptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), sigEvent.lid3_);
        add = add * effZH1.first * effZH2.first;
        if((sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) add = add * effZH3.first;
      }
      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  							       fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        						       TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
      add = add*trigEff*addPU;

      double theWeight = sigEvent.scale1fb_*lumi*add;
      if(useWeightEWKCorr == true) theWeight = theWeight * weightEWKCorr(sigEvent.dilep_.Pt(),0);

      if(fCheckProblem == true && TMath::Abs((sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_)-add)/add>0.05 && sigEvent.sfWeightFR_ > 0 && sigEvent.lid3_ != 0) {
	printf("PROBLEM: %f - %f %f %f %f %f = %f\n",add,sigEvent.sfWeightFR_,sigEvent.sfWeightPU_,sigEvent.sfWeightEff_,sigEvent.sfWeightTrig_,sigEvent.sfWeightHPt_,sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_*sigEvent.sfWeightHPt_);
	printf("  		%f %f %f %f\n",1.0,addPU,leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_)*
  	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_), trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  																			     fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
  																			     TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_)));
      }

      if(passCuts[1][SIGSEL] && theMET > metMin){ // begin making plots
        //printf("FULLSEL %d %6d %d\n",sigEvent.type_,sigEvent.event_,sigEvent.njets_);
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
	else if(thePlot ==12) myVar = sigEvent.trackMet_;
	else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = sigEvent.njets_;
	else if(thePlot ==18) myVar = sigEvent.nvtx_;
	else if(thePlot ==19) myVar = sigEvent.pTrackMet_;
	else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(sigEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(sigEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = collins_Soper(sigEvent.lep1_,sigEvent.lep2_,sigEvent.lq1_);
	else if(thePlot ==25) myVar = fabs(sigEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = sigEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar =DeltaPhi(sigEvent.dilep_.Phi() ,sigEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Min(TMath::Abs(sigEvent.dilep_.Eta()),2.4999);
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
	else if(thePlot ==57) myVar = sigEvent.dR_;
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
      if(passCuts[lType][SIGSEL]) {
        nSigCutSyst[EFFP+lType*nSelTypesSyst]  += theWeight          *addLepEffUp  /addLepEff;
        nSigECutSyst[EFFP+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        nSigCutSyst[EFFM+lType*nSelTypesSyst]  += theWeight          *addLepEffDown/addLepEff;
        nSigECutSyst[EFFM+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }
      if(1){
        if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv                           ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv_CMS_MVALepEffBoundingUp   ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv_CMS_MVALepEffBoundingDown ->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_ZH_hinv_CMS_MVAJESBoundingUp	     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_ZH_hinv_CMS_MVAJESBoundingDown    ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_ZH_hinv_CMS_MVALepResBoundingUp   ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_ZH_hinv_CMS_MVALepResBoundingDown ->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_ZH_hinv_CMS_MVAMETResBoundingUp   ->Fill(MVAVar[5], theWeight);;
      }
    } // if passCuts
  } // Loop over signal
  
  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;
    int charge = (int)(dataEvent.lq1_ + dataEvent.lq2_);

    bool passCuts[3][nSelTypes] = {{false, false, false, false, false, false, false, false, false, false},
                                   {false, false, false, false, false, false, false, false, false, false}};

    double metNew[3]; metChange(dataEvent.met_,dataEvent.metPhi_,dataEvent.nvtx_,year,dataEvent.dstype_ == SmurfTree::data,metNew);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 
    //bool passBtagVeto      = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    //bool passBtagVeto      = dataEvent.nSoftMuons_ == 0 && (dataEvent.jet1_.Pt() <= 20 || dataEvent.jet1Btag_     < 2.100) && (dataEvent.jet2_.Pt() <= 20 || dataEvent.jet2Btag_     < 2.100) && (dataEvent.jet3_.Pt() <= 20 || dataEvent.jet3Btag_     < 2.100) && (dataEvent.jet4_.Pt() <= 20 || dataEvent.jet4Btag_     < 2.100);
    bool passBtagVeto      = dataEvent.nSoftMuons_ == 0 && (dataEvent.jet1_.Pt() <= 20 || dataEvent.jet1ProbBtag_ < 0.244) && (dataEvent.jet2_.Pt() <= 20 || dataEvent.jet2ProbBtag_ < 0.244) && (dataEvent.jet3_.Pt() <= 20 || dataEvent.jet3ProbBtag_ < 0.244) && (dataEvent.jet4_.Pt() <= 20 || dataEvent.jet4ProbBtag_ < 0.244);
    bool pass3rLVeto       = (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passZMass         = fabs(dataEvent.dilep_.M()-91.1876) < 15.;
    bool passZMassLarge    = fabs(dataEvent.dilep_.M()-91.1876) < 30.;
    bool passZMassSB       = (dataEvent.dilep_.M() > 110.0 && dataEvent.dilep_.M() < 200.0);
    bool passMET           = theMET > cutValue[0];
    bool passDPhiZMET      = DeltaPhi(dataEvent.dilep_.Phi() ,theMETPHI) > cutValue[1];
    bool passPTFrac        = fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt() < cutValue[2];
    bool passDPhiLL        = dataEvent.dPhi_*180.0/TMath::Pi() < cutValue[3];
    bool passPTLL          = dataEvent.dilep_.Pt() > 60.;
    int lType = -1;
    if     (dataEvent.type_ == SmurfTree::em || dataEvent.type_ == SmurfTree::me) lType = 0;
    else if(lSel == 0 && dataEvent.type_ == SmurfTree::mm) lType = 1;
    else if(lSel == 1 && dataEvent.type_ == SmurfTree::ee) lType = 1;
    else if(lSel == 2 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) lType = 1;
    else lType = 2;
    bool passBCK1REG       = !passDPhiZMET &&  passPTFrac &&
 			      fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt() < 1.00;
    bool passBCK2REG	   =  passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt() < 1.00;
    bool passBCK12REG	   = !passDPhiZMET && !passPTFrac &&
 			      fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt() < 1.00;

    if(channel > 1000 && channel < 2000){ // ZH->2l+inv selection
      if(
         charge == 0 &&
	 dataEvent.jet1_.Pt() > ptJetMin && dataEvent.jet1_.Pt() < ptJet1st && dataEvent.jet2_.Pt() < ptJet2nd &&
	 theMET > metMin &&
         dataEvent.lep1_.Pt() > 20. &&
         dataEvent.lep2_.Pt() > 20.) {
	 
         passCuts[lType][ZSEL]   = true;	 
	 if(!passBtagVeto &&  pass3rLVeto && !passZMass && passZMassSB                                                                      ) passCuts[lType][WWLOOSESEL] = true;
	 if( passBtagVeto &&  pass3rLVeto && !passZMass && passZMassLarge && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][WWSEL]      = true;
	 if(!passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][BTAGSEL]    = true;
	 if( passBtagVeto && !pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac && dataEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL]    = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGSEL]   = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK1REG		       ) passCuts[lType][BCK1REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK2REG		       ) passCuts[lType][BCK2REG]  = true;
	 if( passBtagVeto &&  pass3rLVeto &&  passZMass                   && passMET && passPTLL && passDPhiLL && passBCK12REG		       ) passCuts[lType][BCK12REG] = true;
         if     (dataEvent.type_ == SmurfTree::mm && pass3rLVeto && passPTLL && passZMassLarge) passCuts[0][ZLLSEL] = true;
         else if(dataEvent.type_ == SmurfTree::ee && pass3rLVeto && passPTLL && passZMassLarge) passCuts[1][ZLLSEL] = true;

      }
    } // ZH->2l+inv selection

    if(passCuts[lType][ZSEL]){

      if(passCuts[1][SIGSEL] && theMET > metMin){ // begin making plots
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
	else if(thePlot ==12) myVar = dataEvent.trackMet_;
	else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = dataEvent.njets_;
	else if(thePlot ==18) myVar = dataEvent.nvtx_;
	else if(thePlot ==19) myVar = dataEvent.pTrackMet_;
	else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(dataEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = collins_Soper(dataEvent.lep1_,dataEvent.lep2_,dataEvent.lq1_);
	else if(thePlot ==25) myVar = fabs(dataEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = dataEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,dataEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Min(TMath::Abs(dataEvent.dilep_.Eta()),2.4999);
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
	else if(thePlot ==57) myVar = dataEvent.dR_;
      	histo5->Fill(myVar,1.0);
      } // end making plots

      double outputVar[15];
      makeSystematicEffects(dataEvent.lid1_, dataEvent.lid2_, dataEvent.lep1_, dataEvent.lep2_, dataEvent.dilep_, 
                            dataEvent.mt_, theMET, theMETPHI, 
                            dataEvent.trackMet_, dataEvent.trackMetPhi_, 
			    dataEvent.njets_, dataEvent.jet1_, dataEvent.jet2_, 
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[8],0,0,0,0,0};
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][SIGSEL]){
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

  char output[200];
  sprintf(output,Form("histo_nice%s.root",ECMsb.Data()));	 
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
    double oldNorm0 = histo0->GetSumOfWeights();
    double oldNorm4 = histo4->GetSumOfWeights();
    for(int i=1; i<=histo0->GetNbinsX(); i++){
      if(histo0->GetBinContent(i) < 0) histo0->SetBinContent(i,0.0);
      if(histo4->GetBinContent(i) < 0) histo4->SetBinContent(i,0.0);
    }
    if(histo0->GetSumOfWeights() > 0) histo0->Scale(oldNorm0/histo0->GetSumOfWeights());
    if(histo4->GetSumOfWeights() > 0) histo4->Scale(oldNorm4/histo4->GetSumOfWeights());
    printf("histo -> s: %8.2f d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f\n",histos->GetSumOfWeights(),histo5->GetSumOfWeights(),histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights());
    double scaleZ = 1.0;
    if(histo0->GetSumOfWeights() > 0) scaleZ = (histo5->GetSumOfWeights()-(histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights()))/histo0->GetSumOfWeights();
    if(scaleZ <= 0) scaleZ = 1.0;
    printf("scaleZ = %f\n",scaleZ);
    histo0->Scale(scaleZ);
    histos->Write();
    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();
    histo_GAMMA_JETS->Write();
    histo_GAMMA_JETS_LOOSE->Write();

    TH1D* histoNoZ = (TH1D*) histo1->Clone("histoNoZ");
    histoNoZ->Add(histo2);
    histoNoZ->Add(histo3);
    histoNoZ->Add(histo4);
    TH1D* histoZ   = (TH1D*) histo0->Clone("histoZ");
    TH1D* histoD   = (TH1D*) histo5->Clone("histoD");  
    printf("%8.2f %8.2f %8.2f\n",histoNoZ->GetSumOfWeights(),histoZ->GetSumOfWeights(),histoD->GetSumOfWeights());
    histoD->Add(histoNoZ,-1.0);
    histoD->Write();
    histoZ->Write();

  outFilePlotsNote->Close();
  
  const unsigned int nBkg = 7;
  double nTot[nSelTypes*2]; double nETot[nSelTypes*2];
  double bgdCombined[nSelTypes*2][nBkg],bgdCombinedE[nSelTypes*2][nBkg];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) {bgdCombined[i][j] = 0.0; bgdCombinedE[i][j] = 0.0;}
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("selection: %s\n",selTypeName[i].Data());
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("data(%2d): %f\n",i,nSelectedData[i]);
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("nSigCut(%2d): %11.3f +/- %8.3f\n",i,nSigCut[i],sqrt(nSigECut[i]));
    nTot[i] = 0.0; nETot[i] = 0.0;
    for(int j=0; j<45; j++){
      if(showSignalOnly == false || i%nSelTypes == SIGSEL) if(bgdDecay[i][j] != 0) printf("bdg(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecay[i][j],sqrt(weiDecay[i][j]));
      nTot[i]  += bgdDecay[i][j];
      nETot[i] += weiDecay[i][j];

      if     (j == 9 || j == 19)               {bgdCombined[i][0] += bgdDecay[i][j]; bgdCombinedE[i][0] += weiDecay[i][j];}
      else if(j == 11)			       {bgdCombined[i][1] += bgdDecay[i][j]; bgdCombinedE[i][1] += weiDecay[i][j];}
      else if(j == 17)			       {bgdCombined[i][2] += bgdDecay[i][j]; bgdCombinedE[i][2] += weiDecay[i][j];}
      else if(j == 18)                         {bgdCombined[i][3] += bgdDecay[i][j]; bgdCombinedE[i][3] += weiDecay[i][j];}
      else if(j == 21 || j == 27 || j == 28 ||
              j == 29 || j == 30 || j ==  5 ||
	      j == 13 || j == 20 || j == 10 ||
	      j == 1  || j == 23)              {bgdCombined[i][4] += bgdDecay[i][j]; bgdCombinedE[i][4] += weiDecay[i][j];}
      //else if(j == 1)                          {bgdCombined[i][5] += bgdDecay[i][j]; bgdCombinedE[i][5] += weiDecay[i][j];}
      //else if(j == 23)                         {bgdCombined[i][6] += bgdDecay[i][j]; bgdCombinedE[i][6] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]-nTot[i])/nTot[i] > 0.00001) {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xxZ) = %11.3f +/- %8.3f\n",bgdCombined[i][0],sqrt(bgdCombinedE[i][0]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(VVV) = %11.3f +/- %8.3f\n",bgdCombined[i][1],sqrt(bgdCombinedE[i][1]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xWZ) = %11.3f +/- %8.3f\n",bgdCombined[i][2],sqrt(bgdCombinedE[i][2]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xZZ) = %11.3f +/- %8.3f\n",bgdCombined[i][3],sqrt(bgdCombinedE[i][3]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xEM) = %11.3f +/- %8.3f\n",bgdCombined[i][4],sqrt(bgdCombinedE[i][4]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(WjE) = %11.3f +/- %8.3f\n",bgdCombined[i][5],sqrt(bgdCombinedE[i][5]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(WjM) = %11.3f +/- %8.3f\n",bgdCombined[i][6],sqrt(bgdCombinedE[i][6]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("*******************************\n");
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

      if     (j == 9 || j == 19)               {bgdCombinedSyst[i][0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][0] += weiDecaySyst[i][j];}
      else if(j == 11)			       {bgdCombinedSyst[i][1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][1] += weiDecaySyst[i][j];}
      else if(j == 17)			       {bgdCombinedSyst[i][2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][2] += weiDecaySyst[i][j];}
      else if(j == 18)                         {bgdCombinedSyst[i][3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][3] += weiDecaySyst[i][j];}
      else if(j == 21 || j == 27 || j == 28 ||
              j == 29 || j == 30 || j ==  5 ||
	      j == 13 || j == 20 || j == 10 ||
	      j == 1  || j == 23)              {bgdCombinedSyst[i][4] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][4] += weiDecaySyst[i][j];}
      //else if(j == 1)                          {bgdCombinedSyst[i][5] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][5] += weiDecaySyst[i][j];}
      //else if(j == 23)                         {bgdCombinedSyst[i][6] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][6] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]-nTotSyst[i])/nTotSyst[i] > 0.00001) {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]);assert(0);}
    if(showSignalOnly == false) printf("------\n");
    if(showSignalOnly == false) printf("bgdSyst(xxZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][0],sqrt(bgdCombinedESyst[i][0]));
    if(showSignalOnly == false) printf("bgdSyst(VVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][1],sqrt(bgdCombinedESyst[i][1]));
    if(showSignalOnly == false) printf("bgdSyst(xWZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][2],sqrt(bgdCombinedESyst[i][2]));
    if(showSignalOnly == false) printf("bgdSyst(xZZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][3],sqrt(bgdCombinedESyst[i][3]));
    if(showSignalOnly == false) printf("bgdSyst(xEM) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][4],sqrt(bgdCombinedESyst[i][4]));
    //if(showSignalOnly == false) printf("bgdSyst(WjE) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][5],sqrt(bgdCombinedESyst[i][5]));
    //if(showSignalOnly == false) printf("bgdSyst(WjM) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][6],sqrt(bgdCombinedESyst[i][6]));
    if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }
  double NFinal[nBkg+1]  = {0.,bgdCombined[SIGSEL+nSelTypes][1],bgdCombined[SIGSEL+nSelTypes][2],bgdCombined[SIGSEL+nSelTypes][3],bgdCombined[SIGSEL+nSelTypes][4],
                               bgdCombined[SIGSEL+nSelTypes][5]+bgdCombined[SIGSEL+nSelTypes][6],0.,
			       nSigCut[SIGSEL+nSelTypes]};
  double NFinalE[nBkg+1] = {0.,1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][1])/bgdCombined[SIGSEL+nSelTypes][1],
                     	       1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][2])/bgdCombined[SIGSEL+nSelTypes][2],1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][3])/bgdCombined[SIGSEL+nSelTypes][3],1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][4])/bgdCombined[SIGSEL+nSelTypes][4],
		               1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][5]+bgdCombinedE[SIGSEL+nSelTypes][6])/(bgdCombined[SIGSEL+nSelTypes][5]+bgdCombined[SIGSEL+nSelTypes][6]),0.0,
		               1.0+sqrt(nSigECut[SIGSEL+nSelTypes])/nSigCut[SIGSEL+nSelTypes]};

  // This is the alternative method
  double kEoverM  = sqrt(nSelectedData[ZLLSEL+nSelTypes]/nSelectedData[ZLLSEL]);
  double NemFact = (kEoverM+1.0/kEoverM)*0.5;
  if     (lSel == 0)  NemFact = (        1.0/kEoverM)*0.5;
  else if(lSel == 1)  NemFact = (kEoverM            )*0.5;

  // This is the default method
  double NemFact_FromMLLSB = 1.0; double NemFact_FromMLLSBE = 1.0;
  if(nSelectedData[WWLOOSESEL+nSelTypes] > 0 && nSelectedData[WWLOOSESEL] > 0) {
    NemFact_FromMLLSB  = nSelectedData[WWLOOSESEL+nSelTypes]/nSelectedData[WWLOOSESEL];
    NemFact_FromMLLSBE = sqrt(1.0/nSelectedData[WWLOOSESEL+nSelTypes]+
                              1.0/nSelectedData[WWLOOSESEL	    ])*NemFact_FromMLLSB;
  }
  else {NemFact_FromMLLSB = NemFact;}

  printf("kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact,NemFact_FromMLLSB,NemFact_FromMLLSBE);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  double EMSystTotal = 1.0; double EMSyst[2] = {0.0, 0.0};
  EMSyst[0] = bgdCombined[SIGSEL+nSelTypes][4]/bgdCombined[SIGSEL][4]/NemFact_FromMLLSB;
  if(EMSyst[0] < 1.0) EMSyst[0] = 1/EMSyst[0]; EMSyst[0] = EMSyst[0] - 1.0;

  EMSyst[1] = NemFact_FromMLLSB/NemFact;
  if(EMSyst[1] < 1.0) EMSyst[1] = 1/EMSyst[1]; EMSyst[1] = EMSyst[1] - 1.0;

  // At least one data event is required
  if(nSelectedData[SIGSEL] == 0) nSelectedData[SIGSEL] = 1;

  if(nSelectedData[SIGSEL] > 0) EMSystTotal = 1.0 + sqrt(EMSyst[0]*EMSyst[0] + EMSyst[1]*EMSyst[1] + 1/nSelectedData[SIGSEL]);
  else                          EMSystTotal = 1.0 + sqrt(EMSyst[0]*EMSyst[0] + EMSyst[1]*EMSyst[1] + 1.0);

  double EMbkg = bgdCombined[SIGSEL][1]+bgdCombined[SIGSEL][2]+bgdCombined[SIGSEL][3]+
                 bgdCombined[SIGSEL][5]+bgdCombined[SIGSEL][6]; // background in EM region
  if(properWjetsBkg == false) {
    EMbkg = bgdCombined[SIGSEL][1]+bgdCombined[SIGSEL][2]+bgdCombined[SIGSEL][3];
    bgdCombined[SIGSEL][5] = 0.0; NFinal[5]  = 0.0;
    bgdCombined[SIGSEL][6] = 0.0; NFinalE[5] = 0.0;
  }
  printf("EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdCombined[SIGSEL+nSelTypes][4],sqrt(bgdCombinedE[SIGSEL+nSelTypes][4]),
	 bgdCombined[SIGSEL][4]*NemFact_FromMLLSB  ,sqrt(bgdCombinedE[SIGSEL][4])*NemFact_FromMLLSB,nSelectedData[SIGSEL],EMbkg,EMSystTotal-1.0,EMSyst[0],EMSyst[1],sqrt((EMSystTotal-1.0)*(EMSystTotal-1.0) - EMSyst[0]*EMSyst[0] - EMSyst[1]*EMSyst[1]));

  if(nSelectedData[SIGSEL]-EMbkg > 0) {
    NFinal[4] = (nSelectedData[SIGSEL]-EMbkg)*NemFact_FromMLLSB;
  }
  NFinalE[4] = EMSystTotal;

  double scaleFactorEM = NFinal[4]/bgdCombined[SIGSEL+nSelTypes][4]; if(scaleFactorEM <= 0) scaleFactorEM = 0.0;
  printf("EM scale Factor: %5.3f +/- %5.3f\n",scaleFactorEM,EMSystTotal-1.0);
  printf("EM yield: %5.3f +/- %5.3f\n",NFinal[4],sqrt(nSelectedData[SIGSEL])*NemFact_FromMLLSB);

  double ZLLMCUnc = sqrt(bgdCombinedE[BCK2REG+nSelTypes][0]+bgdCombinedE[BCK12REG+nSelTypes][0])/
                        ( bgdCombined[BCK2REG+nSelTypes][0]+ bgdCombined[BCK12REG+nSelTypes][0]);

  printf("ZLL MC BCK2REG+BCK12REG = %f+%f = %8.3f +/- %5.3f --> %8.3f +/- %5.3f ==> %8.3f +/- %5.3f\n",
         bgdCombined[BCK2REG+nSelTypes][0],bgdCombined[BCK12REG+nSelTypes][0],
	 bgdCombined[BCK2REG+nSelTypes][0]+bgdCombined[BCK12REG+nSelTypes][0],
	(bgdCombined[BCK2REG+nSelTypes][0]+bgdCombined[BCK12REG+nSelTypes][0])*ZLLMCUnc,
	 bgdCombined[SIGSEL+nSelTypes][0],sqrt(bgdCombinedE[SIGSEL+nSelTypes][0]),
	 bgdCombined[SIGSEL+nSelTypes][0]/(bgdCombined[BCK2REG+nSelTypes][0]+bgdCombined[BCK12REG+nSelTypes][0]),
	 bgdCombined[SIGSEL+nSelTypes][0]/(bgdCombined[BCK2REG+nSelTypes][0]+bgdCombined[BCK12REG+nSelTypes][0])*sqrt(ZLLMCUnc*ZLLMCUnc+bgdCombinedE[SIGSEL+nSelTypes][0]/bgdCombined[SIGSEL+nSelTypes][0]/bgdCombined[SIGSEL+nSelTypes][0]));

  double ZdataBCK2REG  = nSelectedData[BCK2REG+nSelTypes] -bgdCombined[BCK2REG+nSelTypes][1] -bgdCombined[BCK2REG+nSelTypes][2]-
                                                           bgdCombined[BCK2REG+nSelTypes][3] -bgdCombined[BCK2REG+nSelTypes][4];
  double ZdataBCK12REG = nSelectedData[BCK12REG+nSelTypes]-bgdCombined[BCK12REG+nSelTypes][1]-bgdCombined[BCK12REG+nSelTypes][2]-
                                                           bgdCombined[BCK12REG+nSelTypes][3]-bgdCombined[BCK12REG+nSelTypes][4];

  double ZdataSIGSEL  = nSelectedData[SIGSEL+nSelTypes] -bgdCombined[SIGSEL+nSelTypes][1] -bgdCombined[SIGSEL+nSelTypes][2]-
                                                         bgdCombined[SIGSEL+nSelTypes][3] -bgdCombined[SIGSEL+nSelTypes][4];

  double ZdataBCK2REGE2 = (nSelectedData[BCK2REG+nSelTypes]
         +bgdCombinedE[BCK2REG+nSelTypes][1]+bgdCombinedE[BCK2REG+nSelTypes][2]
         +bgdCombinedE[BCK2REG+nSelTypes][3]+bgdCombinedE[BCK2REG+nSelTypes][4]);
  double ZdataBCK12REGE2 = (nSelectedData[BCK12REG+nSelTypes]
         +bgdCombinedE[BCK12REG+nSelTypes][1]+bgdCombinedE[BCK12REG+nSelTypes][2]
         +bgdCombinedE[BCK12REG+nSelTypes][3]+bgdCombinedE[BCK12REG+nSelTypes][4]);

  double ZdataSIGSEL2 = (nSelectedData[SIGSEL+nSelTypes]
         +bgdCombinedE[SIGSEL+nSelTypes][1]+bgdCombinedE[SIGSEL+nSelTypes][2]
         +bgdCombinedE[SIGSEL+nSelTypes][3]+bgdCombinedE[SIGSEL+nSelTypes][4]);

  double ZdataBCK = ZdataBCK2REG+ZdataBCK12REG;
  double ZdataBCKE = sqrt(ZdataBCK2REGE2+ZdataBCK12REGE2)/
                         (ZdataBCK2REG  +ZdataBCK12REG  );

  printf("ZLL DA = (%f +/- %f) + (%f +/- %f) = (%f +/- %f) --> %8.3f +/- %5.3f ==> %8.3f +/- %5.3f\n",
         ZdataBCK2REG ,sqrt(ZdataBCK2REGE2),
         ZdataBCK12REG,sqrt(ZdataBCK12REGE2),
         ZdataBCK     ,ZdataBCKE*ZdataBCK,
	 ZdataSIGSEL  ,sqrt(ZdataSIGSEL2),
	 ZdataSIGSEL/ZdataBCK,
	 ZdataSIGSEL/ZdataBCK*sqrt(ZdataBCKE*ZdataBCKE+ZdataSIGSEL2/ZdataSIGSEL/ZdataSIGSEL));

  NFinal[0] = bgdCombined[SIGSEL+nSelTypes][0]; NFinalE[0] = 2.0;
  if(ZdataBCK > 0.0) NFinal[0] = ZdataBCK*0.20;
  printf("ZLL yield: %8.3f +/- %5.3f\n",ZdataBCK*0.20,ZdataBCKE*ZdataBCK*0.20);

  double QCDscale_VH = 1.072;
  if(nJetsType == 1) QCDscale_VH = 1.105;
  
  double systEffect[nSelTypesSyst][nBkg+1];
  for(unsigned int i=0 ; i<nSelTypesSyst; i++){
    for(unsigned int j=0 ; j<nBkg; j++){
      if(bgdCombinedE[SIGSEL+nSelTypes][j] > 0){
        systEffect[i][j] = bgdCombinedSyst[i+nSelTypesSyst][j]/bgdCombined[SIGSEL+nSelTypes][j];
        if(systEffect[i][j] < 1) systEffect[i][j] = 1.0/systEffect[i][j];
      } else {systEffect[i][j] = 1.0;}
    }
    systEffect[i][nBkg] = nSigCutSyst[i+nSelTypesSyst]/nSigCut[SIGSEL+nSelTypes];    
    if(systEffect[i][nBkg] < 1) systEffect[i][nBkg] = 1.0/systEffect[i][nBkg];
  }
  if(showSignalOnly == false) printf("Syst(sig) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][nBkg]-1,systEffect[JESDOWN][nBkg]-1,systEffect[LEPP][nBkg]-1,systEffect[LEPM][nBkg]-1,systEffect[MET][nBkg]-1,systEffect[EFFP][nBkg]-1,systEffect[EFFM][nBkg]-1);
  if(showSignalOnly == false) printf("Syst(xxZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][0]-1,systEffect[JESDOWN][0]-1,systEffect[LEPP][0]-1,systEffect[LEPM][0]-1,systEffect[MET][0]-1,systEffect[EFFP][0]-1,systEffect[EFFM][0]-1);
  if(showSignalOnly == false) printf("Syst(VVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][1]-1,systEffect[JESDOWN][1]-1,systEffect[LEPP][1]-1,systEffect[LEPM][1]-1,systEffect[MET][1]-1,systEffect[EFFP][1]-1,systEffect[EFFM][1]-1);
  if(showSignalOnly == false) printf("Syst(xWZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][2]-1,systEffect[JESDOWN][2]-1,systEffect[LEPP][2]-1,systEffect[LEPM][2]-1,systEffect[MET][2]-1,systEffect[EFFP][2]-1,systEffect[EFFM][2]-1);
  if(showSignalOnly == false) printf("Syst(xZZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][3]-1,systEffect[JESDOWN][3]-1,systEffect[LEPP][3]-1,systEffect[LEPM][3]-1,systEffect[MET][3]-1,systEffect[EFFP][3]-1,systEffect[EFFM][3]-1);
  if(showSignalOnly == false) printf("Syst(xEM) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][4]-1,systEffect[JESDOWN][4]-1,systEffect[LEPP][4]-1,systEffect[LEPM][4]-1,systEffect[MET][4]-1,systEffect[EFFP][4]-1,systEffect[EFFM][4]-1);

  double WjetsSyst = 1.0;
  if(bgdCombined[SIGSEL+nSelTypes][5]+bgdCombined[SIGSEL+nSelTypes][6] > 0){
    WjetsSyst = sqrt(bgdCombined[SIGSEL+nSelTypes][5]*bgdCombined[SIGSEL+nSelTypes][5]*0.36*0.36+
                     bgdCombined[SIGSEL+nSelTypes][6]*bgdCombined[SIGSEL+nSelTypes][6]*0.36*0.36);
    WjetsSyst = 1.0+WjetsSyst/(bgdCombined[SIGSEL+nSelTypes][5]+bgdCombined[SIGSEL+nSelTypes][6]);
  }
  if(showSignalOnly == false) printf("WjetsSyst: %f %f --> %f\n",bgdCombined[SIGSEL+nSelTypes][5],bgdCombined[SIGSEL+nSelTypes][6],WjetsSyst);
  double pdf_qqbar[3] = {1.055,1.048,1.057};
  if(nJetsType == 1) {pdf_qqbar[0] = 1.060; pdf_qqbar[0] = 1.043; pdf_qqbar[0] = 1.043;}
  double syst_WZ3l = 1.010;
  if(nJetsType == 1) syst_WZ3l = 1.013;
  double syst_btag = 1.000;
  if     (year == 2011 && nJetsType == 0) syst_btag = 1.002;
  else if(year == 2012 && nJetsType == 0) syst_btag = 1.007;
  else if(year == 2011 && nJetsType == 1) syst_btag = 1.027;
  else if(year == 2012 && nJetsType == 1) syst_btag = 1.018;
  else assert(0);

  char outputLimitsCut[200];
  sprintf(outputLimitsCut,"histo_limits_zllhinv%2s_mh%3d_cut_%4s.txt",finalStateName,mH,ECMsb.Data());
  ofstream newcardCut;
  newcardCut.open(outputLimitsCut);
  newcardCut << Form("imax 1 number of channels\n");
  newcardCut << Form("jmax * number of background\n");
  newcardCut << Form("kmax * number of nuisance parameters\n");
  newcardCut << Form("Observation %d\n",(int)nSelectedData[SIGSEL+nSelTypes]);
  newcardCut << Form("bin hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s\n",finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data());
  newcardCut << Form("process ZH_hinv Zjets VVV WZ ZZ EM Wjets\n");
  newcardCut << Form("process 0 1 2 3 4 5 6\n");
  newcardCut << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",NFinal[nBkg],NFinal[0],NFinal[1],NFinal[2],NFinal[3],NFinal[4],TMath::Max(NFinal[5],0.0));
  newcardCut << Form("lumi_%4s                    lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE);		     
  newcardCut << Form("%s                          lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",effName,(systEffect[EFFP][nBkg]+systEffect[EFFM][nBkg])/2.0,(systEffect[EFFP][1]+systEffect[EFFM][1])/2.0,(systEffect[EFFP][2]+systEffect[EFFM][2])/2.0,(systEffect[EFFP][3]+systEffect[EFFM][3])/2.0);		   
  newcardCut << Form("%s                          lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",momName,(systEffect[LEPP][nBkg]+systEffect[LEPM][nBkg])/2.0,(systEffect[LEPP][1]+systEffect[LEPM][1])/2.0,(systEffect[LEPP][2]+systEffect[LEPM][2])/2.0,(systEffect[LEPP][3]+systEffect[LEPM][3])/2.0);		   
  newcardCut << Form("CMS_scale_met		  lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",systEffect[MET][nBkg],systEffect[MET][1],systEffect[MET][2],systEffect[MET][3]);
  newcardCut << Form("CMS_scale_j 		  lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",(systEffect[JESUP][nBkg]+systEffect[JESDOWN][nBkg])/2.0,(systEffect[JESUP][1]+systEffect[JESDOWN][1])/2.0,(systEffect[JESUP][2]+systEffect[JESDOWN][2])/2.0,(systEffect[JESUP][3]+systEffect[JESDOWN][3])/2.0);
  newcardCut << Form("CMS_eff_b                   lnN %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",syst_btag,syst_btag,syst_btag,syst_btag);
  newcardCut << Form("UEPS			  lnN 1.030   -     -     -     -     -     -  \n");
  newcardCut << Form("pdf_qqbar                   lnN %5.3f   -     -   %5.3f %5.3f   -     -  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2]);
  newcardCut << Form("QCDscale_VH		  lnN %5.3f   -     -     -     -     -     -  \n",QCDscale_VH);	
  newcardCut << Form("QCDscale_VV		  lnN   -     -     -   1.107 1.065   -     -  \n");		
  newcardCut << Form("CMS_zllhinv_WZ3l            lnN   -     -     -   %5.3f   -     -     -  \n",syst_WZ3l);		
  if(NFinal[1] > 0)
  newcardCut << Form("QCDscale_VVV		  lnN   -     -   1.500   -     -     -     -  \n");		
  if(NFinal[5] > 0)
  newcardCut << Form("CMS_FakeRate                lnN   -     -     -	 -     -      -    %5.3f\n",WjetsSyst);  
  newcardCut << Form("CMS_zllhinv_ZLL_%4s         lnN   -   %5.3f   -     -     -     -     -  \n",ECMsb.Data(),NFinalE[0]);		
  if(NFinal[4] > 0)
  newcardCut << Form("CMS_zllhinv_EM_%4s          lnN   -     -     -     -     -   %5.3f   -  \n",ECMsb.Data(),NFinalE[4]);		
  newcardCut << Form("CMS_zllhinv%2s_stat_ZH_hinv_%4s  lnN %5.3f   -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data(),NFinalE[nBkg]);  
  if(NFinal[1] > 0)
  newcardCut << Form("CMS_zllhinv%2s_stat_VVV_%4s lnN   -     -	  %5.3f   -     -     -	    -  \n",finalStateName,ECMsb.Data(),NFinalE[1]);  
  newcardCut << Form("CMS_zllhinv%2s_stat_WZ_%4s  lnN   -     -	    -   %5.3f   -     -	    -  \n",finalStateName,ECMsb.Data(),NFinalE[2]);  
  newcardCut << Form("CMS_zllhinv%2s_stat_ZZ_%4s  lnN   -     -	    -     -   %5.3f   -	    -  \n",finalStateName,ECMsb.Data(),NFinalE[3]);  
  if(NFinal[5] > 0)
  newcardCut << Form("CMS_zllhinv%2s_stat_Wjets_%4s lnN	-     -     -     -     -     -	  %5.3f\n",finalStateName,ECMsb.Data(),NFinalE[5]);  
  newcardCut.close();

  double nOld = histo_Wjets->GetSumOfWeights();
  for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
    if(histo_Wjets->GetBinContent(i)                    <= 0) {histo_Wjets                   ->SetBinContent(i,0.000001);histo_Wjets                   ->SetBinError(i,0.000001);}
    if(histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i) <= 0) {histo_Wjets_CMS_MVAWBoundingUp->SetBinContent(i,0.000001);histo_Wjets_CMS_MVAWBoundingUp->SetBinError(i,0.000001);}
  }
  histo_Wjets                   ->Scale(nOld/histo_Wjets                   ->GetSumOfWeights());
  histo_Wjets_CMS_MVAWBoundingUp->Scale(nOld/histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());

  histo_EM->Scale(scaleFactorEM);

  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp     ->SetBinContent(i,TMath::Max(histo_ZH_hinv    ->GetBinContent(i)+factorUp  *histo_ZH_hinv 	->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown   ->SetBinContent(i,TMath::Max(histo_ZH_hinv    ->GetBinContent(i)+factorDown*histo_ZH_hinv 	->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_VVV   ->GetBinContent(i)+factorUp  *histo_VVV   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_VVV   ->GetBinContent(i)+factorDown*histo_VVV   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_WZ    ->GetBinContent(i)+factorUp  *histo_WZ	 ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_WZ    ->GetBinContent(i)+factorDown*histo_WZ	 ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_ZZ    ->GetBinContent(i)+factorUp  *histo_ZZ	 ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_ZZ    ->GetBinContent(i)+factorDown*histo_ZZ	 ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_EM    ->GetBinContent(i)+factorUp  *histo_EM	 ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_EM    ->GetBinContent(i)+factorDown*histo_EM	 ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingUp    ->SetBinContent(i,TMath::Max(histo_Wjets->GetBinContent(i)+factorUp  *histo_Wjets  ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingDown  ->SetBinContent(i,TMath::Max(histo_Wjets->GetBinContent(i)+factorDown*histo_Wjets  ->GetBinError(i),0.000001));

    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]    ->Add(histo_ZH_hinv   ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_ZH_hinv   ->GetBinContent(i)+factorUp  *histo_ZH_hinv	 ->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]  ->Add(histo_ZH_hinv   ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_ZH_hinv   ->GetBinContent(i)+factorDown*histo_ZH_hinv	 ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	     ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorUp  *histo_VVV   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]     ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorDown*histo_VVV   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	     ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorUp  *histo_WZ    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	     ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorDown*histo_WZ    ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	     ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorUp  *histo_ZZ    ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	     ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorDown*histo_ZZ    ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]	     ->Add(histo_EM   ); histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_EM   ->GetBinContent(i)+factorUp  *histo_EM    ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]       ->Add(histo_EM   ); histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_EM   ->GetBinContent(i)+factorDown*histo_EM    ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[i-1]   ->Add(histo_Wjets); histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[i-1]   ->SetBinContent(i,TMath::Max(histo_Wjets->GetBinContent(i)+factorUp  *histo_Wjets ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[i-1] ->Add(histo_Wjets); histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[i-1] ->SetBinContent(i,TMath::Max(histo_Wjets->GetBinContent(i)+factorDown*histo_Wjets ->GetBinError(i),0.000001));
  }
  double mean,up,diff;

  if(showSignalOnly == false) {
    printf("nuisance WZ | ZZ | Wj: %f/%f | %f/%f | %f/%f\n",histo_WZ->GetSumOfWeights(),histo_WZ_CMS_WZNLOBoundingUp->GetSumOfWeights(),
                                                            histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ZZNLOBoundingUp->GetSumOfWeights(),
						 	    histo_Wjets->GetSumOfWeights(),histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());
  }
  histo_Wjets_CMS_MVAWBoundingUp->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());
  histo_WZ_CMS_WZNLOBoundingUp  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingUp  ->GetSumOfWeights());
  histo_ZZ_CMS_ZZNLOBoundingUp  ->Scale(histo_ZZ   ->GetSumOfWeights()/histo_ZZ_CMS_ZZNLOBoundingUp  ->GetSumOfWeights());

  double oldNormGJ       = histo_GAMMA_JETS      ->GetSumOfWeights();
  double oldNormGJ_LOOSE = histo_GAMMA_JETS_LOOSE->GetSumOfWeights();
  for(int i=1; i<=histo_GAMMA_JETS->GetNbinsX(); i++){
    if(histo_GAMMA_JETS      ->GetBinContent(i) < 0) histo_GAMMA_JETS      ->SetBinContent(i,0.0);
    if(histo_GAMMA_JETS_LOOSE->GetBinContent(i) < 0) histo_GAMMA_JETS_LOOSE->SetBinContent(i,0.0);
  }
  printf("renormalize g+jets: %f vs. %f\n",oldNormGJ,histo_GAMMA_JETS->GetSumOfWeights());
  if(histo_GAMMA_JETS->GetSumOfWeights() > 0) histo_GAMMA_JETS->Scale(oldNormGJ/histo_GAMMA_JETS->GetSumOfWeights());
  printf("renormalize g+jets_LOOSE: %f vs. %f\n",oldNormGJ_LOOSE,histo_GAMMA_JETS_LOOSE->GetSumOfWeights());
  if(histo_GAMMA_JETS_LOOSE->GetSumOfWeights() > 0) histo_GAMMA_JETS_LOOSE->Scale(oldNormGJ_LOOSE/histo_GAMMA_JETS_LOOSE->GetSumOfWeights());

  printf("===> g+jets | z+jets (hist) | z+jets (estimation): %f | %f | %f\n",histo_GAMMA_JETS->GetSumOfWeights(),histo_Zjets->GetSumOfWeights(),NFinal[0]);
  if     (histo_GAMMA_JETS->GetSumOfWeights() > 0 && NFinal[0] > 0) {
    histo_GAMMA_JETS->Scale(NFinal[0]/histo_GAMMA_JETS->GetSumOfWeights());
  }
  else if(histo_Zjets->GetSumOfWeights() > 0 && NFinal[0] > 0) {
    histo_Zjets->Scale(NFinal[0]/histo_Zjets->GetSumOfWeights());
  }
  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){
    //histo_Zjets->SetBinContent(i,nZjets/histo_ZH_hinv->GetNbinsX());
    //histo_Zjets->SetBinError  (i,nZjets/histo_ZH_hinv->GetNbinsX());
    if(histo_GAMMA_JETS->GetSumOfWeights() > 0 && NFinal[0] > 0){
      histo_Zjets->SetBinContent(i, histo_GAMMA_JETS->GetBinContent(i));
      histo_Zjets->SetBinError  (i, histo_GAMMA_JETS->GetBinError(i));
    }
    else if(histo_Zjets->GetSumOfWeights() == 0 && histo_GAMMA_JETS->GetSumOfWeights() == 0 && NFinal[0] > 0 && i == 1){
      histo_Zjets->SetBinContent(i, NFinal[0]);
      histo_Zjets->SetBinError  (i, NFinal[0]);
    }

    mean = histo_ZH_hinv                        ->GetBinContent(i);
    up   = histo_ZH_hinv_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZH_hinv_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZH_hinv_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_ZZ			   ->GetBinContent(i);
    up   = histo_ZZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WZ		       ->GetBinContent(i);
    up   = histo_WZ_CMS_WZNLOBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_ZZ		       ->GetBinContent(i);
    up   = histo_ZZ_CMS_ZZNLOBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_CMS_ZZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_CMS_ZZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_Wjets 		         ->GetBinContent(i);
    up   = histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  }
  histo_Wjets_CMS_MVAWBoundingDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingDown->GetSumOfWeights());
  histo_WZ_CMS_WZNLOBoundingDown  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingDown  ->GetSumOfWeights());
  histo_ZZ_CMS_ZZNLOBoundingDown  ->Scale(histo_ZZ   ->GetSumOfWeights()/histo_ZZ_CMS_ZZNLOBoundingDown  ->GetSumOfWeights());

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"zllhinv%2s_%3d_input_%4s.root",finalStateName,mH,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data  ->Write();
  histo_ZH_hinv->Write();
  histo_Zjets ->Write();
  histo_VVV   ->Write();
  histo_WZ    ->Write();
  histo_ZZ    ->Write();
  histo_EM    ->Write();
  histo_Wjets ->Write();

  cout << histo_Data  ->GetSumOfWeights() << " ";
  cout << histo_ZH_hinv->GetSumOfWeights() << " ";
  cout << histo_Zjets ->GetSumOfWeights() << " ";
  cout << histo_VVV   ->GetSumOfWeights() << " ";
  cout << histo_WZ    ->GetSumOfWeights() << " ";
  cout << histo_ZZ    ->GetSumOfWeights() << " ";
  cout << histo_EM    ->GetSumOfWeights() << " ";
  cout << histo_Wjets ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingUp	     ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingDown      ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingUp	     ->GetBinContent(i)/histo_EM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingDown      ->GetBinContent(i)/histo_EM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWjetsStatBoundingUp  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWjetsStatBoundingUp  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWjetsStatBoundingDown->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWjetsStatBoundingDown->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_ZH_hinv_CMS_MVALepEffBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffBoundingUp	     ->GetBinContent(i)/histo_ZH_hinv	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepEffBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffBoundingUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_ZH_hinv_CMS_MVALepResBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepResBoundingUp	     ->GetBinContent(i)/histo_ZH_hinv	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepResBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepResBoundingUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_ZH_hinv_CMS_MVAMETResBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETResBoundingUp	     ->GetBinContent(i)/histo_ZH_hinv	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAMETResBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETResBoundingUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_ZH_hinv_CMS_MVAJESBoundingUp      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAJESBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingDown	     ->GetBinContent(i)/histo_ZH_hinv	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp	  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown	  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties GEN\n");
  histo_Wjets_CMS_MVAWBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWBoundingUp	  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWBoundingDown	  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLOBoundingUp            ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLOBoundingUp	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLOBoundingDown          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLOBoundingDown	  ->GetBinContent(i)/histo_WZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_ZZNLOBoundingUp            ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ZZNLOBoundingUp	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_ZZNLOBoundingDown          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ZZNLOBoundingDown	  ->GetBinContent(i)/histo_ZZ	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=0; nb<nBinMVA; nb++){
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]    ->Write();
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]  ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    ->Write();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	    ->Write();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]       ->Write();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]   ->Write();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb] ->Write();
  }

  char outputLimitsShape[200];
  sprintf(outputLimitsShape,"histo_limits_zllhinv%2s_mh%3d_shape_%4s.txt",finalStateName,mH,ECMsb.Data());
  ofstream newcardShape;
  newcardShape.open(outputLimitsShape);
  newcardShape << Form("imax 1 number of channels\n");
  newcardShape << Form("jmax * number of background\n");
  newcardShape << Form("kmax * number of nuisance parameters\n");
  newcardShape << Form("Observation %d\n",(int)nSelectedData[SIGSEL+nSelTypes]);
  newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
  newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
  newcardShape << Form("bin hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s\n",finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data());
  newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM Wjets\n");
  newcardShape << Form("process 0 1 2 3 4 5 6\n");
  newcardShape << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",NFinal[nBkg],NFinal[0],NFinal[1],NFinal[2],NFinal[3],NFinal[4],TMath::Max(NFinal[5],0.0));
  newcardShape << Form("lumi_%4s                               lnN  %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE); 		     
  newcardShape << Form("%s                                   shape  1.000   -   1.000 1.000 1.000   -	  -  \n",effName);
  newcardShape << Form("%s                                   shape  1.000   -   1.000 1.000 1.000   -	  -  \n",momName);
  newcardShape << Form("CMS_scale_met                        shape  1.000   -   1.000 1.000 1.000   -	  -  \n");
  newcardShape << Form("CMS_scale_j                          shape  1.000   -   1.000 1.000 1.000   -	  -  \n");			 
  newcardShape << Form("UEPS			               lnN  1.030   -     -     -     -     -     -  \n");
  newcardShape << Form("CMS_eff_b                              lnN  %5.3f   -   %5.3f %5.3f %5.3f   -     -  \n",syst_btag,syst_btag,syst_btag,syst_btag);
  newcardShape << Form("pdf_qqbar                              lnN  %5.3f   -     -   %5.3f %5.3f   -     -  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2]);
  newcardShape << Form("QCDscale_VH		               lnN  %5.3f   -     -     -     -     -     -  \n",QCDscale_VH);	
  newcardShape << Form("QCDscale_VV		               lnN    -     -     -   1.107 1.065   -     -  \n");		
  newcardShape << Form("CMS_zllhinv_WZ3l                       lnN    -     -     -   %5.3f   -     -     -  \n",syst_WZ3l);		
  if(NFinal[1] > 0)
  newcardShape << Form("QCDscale_VVV		               lnN    -     -   1.500   -     -     -     -  \n");		
  if(NFinal[5] > 0)
  newcardShape << Form("CMS_FakeRate                           lnN    -     -     -     -     -	    -   %5.3f\n",WjetsSyst);  
  newcardShape << Form("CMS_zllhinv_ZLL_%4s                    lnN    -   %5.3f   -	-     -     -	  -  \n",ECMsb.Data(),NFinalE[0]);	      
  if(NFinal[4] > 0)
  newcardShape << Form("CMS_zllhinv_EM_%4s                     lnN    -     -	  -	-     -   %5.3f   -  \n",ECMsb.Data(),NFinalE[4]);	      
  if(NFinal[5] > 0)
  newcardShape << Form("CMS_zllhinv_MVAWBounding             shape    -     -     -	-     -     -   1.000\n");
  //newcardShape << Form("CMS_zllhinv_WZNLOBounding              shape    -     -	  -   1.000   -     -	  -  \n");
  //newcardShape << Form("CMS_zllhinv_ZZNLOBounding              shape    -     -	  -	-   1.000   -	  -  \n");
  if(useFullStatTemplates == false){
    if(NFinal[nBkg] > 0) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s     shape  1.000   -	 -     -     -     -	 -  \n",finalStateName,ECMsb.Data());
    if(NFinal[1]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s    shape    -     -   1.000   -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[2]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s     shape    -     -     -   1.000   -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[3]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s     shape    -     -     -     -   1.000   -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[4]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding_%s     shape    -     -     -     -	  -   1.000   -  \n",finalStateName,ECMsb.Data());
    if(NFinal[5]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAWjetsStatBounding_%s  shape    -     -     -     -	  -	-   1.000\n",finalStateName,ECMsb.Data());
  } else {
    for(int nb=1; nb<=nBinMVA; nb++){
      if(NFinal[nBkg] > 0 && histo_ZH_hinv->GetBinContent(nb) > 0 && histo_ZH_hinv->GetBinError(nb) / histo_ZH_hinv->GetBinContent(nb) > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%d	  shape  1.000   -	 -	-     -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[1]    > 0 && histo_VVV->GetBinContent(nb)     > 0 && histo_VVV->GetBinError(nb)     / histo_VVV->GetBinContent(nb)     > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%d    shape    -	 -     1.000    -     -	    -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[2]    > 0 && histo_WZ->GetBinContent(nb)      > 0 && histo_WZ->GetBinError(nb)      / histo_WZ->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%d	  shape    -	 -	 -    1.000   -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[3]    > 0 && histo_ZZ->GetBinContent(nb)      > 0 && histo_ZZ->GetBinError(nb)      / histo_ZZ->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%d	  shape    -	 -	 -      -   1.000   -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[4]    > 0 && histo_EM->GetBinContent(nb)      > 0 && histo_EM->GetBinError(nb)      / histo_EM->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%d	  shape    -	 -	 -      -     -   1.000   -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[5]    > 0 && histo_Wjets->GetBinContent(nb)   > 0 && histo_Wjets->GetBinError(nb)   / histo_Wjets->GetBinContent(nb)   > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAWjetsStatBounding_%s_Bin%d  shape    -     -       -      -     -     -   1.000\n",finalStateName,ECMsb.Data(),nb-1);
    }
  }
  newcardShape.close();
  }

  return;
}

void metChange(double met, double metPhi, double nvtx, int year, bool isData, double metNew[2]){

double metx = met * cos(metPhi);
double mety = met * sin(metPhi);
if(year == 999 && isData == true){metx = nvtx; assert(0);} // crazy line to avoid warnings
//if(year == 2012) {
//  if(isData == true){
//    metx += (+0.035912 + 0.326549*nvtx);
//    mety += (-0.120266 - 0.225262*nvtx);
//  }
//  else {
//    metx += (+0.108857 + 0.064336*nvtx);
//    mety += (+0.086738 - 0.069474*nvtx);
//  }
//}

metNew[0] = sqrt(metx*metx+mety*mety);
metNew[1] = atan2(mety,metx);

}

double collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2, const int q1){

TLorentzVector LVlep1(lepton1.Px(), lepton1.Py(), lepton1.Pz(), lepton1.E());
TLorentzVector LVlep2(lepton2.Px(), lepton2.Py(), lepton2.Pz(), lepton2.E());

if(q1 > 0){
  LVlep1.SetPxPyPzE(lepton2.Px(), lepton2.Py(), lepton2.Pz(), lepton2.E());
  LVlep2.SetPxPyPzE(lepton1.Px(), lepton1.Py(), lepton1.Pz(), lepton1.E());
}
TLorentzVector LVZ = LVlep1+LVlep2;

//do transformation to Collins-Soper frame (1 rotation, 2 boosts)

// 1st transormation - rotate to ptZ direction 
double zrot = -LVZ.Phi();
LVlep1.RotateZ(zrot);
LVlep2.RotateZ(zrot);

// 2nd transformation - boost in z direction
double beta_boostz = -LVZ.Pz()/LVZ.E();
LVlep1.Boost(0.,0.,beta_boostz);
LVlep2.Boost(0.,0.,beta_boostz);

// 3rd transformation: boost in transverse direction (x-prime)
double beta_boostx = -(LVlep1.Px()+LVlep2.Px())/(LVlep1.E()+LVlep2.E());
LVlep1.Boost(beta_boostx,0.,0.);
LVlep2.Boost(beta_boostx,0.,0.);

// compute cos(theta*) in Colin-Soper frame 
double cos_theta_CS = LVlep1.CosTheta();
if (LVZ.Pz() < 0) {
       cos_theta_CS *= -1.;
}

if     (cos_theta_CS >= +1) cos_theta_CS = +0.999999999;
else if(cos_theta_CS <= -1) cos_theta_CS = -0.999999999;

return cos_theta_CS;

}
