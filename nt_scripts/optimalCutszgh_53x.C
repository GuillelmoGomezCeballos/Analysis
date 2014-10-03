#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Analysis/nt_scripts/makeSystematicEffects.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"

void metChange(double met, double metPhi, double metNew[2], LorentzVector gamma);

const int verboseLevel =   1;
bool UseDyttDataDriven = false; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 10;
const unsigned int nSelTypesSyst = 7;
const bool showSignalOnly = false;

enum selType {ZLGBSEL=0, ZLGESEL, ZLGFSEL, ZLLGSEL, ZLLGFSEL, WWSEL, BTAGSEL, SIGSEL, SIGFSEL, WWLOOSESEL};
TString selTypeName[nSelTypes*2] = {"ZLGBSEL-EM", "ZLGESEL-EM", "ZLGFSEL-EM", "ZLLGSEL-EM", "ZLLGFSEL-EM", "WWSEL-EM", "BTAGSEL-EM", "SIGSEL-EM", "SIGFSEL-EM", "WWLOOSESEL-EM",
                                    "ZLGBSEL-LL", "ZLGESEL-LL", "ZLGFSEL-LL", "ZLLGSEL-LL", "ZLLGFSEL-LL", "WWSEL-LL", "BTAGSEL-LL", "SIGSEL-LL", "SIGFSEL-LL", "WWLOOSESEL-LL"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-EM", "JESDOWN-EM", "LEPP-EM", "LEPM-EM", "MET-EM", "EFFP-EM", "EFFM-EM",
                                            "JESUP-LL", "JESDOWN-LL", "LEPP-LL", "LEPM-LL", "MET-LL", "EFFP-LL", "EFFM-LL"};

// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void optimalCutszgh_53x
(
 int     mH  	 = 95,
 int thePlot = 29,
 TString bgdInputFile    = "ntuples_53x/lgammachi95.root",
 TString dataInputFile   = "ntuples_53x/data_llg.root",
 int period = 3,
 int lSel = 2, 
 int nJetsType = 0
 ,double var0 = 60, double var1 = 2.7, double var2 = 0.50, double var3 = 2.25
 )
{
  double lumi = 1.0;
  
  //                    MET,   DPhiZMET, PTFrac, dPhi
  double cutValue[4] = {60.0,  2.7,     0.50,    2.25};
  cutValue[0] = var0; cutValue[1] = var1; cutValue[2] = var2; cutValue[3] = var3;
  double ptJetMin = 30.0;

  double useFullStatTemplates = true;
  bool useWeightEWKCorr       = true;

  // lepton to photon scale factor (barrel/endcap)
  double sfLepPho[2] = {1.50,1.17};
  // jet to photon scale factor (barrel/endcap)
  double sfJetPho[2] = {1.50,2.44};

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

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
    UseDyttDataDriven = false;
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

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 400.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 18;  xminPlot =   0.0; xmaxPlot = 900.0;}
  else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 100; xminPlot = 40.0; xmaxPlot = 140.0;}
  else if(thePlot >= 12 && thePlot <= 13) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 14 && thePlot <= 14) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 15 && thePlot <= 15) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot == 16 || thePlot == 27) {nBinPlot = 50; xminPlot =  0.0; xmaxPlot = 2.5;}
  else if(thePlot >= 17 && thePlot <= 25) {nBinPlot = 18; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot ==  7 || thePlot == 26) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
  else if(thePlot >= 28 && thePlot <= 28) {nBinPlot = 60; xminPlot = 0.0; xmaxPlot = 6.0;}

  const int nBinMVA = 2;
  Float_t xbins[nBinMVA+1] = {0, 1.479, 2.5};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D* histos;
  if(thePlot != 29) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
  else              histos = new TH1D("histos", "histos", nBinMVA, xbins);  
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

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_ZH_hinv= (TH1D*) histoMVA->Clone("histo_ZH_hinv");
  TH1D *histo_Zjets  = (TH1D*) histoMVA->Clone("histo_Zjets");
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EM     = (TH1D*) histoMVA->Clone("histo_EM");

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
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp     = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown   = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingUp           = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingUp  ->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingDown         = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinUp[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]    = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]  = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb] ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb] = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	    = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	    = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	  ->Sumw2();
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
    nSigCutSyst[i] = 0.0;
    for(int j=0; j<45; j++) {
      bgdDecaySyst[i][j] = 0.0; weiDecaySyst[i][j] = 0.0; 
    }
  }

  unsigned int patternTopVeto = SmurfTree::TopVeto;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int nRealLeptons = 0; int nZLeptons = 0;
    if(TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) nRealLeptons++;
    if(TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13) nRealLeptons++;
    if(TMath::Abs(bgdEvent.lep3McId_) == 11 || TMath::Abs(bgdEvent.lep3McId_) == 13) nRealLeptons++;
    LorentzVector lep1(0,0,0,0), lep2(0,0,0,0), gamma(0,0,0,0), gammaf(0,0,0,0), dilep(0,0,0,0), leppho(0,0,0,0), llpho(0,0,0,0);
    int charge = 0; int lType = 0; int lid1_ = 0; int lid2_ = 0; int lidPho = 0;
    if     ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) {
      charge += (int)bgdEvent.lq1_;
      if(bgdEvent.lep1MotherMcId_ == 23) nZLeptons++;
      if(lep1.pt() == 0) {lep1 = bgdEvent.lep1_; lid1_ = abs(bgdEvent.lid1_);}
      else               {lep2 = bgdEvent.lep1_; lid2_ = abs(bgdEvent.lid1_);}
      if     (abs(bgdEvent.lid1_) == 13) lType += 1;
      else if(abs(bgdEvent.lid1_) == 11) lType += 10;
    }
    else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2 && TMath::Abs(bgdEvent.lep1_.eta()) < 2.4) {gamma  = bgdEvent.lep1_; lidPho = TMath::Abs(bgdEvent.lep1McId_);}
    else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV1)  == SmurfTree::Lep1LooseMuV1  && TMath::Abs(bgdEvent.lep1_.eta()) < 2.4) {gammaf = bgdEvent.lep1_; lidPho = TMath::Abs(bgdEvent.lep1McId_);}

    if     ((bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) {
      charge += (int)bgdEvent.lq2_;
      if(bgdEvent.lep2MotherMcId_ == 23) nZLeptons++;
      if(lep1.pt() == 0) {lep1 = bgdEvent.lep2_; lid1_ = abs(bgdEvent.lid2_);}
      else               {lep2 = bgdEvent.lep2_; lid2_ = abs(bgdEvent.lid2_);}
      if     (abs(bgdEvent.lid2_) == 13) lType += 1;
      else if(abs(bgdEvent.lid2_) == 11) lType += 10;
    }
    else if((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2 && TMath::Abs(bgdEvent.lep2_.eta()) < 2.4) {gamma  = bgdEvent.lep2_; lidPho = TMath::Abs(bgdEvent.lep2McId_);}
    else if((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV1)  == SmurfTree::Lep2LooseMuV1  && TMath::Abs(bgdEvent.lep2_.eta()) < 2.4) {gammaf = bgdEvent.lep2_; lidPho = TMath::Abs(bgdEvent.lep2McId_);}

    if     ((bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection) {
      charge += (int)bgdEvent.lq3_;
      if(bgdEvent.lep3MotherMcId_ == 23) nZLeptons++;
      if(lep1.pt() == 0) {lep1 = bgdEvent.lep3_; lid1_ = abs(bgdEvent.lid3_);}
      else               {lep2 = bgdEvent.lep3_; lid2_ = abs(bgdEvent.lid3_);}
      if     (abs(bgdEvent.lid3_) == 13) lType += 1;
      else if(abs(bgdEvent.lid3_) == 11) lType += 10;
    }
    else if((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2 && TMath::Abs(bgdEvent.lep3_.eta()) < 2.4) {gamma  = bgdEvent.lep3_; lidPho = TMath::Abs(bgdEvent.lep3McId_);}
    else if((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV1)  == SmurfTree::Lep3LooseMuV1  && TMath::Abs(bgdEvent.lep3_.eta()) < 2.4) {gammaf = bgdEvent.lep3_; lidPho = TMath::Abs(bgdEvent.lep3McId_);}

    dilep = lep1+lep2;
    leppho = lep1+gamma;
    llpho = lep1+lep2+gamma;
    if(gammaf.pt() > 0) leppho = lep1+gammaf;

    if(gamma.pt() > 0 && gammaf.pt() > 0) assert(0);

    //bool isRealLepton = false;
    //if(nRealLeptons >= 2) isRealLepton = true;
    if(nRealLeptons >  3) assert(0);

    if     (             lType == 11) lType = 0;
    else if(lSel == 0 && lType == 2 ) lType = 1;
    else if(lSel == 1 && lType == 20) lType = 1;
    else if(lSel == 2 && (lType == 2 || lType == 20)) lType = 1;
    else if(              lType == 1 ) lType = 0;
    else if(              lType == 10) lType = 1;
    else lType = 2;

    if(lType == 2) continue;

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
            bgdEvent.processId_==122)   fDecay = 42;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 42;
    else if(bgdEvent.processId_==10001) fDecay = 42;
    else if(bgdEvent.processId_==10010) fDecay = 42;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

    if(nZLeptons >= 2) {
      if     (fDecay == 21) fDecay = 11;
      else if(fDecay == 27) fDecay = 17;
      else if(fDecay == 28) fDecay = 18;
      if(nRealLeptons == 3) {
        if     (fDecay == 11) fDecay = 31;
        else if(fDecay == 17) fDecay = 37;
        else if(fDecay == 18) fDecay = 38;
      }
    }
    if((fDecay == 9 || fDecay == 19) && (lidPho == 11 || lidPho == 13)) {
      fDecay = 39;
    }

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false, false, false, false, false, false, false, false, false,},
                                   {false, false, false, false, false, false, false, false, false }};

    double MT = sqrt(2.0*dilep.pt()*bgdEvent.met_*(1.0-cos(DeltaPhi(dilep.phi() ,bgdEvent.metPhi_))));

    // 0      1      2       3     4   5      6        7           8  9            10            11     12  13    14
    // lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtZ,mlljj,mjj;
    double outputVarLepP[15];
    makeSystematicEffects(lid1_, lid2_, lep1, lep2, dilep, 
                          MT, bgdEvent.met_, bgdEvent.metPhi_, 
                          bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			  bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			  year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(lid1_, lid2_, lep1, lep2, dilep, 
                          MT, bgdEvent.met_, bgdEvent.metPhi_, 
                          bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			  bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			  year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(lid1_, lid2_, lep1, lep2, dilep, 
                          MT, bgdEvent.met_, bgdEvent.metPhi_, 
                          bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			  bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			  year, 2, outputVarMET);

    double metNew[2]; metChange(bgdEvent.met_,bgdEvent.metPhi_,metNew,gamma);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 

    double metNewSyst[2]; metChange(outputVarMET[4],outputVarMET[5],metNewSyst,gamma);

    bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;

    bool passZMass         = fabs(dilep.mass()-91.1876) < 15.;
    bool passZMassSystP    = fabs(outputVarLepP[2]-91.1876) < 15.;
    bool passZMassSystM    = fabs(outputVarLepM[2]-91.1876) < 15.;

    bool passZMassLarge    = fabs(dilep.mass()-91.1876) < 30.;
    bool passZMassSB       = (dilep.mass() > 110.0 && dilep.mass() < 200.0);

    bool passMET           = bgdEvent.met_ > cutValue[0];
    bool passMETSystMET    = outputVarMET[4] > cutValue[0];

    bool passDPhiZMET        = DeltaPhi(dilep.phi() ,theMETPHI) > cutValue[1];
    bool passDPhiZMETSystMET = DeltaPhi(dilep.phi() ,metNewSyst[1]) > cutValue[1];

    bool passPTFrac        = fabs(theMET-dilep.pt())/dilep.pt() < cutValue[2];
    bool passPTFracSystP   = fabs(theMET-outputVarLepP[3])/outputVarLepP[3] < cutValue[2];
    bool passPTFracSystM   = fabs(theMET-outputVarLepM[3])/outputVarLepM[3] < cutValue[2];
    bool passPTFracSystMET = fabs(metNewSyst[0]-dilep.pt())/dilep.pt() < cutValue[2];

    bool passDPhiLL        = DeltaPhi(lep1.phi() ,lep2.phi()) < cutValue[3];

    bool passPTLL          = dilep.pt() > 60.;
    bool passPTLLSystP     = outputVarLepP[3] > 60.;
    bool passPTLLSystM     = outputVarLepM[3] > 60.;

    bool passLLG           = charge == 0 && lep1.pt() > 20.        && lep2.pt() > 20.        && gamma.pt() > 20;
    bool passLLGSystP      = charge == 0 && outputVarLepP[0] > 20. && outputVarLepP[1] > 20. && gamma.pt() > 20;
    bool passLLGSystM	   = charge == 0 && outputVarLepM[0] > 20. && outputVarLepM[1] > 20. && gamma.pt() > 20;

    bool passLLGF          = charge == 0 && lep1.pt() > 20. && lep2.pt() > 20. && gammaf.pt() > 20;
    bool passLG            = lep1.pt() > 30. && lep2.pt() <= 0. && gamma.pt() > 20;
    bool passLGF           = lep1.pt() > 30. && lep2.pt() <= 0. && gammaf.pt() > 20;

    // jet1McId and jet2McId actually provide lepton rejection information
    // hasZCand = reject events with |mll-mZ|<15
    // trackSel[0] == reject events with |mll-mZ|<10 for pf candidates with |eta|<3.0
    // trackSel[1] == reject events with |mll-mZ|<10 for pf candidates with |eta|<4.7
    // trackSel[2] == reject events with isolated reconstructed leptons with pt>10 and iso/pt<0.1
    // trackSel[3] == reject events with isolated tracks with pt>10 and iso/pt<0.1 (not used by default)
    int trackSel[4] = {int((bgdEvent.jet2McId_%100-bgdEvent.jet2McId_%10)/10),int((bgdEvent.jet2McId_%1000-bgdEvent.jet2McId_%100)/100),int((bgdEvent.jet2McId_%10000-bgdEvent.jet2McId_%1000)/1000),int(bgdEvent.jet2McId_/10000)};

    double MVAVar = 0.5; if(TMath::Abs(gamma.Eta())>1.479) MVAVar = 2.0;

    double addLepEff = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff =                                leptonEfficiency(lep1.pt(), lep1.eta(), fhDEffMu, fhDEffEl, lid1_, 0);
    if(lep2.pt() > 0) addLepEff  = addLepEff * leptonEfficiency(lep2.pt(), lep2.eta(), fhDEffMu, fhDEffEl, lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp =                                  leptonEfficiency(lep1.pt(), lep1.eta(), fhDEffMu, fhDEffEl, lid1_, 1);
      if(lep2.pt() > 0) addLepEffUp  = addLepEffUp * leptonEfficiency(lep2.pt(), lep2.eta(), fhDEffMu, fhDEffEl, lid2_, 1);
      addLepEffDown =                                    leptonEfficiency(lep1.pt(), lep1.eta(), fhDEffMu, fhDEffEl, lid1_, -1);
      if(lep2.pt() > 0) addLepEffDown  = addLepEffDown * leptonEfficiency(lep2.pt(), lep2.eta(), fhDEffMu, fhDEffEl, lid2_, -1);
    } else {addLepEff = 1.0;}

    int NjetSyst[3] = {0,0,0};
    if(bgdEvent.jet1_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[2]++;
    if(bgdEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[2]++;
    if(bgdEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[2]++;
    if(bgdEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[2]++;
    //NjetSyst[0] = nJetsType;NjetSyst[1] = nJetsType;NjetSyst[2] = nJetsType;

    if(1) {

       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLG  && ((fabs(leppho.mass()-91.1876) < 15. && lid1_ == 11 && lType == 1) || (leppho.mass() > 100 && lid1_ == 13 && lType == 0)) && TMath::Abs(gamma.Eta()) <= 1.479) passCuts[lType][ZLGBSEL] = true;
       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLG  && ((fabs(leppho.mass()-91.1876) < 15. && lid1_ == 11 && lType == 1) || (leppho.mass() > 100 && lid1_ == 13 && lType == 0)) && TMath::Abs(gamma.Eta()) >  1.479) passCuts[lType][ZLGESEL] = true;
       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLGF && ((fabs(leppho.mass()-91.1876) < 45. && lType == 1) || (leppho.mass() > 100 && lType == 0))) passCuts[lType][ZLGFSEL] = true;
       if(passLLGF && passZMass) passCuts[lType][ZLLGFSEL] = true;
       if(NjetSyst[0] >= 1	   && trackSel[2]+trackSel[3] == 0 && passLLG  && !passBtagVeto && !passZMass && passZMassSB								          ) passCuts[lType][WWLOOSESEL]  = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  &&  passBtagVeto && !passZMass && passZMassLarge && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][WWSEL]      = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  && !passBtagVeto &&  passZMass		       && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac ) passCuts[lType][BTAGSEL]     = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLGF &&  passBtagVeto &&  passZMass		       && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac ) passCuts[lType][SIGFSEL]     = true;

       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG       &&  passBtagVeto &&  passZMass             && passMET        && passPTLL      && passDPhiLL && passDPhiZMET        && passPTFrac       ) passCuts[lType][SIGSEL]      = true;
       if(NjetSyst[1] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG       &&  passBtagVeto &&  passZMass	       && passMET        && passPTLL      && passDPhiLL && passDPhiZMET        && passPTFrac       ) passSystCuts[lType][JESUP]   = true;
       if(NjetSyst[2] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG       &&  passBtagVeto &&  passZMass	       && passMET        && passPTLL      && passDPhiLL && passDPhiZMET        && passPTFrac       ) passSystCuts[lType][JESDOWN] = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLGSystP  &&  passBtagVeto &&  passZMassSystP	       && passMET        && passPTLLSystP && passDPhiLL && passDPhiZMET        && passPTFracSystP  ) passSystCuts[lType][LEPP]    = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLGSystM  &&  passBtagVeto &&  passZMassSystM	       && passMET        && passPTLLSystM && passDPhiLL && passDPhiZMET        && passPTFracSystM  ) passSystCuts[lType][LEPM]    = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG       &&  passBtagVeto &&  passZMass             && passMETSystMET && passPTLL      && passDPhiLL && passDPhiZMETSystMET && passPTFracSystMET) passSystCuts[lType][MET]     = true;

       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  &&  fabs(llpho.mass()-91.1876) < 100.) passCuts[lType][ZLLGSEL] = true;

       if(fDecay == 9 && (lidPho == 8 || lidPho == 22)) passCuts[lType][ZLLGSEL] = false;

      //if(isRealLepton == false &&
      //   (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
      //    bgdEvent.dstype_ == SmurfTree::qqww   || bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz	|| bgdEvent.dstype_ == SmurfTree::zz   ||
      //    bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) 
      //  {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false; passCuts[0][ZLLSEL] = false; passCuts[1][ZLLSEL] = false;}
    }

    if(1){
      double trigEff = trigLookup.GetExpectedTriggerEfficiency(lep1.eta(), lep1.pt(), 
	    						       lep2.eta(), lep2.pt(), 
        						       lid1_, lid2_);
      double puWeight = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);

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
 
      if(nFake >= 1){
        theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
        theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,bgdEvent.met_);

        if(fDecay == 31 || fDecay == 37 || fDecay == 38 || fDecay == 39) {
	  if(TMath::Abs(gamma.Eta()) <= 1.479) add = add*sfLepPho[0];
	  else                                 add = add*sfLepPho[1];
	}
        if(fDecay == 11 || fDecay == 17 || fDecay == 18 || fDecay == 9) {
	  if(TMath::Abs(gamma.Eta()) <= 1.479) add = add*sfJetPho[0];
          else                                 add = add*sfJetPho[1];
	}

	theWeight = bgdEvent.scale1fb_*lumi*add*puWeight*trigEff*addLepEff;
      }

      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::wz) theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,1);
      if(useWeightEWKCorr == true && bgdEvent.dstype_ == SmurfTree::zz) theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,2);
      if(useWeightEWKCorr == true && fDecay == 42)                      theWeight = theWeight * weightEWKCorr(bgdEvent.higgsPt_,0);

      if(passCuts[1][SIGSEL]){ // begin making plots
	double myVar = bgdEvent.met_; // var0
	if     (thePlot == 1) myVar = lep1.pt();
	else if(thePlot == 2) myVar = lep2.pt();
        else if(thePlot == 3) myVar = TMath::Max(gamma.pt(),gammaf.pt());
	else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = dilep.mass();
	else if(thePlot == 8) myVar = MT;
	else if(thePlot == 9) myVar = leppho.mass();
	else if(thePlot ==10) myVar = dilep.pt();
	else if(thePlot ==11) myVar = fabs(dilep.mass()-91.1876);
	else if(thePlot ==12) myVar = fabs(theMET-dilep.pt())/dilep.pt(); // var2
	else if(thePlot ==13) myVar = lep2.pt()/lep1.pt();
	else if(thePlot ==14) myVar = bgdEvent.njets_;
	else if(thePlot ==15) myVar = bgdEvent.nvtx_;
	else if(thePlot ==16) myVar = TMath::Min(TMath::Abs(dilep.eta()),2.4999);
	else if(thePlot ==17) myVar = DeltaPhi(lep1.phi() ,lep2.phi())*180.0/TMath::Pi();
	else if(thePlot ==18) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==19) myVar = DeltaPhi(dilep.phi() ,theMETPHI)*180.0/TMath::Pi(); // var1
	else if(thePlot ==20) myVar = DeltaPhi(bgdEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = DeltaPhi(dilep.Phi() ,bgdEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(gamma.phi(),bgdEvent.metPhi_)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(gamma.phi(),theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = DeltaPhi(dilep.Phi() ,bgdEvent.metPhi_)*180.0/TMath::Pi();
	else if(thePlot ==25) myVar = DeltaPhi(dilep.Phi() ,gamma.phi())*180.0/TMath::Pi();
	else if(thePlot ==26) myVar = llpho.mass();
	else if(thePlot ==27) myVar = TMath::Abs(gamma.Eta());
	else if(thePlot ==28) myVar = TMath::Min(DeltaR(gamma.Phi(),gamma.Eta(),lep1.Phi(),lep1.Eta()),DeltaR(gamma.Phi(),gamma.Eta(),lep2.Phi(),lep2.Eta()));
	else if(thePlot ==29) myVar = MVAVar;

      	if     (fDecay == 9 || fDecay == 19){
      	  histo0->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 11 || fDecay == 17 || fDecay == 18){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 31 || fDecay == 37 || fDecay == 38){
      	  histo2->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 39){
      	  histo3->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 21 || fDecay == 27 || fDecay == 28 ||
	        fDecay ==  1 || fDecay ==  3 || fDecay == 23 || 
                fDecay == 29 || fDecay == 30 || fDecay ==  5 ||
	        fDecay == 13 || fDecay == 20 || fDecay == 10){
      	  histo4->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 41 || fDecay == 42 || fDecay == 43){
	  histos->Fill(myVar,theWeight);
      	}
      	else {
      	  printf("NOOOOOOOOOOOOOOOOOOOO\n");
      	}
      } // end making plots
      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i] && fDecay != 42) {
            bgdDecay[i+j*nSelTypes][(int)fDecay] += theWeight;
            weiDecay[i+j*nSelTypes][(int)fDecay] += theWeight*theWeight;
          }
	  else if(passCuts[j][i]) {
            nSigCut[i+j*nSelTypes]  += theWeight;
            nSigECut[i+j*nSelTypes] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if     (passSystCuts[j][i] && fDecay != 42) {
            bgdDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight;
            weiDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight*theWeight;
          }
          else if(passSystCuts[j][i]) {
            nSigCutSyst [i+j*nSelTypesSyst] += theWeight;
            nSigECutSyst[i+j*nSelTypesSyst] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][SIGSEL]) {
        if(fDecay != 42) {
          bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
          weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
          bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
          weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
        } else {
          nSigCutSyst [EFFP+lType*nSelTypesSyst] += theWeight	       *addLepEffUp  /addLepEff;
          nSigECutSyst[EFFP+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
          nSigCutSyst [EFFM+lType*nSelTypesSyst] += theWeight	       *addLepEffDown/addLepEff;
          nSigECutSyst[EFFM+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
	}

        if     (fDecay == 9 || fDecay == 19 || fDecay == 39){
	  if(passCuts[1][SIGSEL])              histo_Zjets                        ->Fill(MVAVar, theWeight);
        }
	else if(fDecay == 11 || fDecay == 31){
	  if(passCuts[1][SIGSEL])	       histo_VVV			  ->Fill(MVAVar, theWeight);
          if(passCuts[1][SIGSEL])              histo_VVV_CMS_MVALepEffBoundingUp  ->Fill(MVAVar, theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])              histo_VVV_CMS_MVALepEffBoundingDown->Fill(MVAVar, theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_VVV_CMS_MVAJESBoundingUp     ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_VVV_CMS_MVAJESBoundingDown   ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_VVV_CMS_MVALepResBoundingUp  ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_VVV_CMS_MVALepResBoundingDown->Fill(MVAVar, theWeight);
          if(passSystCuts[1][MET]     == true) histo_VVV_CMS_MVAMETResBoundingUp  ->Fill(MVAVar, theWeight);;
        }
	else if(fDecay == 17 || fDecay == 37){
	  if(passCuts[1][SIGSEL])	       histo_WZ 			 ->Fill(MVAVar, theWeight);
          if(passCuts[1][SIGSEL])	       histo_WZ_CMS_MVALepEffBoundingUp  ->Fill(MVAVar, theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])	       histo_WZ_CMS_MVALepEffBoundingDown->Fill(MVAVar, theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_WZ_CMS_MVAJESBoundingUp     ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_WZ_CMS_MVAJESBoundingDown   ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_WZ_CMS_MVALepResBoundingUp  ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_WZ_CMS_MVALepResBoundingDown->Fill(MVAVar, theWeight);
          if(passSystCuts[1][MET]     == true) histo_WZ_CMS_MVAMETResBoundingUp  ->Fill(MVAVar, theWeight);;
        }
	else if(fDecay == 18 || fDecay == 38){
	  if(passCuts[1][SIGSEL])	       histo_ZZ 			 ->Fill(MVAVar, theWeight);
          if(passCuts[1][SIGSEL])	       histo_ZZ_CMS_MVALepEffBoundingUp  ->Fill(MVAVar, theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL])	       histo_ZZ_CMS_MVALepEffBoundingDown->Fill(MVAVar, theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_ZZ_CMS_MVAJESBoundingUp     ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_ZZ_CMS_MVAJESBoundingDown   ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_ZZ_CMS_MVALepResBoundingUp  ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_ZZ_CMS_MVALepResBoundingDown->Fill(MVAVar, theWeight);
          if(passSystCuts[1][MET]     == true) histo_ZZ_CMS_MVAMETResBoundingUp  ->Fill(MVAVar, theWeight);;
        }
	else if(fDecay == 21 || fDecay == 27 || fDecay == 28 ||
	        fDecay == 29 || fDecay == 30 || fDecay ==  5 ||
	        fDecay == 13 || fDecay == 20 || fDecay == 10 ||
		fDecay == 1  || fDecay == 23){
	  if(passCuts[1][SIGSEL])              histo_EM                          ->Fill(MVAVar, theWeight);
        }
        else if(fDecay == 42){
          if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv                            ->Fill(MVAVar, theWeight);
          if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv_CMS_MVALepEffBoundingUp    ->Fill(MVAVar, theWeight*addLepEffUp  /addLepEff);
          if(passCuts[1][SIGSEL]) 	     histo_ZH_hinv_CMS_MVALepEffBoundingDown  ->Fill(MVAVar, theWeight*addLepEffDown/addLepEff);
          if(passSystCuts[1][JESUP  ] == true) histo_ZH_hinv_CMS_MVAJESBoundingUp     ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][JESDOWN] == true) histo_ZH_hinv_CMS_MVAJESBoundingDown   ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPP]    == true) histo_ZH_hinv_CMS_MVALepResBoundingUp  ->Fill(MVAVar, theWeight);
          if(passSystCuts[1][LEPM]    == true) histo_ZH_hinv_CMS_MVALepResBoundingDown->Fill(MVAVar, theWeight);
          if(passSystCuts[1][MET]     == true) histo_ZH_hinv_CMS_MVAMETResBoundingUp  ->Fill(MVAVar, theWeight);;
	}
	else {
	  printf("%d\n",fDecay);assert(0);
	}
      }
    } // if passCuts
  } // end background loop
  
  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    if(dataEvent.dstype_ == SmurfTree::data &&
      (dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    LorentzVector lep1(0,0,0,0), lep2(0,0,0,0), gamma(0,0,0,0), gammaf(0,0,0,0), dilep(0,0,0,0), leppho(0,0,0,0), llpho(0,0,0,0);
    int charge = 0; int lType = 0; int lid1_ = 0; int lid2_ = 0;
    if     ((dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) {
      charge += (int)dataEvent.lq1_;
      if(lep1.pt() == 0) {lep1 = dataEvent.lep1_; lid1_ = abs(dataEvent.lid1_);}
      else               {lep2 = dataEvent.lep1_; lid2_ = abs(dataEvent.lid1_);}
      if     (abs(dataEvent.lid1_) == 13) lType += 1;
      else if(abs(dataEvent.lid1_) == 11) lType += 10;
    }
    else if((dataEvent.cuts_ & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2 && TMath::Abs(dataEvent.lep1_.eta()) < 2.4) gamma  = dataEvent.lep1_;
    else if((dataEvent.cuts_ & SmurfTree::Lep1LooseMuV1 ) == SmurfTree::Lep1LooseMuV1  && TMath::Abs(dataEvent.lep1_.eta()) < 2.4) gammaf = dataEvent.lep1_;

    if     ((dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) {
      charge += (int)dataEvent.lq2_;
      if(lep1.pt() == 0) {lep1 = dataEvent.lep2_; lid1_ = abs(dataEvent.lid2_);}
      else               {lep2 = dataEvent.lep2_; lid2_ = abs(dataEvent.lid2_);}
      if     (abs(dataEvent.lid2_) == 13) lType += 1;
      else if(abs(dataEvent.lid2_) == 11) lType += 10;
    }
    else if((dataEvent.cuts_ & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2 && TMath::Abs(dataEvent.lep2_.eta()) < 2.4) gamma  = dataEvent.lep2_;
    else if((dataEvent.cuts_ & SmurfTree::Lep2LooseMuV1)  == SmurfTree::Lep2LooseMuV1  && TMath::Abs(dataEvent.lep2_.eta()) < 2.4) gammaf = dataEvent.lep2_;

    if     ((dataEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection) {
      charge += (int)dataEvent.lq3_;
      if(lep1.pt() == 0) {lep1 = dataEvent.lep3_; lid1_ = abs(dataEvent.lid3_);}
      else               {lep2 = dataEvent.lep3_; lid2_ = abs(dataEvent.lid3_);}
      if     (abs(dataEvent.lid3_) == 13) lType += 1;
      else if(abs(dataEvent.lid3_) == 11) lType += 10;
    }
    else if((dataEvent.cuts_ & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2 && TMath::Abs(dataEvent.lep3_.eta()) < 2.4) gamma  = dataEvent.lep3_;
    else if((dataEvent.cuts_ & SmurfTree::Lep3LooseMuV1)  == SmurfTree::Lep3LooseMuV1  && TMath::Abs(dataEvent.lep3_.eta()) < 2.4) gammaf = dataEvent.lep3_;

    dilep = lep1+lep2;
    leppho = lep1+gamma;
    llpho = lep1+lep2+gamma;
    if(gammaf.pt() > 0) leppho = lep1+gammaf;

    if     (             lType == 11) lType = 0;
    else if(lSel == 0 && lType == 2 ) lType = 1;
    else if(lSel == 1 && lType == 20) lType = 1;
    else if(lSel == 2 && (lType == 2 || lType == 20)) lType = 1;
    else if(              lType == 1 ) lType = 0;
    else if(              lType == 10) lType = 1;
    else lType = 2;

    if(lType == 2) continue;
    if(lid2_ == 9) assert(0); // should never happen

    bool passCuts[2][nSelTypes] = {{false, false, false, false, false, false, false, false, false},
                                   {false, false, false, false, false, false, false, false, false}};

    double MT = sqrt(2.0*dilep.pt()*dataEvent.met_*(1.0-cos(DeltaPhi(dilep.phi() ,dataEvent.metPhi_))));

    double metNew[2]; metChange(dataEvent.met_,dataEvent.metPhi_,metNew,gamma);
    double theMET = metNew[0]; double theMETPHI = metNew[1]; 
    bool passBtagVeto      = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool passZMass         = fabs(dilep.mass()-91.1876) < 15.;
    bool passZMassLarge    = fabs(dilep.mass()-91.1876) < 30.;
    bool passZMassSB       = (dilep.mass() > 110.0 && dilep.mass() < 200.0);
    bool passMET           = dataEvent.met_ > cutValue[0];
    bool passDPhiZMET      = DeltaPhi(dilep.phi() ,theMETPHI) > cutValue[1];
    bool passPTFrac        = fabs(theMET-dilep.pt())/dilep.pt() < cutValue[2];
    bool passDPhiLL        = DeltaPhi(lep1.phi() ,lep2.phi()) < cutValue[3];
    bool passPTLL          = dilep.pt() > 60.;
    bool passLLG           = charge == 0 && lep1.pt() > 20. && lep2.pt() > 20. && gamma.pt()  > 20;
    bool passLLGF          = charge == 0 && lep1.pt() > 20. && lep2.pt() > 20. && gammaf.pt() > 20;
    bool passLG            = lep1.pt() > 30. && lep2.pt() <= 0. && gamma.pt() > 20;
    bool passLGF           = lep1.pt() > 30. && lep2.pt() <= 0. && gammaf.pt() > 20;

    int trackSel[4] = {int((dataEvent.jet2McId_%100-dataEvent.jet2McId_%10)/10),int((dataEvent.jet2McId_%1000-dataEvent.jet2McId_%100)/100),int((dataEvent.jet2McId_%10000-dataEvent.jet2McId_%1000)/1000),int(dataEvent.jet2McId_/10000)};

    int NjetSyst[3] = {0,0,0};
    if(dataEvent.jet1_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(dataEvent.jet2_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(dataEvent.jet3_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    if(dataEvent.jet4_.Pt()*1.00 > ptJetMin) NjetSyst[0]++;
    //NjetSyst[0] = nJetsType;NjetSyst[1] = nJetsType;NjetSyst[2] = nJetsType;

    if(1) {

       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLG  && ((fabs(leppho.mass()-91.1876) < 15. && lid1_ == 11 && lType == 1) || (leppho.mass() > 100 && lid1_ == 13 && lType == 0)) && TMath::Abs(gamma.Eta()) <= 1.479) passCuts[lType][ZLGBSEL] = true;
       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLG  && ((fabs(leppho.mass()-91.1876) < 15. && lid1_ == 11 && lType == 1) || (leppho.mass() > 100 && lid1_ == 13 && lType == 0)) && TMath::Abs(gamma.Eta()) >  1.479) passCuts[lType][ZLGESEL] = true;
       if(NjetSyst[0] == nJetsType &&  passBtagVeto && passLGF && ((fabs(leppho.mass()-91.1876) < 45. && lType == 1) || (leppho.mass() > 100 && lType == 0))) passCuts[lType][ZLGFSEL] = true;
       if(passLLGF && passZMass) passCuts[lType][ZLLGFSEL] = true;
       if(NjetSyst[0] >= 1	   && trackSel[2]+trackSel[3] == 0 && passLLG  && !passBtagVeto && !passZMass && passZMassSB								          ) passCuts[lType][WWLOOSESEL]  = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  &&  passBtagVeto && !passZMass && passZMassLarge && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][WWSEL]       = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  && !passBtagVeto &&  passZMass		        && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][BTAGSEL]     = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  &&  passBtagVeto &&  passZMass		        && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGSEL]      = true;
       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLGF &&  passBtagVeto &&  passZMass		        && passMET && passPTLL && passDPhiLL && passDPhiZMET && passPTFrac) passCuts[lType][SIGFSEL]     = true;

       if(NjetSyst[0] == nJetsType && trackSel[2]+trackSel[3] == 0 && passLLG  &&  fabs(llpho.mass()-91.1876) < 100.) passCuts[lType][ZLLGSEL]      = true;

       // blinded!
       //passCuts[1][SIGSEL] = false;

    }

    double MVAVar = 0.5; if(TMath::Abs(gamma.Eta())>1.479) MVAVar = 2.0;

    if(passCuts[1][SIGSEL]){
      histo_Data->Fill(MVAVar, 1.0);
    }

    for(unsigned int i=0; i<nSelTypes; i++) {
      for(int j=0; j<2; j++){
    	if(passCuts[j][i]) {
    	  nSelectedData[i+j*nSelTypes]  += 1.0;
    	}
      }
    }

    if(passCuts[1][SIGSEL]){ // begin making plots
      double myVar = dataEvent.met_; // var0
      if     (thePlot == 1) myVar = lep1.pt();
      else if(thePlot == 2) myVar = lep2.pt();
      else if(thePlot == 3) myVar = TMath::Max(gamma.pt(),gammaf.pt());
      else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
      else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
      else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
      else if(thePlot == 7) myVar = dilep.mass();
      else if(thePlot == 8) myVar = MT;
      else if(thePlot == 9) myVar = leppho.mass();
      else if(thePlot ==10) myVar = dilep.pt();
      else if(thePlot ==11) myVar = fabs(dilep.mass()-91.1876);
      else if(thePlot ==12) myVar = fabs(theMET-dilep.pt())/dilep.pt(); // var2
      else if(thePlot ==13) myVar = lep2.pt()/lep1.pt();
      else if(thePlot ==14) myVar = dataEvent.njets_;
      else if(thePlot ==15) myVar = dataEvent.nvtx_;
      else if(thePlot ==16) myVar = TMath::Min(TMath::Abs(dilep.eta()),2.4999);
      else if(thePlot ==17) myVar = DeltaPhi(lep1.phi() ,lep2.phi())*180.0/TMath::Pi();
      else if(thePlot ==18) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
      else if(thePlot ==19) myVar = DeltaPhi(dilep.phi() ,theMETPHI)*180.0/TMath::Pi(); // var1
      else if(thePlot ==20) myVar = DeltaPhi(dataEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
      else if(thePlot ==21) myVar = DeltaPhi(dilep.Phi() ,dataEvent.jet1_.Phi())*180.0/TMath::Pi();
      else if(thePlot ==22) myVar = DeltaPhi(gamma.phi(),dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==23) myVar = DeltaPhi(gamma.phi(),theMETPHI)*180.0/TMath::Pi();
      else if(thePlot ==24) myVar = DeltaPhi(dilep.Phi() ,dataEvent.metPhi_)*180.0/TMath::Pi();
      else if(thePlot ==25) myVar = DeltaPhi(dilep.Phi() ,gamma.phi())*180.0/TMath::Pi();
      else if(thePlot ==26) myVar = llpho.mass();
      else if(thePlot ==27) myVar = TMath::Abs(gamma.Eta());
      else if(thePlot ==28) myVar = TMath::Min(DeltaR(gamma.Phi(),gamma.Eta(),lep1.Phi(),lep1.Eta()),DeltaR(gamma.Phi(),gamma.Eta(),lep2.Phi(),lep2.Eta()));
      else if(thePlot ==29) myVar = MVAVar;

      histo5->Fill(myVar,1.0);
    } // end making plots
  } // End loop data

  const unsigned int nBkg = 9;
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
      else if(j == 39)                         {bgdCombined[i][4] += bgdDecay[i][j]; bgdCombinedE[i][4] += weiDecay[i][j];}
      else if(j == 31)			       {bgdCombined[i][5] += bgdDecay[i][j]; bgdCombinedE[i][5] += weiDecay[i][j];}
      else if(j == 37)			       {bgdCombined[i][6] += bgdDecay[i][j]; bgdCombinedE[i][6] += weiDecay[i][j];}
      else if(j == 38)                         {bgdCombined[i][7] += bgdDecay[i][j]; bgdCombinedE[i][7] += weiDecay[i][j];}
      else if(j == 21 || j == 27 || j == 28 ||
	      j == 29 || j == 30 || j ==  5 ||
	      j == 13 || j == 20 || j == 10 ||
	      j == 1  || j == 23)              {bgdCombined[i][8] += bgdDecay[i][j]; bgdCombinedE[i][8] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+
                                   bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]+bgdCombined[i][7]+
				   bgdCombined[i][8]-nTot[i])/nTot[i] > 0.00001) {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]+bgdCombined[i][6]+bgdCombined[i][7]+bgdCombined[i][8]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xxZ-1l) = %11.3f +/- %8.3f\n",bgdCombined[i][0],sqrt(bgdCombinedE[i][0]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(VVV-2l) = %11.3f +/- %8.3f\n",bgdCombined[i][1],sqrt(bgdCombinedE[i][1]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xWZ-2l) = %11.3f +/- %8.3f\n",bgdCombined[i][2],sqrt(bgdCombinedE[i][2]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xZZ-2l) = %11.3f +/- %8.3f\n",bgdCombined[i][3],sqrt(bgdCombinedE[i][3]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xxZ-2l) = %11.3f +/- %8.3f\n",bgdCombined[i][4],sqrt(bgdCombinedE[i][4]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(VVV-3l) = %11.3f +/- %8.3f\n",bgdCombined[i][5],sqrt(bgdCombinedE[i][5]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xWZ-3l) = %11.3f +/- %8.3f\n",bgdCombined[i][6],sqrt(bgdCombinedE[i][6]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xZZ-3l) = %11.3f +/- %8.3f\n",bgdCombined[i][7],sqrt(bgdCombinedE[i][7]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("bgd(xxxxEM) = %11.3f +/- %8.3f\n",bgdCombined[i][8],sqrt(bgdCombinedE[i][8]));
    if(showSignalOnly == false || i%nSelTypes == SIGSEL) printf("*******************************\n");
  }

  //if(showSignalOnly == false) printf("+++++++++++++++++++++++++++++++\n");
  double nTotSyst[nSelTypesSyst*2]; double nETotSyst[nSelTypesSyst*2];
  double bgdCombinedSyst[nSelTypesSyst*2][nBkg],bgdCombinedESyst[nSelTypesSyst*2][nBkg];
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    if(showSignalOnly == false && nSigCutSyst[i] > 0) printf("selectionSyst: %s\n",selTypeNameSyst[i].Data());
    for(unsigned int j=0; j<nBkg; j++) {bgdCombinedSyst[i][j] = 0.0; bgdCombinedESyst[i][j] = 0.0;}
    if(showSignalOnly == false && nSigCutSyst[i] > 0) printf("nSigCutSyst(%2d): %11.3f +/- %8.3f\n",i,nSigCutSyst[i],sqrt(nSigECutSyst[i]));
    nTotSyst[i] = 0.0; nETotSyst[i] = 0.0;
    for(int j=0; j<45; j++){
      if(showSignalOnly == false) if(bgdDecaySyst[i][j] != 0) printf("bdgSyst(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecaySyst[i][j],sqrt(weiDecaySyst[i][j]));
      nTotSyst[i]  += bgdDecaySyst[i][j];
      nETotSyst[i] += weiDecaySyst[i][j];

      if     (j == 9 || j == 19)               {bgdCombinedSyst[i][0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][0] += weiDecaySyst[i][j];}
      else if(j == 11)			       {bgdCombinedSyst[i][1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][1] += weiDecaySyst[i][j];}
      else if(j == 17)			       {bgdCombinedSyst[i][2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][2] += weiDecaySyst[i][j];}
      else if(j == 18)                         {bgdCombinedSyst[i][3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][3] += weiDecaySyst[i][j];}
      else if(j == 39)                         {bgdCombinedSyst[i][4] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][4] += weiDecaySyst[i][j];}
      else if(j == 31)			       {bgdCombinedSyst[i][5] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][5] += weiDecaySyst[i][j];}
      else if(j == 37)			       {bgdCombinedSyst[i][6] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][6] += weiDecaySyst[i][j];}
      else if(j == 38)                         {bgdCombinedSyst[i][7] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][7] += weiDecaySyst[i][j];}
      else if(j == 21 || j == 27 || j == 28 ||
	      j == 29 || j == 30 || j ==  5 ||
	      j == 13 || j == 20 || j == 10 ||
	      j == 1  || j == 23)              {bgdCombinedSyst[i][8] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][8] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false && nTotSyst[i] > 0) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+
                                   bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]+bgdCombinedSyst[i][7]+
				   bgdCombinedSyst[i][8]+-nTotSyst[i])/nTotSyst[i] > 0.00001) {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]+bgdCombinedSyst[i][6]+bgdCombinedSyst[i][7]+bgdCombinedSyst[i][8]);assert(0);}
    //if(showSignalOnly == false) printf("------\n");
    //if(showSignalOnly == false) printf("bgdSyst(xxZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][0],sqrt(bgdCombinedESyst[i][0]));
    //if(showSignalOnly == false) printf("bgdSyst(VVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][1],sqrt(bgdCombinedESyst[i][1]));
    //if(showSignalOnly == false) printf("bgdSyst(xWZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][2],sqrt(bgdCombinedESyst[i][2]));
    //if(showSignalOnly == false) printf("bgdSyst(xZZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][3],sqrt(bgdCombinedESyst[i][3]));
    //if(showSignalOnly == false) printf("bgdSyst(xEM) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][4],sqrt(bgdCombinedESyst[i][4]));
    //if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }
  double NFinal[nBkg+1]  = {bgdCombined[SIGSEL+nSelTypes][0],bgdCombined[SIGSEL+nSelTypes][1],bgdCombined[SIGSEL+nSelTypes][2],bgdCombined[SIGSEL+nSelTypes][3],
                            bgdCombined[SIGSEL+nSelTypes][4],bgdCombined[SIGSEL+nSelTypes][5],bgdCombined[SIGSEL+nSelTypes][6],bgdCombined[SIGSEL+nSelTypes][7],
			    bgdCombined[SIGSEL+nSelTypes][8],nSigCut[SIGSEL+nSelTypes]};
  double NFinalE[nBkg+1] = {1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][0])/bgdCombined[SIGSEL+nSelTypes][0],
                            1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][1])/bgdCombined[SIGSEL+nSelTypes][1],
                     	    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][2])/bgdCombined[SIGSEL+nSelTypes][2],
			    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][3])/bgdCombined[SIGSEL+nSelTypes][3],
			    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][4])/bgdCombined[SIGSEL+nSelTypes][4],
                            1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][5])/bgdCombined[SIGSEL+nSelTypes][5],
                     	    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][6])/bgdCombined[SIGSEL+nSelTypes][6],
			    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][7])/bgdCombined[SIGSEL+nSelTypes][7],
			    1.0+sqrt(bgdCombinedE[SIGSEL+nSelTypes][8])/bgdCombined[SIGSEL+nSelTypes][8],
		            1.0+sqrt(nSigECut[SIGSEL+nSelTypes])/nSigCut[SIGSEL+nSelTypes]};

  // This is the default method
  double NemFact_FromMLLSB = 1.0; double NemFact_FromMLLSBE = 1.0;
  if(nSelectedData[WWLOOSESEL+nSelTypes] > 0 && nSelectedData[WWLOOSESEL] > 0) {
    NemFact_FromMLLSB  = (nSelectedData[WWLOOSESEL+nSelTypes]-bgdCombined[WWLOOSESEL+nSelTypes][0]-bgdCombined[WWLOOSESEL+nSelTypes][1]-bgdCombined[WWLOOSESEL+nSelTypes][2]-bgdCombined[WWLOOSESEL+nSelTypes][3])/
                         (nSelectedData[WWLOOSESEL]          -bgdCombined[WWLOOSESEL][0]          -bgdCombined[WWLOOSESEL][1]          -bgdCombined[WWLOOSESEL][2]          -bgdCombined[WWLOOSESEL][3]);
    NemFact_FromMLLSBE = sqrt(1.0/nSelectedData[WWLOOSESEL+nSelTypes]+
                              1.0/nSelectedData[WWLOOSESEL	    ])*NemFact_FromMLLSB;
  }
  else {NemFact_FromMLLSB = 1.0;}

  printf("NemFact_FromMLLSB = %f +/- %f\n", NemFact_FromMLLSB,NemFact_FromMLLSBE);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  double EMSystTotal = 1.0; double EMSyst[1] = {0.0};
  EMSyst[0] = bgdCombined[SIGSEL+nSelTypes][8]/bgdCombined[SIGSEL][8]/NemFact_FromMLLSB;
  if(EMSyst[0] < 1.0) EMSyst[0] = 1/EMSyst[0]; EMSyst[0] = EMSyst[0] - 1.0;

  // At least one data event is required
  if(nSelectedData[SIGSEL] == 0) nSelectedData[SIGSEL] = 1;

  if(nSelectedData[SIGSEL] > 0) EMSystTotal = 1.0 + sqrt(EMSyst[0]*EMSyst[0] + 1/nSelectedData[SIGSEL]);
  else                          EMSystTotal = 1.0 + sqrt(EMSyst[0]*EMSyst[0] + 1.0);

  double EMbkg = bgdCombined[SIGSEL][0]+bgdCombined[SIGSEL][1]+bgdCombined[SIGSEL][2]+bgdCombined[SIGSEL][3]+
                 bgdCombined[SIGSEL][4]+bgdCombined[SIGSEL][5]+bgdCombined[SIGSEL][6]+bgdCombined[SIGSEL][7];

  printf("EM MC: %8.3f +/- %5.3f --> EM MC Prediction: %8.3f +/- %5.3f, EM data/bkg: %d/%f --> syst: %f (%f,%f)\n",
         bgdCombined[SIGSEL+nSelTypes][8],sqrt(bgdCombinedE[SIGSEL+nSelTypes][8]),
	 bgdCombined[SIGSEL][8]*NemFact_FromMLLSB  ,sqrt(bgdCombinedE[SIGSEL][8])*NemFact_FromMLLSB,
	 (int)nSelectedData[SIGSEL],EMbkg,EMSystTotal-1.0,EMSyst[0],sqrt((EMSystTotal-1.0)*(EMSystTotal-1.0) - EMSyst[0]*EMSyst[0]));

  if(nSelectedData[SIGSEL]-EMbkg > 0) {
    NFinal[8] = (nSelectedData[SIGSEL]-EMbkg)*NemFact_FromMLLSB;
  }
  NFinalE[8] = EMSystTotal;

  double scaleFactorEM = NFinal[8]/bgdCombined[SIGSEL+nSelTypes][8]; if(scaleFactorEM <= 0) scaleFactorEM = 0.0;
  printf("EM scale Factor: %5.3f +/- %5.3f\n",scaleFactorEM,EMSystTotal-1.0);
  printf("EM yield: %5.3f +/- %5.3f\n",NFinal[8],sqrt(nSelectedData[SIGSEL])*NemFact_FromMLLSB);

  //NFinal[0] = bgdCombined[SIGSEL+nSelTypes][0]; NFinalE[0] = 1.5;

  histo_EM->Scale(scaleFactorEM);
  histo4->Scale(scaleFactorEM);

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
    //histo0->Scale(scaleZ);
    histos->Write();
    histo0->Write();
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();
  outFilePlotsNote->Close();

  // Estimation of electron->photon and jet->photon rates (a-posteriori)
  double bck_fakeJetB = bgdCombined[ZLGBSEL][4]+bgdCombined[ZLGBSEL][5]+bgdCombined[ZLGBSEL][6]+bgdCombined[ZLGBSEL][7];
  double sig_fakeJetB = bgdCombined[ZLGBSEL][0]+bgdCombined[ZLGBSEL][1]+bgdCombined[ZLGBSEL][2]+bgdCombined[ZLGBSEL][3]+bgdCombined[ZLGBSEL][8];
  double cor_fakeJetB = (nSelectedData[ZLGBSEL]-bck_fakeJetB)/sig_fakeJetB;

  double bck_fakeJetE = bgdCombined[ZLGESEL][4]+bgdCombined[ZLGESEL][5]+bgdCombined[ZLGESEL][6]+bgdCombined[ZLGESEL][7];
  double sig_fakeJetE = bgdCombined[ZLGESEL][0]+bgdCombined[ZLGESEL][1]+bgdCombined[ZLGESEL][2]+bgdCombined[ZLGESEL][3]+bgdCombined[ZLGESEL][8];
  double cor_fakeJetE = (nSelectedData[ZLGESEL]-bck_fakeJetE)/sig_fakeJetE;

  double bck_fakeEleB = bgdCombined[ZLGBSEL+nSelTypes][0]+bgdCombined[ZLGBSEL+nSelTypes][1]+bgdCombined[ZLGBSEL+nSelTypes][2]+bgdCombined[ZLGBSEL+nSelTypes][3]+bgdCombined[ZLGBSEL+nSelTypes][8];
  double sig_fakeEleB = bgdCombined[ZLGBSEL+nSelTypes][4]+bgdCombined[ZLGBSEL+nSelTypes][5]+bgdCombined[ZLGBSEL+nSelTypes][6]+bgdCombined[ZLGBSEL+nSelTypes][7];
  double cor_fakeEleB = (nSelectedData[ZLGBSEL+nSelTypes]-bck_fakeEleB*cor_fakeJetB)/sig_fakeEleB;

  double bck_fakeEleE = bgdCombined[ZLGESEL+nSelTypes][0]+bgdCombined[ZLGESEL+nSelTypes][1]+bgdCombined[ZLGESEL+nSelTypes][2]+bgdCombined[ZLGESEL+nSelTypes][3]+bgdCombined[ZLGESEL+nSelTypes][8];
  double sig_fakeEleE = bgdCombined[ZLGESEL+nSelTypes][4]+bgdCombined[ZLGESEL+nSelTypes][5]+bgdCombined[ZLGESEL+nSelTypes][6]+bgdCombined[ZLGESEL+nSelTypes][7];
  double cor_fakeEleE = (nSelectedData[ZLGESEL+nSelTypes]-bck_fakeEleE*cor_fakeJetE)/sig_fakeEleE;

  printf("cor_fakeJetB: %f cor_fakeJetE: %f cor_fakeEleB: %f cor_fakeEleE: %f\n",cor_fakeJetB,cor_fakeJetE,cor_fakeEleB,cor_fakeEleE); 

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
    if(nSigCut[SIGSEL+nSelTypes] > 0) systEffect[i][nBkg] = nSigCutSyst[i+nSelTypesSyst]/nSigCut[SIGSEL+nSelTypes];    
    if(systEffect[i][nBkg] < 1) systEffect[i][nBkg] = 1.0/systEffect[i][nBkg];
  }
  //if(showSignalOnly == false) printf("Syst(sig) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][nBkg]-1,systEffect[JESDOWN][nBkg]-1,systEffect[LEPP][nBkg]-1,systEffect[LEPM][nBkg]-1,systEffect[MET][nBkg]-1,systEffect[EFFP][nBkg]-1,systEffect[EFFM][nBkg]-1);
  //if(showSignalOnly == false) printf("Syst(xxZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][0]-1,systEffect[JESDOWN][0]-1,systEffect[LEPP][0]-1,systEffect[LEPM][0]-1,systEffect[MET][0]-1,systEffect[EFFP][0]-1,systEffect[EFFM][0]-1);
  //if(showSignalOnly == false) printf("Syst(VVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][1]-1,systEffect[JESDOWN][1]-1,systEffect[LEPP][1]-1,systEffect[LEPM][1]-1,systEffect[MET][1]-1,systEffect[EFFP][1]-1,systEffect[EFFM][1]-1);
  //if(showSignalOnly == false) printf("Syst(xWZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][2]-1,systEffect[JESDOWN][2]-1,systEffect[LEPP][2]-1,systEffect[LEPM][2]-1,systEffect[MET][2]-1,systEffect[EFFP][2]-1,systEffect[EFFM][2]-1);
  //if(showSignalOnly == false) printf("Syst(xZZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][3]-1,systEffect[JESDOWN][3]-1,systEffect[LEPP][3]-1,systEffect[LEPM][3]-1,systEffect[MET][3]-1,systEffect[EFFP][3]-1,systEffect[EFFM][3]-1);
  //if(showSignalOnly == false) printf("Syst(xEM) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][4]-1,systEffect[JESDOWN][4]-1,systEffect[LEPP][4]-1,systEffect[LEPM][4]-1,systEffect[MET][4]-1,systEffect[EFFP][4]-1,systEffect[EFFM][4]-1);

  double pdf_qqbar[3] = {1.055,1.048,1.057};
  if(nJetsType == 1) {pdf_qqbar[0] = 1.060; pdf_qqbar[0] = 1.043; pdf_qqbar[0] = 1.043;}
  double syst_WZ3l = 1.100;
  double syst_btag = 1.000;
  if     (year == 2011 && nJetsType == 0) syst_btag = 1.002;
  else if(year == 2012 && nJetsType == 0) syst_btag = 1.007;
  else if(year == 2011 && nJetsType == 1) syst_btag = 1.027;
  else if(year == 2012 && nJetsType == 1) syst_btag = 1.018;
  else assert(0);

  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp     ->SetBinContent(i,TMath::Max(histo_ZH_hinv ->GetBinContent(i)+factorUp  *histo_ZH_hinv	     ->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown   ->SetBinContent(i,TMath::Max(histo_ZH_hinv ->GetBinContent(i)+factorDown*histo_ZH_hinv	     ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingUp    ->SetBinContent(i,TMath::Max(histo_Zjets   ->GetBinContent(i)+factorUp  *histo_Zjets   ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingDown  ->SetBinContent(i,TMath::Max(histo_Zjets   ->GetBinContent(i)+factorDown*histo_Zjets   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_VVV     ->GetBinContent(i)+factorUp  *histo_VVV   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_VVV     ->GetBinContent(i)+factorDown*histo_VVV   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_WZ      ->GetBinContent(i)+factorUp  *histo_WZ	 ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_WZ      ->GetBinContent(i)+factorDown*histo_WZ	 ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_ZZ      ->GetBinContent(i)+factorUp  *histo_ZZ	 ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_ZZ      ->GetBinContent(i)+factorDown*histo_ZZ	 ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingUp	      ->SetBinContent(i,TMath::Max(histo_EM      ->GetBinContent(i)+factorUp  *histo_EM	 ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_EM      ->GetBinContent(i)+factorDown*histo_EM	 ->GetBinError(i),0.000001));

    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]    ->Add(histo_ZH_hinv   ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_ZH_hinv   ->GetBinContent(i)+factorUp  *histo_ZH_hinv	 ->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]  ->Add(histo_ZH_hinv   ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_ZH_hinv   ->GetBinContent(i)+factorDown*histo_ZH_hinv	 ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[i-1]   ->Add(histo_Zjets  ); histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_Zjets  ->GetBinContent(i)+factorUp  *histo_Zjets   ->GetBinError(i),0.000001));
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[i-1] ->Add(histo_Zjets  ); histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_Zjets  ->GetBinContent(i)+factorDown*histo_Zjets   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	     ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorUp  *histo_VVV   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]     ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorDown*histo_VVV   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	     ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]            ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorUp  *histo_WZ    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	     ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]          ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorDown*histo_WZ    ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	     ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]            ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorUp  *histo_ZZ    ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	     ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]          ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorDown*histo_ZZ    ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]	     ->Add(histo_EM   ); histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]            ->SetBinContent(i,TMath::Max(histo_EM   ->GetBinContent(i)+factorUp  *histo_EM    ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]       ->Add(histo_EM   ); histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]          ->SetBinContent(i,TMath::Max(histo_EM   ->GetBinContent(i)+factorDown*histo_EM    ->GetBinError(i),0.000001));
  }
  double mean,up,diff;

  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){

    mean = histo_ZH_hinv                        ->GetBinContent(i);
    up   = histo_ZH_hinv_CMS_MVALepResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZH_hinv_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZH_hinv_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_CMS_MVALepResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_CMS_MVALepResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_ZZ			   ->GetBinContent(i);
    up   = histo_ZZ_CMS_MVALepResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_ZZ_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_ZZ_CMS_MVALepResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

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

  }

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"zllhinv%2s_%d_input_%4s.root",finalStateName,mH,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  histo_ZH_hinv->Write();
  histo_Zjets  ->Write();
  histo_VVV    ->Write();
  histo_WZ     ->Write();
  histo_ZZ     ->Write();
  histo_EM     ->Write();

  cout << histo_Data  ->GetSumOfWeights() << " ";
  cout << histo_ZH_hinv->GetSumOfWeights() << " ";
  cout << histo_Zjets ->GetSumOfWeights() << " ";
  cout << histo_VVV   ->GetSumOfWeights() << " ";
  cout << histo_WZ    ->GetSumOfWeights() << " ";
  cout << histo_ZZ    ->GetSumOfWeights() << " ";
  cout << histo_EM    ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv       ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingUp	   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv       ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingDown	->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Zjets  ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAZjetsStatBoundingUp        ->GetBinContent(i)/histo_Zjets  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_MVAZjetsStatBoundingDown->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Zjets  ->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_MVAZjetsStatBoundingDown	->GetBinContent(i)/histo_Zjets  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp	  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown	   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown	   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_EM    ->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingUp	   ->GetBinContent(i)/histo_EM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_EM    ->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingDown	   ->GetBinContent(i)/histo_EM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_ZH_hinv_CMS_MVALepEffBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	 ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffBoundingUp	   ->GetBinContent(i)/histo_ZH_hinv	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepEffBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	 ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffBoundingDown	->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingUp	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingDown	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffBoundingUp	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffBoundingDown	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_ZH_hinv_CMS_MVALepResBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepResBoundingUp	 ->GetBinContent(i)/histo_ZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepResBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepResBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepResBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_ZH_hinv_CMS_MVAMETResBoundingUp   ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETResBoundingUp	 ->GetBinContent(i)/histo_ZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAMETResBoundingDown ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingUp       ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETResBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETResBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_ZH_hinv_CMS_MVAJESBoundingUp      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	 ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ZH_hinv	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAJESBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv 	 ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingDown	   ->GetBinContent(i)/histo_ZH_hinv	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown        ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown	->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=0; nb<nBinMVA; nb++){
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]    ->Write();
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]  ->Write();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]   ->Write();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb] ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    ->Write();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	    ->Write();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]       ->Write();
  }

  char outputLimitsShape[200];
  sprintf(outputLimitsShape,"histo_limits_zllhinv%2s_mh%d_shape_%4s.txt",finalStateName,mH,ECMsb.Data());
  ofstream newcardShape;
  newcardShape.open(outputLimitsShape);
  newcardShape << Form("imax 1 number of channels\n");
  newcardShape << Form("jmax * number of background\n");
  newcardShape << Form("kmax * number of nuisance parameters\n");
  newcardShape << Form("Observation %d\n",(int)nSelectedData[SIGSEL+nSelTypes]);
  newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
  newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
  newcardShape << Form("bin hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s hinv%2s%4s\n",finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data());
  newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM\n");
  newcardShape << Form("process 0 1 2 3 4 5\n");
  newcardShape << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",NFinal[nBkg],NFinal[0]+NFinal[4],NFinal[1]+NFinal[5],NFinal[2]+NFinal[6],NFinal[3]+NFinal[7],NFinal[8]);
  newcardShape << Form("lumi_%4s                               lnN  %5.3f   -   %5.3f %5.3f %5.3f   -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE);		       
  newcardShape << Form("%s                                   shape  1.000   -   1.000 1.000 1.000   -  \n",effName);
  newcardShape << Form("%s                                   shape  1.000   -   1.000 1.000 1.000   -  \n",momName);
  newcardShape << Form("CMS_eff_photon                         lnN  1.030   -   1.030 1.030 1.030   -  \n");
  newcardShape << Form("CMS_scale_met                        shape  1.000   -   1.000 1.000 1.000   -  \n");
  newcardShape << Form("CMS_scale_j                          shape  1.000   -   1.000 1.000 1.000   -  \n");			   
  newcardShape << Form("UEPS			               lnN  1.030   -     -     -     -     -  \n");
  newcardShape << Form("CMS_eff_b                              lnN  %5.3f   -   %5.3f %5.3f %5.3f   -  \n",syst_btag,syst_btag,syst_btag,syst_btag);
  newcardShape << Form("pdf_qqbar                              lnN  %5.3f   -     -   %5.3f %5.3f   -  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2]);
  newcardShape << Form("QCDscale_VH		               lnN  %5.3f   -     -     -     -     -  \n",QCDscale_VH);  
  newcardShape << Form("QCDscale_VV		               lnN    -     -     -   1.107 1.065   -  \n");		  
  newcardShape << Form("CMS_zllghinv_3l                        lnN    -     -   %5.3f %5.3f %5.3f   -  \n",syst_WZ3l,syst_WZ3l,syst_WZ3l);  	  
  if(NFinal[0]+NFinal[4] > 0)
  newcardShape << Form("CMS_zllhinv_ZLL_%4s                    lnN    -   %5.3f   -	-     -     -  \n",ECMsb.Data(),1.5);		
  if(NFinal[1]+NFinal[5] > 0)
  newcardShape << Form("QCDscale_VVV		               lnN    -     -   1.500   -     -     -  \n");		  
  if(NFinal[8] > 0)
  newcardShape << Form("CMS_zllhinv_EM_%4s                     lnN    -     -	  -	-     -   %5.3f\n",ECMsb.Data(),NFinalE[8]);		
  if(useFullStatTemplates == false){
    if(NFinal[nBkg]           > 0) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s     shape  1.000   -     -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[0]+NFinal[4]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding_%s  shape    -   1.000   -     -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[1]+NFinal[5]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s    shape    -     -   1.000   -     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[2]+NFinal[6]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s     shape    -     -     -   1.000   -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[3]+NFinal[7]    > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s     shape    -     -     -     -   1.000   -  \n",finalStateName,ECMsb.Data());
    if(NFinal[8]              > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding_%s     shape    -     -     -     -     -   1.000\n",finalStateName,ECMsb.Data());
  } else {
    for(int nb=1; nb<=nBinMVA; nb++){
      if(NFinal[nBkg]           > 0 && histo_ZH_hinv->GetBinContent(nb) > 0 && histo_ZH_hinv->GetBinError(nb) / histo_ZH_hinv->GetBinContent(nb) > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%d	   shape  1.000   -	-     -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[0]+NFinal[4]    > 0 && histo_Zjets->GetBinContent(nb)   > 0 && histo_Zjets->GetBinError(nb)   / histo_Zjets->GetBinContent(nb)   > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%d shape    -	1.000	-     -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[1]+NFinal[5]    > 0 && histo_VVV->GetBinContent(nb)     > 0 && histo_VVV->GetBinError(nb)     / histo_VVV->GetBinContent(nb)     > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%d   shape    -	  -   1.000   -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[2]+NFinal[6]    > 0 && histo_WZ->GetBinContent(nb)      > 0 && histo_WZ->GetBinError(nb)      / histo_WZ->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%d	   shape    -	  -	-   1.000   -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[3]+NFinal[7]    > 0 && histo_ZZ->GetBinContent(nb)      > 0 && histo_ZZ->GetBinError(nb)      / histo_ZZ->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%d	   shape    -	  -	-     -   1.000   -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[8]              > 0 && histo_EM->GetBinContent(nb)      > 0 && histo_EM->GetBinError(nb)      / histo_EM->GetBinContent(nb)      > 0.05) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%d	   shape    -	  -	-     -     -	1.000\n",finalStateName,ECMsb.Data(),nb-1);
    }
  }
  newcardShape.close();
  }

  return;
}

void metChange(double met, double metPhi, double metNew[2], LorentzVector gamma){

double metx = met * cos(metPhi) + gamma.px();
double mety = met * sin(metPhi) + gamma.py();

metNew[0] = sqrt(metx*metx+mety*mety);
metNew[1] = atan2(mety,metx);

}
