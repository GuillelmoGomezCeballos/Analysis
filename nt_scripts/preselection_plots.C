#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Analysis/HWWlvlv/factors.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"

const int verboseLevel =   1;
const int nDecays      =  25;

// GF  == 10010, WBF == 10001, WH == 24, ZH == 26, ttH=121/122
void preselection_plots
(
 int channel = 3,
 int thePlot = 12,
 unsigned int option  =  0,
 //TString bgdInputFile    = "ntuples_42x_pu/backgroundD.root",
 //TString bgdInputFile    = "/home/ceballos/releases/CMSSW_5_3_14/src/test/histo_ww_new_smurf0_all_noskim.root",
 TString bgdInputFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_noweights/data.root",
 TString signalInputFile = "ntuples_42x_pu/hww160.root",
 TString dataInputFile   = "ntuples_42x_pu/data.root"
 )
{
  bool makeGoodPlots = true;
  double lumi = 1.0;
  bool is2011Ana = false;
  double iniEvents = 1; //19332.35;

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile, -1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  unsigned int patternTopVeto         = SmurfTree::TopVeto;

  int nBin    = 200;
  double xmin = 0.0;
  double xmax = 1.0;

  // MET selection
  TH1D* hDRin   = new TH1D("hDRin"   , "hDRin"   , nBin, xmin, xmax);
  TH1D* hDRou   = new TH1D("hDRou"   , "hDRou"   , nBin, xmin, xmax);
  TH1D* hDRouH  = new TH1D("hDRouH"  , "hDRouH"  , nBin, xmin, xmax);
  TH1D* hDRinmm = new TH1D("hDRinmm" , "hDRinmm" , nBin, xmin, xmax);
  TH1D* hDRoumm = new TH1D("hDRoumm" , "hDRoumm" , nBin, xmin, xmax);
  TH1D* hDRoummH= new TH1D("hDRoummH", "hDRoummH", nBin, xmin, xmax);
  TH1D* hDRinee = new TH1D("hDRinee" , "hDRinee" , nBin, xmin, xmax);
  TH1D* hDRouee = new TH1D("hDRouee" , "hDRouee" , nBin, xmin, xmax);
  TH1D* hDRoueeH= new TH1D("hDRoueeH", "hDRoueeH", nBin, xmin, xmax);
  hDRin   ->Sumw2();
  hDRou   ->Sumw2();
  hDRouH  ->Sumw2();
  hDRinmm ->Sumw2();
  hDRoumm ->Sumw2();
  hDRoummH->Sumw2();
  hDRinee ->Sumw2();
  hDRouee ->Sumw2();
  hDRoueeH->Sumw2();
  TH1D* hDVarin   = new TH1D("hDVarin"   , "hDVarin"   , nBin, xmin, xmax);
  TH1D* hDVarou   = new TH1D("hDVarou"   , "hDVarou"   , nBin, xmin, xmax);
  TH1D* hDVarouH  = new TH1D("hDVarouH"  , "hDVarouH"  , nBin, xmin, xmax);
  TH1D* hDVarinmm = new TH1D("hDVarinmm" , "hDVarinmm" , nBin, xmin, xmax);
  TH1D* hDVaroumm = new TH1D("hDVaroumm" , "hDVaroumm" , nBin, xmin, xmax);
  TH1D* hDVaroummH= new TH1D("hDVaroummH", "hDVaroummH", nBin, xmin, xmax);
  TH1D* hDVarinee = new TH1D("hDVarinee" , "hDVarinee" , nBin, xmin, xmax);
  TH1D* hDVarouee = new TH1D("hDVarouee" , "hDVarouee" , nBin, xmin, xmax);
  TH1D* hDVaroueeH= new TH1D("hDVaroueeH", "hDVaroueeH", nBin, xmin, xmax);
  hDVarin   ->Sumw2();
  hDVarou   ->Sumw2();
  hDVarouH  ->Sumw2();
  hDVarinmm ->Sumw2();
  hDVaroumm ->Sumw2();
  hDVaroummH->Sumw2();
  hDVarinee ->Sumw2();
  hDVarouee ->Sumw2();
  hDVaroueeH->Sumw2();
  double sumEvol[3] = {0., 0., 0.};
  TH1D* hDWWLL0 = new TH1D("hDWWLL0", "hDWWLL0", 12, -0.5, 11.5);
  TH1D* hDWWMM0 = new TH1D("hDWWMM0", "hDWWMM0", 12, -0.5, 11.5);
  TH1D* hDWWEM0 = new TH1D("hDWWEM0", "hDWWEM0", 12, -0.5, 11.5);
  TH1D* hDWWME0 = new TH1D("hDWWME0", "hDWWME0", 12, -0.5, 11.5);
  TH1D* hDWWEE0 = new TH1D("hDWWEE0", "hDWWEE0", 12, -0.5, 11.5);
  TH1D* hDWWLL1 = new TH1D("hDWWLL1", "hDWWLL1", 12, -0.5, 11.5);
  TH1D* hDWWMM1 = new TH1D("hDWWMM1", "hDWWMM1", 12, -0.5, 11.5);
  TH1D* hDWWEM1 = new TH1D("hDWWEM1", "hDWWEM1", 12, -0.5, 11.5);
  TH1D* hDWWME1 = new TH1D("hDWWME1", "hDWWME1", 12, -0.5, 11.5);
  TH1D* hDWWEE1 = new TH1D("hDWWEE1", "hDWWEE1", 12, -0.5, 11.5);
  TH1D* hDWWLL2 = new TH1D("hDWWLL2", "hDWWLL2", 12, -0.5, 11.5);
  TH1D* hDWWMM2 = new TH1D("hDWWMM2", "hDWWMM2", 12, -0.5, 11.5);
  TH1D* hDWWEM2 = new TH1D("hDWWEM2", "hDWWEM2", 12, -0.5, 11.5);
  TH1D* hDWWME2 = new TH1D("hDWWME2", "hDWWME2", 12, -0.5, 11.5);
  TH1D* hDWWEE2 = new TH1D("hDWWEE2", "hDWWEE2", 12, -0.5, 11.5);

  TH1D* hDMetSig   = new TH1D("hDMetSig" , "hDMetSig" , nBin, xmin, xmax);
  TH1D* hDMetBck   = new TH1D("hDMetBck" , "hDMetBck" , nBin, xmin, xmax);

  char sb[200];
  TH1D* hDbtag[10],*hDpttag[10],*hDetatag[10];
  if(channel == 1){ // btag selection
    nBin = 7;
    xmin = -0.5;
    xmax =  6.5;
    for(int n=0; n<10; n++){
      sprintf(sb,"hDbtag_%d",n);    hDbtag[n]    = new TH1D(sb, sb, nBin, xmin, xmax);
      sprintf(sb,"hDpttag_%d",n);   hDpttag[n]   = new TH1D(sb, sb,   10,    0,  200);
      sprintf(sb,"hDetatag_%d",n);  hDetatag[n]  = new TH1D(sb, sb,    6,   -3,  3.0);
      hDbtag[n]  ->Sumw2();
      hDpttag[n] ->Sumw2();
      hDetatag[n]->Sumw2();
    }
  }
  double bgdDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double weiDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Char_t xTitle[]="myX"; Char_t yTitle[]="Fraction";
                         
  int nBgd=bgdEvent.tree_->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.tree_->GetEntry(i);

    if(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection || 
        (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection)) continue;

    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay = 3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay = 5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay = 9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::qcd             ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.dstype_ == SmurfTree::ggzz            ) fDecay = 29;
    else if(bgdEvent.dstype_ == 3                          ) fDecay = 29;
    else if(bgdEvent.dstype_ == 4                          ) fDecay = 29;
    else if(bgdEvent.dstype_ == 23                         ) fDecay = 29;
    else if(bgdEvent.dstype_ == 24                         ) fDecay = 29;
    else {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;assert(0);}
    if(fDecay == -1 || fDecay > 100) fDecay = 0;//44;

    if(fDecay == 27 || fDecay == 28) {
      if(bgdEvent.lep1McId_ == 21 && bgdEvent.lep2McId_ == 21) {
        fDecay = 31;
      }
    }
    int charge = (int)(bgdEvent.lq1_ + bgdEvent.lq2_);
    double add = 1.;
    //add = scaleFactor(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), bgdEvent.nvtx_, 2);

    unsigned int Njet3 = bgdEvent.njets_;
    if(bgdEvent.njets_ >= 20){ // nJetsType = 0/1/2-jet selection
      if(bgdEvent.jet3_.pt() <= 30)					                                     Njet3 = 2;
      else if(bgdEvent.jet3_.pt() > 30 && (
    	(bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() < 0) ||
    	(bgdEvent.jet2_.eta()-bgdEvent.jet3_.eta() > 0 && bgdEvent.jet1_.eta()-bgdEvent.jet3_.eta() < 0)))   Njet3 = 3;
      else							                                             Njet3 = 2;
      if(bgdEvent.njets_ < 2 || bgdEvent.njets_ > 3)                                                         Njet3 = 3;
      //if(TMath::Abs(bgdEvent.jet1_.eta()) >= 4.5 || TMath::Abs(bgdEvent.jet2_.eta()) >= 4.5)                 Njet3 = 3;
    }
    bool passNewCuts = true;
    double usedMet = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_);
    bool   passMET = true;
    if(is2011Ana == true){
      if(bgdEvent.lep2_.Pt() <= 15 && (bgdEvent.type_ == SmurfTree::mm||bgdEvent.type_ == SmurfTree::ee)) passNewCuts = false;
      if(bgdEvent.dilep_.Pt() <= 45) passNewCuts = false;
      passMET = usedMet > 20. &&
               (usedMet > 37.+bgdEvent.nvtx_/2.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    } else {
      if(bgdEvent.lep2_.Pt() <= 10 && (bgdEvent.type_ == SmurfTree::mm||bgdEvent.type_ == SmurfTree::ee)) passNewCuts = false;
      if(bgdEvent.dilep_.Pt() <= 30) passNewCuts = false;
      if(bgdEvent.dilep_.Pt() <= 45 && (bgdEvent.type_ == SmurfTree::mm||bgdEvent.type_ == SmurfTree::ee)) passNewCuts = false;
      passMET = usedMet > 20.;
      if     (bgdEvent.njets_ == 0) passMET = passMET && (bgdEvent.dymva_ >  0.88 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else if(bgdEvent.njets_ == 1) passMET = passMET && (bgdEvent.dymva_ >  0.84 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      else                          passMET = passMET && (bgdEvent.met_ > 45.0 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      //if     (bgdEvent.njets_ == 0) passMET = passMET && (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      //else if(bgdEvent.njets_ == 1) passMET = passMET && (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
      //else			    passMET = passMET && (usedMet > 40. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me);
    }

    bool dPhiDiLepJetCut = true;
    if(is2011Ana == true) {
      if(bgdEvent.njets_ <= 1) dPhiDiLepJetCut = bgdEvent.jet1_.Pt() <= 15. || bgdEvent.dPhiDiLepJet1_*180.0/TMath::Pi() < 165. || 
    	                                         bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
      else                     dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. || 
   	                                         bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    }
    if(bgdEvent.njets_ >= 2) dPhiDiLepJetCut = DeltaPhi((bgdEvent.jet1_+bgdEvent.jet2_).Phi(),bgdEvent.dilep_.Phi())*180.0/TMath::Pi() < 165. ||
                                                         bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me;
    float pCHSMet = 0.0; float pNHSMet = 0.0; float pMet3 = bgdEvent.pmet_; float met3 = bgdEvent.met_;
    if(thePlot == 18 || thePlot == 19 || thePlot == 20 || thePlot == 21 || thePlot == 22){
      float deltaPhiCHS = TMath::Min(DeltaPhi(bgdEvent.lep1_.Phi(),bgdEvent.CHSMetPhi_),
                                     DeltaPhi(bgdEvent.lep2_.Phi(),bgdEvent.CHSMetPhi_));
      pCHSMet = bgdEvent.CHSMet_;
      if(deltaPhiCHS < TMath::Pi()/2) pCHSMet = pCHSMet * sin(deltaPhiCHS);
      
      float deltaPhiNHS = TMath::Min(DeltaPhi(bgdEvent.lep1_.Phi(),bgdEvent.NHSMetPhi_),
                                     DeltaPhi(bgdEvent.lep2_.Phi(),bgdEvent.NHSMetPhi_));
      pNHSMet = bgdEvent.NHSMet_;
      if(deltaPhiNHS < TMath::Pi()/2) pNHSMet = pNHSMet * sin(deltaPhiNHS);
      
      pMet3 = TMath::Min(pMet3,pCHSMet);
      pMet3 = TMath::Min(pMet3,pNHSMet);
      met3  = TMath::Min(met3,bgdEvent.CHSMet_);
      met3  = TMath::Min(met3,bgdEvent.NHSMet_);
    }

    if(channel == 0){ // MET selection
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
	 bgdEvent.njets_ == option &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
       (!(fDecay == 1 || fDecay == 9) || fabs(bgdEvent.dilep_.M()-91.1876) < 25.) && 
        (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto &&
	 dPhiDiLepJetCut == true &&
         bgdEvent.dilep_.Pt() > 45.0 &&
         TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > 20.0 &&
	 //bgdEvent.mt_ > 80.0 &&
	 1 == 1
	){

        if     (fDecay == 1 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) add =  1.0;
	else if(fDecay == 1)                                                                         add = -1.0;
        else if(fDecay == 9 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) add =  1.0;
	else if(fDecay == 9)                                                                         add = -1.0;
	else                                                                                         add =  1.0;
	double weight = bgdEvent.scale1fb_*lumi*add;
	bgdDecay[(int)fDecay] += weight;
	weiDecay[(int)fDecay] += weight*weight;

	double myVar = bgdEvent.met_/100.;
        if     (thePlot == 1) myVar = bgdEvent.lep1_.Pt()/100.;
        else if(thePlot == 2) myVar = bgdEvent.lep2_.Pt()/100.;
        else if(thePlot == 3) myVar = bgdEvent.dPhi_/TMath::Pi();
        else if(thePlot == 4) myVar = bgdEvent.mt_/500.;
	else if(thePlot == 5) myVar = bgdEvent.jet1_.Pt()/100;
        else if(thePlot == 6) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
	else if(thePlot == 7) myVar = bgdEvent.met_/100.;
	else if(thePlot == 8) myVar = bgdEvent.metMVA_/100.;
	else if(thePlot == 9) myVar = bgdEvent.trackMet_/100.;
	else if(thePlot ==10) myVar = TMath::Min(bgdEvent.met_,bgdEvent.trackMet_)/100;
	else if(thePlot ==11) myVar = TMath::Min(TMath::Min(bgdEvent.met_,bgdEvent.trackMet_),bgdEvent.metMVA_)/100;
	else if(thePlot ==12) myVar = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_)/100;
	else if(thePlot ==13) myVar = TMath::Min(bgdEvent.pmetMVA_,bgdEvent.pTrackMet_)/100;
	else if(thePlot ==14) myVar = TMath::Min(TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_),bgdEvent.pmetMVA_)/100;
	else if(thePlot ==15) myVar = bgdEvent.pmet_/100;
	else if(thePlot ==16) myVar = bgdEvent.pmetMVA_/100;
	else if(thePlot ==17) myVar = bgdEvent.pTrackMet_/100;
	else if(thePlot ==18) myVar = pCHSMet/100;
	else if(thePlot ==19) myVar = pNHSMet/100;
	else if(thePlot ==20) myVar = pMet3/100;
	else if(thePlot ==21) myVar = met3/100;
	else if(thePlot ==22) myVar = (bgdEvent.dymva_+1.0)/2.0;

	if     (myVar < 0) myVar = 0.001;
	else if(myVar > 1) myVar = 0.999;
	if(fDecay == 1 ||fDecay == 9){
	  hDMetBck->Fill(myVar,weight);
        }
	if((fDecay == 29 || fDecay == 30)){
	  hDMetSig->Fill(myVar,weight);
        }
	
	if((fDecay == 9) &&
	   (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)){
	  if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.) hDVarin ->Fill(myVar,weight);
	  if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.) hDVarou ->Fill(myVar,weight);
	  if(bgdEvent.dilep_.M() < 50                ) hDVarouH->Fill(myVar,weight);
	  if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.&& bgdEvent.type_ == SmurfTree::mm) hDVarinmm ->Fill(myVar,weight);
	  if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.&& bgdEvent.type_ == SmurfTree::mm) hDVaroumm ->Fill(myVar,weight);
	  if(bgdEvent.dilep_.M() < 50                && bgdEvent.type_ == SmurfTree::mm) hDVaroummH->Fill(myVar,weight);
	  if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.&& bgdEvent.type_ == SmurfTree::ee) hDVarinee ->Fill(myVar,weight);
	  if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.&& bgdEvent.type_ == SmurfTree::ee) hDVarouee ->Fill(myVar,weight);
	  if(bgdEvent.dilep_.M() < 50                && bgdEvent.type_ == SmurfTree::ee) hDVaroueeH->Fill(myVar,weight);
	  for(int n1=0; n1<nBin; n1++){
	    if(myVar > 1.0*n1/nBin){
	      if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.) hDRin ->SetBinContent(n1,hDRin ->GetBinContent(n1)+weight);
	      if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.) hDRou ->SetBinContent(n1,hDRou ->GetBinContent(n1)+weight);
   	      if(bgdEvent.dilep_.M() < 50                ) hDRouH->SetBinContent(n1,hDRouH->GetBinContent(n1)+weight);
	      if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.&& bgdEvent.type_ == SmurfTree::mm) hDRinmm ->SetBinContent(n1,hDRinmm ->GetBinContent(n1)+weight);
	      if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.&& bgdEvent.type_ == SmurfTree::mm) hDRoumm ->SetBinContent(n1,hDRoumm ->GetBinContent(n1)+weight);
   	      if(bgdEvent.dilep_.M() < 50                && bgdEvent.type_ == SmurfTree::mm) hDRoummH->SetBinContent(n1,hDRoummH->GetBinContent(n1)+weight);
	      if(fabs(bgdEvent.dilep_.M()-91.1876) <= 15.&& bgdEvent.type_ == SmurfTree::ee) hDRinee ->SetBinContent(n1,hDRinee ->GetBinContent(n1)+weight);
	      if(fabs(bgdEvent.dilep_.M()-91.1876)  > 15.&& bgdEvent.type_ == SmurfTree::ee) hDRouee ->SetBinContent(n1,hDRouee ->GetBinContent(n1)+weight);
   	      if(bgdEvent.dilep_.M() < 50                && bgdEvent.type_ == SmurfTree::ee) hDRoueeH->SetBinContent(n1,hDRoueeH->GetBinContent(n1)+weight);
	    }
	  }
	}
      }
    } // MET selection
    else if(channel == 1){ // btag selection
      bool whichDecay = (fDecay == 5 || fDecay == 13);
      if     (option == 1) whichDecay = (fDecay == 5);
      else if(option == 2) whichDecay = (fDecay == 13);
      else if(option == 3) whichDecay = (fDecay == 5 || fDecay == 13) && 
                                         TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 4.0 &&
					 bgdEvent.jet1_.Eta()*bgdEvent.jet2_.Eta() < 0 &&
					 (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 &&
					 bgdEvent.njets_ >= 2;
      else if(option == 4) whichDecay = fDecay != 999;
      if(
         bgdEvent.dilep_.M()   > 12 &&
         bgdEvent.lid3_ == 0 &&
         charge == 0 &&
         bgdEvent.lep1_.Pt() > 20. &&
         bgdEvent.lep2_.Pt() > 10. &&
        (bgdEvent.lep2_.Pt() > 15. || bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::em) &&
         bgdEvent.pmet_ > 20. &&
        (bgdEvent.pmet_ > 35. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) &&
        (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
         whichDecay == true &&
	 1 == 1
	){
	//
	double weight = bgdEvent.scale1fb_*lumi*add; 
	if(fDecay != 0) {
	  bgdDecay[(int)fDecay] += weight; 
	  weiDecay[(int)fDecay] += weight*weight;
        }
        bool combinedBtag = bgdEvent.nSoftMuons_ > 0. || bgdEvent.jetLowBtag_ >= 2.1 || bgdEvent.jet1Btag_ >= 2.1 ||  bgdEvent.jet2Btag_ >= 2.1 ||  bgdEvent.jet3Btag_ >= 2.1;
	hDbtag[0]->Fill(1.0*bgdEvent.njets_,weight);
	if(combinedBtag == true)                                      hDbtag[1]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.nSoftMuons_ > 0. || bgdEvent.jetLowBtag_ >= 2.1)  hDbtag[2]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.jet1Btag_ >= 2.1)                                 hDbtag[3]->Fill(1.0*bgdEvent.njets_,weight);
	if( bgdEvent.nSoftMuons_ == 0. && bgdEvent.jetLowBtag_ < 2.1) hDbtag[4]->Fill(1.0*bgdEvent.njets_,weight);
	if((bgdEvent.nSoftMuons_ == 0. && bgdEvent.jetLowBtag_ < 2.1)
	   && bgdEvent.jet1Btag_ >= 2.1)                              hDbtag[5]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.jet2Btag_ < 2.1)                                  hDbtag[6]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.jet2Btag_ < 2.1 && bgdEvent.jet1Btag_ >= 2.1)     hDbtag[7]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.nSoftMuons_ > 0.)                                 hDbtag[8]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.jetLowBtag_ >= 2.1)                               hDbtag[9]->Fill(1.0*bgdEvent.njets_,weight);
	if(bgdEvent.njets_ >= 2){
	  hDpttag[0] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDpttag[1] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  hDetatag[0]->Fill(bgdEvent.jet1_.Eta(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDetatag[1]->Fill(bgdEvent.jet1_.Eta(),weight);
	}
	if(bgdEvent.njets_ >= 2 && bgdEvent.nSoftMuons_ == 0. && bgdEvent.jetLowBtag_ < 2.1){
	  hDpttag[2] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDpttag[3] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  hDetatag[2]->Fill(bgdEvent.jet1_.Eta(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDetatag[3]->Fill(bgdEvent.jet1_.Eta(),weight);
	}
	if(bgdEvent.njets_ >= 2 && bgdEvent.nSoftMuons_ == 0. && bgdEvent.jetLowBtag_ < 2.1 && bgdEvent.jet2Btag_ >= 2.1){
	  hDpttag[4] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDpttag[5] ->Fill(bgdEvent.jet1_.Pt(),weight);
	  hDetatag[4]->Fill(bgdEvent.jet1_.Eta(),weight);
	  if(bgdEvent.jet1Btag_ >= 2.1) hDetatag[5]->Fill(bgdEvent.jet1_.Eta(),weight);
	}
      }
    } // btag selection
     //WW selection
    if(channel == 2) {
      bool passCuts[2] = {false, false};
      unsigned int pattern = SmurfTree::BaseLine|SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|SmurfTree::ChargeMatch|SmurfTree::FullMET|SmurfTree::ZVeto|SmurfTree::ExtraLeptonVeto|SmurfTree::TopVeto;
      //unsigned int pattern = SmurfTree::BaseLine|SmurfTree::Lep1FullSelection|SmurfTree::Lep2FullSelection|SmurfTree::ChargeMatch|SmurfTree::ExtraLeptonVeto;
      if((bgdEvent.cuts_ & pattern) == pattern &&
         // (bgdEvent.cuts_ & SmurfTree::TopTag) == 0 &&
	  (bgdEvent.njets_ == option || (int)option == -1)){
	passCuts[0] = true; 
	bgdDecay[(int)fDecay] += 1;	
	weiDecay[(int)fDecay] += 1;
      }
      if(
         bgdEvent.lid3_ == 0 &&
	 charge == 0 &&
	 bgdEvent.dilep_.M()   > 12 && bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 10. && (bgdEvent.lep2_.Pt() > 15. || bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::em) &&
         passMET == true &&
         (fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) && 
          (bgdEvent.nSoftMuons_ == 0. && bgdEvent.jetLowBtag_ < 2.1 &&
	  (bgdEvent.jet1_.Pt() <= 7 || bgdEvent.jet1Btag_ < 2.1) &&
	  (bgdEvent.jet2_.Pt() <= 7 || bgdEvent.jet2Btag_ < 2.1) && 
	  (bgdEvent.jet3_.Pt() <= 7 || bgdEvent.jet3Btag_ < 2.1)) &&
	 (bgdEvent.njets_ == option || (int)option == -1) &&
        1 & 1) {
	passCuts[1] = true; 
      }
      if(passCuts[0] != passCuts[1]) cout << passCuts[0] << " " << passCuts[1] << endl;
      if(passCuts[0] != passCuts[1]) cout << bgdEvent.lep1_.Pt()  << " " << bgdEvent.lep2_.Pt()<<" " <<bgdEvent.lep3_.Pt()<<endl;
    }
    // cut evolution
    if(channel == 3) {
      bool goodEvent = true;
      sumEvol[0] = 0.0; sumEvol[1] = 0.0; sumEvol[2] = 0.0;
      if((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
         (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection &&
          bgdEvent.lep1_.Pt() > 20. &&
          bgdEvent.lep2_.Pt() > 20. &&
	  charge == 0 && goodEvent) {
	 sumEvol[0] = 1.0; sumEvol[1] = 1.0; sumEvol[2] = 1.0;
	 //cout << "SSS " << bgdEvent.event_ << " " << bgdEvent.type_ << " " << Njet3 << endl;
         //     printf("SSS %d %d %d %f %f %f %f %f\n",bgdEvent.event_,bgdEvent.type_,bgdEvent.njets_,bgdEvent.met_,bgdEvent.pmet_,bgdEvent.dymva_,usedMet,bgdEvent.dilep_.M());
        if(bgdEvent.lid3_ == 0) {
          sumEvol[0] = 2.0; sumEvol[1] = 2.0; sumEvol[2] = 2.0;
      	  if(bgdEvent.met_ > 20.) {
	     sumEvol[0] = 3.0; sumEvol[1] = 3.0; sumEvol[2] = 3.0;
            if(bgdEvent.dilep_.M() > 12 && 
	      (bgdEvent.dilep_.M() > 20 || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me || is2011Ana == false)) {
              sumEvol[0] = 4.0; sumEvol[1] = 4.0; sumEvol[2] = 4.0;
              if(fabs(bgdEvent.dilep_.M()-91.1876) > 15. || bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me) {
        	sumEvol[0] = 5.0; sumEvol[1] = 5.0; sumEvol[2] = 5.0;
        	//printf("SSS %d %d %d %f %f %f %f %f\n",bgdEvent.event_,bgdEvent.type_,bgdEvent.njets_,
		//bgdEvent.met_,bgdEvent.pmet_,bgdEvent.dymva_,usedMet,bgdEvent.dilep_.M());
        	if(passMET == true) {
	  	  sumEvol[0] = 6.0; sumEvol[1] = 6.0; sumEvol[2] = 6.0;
          	  //printf("SSS %d %d %d %f %f %f %f %f\n",bgdEvent.event_,bgdEvent.type_,bgdEvent.njets_,bgdEvent.met_,bgdEvent.pmet_,bgdEvent.dymva_,usedMet,bgdEvent.dilep_.M());
          	  if(dPhiDiLepJetCut == true) {
          	    sumEvol[0] = 7.0; sumEvol[1] = 7.0; sumEvol[2] = 7.0;
          	    if(bgdEvent.nSoftMuons_ == 0.) {
              //printf("SSS %d %d %d %f %f %f %f %f\n",bgdEvent.event_,bgdEvent.type_,Njet3,bgdEvent.met_,bgdEvent.pmet_,bgdEvent.dymva_,usedMet,bgdEvent.dilep_.M());
          	      sumEvol[0] = 8.0; sumEvol[1] = 8.0; sumEvol[2] = 8.0;
                      if((bgdEvent.cuts_ & patternTopVeto) == patternTopVeto) {
                        sumEvol[0] = 9.0; sumEvol[1] = 9.0; sumEvol[2] = 9.0;
                        if(passNewCuts == true) {
              //printf("SSS %d %d %f %f %f %f %f %f %f %f\n",bgdEvent.event_,bgdEvent.type_,bgdEvent.jet1_.pt(),bgdEvent.jet1_.eta(),bgdEvent.jet2_.pt(),bgdEvent.jet2_.eta(),bgdEvent.jet3_.pt(),bgdEvent.jet3_.eta(),bgdEvent.jet4_.pt(),bgdEvent.jet4_.eta());
                          sumEvol[0] = 10.0; sumEvol[1] = 10.0; sumEvol[2] = 10.0;
                          if	   (Njet3 == 0) {
                            sumEvol[0] = 11.0;
		          } else if(Njet3 == 1) {
                            sumEvol[1] = 11.0;
		          } else if(Njet3 == 2) {
                            sumEvol[2] = 11.0;
		          }
                        } // 10
                      } // 9
                    } // 8
                  } // 7
                } // 6
              } // 5
            } // 4
          } // 3
        } // 2
      } // 1
    } // 0
    for(int nl = 0; nl <=(int)sumEvol[0]; nl++) hDWWLL0->Fill((double)nl);
    for(int nl = 0; nl <=(int)sumEvol[1]; nl++) hDWWLL1->Fill((double)nl);
    for(int nl = 0; nl <=(int)sumEvol[2]; nl++) hDWWLL2->Fill((double)nl);
    if(bgdEvent.type_ == SmurfTree::mm){
      for(int nl = 0; nl <=(int)sumEvol[0]; nl++) hDWWMM0->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[1]; nl++) hDWWMM1->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[2]; nl++) hDWWMM2->Fill((double)nl);    
    }
    if(bgdEvent.type_ == SmurfTree::em){
      for(int nl = 0; nl <=(int)sumEvol[0]; nl++) hDWWEM0->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[1]; nl++) hDWWEM1->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[2]; nl++) hDWWEM2->Fill((double)nl);    
    }
    if(bgdEvent.type_ == SmurfTree::me){
      for(int nl = 0; nl <=(int)sumEvol[0]; nl++) hDWWME0->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[1]; nl++) hDWWME1->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[2]; nl++) hDWWME2->Fill((double)nl);	
    }
    if(bgdEvent.type_ == SmurfTree::ee){
      for(int nl = 0; nl <=(int)sumEvol[0]; nl++) hDWWEE0->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[1]; nl++) hDWWEE1->Fill((double)nl);
      for(int nl = 0; nl <=(int)sumEvol[2]; nl++) hDWWEE2->Fill((double)nl);    
    }
  } // big loop

  for(int i=0; i<45; i++){
    if(bgdDecay[i] > 0) printf("bdg(%2d) = %f +/- %f\n",i,bgdDecay[i],sqrt(weiDecay[i]));
  }

  if(channel == 3){ // cut evolution
      char output[200];
      sprintf(output,"cutEvolution.root");     
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      printf("    MM      EE      EM      ME      LL\n");
      if(iniEvents == 1){
        for(int nl =  2; nl <=11; nl++) printf("%6.0f  %6.0f  %6.0f  %6.0f  %6.0f\n",hDWWMM0->GetBinContent(nl),hDWWEE0->GetBinContent(nl),hDWWEM0->GetBinContent(nl),hDWWME0->GetBinContent(nl),hDWWLL0->GetBinContent(nl));
        printf("----0j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.0f  %6.0f  %6.0f  %6.0f  %6.0f\n",hDWWMM0->GetBinContent(nl),hDWWEE0->GetBinContent(nl),hDWWEM0->GetBinContent(nl),hDWWME0->GetBinContent(nl),hDWWLL0->GetBinContent(nl));
        printf("----1j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.0f  %6.0f  %6.0f  %6.0f  %6.0f\n",hDWWMM1->GetBinContent(nl),hDWWEE1->GetBinContent(nl),hDWWEM1->GetBinContent(nl),hDWWME1->GetBinContent(nl),hDWWLL1->GetBinContent(nl));
        printf("----2j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.0f  %6.0f  %6.0f  %6.0f  %6.0f\n",hDWWMM2->GetBinContent(nl),hDWWEE2->GetBinContent(nl),hDWWEM2->GetBinContent(nl),hDWWME2->GetBinContent(nl),hDWWLL2->GetBinContent(nl));
      } else {
        for(int nl =  2; nl <=11; nl++) printf("%6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",hDWWMM0->GetBinContent(nl)/iniEvents,hDWWEE0->GetBinContent(nl)/iniEvents,hDWWEM0->GetBinContent(nl)/iniEvents,hDWWME0->GetBinContent(nl)/iniEvents,hDWWLL0->GetBinContent(nl)/iniEvents);
        printf("----0j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",hDWWMM0->GetBinContent(nl)/iniEvents,hDWWEE0->GetBinContent(nl)/iniEvents,hDWWEM0->GetBinContent(nl)/iniEvents,hDWWME0->GetBinContent(nl)/iniEvents,hDWWLL0->GetBinContent(nl)/iniEvents);
        printf("----1j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",hDWWMM1->GetBinContent(nl)/iniEvents,hDWWEE1->GetBinContent(nl)/iniEvents,hDWWEM1->GetBinContent(nl)/iniEvents,hDWWME1->GetBinContent(nl)/iniEvents,hDWWLL1->GetBinContent(nl)/iniEvents);
        printf("----2j-----\n");
        for(int nl = 12; nl <=12; nl++) printf("%6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",hDWWMM2->GetBinContent(nl)/iniEvents,hDWWEE2->GetBinContent(nl)/iniEvents,hDWWEM2->GetBinContent(nl)/iniEvents,hDWWME2->GetBinContent(nl)/iniEvents,hDWWLL2->GetBinContent(nl)/iniEvents);
      }
      outFilePlotsNote->cd();
   	hDWWLL0->Write();
   	hDWWLL1->Write();
   	hDWWLL2->Write();
   	hDWWMM0->Write();
   	hDWWMM1->Write();
   	hDWWMM2->Write();
   	hDWWEM0->Write();
   	hDWWEM1->Write();
   	hDWWEM2->Write();
   	hDWWME0->Write();
   	hDWWME1->Write();
   	hDWWME2->Write();
   	hDWWEE0->Write();
   	hDWWEE1->Write();
   	hDWWEE2->Write();
      outFilePlotsNote->Close();
      delete outFilePlotsNote;
  }

  if(channel == 0){ // MET selection
    char fAsciiFileName[200];
    sprintf(fAsciiFileName,"metStudy_%d_%1dj.txt",thePlot,(int)option);
    cout << "Writing to... " << fAsciiFileName << endl;
    std::ofstream *fOFile0 = new std::ofstream(fAsciiFileName,ios::out);
    double S = 0.0;
    double B = 0.0;
    double total[2] = {hDMetSig->GetSumOfWeights(), hDMetBck->GetSumOfWeights()};
    for(int i=hDMetSig->GetNbinsX(); i>=1; i--){
    //for(int i=1; i<=hDMetSig->GetNbinsX(); i++){
      S = S + hDMetSig->GetBinContent(i);
      B = B + hDMetBck->GetBinContent(i);
      double effS = S /  total[0];
      double effB = B /  total[1];
      *fOFile0
      <<        setw(10) << effS
      << " " << setw(10) << effB
      << endl;
    }

    hDRou->Divide(hDRin);
    hDRoumm->Divide(hDRinmm);
    hDRouee->Divide(hDRinee);
    hDRouH->Divide(hDRin);
    hDRoummH->Divide(hDRinmm);
    hDRoueeH->Divide(hDRinee);

    atributes(hDRou   ,xTitle,1,yTitle);
    atributes(hDRoumm ,xTitle,2,yTitle);
    atributes(hDRouee ,xTitle,4,yTitle);
    atributes(hDRouH  ,xTitle,1,yTitle);
    atributes(hDRoummH,xTitle,2,yTitle);
    atributes(hDRoueeH,xTitle,4,yTitle);

    //hDRou->Draw("hist,e");
    //hDRoumm->Draw("same,hist,e");
    //hDRouee->Draw("same,hist,e");
    hDRouH->Draw("hist,e");
    hDRoummH->Draw("same,hist,e");
    hDRoueeH->Draw("same,hist,e");
    atributes(hDVarou	,xTitle,1,yTitle);
    atributes(hDVaroumm ,xTitle,2,yTitle);
    atributes(hDVarouee ,xTitle,4,yTitle);
    atributes(hDVarouH  ,xTitle,1,yTitle);
    atributes(hDVaroummH,xTitle,2,yTitle);
    atributes(hDVaroueeH,xTitle,4,yTitle);
    atributes(hDVarin	,xTitle,1,yTitle);
    atributes(hDVarinmm ,xTitle,2,yTitle);
    atributes(hDVarinee ,xTitle,4,yTitle);
    hDVarou   ->Scale(1./hDVarou   ->GetSumOfWeights());
    hDVaroumm ->Scale(1./hDVaroumm ->GetSumOfWeights());
    hDVarouee ->Scale(1./hDVarouee ->GetSumOfWeights());
    hDVarouH  ->Scale(1./hDVarouH  ->GetSumOfWeights());
    hDVaroummH->Scale(1./hDVaroummH->GetSumOfWeights());
    hDVaroueeH->Scale(1./hDVaroueeH->GetSumOfWeights());
    hDVarin   ->Scale(1./hDVarin   ->GetSumOfWeights());
    hDVarinmm ->Scale(1./hDVarinmm ->GetSumOfWeights());
    hDVarinee ->Scale(1./hDVarinee ->GetSumOfWeights());

    if(makeGoodPlots){
      char output[200];
      sprintf(output,"histo_nice.root");     
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
   	hDRou	  ->Write();
   	hDRoumm   ->Write();
   	hDRouee   ->Write();
   	hDRouH    ->Write();
   	hDRoummH  ->Write();
   	hDRoueeH  ->Write();
   	hDVarou   ->Write();
   	hDVaroumm ->Write();
   	hDVarouee ->Write();
   	hDVarouH  ->Write();
   	hDVaroummH->Write();
   	hDVaroueeH->Write();
   	hDVarin   ->Write();
   	hDVarinmm ->Write();
   	hDVarinee ->Write();
      outFilePlotsNote->Close();
    }
  }
  else if(channel == 1){ // btag selection
    if(thePlot == 0){
      sprintf(xTitle,"N_{jets}");
      sprintf(yTitle,"efficiency");
      hDbtag[1]->Divide(hDbtag[0]);
      hDbtag[2]->Divide(hDbtag[0]);
      hDbtag[3]->Divide(hDbtag[0]);
      hDbtag[5]->Divide(hDbtag[4]);
      hDbtag[7]->Divide(hDbtag[6]);
      hDbtag[8]->Divide(hDbtag[0]);
      hDbtag[9]->Divide(hDbtag[0]);
      hDbtag[1]->SetMinimum(0.0);
      hDbtag[1]->SetMaximum(1.5);
      atributes(hDbtag[1],xTitle,1,yTitle);
      atributes(hDbtag[2],xTitle,2,yTitle);
      atributes(hDbtag[3],xTitle,4,yTitle);
      atributes(hDbtag[5],xTitle,6,yTitle);
      atributes(hDbtag[7],xTitle,13,yTitle);
      atributes(hDbtag[8],xTitle,8,yTitle);
      atributes(hDbtag[9],xTitle,11,yTitle);
      hDbtag[1]->Draw();
      hDbtag[2]->Draw("same");
      hDbtag[3]->Draw("same");
      hDbtag[5]->Draw("same");
      hDbtag[7]->Draw("same");
      hDbtag[8]->Draw("same");
      hDbtag[9]->Draw("same");
      TLegend* leg = new TLegend(0.2,0.7,0.85,0.9);			        		
      leg->SetFillColor(10);						        		
      leg->SetTextSize(0.025);
      leg->AddEntry(hDbtag[1] ,"total jet tagging","l");		    
      leg->AddEntry(hDbtag[2] ,"low pt jet and soft muon tagging","l");	          
      leg->AddEntry(hDbtag[3] ,"highest pt jet tagging","l");
      leg->AddEntry(hDbtag[5] ,"highest pt jet tagging with no low jet pt tagged","l");    	  
      leg->AddEntry(hDbtag[7] ,"highest pt jet tagging with no 2nd jet tagged","l");          
      leg->AddEntry(hDbtag[8] ,"soft muon tagging","l");	    
      leg->AddEntry(hDbtag[9] ,"low pt jet tagging","l");	    
      leg->Draw();	    
    }
    else if(thePlot == 1){
      sprintf(xTitle,"N_{jets}");
      sprintf(yTitle,"efficiency");
      hDbtag[2]->Divide(hDbtag[0]);
      hDbtag[8]->Divide(hDbtag[0]);
      hDbtag[9]->Divide(hDbtag[0]);
      hDbtag[2]->SetMinimum(0.0);
      hDbtag[2]->SetMaximum(1.5);
      atributes(hDbtag[2],xTitle,2,yTitle);
      atributes(hDbtag[8],xTitle,8,yTitle);
      atributes(hDbtag[9],xTitle,11,yTitle);
      hDbtag[2]->Draw();
      hDbtag[8]->Draw("same");
      hDbtag[9]->Draw("same");
      TLegend* leg = new TLegend(0.2,0.75,0.85,0.9);			        		
      leg->SetFillColor(10);						        		
      leg->SetTextSize(0.03);
      leg->AddEntry(hDbtag[2] ,"low pt jet and soft muon tagging","l");	          
      leg->AddEntry(hDbtag[8] ,"soft muon tagging","l");	    
      leg->AddEntry(hDbtag[9] ,"low pt jet tagging","l");	    
      leg->Draw();	    
    }
    else if(thePlot == 2){
      sprintf(xTitle,"N_{jets}");
      sprintf(yTitle,"efficiency");
      hDbtag[1]->Divide(hDbtag[0]);
      hDbtag[3]->Divide(hDbtag[0]);
      hDbtag[5]->Divide(hDbtag[4]);
      hDbtag[7]->Divide(hDbtag[6]);
      hDbtag[1]->SetMinimum(0.0);
      hDbtag[1]->SetMaximum(1.5);
      atributes(hDbtag[1],xTitle,1,yTitle);
      atributes(hDbtag[3],xTitle,4,yTitle);
      atributes(hDbtag[5],xTitle,6,yTitle);
      atributes(hDbtag[7],xTitle,13,yTitle);
      hDbtag[1]->Draw();
      hDbtag[3]->Draw("same");
      hDbtag[5]->Draw("same");
      hDbtag[7]->Draw("same");
      TLegend* leg = new TLegend(0.2,0.7,0.85,0.9);			        		
      leg->SetFillColor(10);						        		
      leg->SetTextSize(0.025);
      leg->AddEntry(hDbtag[1] ,"total jet tagging","l");		    
      leg->AddEntry(hDbtag[3] ,"highest pt jet tagging","l");
      leg->AddEntry(hDbtag[5] ,"highest pt jet tagging with no low jet pt tagged","l");    	  
      leg->AddEntry(hDbtag[7] ,"highest pt jet tagging with no 2nd jet tagged","l");          
      leg->Draw();	    
    }
    else if(thePlot == 3){
      sprintf(yTitle,"efficiency");
      hDpttag[1]->SetMinimum(0.0);
      hDpttag[1]->SetMaximum(1.2);
      hDpttag[1]->Divide(hDpttag[0]);
      hDpttag[3]->Divide(hDpttag[2]);
      hDpttag[5]->Divide(hDpttag[4]);
      atributes(hDpttag[1] ,"p_{T} [GeV]",1,yTitle);
      atributes(hDpttag[3] ,"p_{T} [GeV]",2,yTitle);
      atributes(hDpttag[5] ,"p_{T} [GeV]",4,yTitle);
      hDpttag[1] ->Draw();
      hDpttag[3] ->Draw("same");
      hDpttag[5] ->Draw("same");
      TLegend* leg = new TLegend(0.2,0.8,0.9,0.9);			        		
      leg->SetFillColor(10);						        		
      leg->SetTextSize(0.02);
      leg->AddEntry(hDpttag[1] ,"1rst jet","l");	    
      leg->AddEntry(hDpttag[3] ,"1rst jet with no low jet pt tagged","l");	  
      leg->AddEntry(hDpttag[5] ,"1rst jet with no low jet pt tagged and btagged 2nd jet","l");         
      leg->Draw();	    
    }
    else if(thePlot == 4){
      sprintf(yTitle,"efficiency");
      hDetatag[1]->SetMinimum(0.0);
      hDetatag[1]->SetMaximum(1.2);
      hDetatag[1]->Divide(hDetatag[0]);
      hDetatag[3]->Divide(hDetatag[2]);
      hDetatag[5]->Divide(hDetatag[4]);
      atributes(hDetatag[1] ,"#eta",1,yTitle);
      atributes(hDetatag[3] ,"#eta",2,yTitle);
      atributes(hDetatag[5] ,"#eta",4,yTitle);
      hDetatag[1] ->Draw();
      hDetatag[3] ->Draw("same");
      hDetatag[5] ->Draw("same");
      TLegend* leg = new TLegend(0.2,0.8,0.9,0.9);			        		
      leg->SetFillColor(10);						        		
      leg->SetTextSize(0.02);
      leg->AddEntry(hDetatag[1] ,"1rst jet","l");	    
      leg->AddEntry(hDetatag[3] ,"1rst jet with no low jet pt tagged","l");	  
      leg->AddEntry(hDetatag[5] ,"1rst jet with no low jet pt tagged and btagged 2nd jet","l");         
      leg->Draw();	    
    }
  }

  return;

}
