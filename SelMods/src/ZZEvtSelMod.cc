// $Id: ZZEvtSelMod.cc,v 1.11 2012/07/05 07:22:12 ceballos Exp $

#include "Analysis/SelMods/interface/ZZEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/ChargedParticleCol.h"

using namespace mithep;
ClassImp(mithep::ZZEvtSelMod)

//-----------------------------------------------------------------------------
ZZEvtSelMod::ZZEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fMuonName("random"),
  fElectronName("random"),
  fLeptonName("random"),
  fJetName("random"),
  fMuonFakeName("random"),
  fElectronFakeName("random"),
  fLeptonFakeName("random"),
  fPFMetName("PFMet"),
  fIsData(kFALSE),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0)
{
  // Constructor.
}

//-----------------------------------------------------------------------------
void ZZEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//-----------------------------------------------------------------------------
void ZZEvtSelMod::Process()
{
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  MuonOArr     *CleanMuons             = GetObjThisEvt<MuonOArr>    (fMuonName);
  ElectronOArr *CleanElectrons         = GetObjThisEvt<ElectronOArr>(fElectronName);
  ParticleOArr *Leptons                = GetObjThisEvt<ParticleOArr>(fLeptonName);
  JetOArr      *CleanJets              = GetObjThisEvt<JetOArr>     (fJetName);

  //MuonOArr     *CleanMuonsFakeable     = GetObjThisEvt<MuonOArr>    (fMuonFakeName);
  //ElectronOArr *CleanElectronsFakeable = GetObjThisEvt<ElectronOArr>(fElectronFakeName);
  ParticleOArr *LeptonsFakeable        = GetObjThisEvt<ParticleOArr>(fLeptonFakeName);

  LoadBranch(fPFMetName);
  LoadBranch(fEvtHdrName);

  const PFMet *pfMet                   = fPFMet->At(0);
  if(fIsData == kFALSE){
    MCParticleOArr *GenLeptons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
    vector<unsigned int> genLeptonsAllZZUsed;
    vector<unsigned int> genLeptonsHighPtZZUsed;
    for(unsigned int i=0; i<GenLeptons->GetEntries(); i++) genLeptonsAllZZUsed.push_back(0);
    for(unsigned int i=0; i<GenLeptons->GetEntries(); i++) genLeptonsHighPtZZUsed.push_back(0);
    vector<CompositeParticle> diGenLeptonsAllZZ;
    vector<CompositeParticle> diGenLeptonsHighPtZZ;
    if(GenLeptons && GenLeptons->GetEntries() >= 2) {
      for(unsigned int i=0; i<GenLeptons->GetEntries()-1; i++){
	MCParticle *mc1 = GenLeptons->At(i);
	if(mc1->AbsEta() >= 2.5) continue;
	if(mc1->Pt()     <= 5.0) continue;
	for(unsigned int j=i+1; j<GenLeptons->GetEntries(); j++){
	  MCParticle *mc2 = GenLeptons->At(j);
          if(mc2->AbsEta() >= 2.5) continue;
          if(mc2->Pt()     <= 5.0) continue;
          CompositeParticle dilepton;
          dilepton.AddDaughter(mc1);
          dilepton.AddDaughter(mc2);
	  if(TMath::Abs(dilepton.Mass()-91.1876) < 15.0 &&
	     genLeptonsAllZZUsed[i] == 0 && genLeptonsAllZZUsed[j] == 0) {
	    genLeptonsAllZZUsed[i] = 1;
	    genLeptonsAllZZUsed[j] = 1;
	    diGenLeptonsAllZZ.push_back(dilepton);
	  } // Z candidate, ptl > 5
	  if(TMath::Abs(dilepton.Mass()-91.1876) < 15.0 &&
	     genLeptonsHighPtZZUsed[i] == 0 && genLeptonsHighPtZZUsed[j] == 0 &&
	     mc1->Pt() > 10.0 && mc2->Pt() > 10.0) {
	    genLeptonsHighPtZZUsed[i] = 1;
	    genLeptonsHighPtZZUsed[j] = 1;
	    diGenLeptonsHighPtZZ.push_back(dilepton);
	  } // Z candidate, ptl > 10
	}
      }
    }
    hDZZXGenSel[0]->Fill(TMath::Min((double)diGenLeptonsAllZZ.size(),4.499),NNLOWeight->GetVal());
    hDZZXGenSel[1]->Fill(TMath::Min((double)diGenLeptonsHighPtZZ.size(),4.499),NNLOWeight->GetVal());
    if(diGenLeptonsAllZZ.size() == 2){
      vector<MCParticle> AdditionalGenLeptons;
      vector<MCParticle> AdditionalGenLeptonsHighPt;
      for(unsigned int i=0; i<GenLeptons->GetEntries(); i++){
        MCParticle *mc = GenLeptons->At(i);
        if(mc->AbsEta() >= 2.5) continue;
        if(mc->Pt()     <= 5.0) continue;
        if(genLeptonsAllZZUsed[i] == 0) AdditionalGenLeptons.push_back(*mc);
        if(genLeptonsHighPtZZUsed[i] == 0 && mc->Pt() > 10.0) AdditionalGenLeptonsHighPt.push_back(*mc);
      }
      hDZZXGenSel[2]->Fill(TMath::Min((double)AdditionalGenLeptons.size(),4.499),NNLOWeight->GetVal());
      hDZZXGenSel[3]->Fill(TMath::Min((double)AdditionalGenLeptonsHighPt.size(),4.499),NNLOWeight->GetVal());
    }
  }

  hDZZXSel[0]->Fill(TMath::Min((double)Leptons->GetEntries(),9.499),NNLOWeight->GetVal());

  if(Leptons->GetEntries() < 4) {return;}
  if(Leptons->GetEntries()-CleanMuons->GetEntries()-CleanElectrons->GetEntries()!=0) {assert(0);}

  vector<ChargedParticle> leptonsZZ;
  vector<CompositeParticle> diLeptonsZZ;
  CompositeParticle diZZ;
  vector<unsigned int> muonsZZUsed;
  vector<unsigned int> electronsZZUsed;
  for(unsigned int i=0; i<CleanMuons    ->GetEntries(); i++) muonsZZUsed.push_back(0);
  for(unsigned int i=0; i<CleanElectrons->GetEntries(); i++) electronsZZUsed.push_back(0);

  if(CleanMuons->GetEntries() >= 2){
    for(unsigned int i=0; i<CleanMuons->GetEntries()-1; i++){
      Muon *mu1 = CleanMuons->At(i);
      for(unsigned int j=i+1; j<CleanMuons->GetEntries(); j++){
	Muon *mu2 = CleanMuons->At(j);
	if(mu1->Charge() != mu2->Charge()){
          CompositeParticle dilepton;
          dilepton.AddDaughter(mu1);
          dilepton.AddDaughter(mu2);
	  if(TMath::Abs(dilepton.Mass()-91.1876) < 15.0 &&
	     muonsZZUsed[i] == 0 && muonsZZUsed[j] == 0) {
	    muonsZZUsed[i] = 1;
	    muonsZZUsed[j] = 1;
            leptonsZZ.push_back(*mu1);
            leptonsZZ.push_back(*mu2);
	    diLeptonsZZ.push_back(dilepton);
            diZZ.AddDaughter(mu1);
            diZZ.AddDaughter(mu2);	   
	  } // Z candidate
	} // q1 != q2
      } // j loop
    } // i loop
  }
  if(CleanElectrons->GetEntries() >= 2){
    for(unsigned int i=0; i<CleanElectrons->GetEntries()-1; i++){
      Electron *el1 = CleanElectrons->At(i);
      for(unsigned int j=i+1; j<CleanElectrons->GetEntries(); j++){
	Electron *el2 = CleanElectrons->At(j);
	if(el1->Charge() != el2->Charge()){ 
          CompositeParticle dilepton;
          dilepton.AddDaughter(el1);
          dilepton.AddDaughter(el2);
	  if(TMath::Abs(dilepton.Mass()-91.1876) < 15.0 &&
	     electronsZZUsed[i] == 0 && electronsZZUsed[j] == 0) {
	    electronsZZUsed[i] = 1;
	    electronsZZUsed[j] = 1;
            leptonsZZ.push_back(*el1);
            leptonsZZ.push_back(*el2);
	    diLeptonsZZ.push_back(dilepton);
            diZZ.AddDaughter(el1);
            diZZ.AddDaughter(el2);	   
	  } // Z candidate
	} // q1 != q2
      } // j loop
    } // i loop
  }

  hDZZXSel[1]->Fill(TMath::Min((double)diLeptonsZZ.size(),4.499),NNLOWeight->GetVal());

  if(diLeptonsZZ.size() < 2) {return;}

  unsigned int tightLeptons = 0;
  for(unsigned int i=0; i<Leptons->GetEntries(); i++){
    for(unsigned int j=0; j<LeptonsFakeable->GetEntries(); j++){
      if(Leptons->At(i) == LeptonsFakeable->At(j)) {tightLeptons++;continue;}
    }
  }
  hDZZXSel[2]->Fill(TMath::Min((double)(Leptons->GetEntries()-tightLeptons),9.499),NNLOWeight->GetVal());

  if(tightLeptons < 1) {return;}

  for(unsigned int i=0; i<diLeptonsZZ.size(); i++) hDZZXSel[3]->Fill(TMath::Min(TMath::Abs(diLeptonsZZ[i].Mass()-91.1876),14.999),NNLOWeight->GetVal());

  vector<Jet> sortedJets;
  CompositeParticle dijet;
  for(unsigned int i=0; i<CleanJets->GetEntries(); i++) {
    Jet *jet = CleanJets->At(i);
    if(jet->Pt() > 30) {
      sortedJets.push_back(*jet);
      if(sortedJets.size() <= 2) dijet.AddDaughter(jet);
    }
  }
  hDZZXSel[4]->Fill(TMath::Min((double)sortedJets.size(),4.499),NNLOWeight->GetVal());

  if(Leptons->GetEntries() == 4){
    if(sortedJets.size() >= 2) {
      hDZZXSel[5]->Fill(TMath::Min(dijet.Mass(),199.999),NNLOWeight->GetVal());
    }
    hDZZXSel[6]->Fill(TMath::Min(diZZ.Mass(),999.999),NNLOWeight->GetVal());
  } else {
    double mTW = -1.0;
    hDZZXSel[7]->Fill(TMath::Min(pfMet->Pt(),199.999),NNLOWeight->GetVal());
    hDZZXSel[8]->Fill(TMath::Min((double)Leptons->GetEntries()-4.0,4.499),NNLOWeight->GetVal());
    for(unsigned int i=0; i<CleanMuons->GetEntries(); i++){
      Muon *mu = CleanMuons->At(i);
      if(muonsZZUsed[i]     == 0) mTW = TMath::Sqrt(2.0*mu->Pt()*pfMet->Pt()*(1.0 - cos(MathUtils::DeltaPhi(pfMet->Phi(),mu->Phi()))));
    }
    for(unsigned int i=0; i<CleanElectrons->GetEntries(); i++){
      Electron *el = CleanElectrons->At(i);
      if(electronsZZUsed[i] == 0) mTW = TMath::Sqrt(2.0*el->Pt()*pfMet->Pt()*(1.0 - cos(MathUtils::DeltaPhi(pfMet->Phi(),el->Phi()))));
    }
    if(Leptons->GetEntries()-4 == 1) hDZZXSel[9]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
    if(fIsData){
      printf("ZZXSel_data: %d %d %d %d\n",Leptons->GetEntries(),fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
    }
  }

}

//--------------------------------------------------------------------------------------------------
void ZZEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fPFMetName,fPFMet);
  ReqBranch(fEvtHdrName,fEventHeader);

  char sb[200];

  sprintf(sb,"hDZZXSel_0"); hDZZXSel[0] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZZXSel_1"); hDZZXSel[1] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXSel_2"); hDZZXSel[2] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZZXSel_3"); hDZZXSel[3] = new TH1D(sb,sb,30,0.0,15.0);
  sprintf(sb,"hDZZXSel_4"); hDZZXSel[4] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXSel_5"); hDZZXSel[5] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZZXSel_6"); hDZZXSel[6] = new TH1D(sb,sb,500,0.0,1000.0);
  sprintf(sb,"hDZZXSel_7"); hDZZXSel[7] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZZXSel_8"); hDZZXSel[8] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXSel_9"); hDZZXSel[9] = new TH1D(sb,sb,200,0.0,200.0);
  for(int i=0; i<10; i++) AddOutput(hDZZXSel[i]);

  sprintf(sb,"hDZZXGenSel_0"); hDZZXGenSel[0] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXGenSel_1"); hDZZXGenSel[1] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXGenSel_2"); hDZZXGenSel[2] = new TH1D(sb,sb, 5,-0.5,4.5);
  sprintf(sb,"hDZZXGenSel_3"); hDZZXGenSel[3] = new TH1D(sb,sb, 5,-0.5,4.5);
  for(int i=0; i<4; i++) AddOutput(hDZZXGenSel[i]);

}
//--------------------------------------------------------------------------------------------------
void ZZEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ZZEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
