// $Id: FRStudy.cc,v 1.4 2012/01/26 13:17:56 ceballos Exp $

#include "Ana/SelMods/interface/FRStudy.h"
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"

using namespace mithep;
ClassImp(mithep::FRStudy)

//--------------------------------------------------------------------------------------------------
FRStudy::FRStudy(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsFastSim(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("randomName1"),
  fCleanJetsName("randomName2"),
  fVertexName(ModNames::gkGoodVertexesName),
  fVertices(0),
  fIsData(kFALSE),
  fSelectGenLeptons(kTRUE),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void FRStudy::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void FRStudy::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //MuonOArr  *CleanMuonsFakeable        = GetObjThisEvt<MuonOArr>("CleanMuonsFakeable");
  //MuonOArr  *CleanMuons                = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  //ElectronOArr *CleanElectronsFakeable = GetObjThisEvt<ElectronOArr>("CleanElectronsFakeable");
  //ElectronOArr *CleanElectrons         = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  ParticleOArr *leptonsFakeable        = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  ParticleOArr *leptons                = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MCParticleOArr *GenLeptons           = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);

  MetOArr *CleanMet	               = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet	               = CleanMet->At(0);
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  const TriggerObjectCol *objs = GetHLTObjects("myhltobjsFR");

  // No need to go further if no fakeable lepton is found in the event
  if(leptonsFakeable->GetEntries() == 0) return;

  for(UInt_t i=0; i<leptonsFakeable->GetEntries(); i++){
    if(leptonsFakeable->At(i)->Pt() <= 10) continue;

    Int_t ents=objs->GetEntries();
    double deltaRleptonHLT = 999.9;
    for(Int_t nt=0;nt<ents;++nt) {
      const TriggerObject* to = objs->At(nt);
      if     (leptonsFakeable->At(i)->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || 
    						             to->Type() == TriggerObject::TriggerCluster)){
        double DeltaRAux = MathUtils::DeltaR(leptonsFakeable->At(i)->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
      else if(leptonsFakeable->At(i)->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
    							         to->Type() == TriggerObject::TriggerCluster  ||
	  						         to->Type() == TriggerObject::TriggerPhoton)){
        double DeltaRAux = MathUtils::DeltaR(leptonsFakeable->At(i)->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
    }
    int lType = 0;
    if(leptonsFakeable->At(i)->ObjType() == kElectron) lType = 1;
    hDFRStudyPresel[0+lType]->Fill(TMath::Min(deltaRleptonHLT,0.999),NNLOWeight->GetVal());

    if(deltaRleptonHLT > 0.2) return;

    double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptonsFakeable->At(i)->Phi()));
    double mTW = TMath::Sqrt(2.0*leptonsFakeable->At(i)->Pt()*caloMet->Pt()*
    				(1.0 - cos(deltaPhiMetLepton)));

    Bool_t isOnlyFake = kTRUE;
    for(UInt_t j=0; j<leptons->GetEntries(); j++) {
      if(leptons->At(j) == leptonsFakeable->At(i)) {
    	isOnlyFake = kFALSE;
        break;
      }
    }
    if(mTW > 30 && caloMet->Pt() > 30.0 && isOnlyFake == kFALSE) {
      hDFRStudyPresel[2+lType]->Fill(TMath::Min(leptonsFakeable->At(i)->Pt(),99.999),NNLOWeight->GetVal());
      hDFRStudyPresel[4+lType]->Fill(TMath::Min(leptonsFakeable->At(i)->AbsEta(),2.499),NNLOWeight->GetVal());
    }
    if(            caloMet->Pt() > 30.0 && isOnlyFake == kFALSE) {
      hDFRStudyPresel[6+lType]->Fill(TMath::Min(mTW,99.999),NNLOWeight->GetVal());
    }

    bool isGenLepton = kFALSE;
    if(fIsData == kFALSE){
      for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
        MCParticle *gen = GenLeptons->At(j);
        if(MathUtils::DeltaR(gen->Mom(), leptonsFakeable->At(i)->Mom()) < 0.1) {
	  isGenLepton = kTRUE;
	  break;
	}
      }
    }
    if(fIsData == kFALSE && (isGenLepton == kFALSE && fSelectGenLeptons == kTRUE )) continue;
    if(fIsData == kFALSE && (isGenLepton == kTRUE  && fSelectGenLeptons == kFALSE)) continue;

    // No need to go further if MET > 20
    bool passMETCut = caloMet->Pt() < 20.0;
    if(fIsData == kFALSE && fSelectGenLeptons == kFALSE) passMETCut = kTRUE;
    if(passMETCut == kFALSE) return;

    bool passmTWCut = mTW < 15.0;
    if(fIsData == kFALSE && fSelectGenLeptons == kFALSE) passmTWCut = kTRUE;
    if(passmTWCut == kFALSE) continue;

    bool isZ = kFALSE;
    if(fIsData == kTRUE || fSelectGenLeptons == kTRUE){
      for(UInt_t j=0; j<leptonsFakeable->GetEntries(); j++){
	if(i == j) continue;
	if(leptonsFakeable->At(i)->ObjType() != leptonsFakeable->At(j)->ObjType()) continue;
	CompositeParticle dilepton;
	dilepton.AddDaughter(leptonsFakeable->At(i));
	dilepton.AddDaughter(leptonsFakeable->At(j));
	hDFRStudyPresel[8+lType]->Fill(dilepton.Mass(),NNLOWeight->GetVal());
	if(TMath::Abs(dilepton.Mass() - 91.1876) < 15.0) {
          isZ = kTRUE;
          break;
	}
      }
    }
    if(isZ == kTRUE) continue;

    int nCentralJets = 0;
    if(fIsData == kTRUE || fSelectGenLeptons == kTRUE){
      for(UInt_t n=0; n<CleanJets->GetEntries(); n++){
	const Jet *jet = CleanJets->At(n);	
	if(TMath::Abs(jet->Eta()) < fEtaJetCut &&
           jet->Pt() > fPtJetCut){
          double dR = MathUtils::DeltaR(leptonsFakeable->At(i)->Mom(), jet->Mom());
          hDFRStudyPresel[10+lType]->Fill(dR,NNLOWeight->GetVal());
          if(dR >= 0.5) nCentralJets++;	
	}
      }
    }
    else {
      nCentralJets = 1;
    }
    hDFRStudyPresel[12+lType]->Fill(nCentralJets,NNLOWeight->GetVal());
    if(nCentralJets == 0) continue;

    int pairType = -1;
    if     (leptonsFakeable->At(i)->ObjType() == kMuon    ) pairType = 0;
    else if(leptonsFakeable->At(i)->ObjType() == kElectron) pairType = 1;

    hDFRStudySel[0+20*pairType]->Fill(leptonsFakeable->At(i)->Pt(),NNLOWeight->GetVal());
    hDFRStudySel[1+20*pairType]->Fill(leptonsFakeable->At(i)->AbsEta(),NNLOWeight->GetVal());
    int npt = -1;
    if     (leptonsFakeable->At(i)->Pt() < 15) npt = 0;
    else if(leptonsFakeable->At(i)->Pt() < 20) npt = 1;
    else if(leptonsFakeable->At(i)->Pt() < 25) npt = 2;
    else if(leptonsFakeable->At(i)->Pt() < 30) npt = 3;
    else if(leptonsFakeable->At(i)->Pt() < 35) npt = 4;
    else                                       npt = 5;
    hDFRStudySel[2+npt+20*pairType]->Fill(leptonsFakeable->At(i)->AbsEta(),NNLOWeight->GetVal());

    int neta = -1;
    if     (leptonsFakeable->At(i)->AbsEta() < 1.000) neta = 0;
    else if(leptonsFakeable->At(i)->AbsEta() < 1.479) neta = 1;
    else if(leptonsFakeable->At(i)->AbsEta() < 2.000) neta = 2;
    else                                              neta = 3;
    hDFRStudySel[8+neta+20*pairType]->Fill(leptonsFakeable->At(i)->Pt(),NNLOWeight->GetVal());

    if(isOnlyFake == kFALSE){
      pairType = pairType + 2;
      hDFRStudySel[0+20*pairType]->Fill(leptonsFakeable->At(i)->Pt(),NNLOWeight->GetVal());
      hDFRStudySel[1+20*pairType]->Fill(leptonsFakeable->At(i)->AbsEta(),NNLOWeight->GetVal());
      hDFRStudySel[2+npt+20*pairType]->Fill(leptonsFakeable->At(i)->AbsEta(),NNLOWeight->GetVal());
      hDFRStudySel[8+neta+20*pairType]->Fill(leptonsFakeable->At(i)->Pt(),NNLOWeight->GetVal());
    }
  } // Loop over fakeable objects
}
//--------------------------------------------------------------------------------------------------
void FRStudy::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  char sb[200];
  sprintf(sb,"hDFRStudyPresel_%d",0); hDFRStudyPresel[0]  = new TH1D(sb,sb,100,0.0,1.0); 
  sprintf(sb,"hDFRStudyPresel_%d",1); hDFRStudyPresel[1]  = new TH1D(sb,sb,100,0.0,1.0); 
  sprintf(sb,"hDFRStudyPresel_%d",2); hDFRStudyPresel[2]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDFRStudyPresel_%d",3); hDFRStudyPresel[3]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDFRStudyPresel_%d",4); hDFRStudyPresel[4]  = new TH1D(sb,sb,100,0.0,2.5); 
  sprintf(sb,"hDFRStudyPresel_%d",5); hDFRStudyPresel[5]  = new TH1D(sb,sb,100,0.0,2.5); 
  sprintf(sb,"hDFRStudyPresel_%d",6); hDFRStudyPresel[6]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDFRStudyPresel_%d",7); hDFRStudyPresel[7]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDFRStudyPresel_%d",8); hDFRStudyPresel[8]  = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDFRStudyPresel_%d",9); hDFRStudyPresel[9]  = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDFRStudyPresel_%d",10);hDFRStudyPresel[10] = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDFRStudyPresel_%d",11);hDFRStudyPresel[11] = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDFRStudyPresel_%d",12);hDFRStudyPresel[12] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDFRStudyPresel_%d",13);hDFRStudyPresel[13] = new TH1D(sb,sb,10,-0.5,9.5);

  for(int i=0; i<14; i++){
    AddOutput(hDFRStudyPresel[i]);
  }

  for(int j=0; j<4; j++){
    int ind = 20 * j;
    sprintf(sb,"hDFRStudySel_%d",ind+ 0); hDFRStudySel[ind+ 0] = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDFRStudySel_%d",ind+ 1); hDFRStudySel[ind+ 1] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 2); hDFRStudySel[ind+ 2] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 3); hDFRStudySel[ind+ 3] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 4); hDFRStudySel[ind+ 4] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 5); hDFRStudySel[ind+ 5] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 6); hDFRStudySel[ind+ 6] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 7); hDFRStudySel[ind+ 7] = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDFRStudySel_%d",ind+ 8); hDFRStudySel[ind+ 8] = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDFRStudySel_%d",ind+ 9); hDFRStudySel[ind+ 9] = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDFRStudySel_%d",ind+10); hDFRStudySel[ind+10] = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDFRStudySel_%d",ind+11); hDFRStudySel[ind+11] = new TH1D(sb,sb,100,0.0,100.);
  }

  for(int j=0; j<4; j++){
    for(int i=0; i<12; i++){
      AddOutput(hDFRStudySel[i+j*20]);
    }
  }

}

//--------------------------------------------------------------------------------------------------
void FRStudy::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void FRStudy::Terminate()
{
  // Run finishing code on the client computer
}
