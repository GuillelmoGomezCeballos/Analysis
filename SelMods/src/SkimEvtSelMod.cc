// $Id: SkimEvtSelMod.cc,v 1.6 2012/04/24 11:26:27 ceballos Exp $

#include "Ana/SelMods/interface/SkimEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

using namespace mithep;
ClassImp(mithep::SkimEvtSelMod)

//--------------------------------------------------------------------------------------------------
SkimEvtSelMod::SkimEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void SkimEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void SkimEvtSelMod::Process()
{
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    printf("SkimEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  ParticleOArr *leptonsFakeable = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  ParticleOArr *leptons         = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);

  hDskimPresel[0]->Fill(TMath::Min((double)leptonsFakeable->GetEntries(),9.499),NNLOWeight->GetVal());
  hDskimPresel[1]->Fill(TMath::Min((double)leptons->GetEntries(),9.499),NNLOWeight->GetVal());

  for (UInt_t i=0; i<leptonsFakeable->GetEntries(); ++i) {
    const Particle *p = leptonsFakeable->At(i);
    int leptonType = 0;
    if(p->ObjType() == kElectron) leptonType = 1;
    double deltaRleptonHLT = 999.;
    for(UInt_t i=0;i<objs->GetEntries();++i) {
      const TriggerObject* to = objs->At(i);
      if     (p->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || 
    			                to->Type() == TriggerObject::TriggerCluster)){
    	double DeltaRAux = MathUtils::DeltaR(p->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
      else if(p->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
    					    to->Type() == TriggerObject::TriggerCluster  ||
    					    to->Type() == TriggerObject::TriggerPhoton)){
    	double DeltaRAux = MathUtils::DeltaR(p->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
    } // Loop over HLT objects
    hDskimSel[0+10*leptonType]->Fill(TMath::Min(deltaRleptonHLT,0.999),NNLOWeight->GetVal());
    hDskimSel[1+10*leptonType]->Fill(TMath::Min(p->Pt(),99.999),NNLOWeight->GetVal());
    hDskimSel[2+10*leptonType]->Fill(TMath::Min(p->AbsEta(),2.5),NNLOWeight->GetVal());
    if(deltaRleptonHLT < 0.1){
      hDskimSel[3+10*leptonType]->Fill(TMath::Min(p->Pt(),99.999),NNLOWeight->GetVal());
      hDskimSel[4+10*leptonType]->Fill(TMath::Min(p->AbsEta(),2.5),NNLOWeight->GetVal());
    }
  }

  for (UInt_t i=0; i<leptons->GetEntries(); ++i) {
    const Particle *p = leptons->At(i);
    int leptonType = 2;
    if(p->ObjType() == kElectron) leptonType = 3;
    double deltaRleptonHLT = 999.;
    for(UInt_t i=0;i<objs->GetEntries();++i) {
      const TriggerObject* to = objs->At(i);
      if     (p->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || 
    			                to->Type() == TriggerObject::TriggerCluster)){
    	double DeltaRAux = MathUtils::DeltaR(p->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
      else if(p->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
    					    to->Type() == TriggerObject::TriggerCluster  ||
    					    to->Type() == TriggerObject::TriggerPhoton)){
    	double DeltaRAux = MathUtils::DeltaR(p->Mom(), to->Mom());
    	if(DeltaRAux < deltaRleptonHLT) deltaRleptonHLT = DeltaRAux;
      }
    } // Loop over HLT objects
    hDskimSel[0+10*leptonType]->Fill(TMath::Min(deltaRleptonHLT,0.999),NNLOWeight->GetVal());
    hDskimSel[1+10*leptonType]->Fill(TMath::Min(p->Pt(),99.999),NNLOWeight->GetVal());
    hDskimSel[2+10*leptonType]->Fill(TMath::Min(p->AbsEta(),2.5),NNLOWeight->GetVal());
    if(deltaRleptonHLT < 0.1){
      hDskimSel[3+10*leptonType]->Fill(TMath::Min(p->Pt(),99.999),NNLOWeight->GetVal());
      hDskimSel[4+10*leptonType]->Fill(TMath::Min(p->AbsEta(),2.5),NNLOWeight->GetVal());
    }
  }
}
//--------------------------------------------------------------------------------------------------
void SkimEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  char sb[200];
  for(int j=0; j<4; j++){
    int ind = 10 * j;
    sprintf(sb,"hDskimSel_%d",ind+0);  hDskimSel[ind+0]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDskimSel_%d",ind+1);  hDskimSel[ind+1]  = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDskimSel_%d",ind+2);  hDskimSel[ind+2]  = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDskimSel_%d",ind+3);  hDskimSel[ind+3]  = new TH1D(sb,sb,100,0.0,100.);
    sprintf(sb,"hDskimSel_%d",ind+4);  hDskimSel[ind+4]  = new TH1D(sb,sb,100,0.0,2.5);
  }

  for(int i=0; i<5; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDskimSel[i+j*10]);
    }
  }

  sprintf(sb,"hDskimPresel_%d",0);  hDskimPresel[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDskimPresel_%d",1);  hDskimPresel[1]  = new TH1D(sb,sb,10,-0.5,9.5); 

  for(int i=0; i<2; i++){
    AddOutput(hDskimPresel[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void SkimEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void SkimEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
