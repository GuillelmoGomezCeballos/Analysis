// $Id: LowEvtSelMod.cc,v 1.8 2012/05/28 06:47:07 ceballos Exp $

#include "Analysis/SelMods/interface/LowEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

using namespace mithep;
ClassImp(mithep::LowEvtSelMod)

//--------------------------------------------------------------------------------------------------
LowEvtSelMod::LowEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsFastSim(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fMuonName(Names::gkMuonBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMuons(0),
  fCaloMetName(Names::gkCaloMetBrn),
  fCaloMet(0),
  fJetScaleSyst(0.0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void LowEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void LowEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  //Obtain all the good objects from the event cleaning module
  JetOArr *CleanJets            = GetObjThisEvt<JetOArr>(fCleanJetsName);
  MetOArr *CleanMet             = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet            = CleanMet->At(0);

  LoadBranch(fMuonName);
  ObjArray<Muon> *SoftMuons = new ObjArray<Muon>;
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  int NMuons[4] = {0, 0, 0, 0};
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);

    bool isCleanMuon = MuonTools::PassSoftMuonCut(mu, fVertices, 0.2);

    if     (mu->HasGlobalTrk()){
      hDLowPresel[0]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[4]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
      if(isCleanMuon){
    	hDLowPresel[8] ->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
    	hDLowPresel[12]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
    	NMuons[0]++;
      }
    }
    else if(mu->HasTrackerTrk()){
      hDLowPresel[2]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[6]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
      if(isCleanMuon){
    	hDLowPresel[10]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
    	hDLowPresel[14]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
    	NMuons[1]++;
      }
    }
    else if(mu->HasStandaloneTrk()){
      hDLowPresel[1]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[5]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
      if(isCleanMuon){
    	hDLowPresel[9] ->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
    	hDLowPresel[13]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
    	NMuons[2]++;
      }
    }
    else if(mu->IsCaloMuon()){
      hDLowPresel[3]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[7]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
      if(isCleanMuon){
    	hDLowPresel[11]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
    	hDLowPresel[15]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
    	NMuons[3]++;
      }
    }
    else {
      printf("LowSel, it can not be!\n");
    }

    if(isCleanMuon){
      hDLowPresel[16]->Fill(TMath::Min(mu->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[17]->Fill(TMath::Min(mu->AbsEta(),2.499),NNLOWeight->GetVal());
    }

    if(isCleanMuon) SoftMuons->Add(mu);
  }

  SoftMuons->Sort();

  if(NMuons[0]+NMuons[1]+NMuons[2]+NMuons[3] >= 2){
    hDLowPresel[18]->Fill(TMath::Min((double)NMuons[0],4.499),NNLOWeight->GetVal());
    hDLowPresel[19]->Fill(TMath::Min((double)NMuons[1],4.499),NNLOWeight->GetVal());
    hDLowPresel[20]->Fill(TMath::Min((double)NMuons[2],4.499),NNLOWeight->GetVal());
    hDLowPresel[21]->Fill(TMath::Min((double)NMuons[3],4.499),NNLOWeight->GetVal());
  }

  hDLowPresel[22]->Fill(TMath::Min((double)SoftMuons->GetEntries(),9.499),NNLOWeight->GetVal());
  if(SoftMuons->GetEntries() >= 2 &&
     SoftMuons->At(0)->BestTrk()  && SoftMuons->At(1)->BestTrk()){
    // Sort and count the number of central Jets for vetoing
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	 CleanJets->At(i)->Pt()*(1.0+fJetScaleSyst) > fPtJetCut){
        nCentralJets++;
      }
    }
    CompositeParticle dilepton;
    dilepton.AddDaughter(SoftMuons->At(0));
    dilepton.AddDaughter(SoftMuons->At(1));

    CompositeParticle dileptonTrk;
    MCParticle *pm0 = new MCParticle();
    MCParticle *pm1 = new MCParticle();
    
    pm0->SetMom(SoftMuons->At(0)->BestTrk()->Px(),
                SoftMuons->At(0)->BestTrk()->Py(),
                SoftMuons->At(0)->BestTrk()->Pz(),
                sqrt(SoftMuons->At(0)->BestTrk()->P()*SoftMuons->At(0)->BestTrk()->P()+
		     0.105658369*0.105658369)); 			     
    pm1->SetMom(SoftMuons->At(1)->BestTrk()->Px(),
                SoftMuons->At(1)->BestTrk()->Py(),
                SoftMuons->At(1)->BestTrk()->Pz(),
                sqrt(SoftMuons->At(1)->BestTrk()->P()*SoftMuons->At(1)->BestTrk()->P()+
		     0.105658369*0.105658369)); 			      
    dileptonTrk.AddDaughter(pm0);
    dileptonTrk.AddDaughter(pm1);

    if(dilepton.Charge() == 0){
      hDLowPresel[23]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      hDLowPresel[24]->Fill(TMath::Min((double)nCentralJets,9.499),NNLOWeight->GetVal());  
 
      double deltaPhiLeptons = fabs(MathUtils::DeltaPhi(SoftMuons->At(0)->Phi(), 
        					   SoftMuons->At(1)->Phi()))* 180./TMath::Pi();
      int pairType = -1;
      if     (SoftMuons->At(0)->HasGlobalTrk() && SoftMuons->At(1)->HasGlobalTrk()){
        pairType = 0;
      }
      else if((SoftMuons->At(0)->HasGlobalTrk()  && SoftMuons->At(1)->HasTrackerTrk()) ||
              (SoftMuons->At(0)->HasTrackerTrk() && SoftMuons->At(1)->HasGlobalTrk())){
        pairType = 1;
      }
      else {
        pairType = 2;
      }
      if(pairType != -1){
        hDLowSel[ 0+100*pairType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
        hDLowSel[ 1+100*pairType]->Fill(dilepton.Mass(),NNLOWeight->GetVal());
        hDLowSel[ 2+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	if     (TMath::Max(SoftMuons->At(0)->AbsEta(),SoftMuons->At(1)->AbsEta()) < 0.5){
	  hDLowSel[ 3+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	}
	else if(TMath::Max(SoftMuons->At(0)->AbsEta(),SoftMuons->At(1)->AbsEta()) < 1.0){
	  hDLowSel[ 4+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	}
	else if(TMath::Max(SoftMuons->At(0)->AbsEta(),SoftMuons->At(1)->AbsEta()) < 1.5){
	  hDLowSel[ 5+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	}
	else if(TMath::Max(SoftMuons->At(0)->AbsEta(),SoftMuons->At(1)->AbsEta()) < 2.0){
	  hDLowSel[ 6+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	}
	else {
	  hDLowSel[ 7+100*pairType]->Fill(dileptonTrk.Mass(),NNLOWeight->GetVal());
	}
        hDLowSel[ 8+100*pairType]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
        hDLowSel[ 9+100*pairType]->Fill(TMath::Min(SoftMuons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        hDLowSel[10+100*pairType]->Fill(TMath::Min(SoftMuons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
        hDLowSel[11+100*pairType]->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
      }
    }
    delete pm0;
    delete pm1;
  }

  delete SoftMuons;
}
//--------------------------------------------------------------------------------------------------
void LowEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMuonName,    fMuons);

  char sb[200];

  sprintf(sb,"hDLowPresel_%d",0);  hDLowPresel[0]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",1);  hDLowPresel[1]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",2);  hDLowPresel[2]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",3);  hDLowPresel[3]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",4);  hDLowPresel[4]  = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",5);  hDLowPresel[5]  = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",6);  hDLowPresel[6]  = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",7);  hDLowPresel[7]  = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",8);  hDLowPresel[8]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",9);  hDLowPresel[9]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",10); hDLowPresel[10] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",11); hDLowPresel[11] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",12); hDLowPresel[12] = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",13); hDLowPresel[13] = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",14); hDLowPresel[14] = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",15); hDLowPresel[15] = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",16); hDLowPresel[16] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",17); hDLowPresel[17] = new TH1D(sb,sb,50,0.0,2.5);
  sprintf(sb,"hDLowPresel_%d",18); hDLowPresel[18] = new TH1D(sb,sb,5,-0.5,4.5);
  sprintf(sb,"hDLowPresel_%d",19); hDLowPresel[19] = new TH1D(sb,sb,5,-0.5,4.5);
  sprintf(sb,"hDLowPresel_%d",20); hDLowPresel[20] = new TH1D(sb,sb,5,-0.5,4.5);
  sprintf(sb,"hDLowPresel_%d",21); hDLowPresel[21] = new TH1D(sb,sb,5,-0.5,4.5);
  sprintf(sb,"hDLowPresel_%d",22); hDLowPresel[22] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLowPresel_%d",23); hDLowPresel[23] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDLowPresel_%d",24); hDLowPresel[24] = new TH1D(sb,sb,10,-0.5,9.5);

  for(int i=0; i<25; i++){
    AddOutput(hDLowPresel[i]);
  }
  
  for(int i=0; i<3; i++){
    int ind = 100*i;
    sprintf(sb,"hDLowSel_%d", 0+ind); hDLowSel[ 0+ind] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDLowSel_%d", 1+ind); hDLowSel[ 1+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 2+ind); hDLowSel[ 2+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 3+ind); hDLowSel[ 3+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 4+ind); hDLowSel[ 4+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 5+ind); hDLowSel[ 5+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 6+ind); hDLowSel[ 6+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 7+ind); hDLowSel[ 7+ind] = new TH1D(sb,sb,2000,0.0,20.0);
    sprintf(sb,"hDLowSel_%d", 8+ind); hDLowSel[ 8+ind] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDLowSel_%d", 9+ind); hDLowSel[ 9+ind] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDLowSel_%d",10+ind); hDLowSel[10+ind] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDLowSel_%d",11+ind); hDLowSel[11+ind] = new TH1D(sb,sb,90,0.0,180.0);
  }

  for(int j=0; j<12; j++){
    for(int i=0; i<3; i++){
      AddOutput(hDLowSel[j+100*i]);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void LowEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void LowEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
