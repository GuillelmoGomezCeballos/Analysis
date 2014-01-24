// $Id: ZttEvtSelMod.cc,v 1.10 2012/04/24 11:26:27 ceballos Exp $

#include "Ana/SelMods/interface/ZttEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"

using namespace mithep;
ClassImp(mithep::ZttEvtSelMod)

//--------------------------------------------------------------------------------------------------
ZttEvtSelMod::ZttEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fJetScaleSyst(0.0),
  fMetName(Names::gkCaloMetBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fCaloJetName0("AKt5Jets"),
  fCaloJet0(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ZttEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ZttEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  PFTauOArr  *CleanPFTaus      = GetObjThisEvt<PFTauOArr>(ModNames::gkCleanPFTausName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    printf("ZttEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);

  MetOArr *GenMetColl          = GetObjThisEvt<MetOArr>(ModNames::gkMCMETName);
  MCParticleOArr *GenLeptons   = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
  MCParticleOArr *GenTaus      = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCTausName);
  if(GenMetColl->GetEntries() > 0 && GenMetColl->At(0)->Pt() > 0.0){
    const Met *GenMet = GenMetColl->At(0);
    hDZttMCSel[ 0]->Fill(TMath::Min((double)GenLeptons->GetEntries(),9.499),NNLOWeight->GetVal());
    hDZttMCSel[ 1]->Fill(TMath::Min((double)GenTaus->GetEntries(),9.499),NNLOWeight->GetVal());
    hDZttMCSel[ 2]->Fill(TMath::Min((double)(GenLeptons->GetEntries()+GenTaus->GetEntries()),9.499),NNLOWeight->GetVal());

    if     (GenLeptons->GetEntries() >= 2 && 
            GenLeptons->At(0)->Pt() > 10.0 &&  GenLeptons->At(1)->Pt() > 10.0){
      hDZttMCSel[ 3]->Fill(0.0,NNLOWeight->GetVal());
      hDZttMCSel[ 4]->Fill(TMath::Min(TMath::Max((caloMet->Pt()-GenMet->Pt())/GenMet->Pt(),-1.999),1.999),NNLOWeight->GetVal());
      hDZttMCSel[ 5]->Fill(fabs(MathUtils::DeltaPhi(caloMet->Phi(), GenMet->Phi())) * 180./TMath::Pi(),NNLOWeight->GetVal());
      Tau *tau0 = new Tau(); tau0->SetMom(GenLeptons->At(0)->Px(),GenLeptons->At(0)->Py(),GenLeptons->At(0)->Pz(),GenLeptons->At(0)->E());
      Tau *tau1 = new Tau(); tau1->SetMom(GenLeptons->At(1)->Px(),GenLeptons->At(1)->Py(),GenLeptons->At(1)->Pz(),GenLeptons->At(1)->E());
      ParticleOArr *fColOut = new ParticleOArr();
      fColOut->Add(tau0);
      fColOut->Add(tau1);
      DiTauSystem ditau(fColOut->At(0), fColOut->At(1), GenMet);
      hDZttMCSel[6]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZttMCSel[6]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
        hDZttMCSel[7]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
      }
      delete tau0;
      delete tau1;
      delete fColOut;
    }
    else if(GenLeptons->GetEntries() == 1 && GenTaus->GetEntries() >= 1 && 
            GenLeptons->At(0)->Pt() > 10.0 &&  GenTaus->At(0)->Pt() > 10.0){
      hDZttMCSel[ 3]->Fill(1.0,NNLOWeight->GetVal());
      hDZttMCSel[ 8]->Fill(TMath::Min(TMath::Max((caloMet->Pt()-GenMet->Pt())/GenMet->Pt(),-1.999),1.999),NNLOWeight->GetVal());
      hDZttMCSel[ 9]->Fill(fabs(MathUtils::DeltaPhi(caloMet->Phi(), GenMet->Phi())) * 180./TMath::Pi(),NNLOWeight->GetVal());
      Tau *tau0 = new Tau(); tau0->SetMom(GenLeptons->At(0)->Px(),GenLeptons->At(0)->Py(),GenLeptons->At(0)->Pz(),GenLeptons->At(0)->E());
      Tau *tau1 = new Tau(); tau1->SetMom(GenTaus->At(0)->Px(),   GenTaus->At(0)->Py(),   GenTaus->At(0)->Pz(),   GenTaus->At(0)->E());
      ParticleOArr *fColOut = new ParticleOArr();
      fColOut->Add(tau0);
      fColOut->Add(tau1);
      DiTauSystem ditau(fColOut->At(0), fColOut->At(1), GenMet);
      hDZttMCSel[10]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZttMCSel[10]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
        hDZttMCSel[11]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
      }
      delete tau0;
      delete tau1;
      delete fColOut;
    }
    else if(GenLeptons->GetEntries() == 0 && GenTaus->GetEntries() >= 2 && 
            GenTaus->At(0)->Pt() > 10.0 &&  GenTaus->At(1)->Pt() > 10.0){
      hDZttMCSel[ 3]->Fill(2.0,NNLOWeight->GetVal());
      hDZttMCSel[12]->Fill(TMath::Min(TMath::Max((caloMet->Pt()-GenMet->Pt())/GenMet->Pt(),-1.999),1.999),NNLOWeight->GetVal());
      hDZttMCSel[13]->Fill(fabs(MathUtils::DeltaPhi(caloMet->Phi(), GenMet->Phi())) * 180./TMath::Pi(),NNLOWeight->GetVal());
      Tau *tau0 = new Tau(); tau0->SetMom(GenTaus->At(0)->Px(),   GenTaus->At(0)->Py(),   GenTaus->At(0)->Pz(),   GenTaus->At(0)->E());
      Tau *tau1 = new Tau(); tau1->SetMom(GenTaus->At(1)->Px(),   GenTaus->At(1)->Py(),   GenTaus->At(1)->Pz(),   GenTaus->At(1)->E());
      ParticleOArr *fColOut = new ParticleOArr();
      fColOut->Add(tau0);
      fColOut->Add(tau1);
      DiTauSystem ditau(fColOut->At(0), fColOut->At(1), GenMet);
      hDZttMCSel[14]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZttMCSel[14]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
        hDZttMCSel[15]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
      }
      delete tau0;
      delete tau1;
      delete fColOut;
    }
    else {
      hDZttMCSel[ 3]->Fill(3.0,NNLOWeight->GetVal());
    }
  }


  double zAverage = 0.0;

  // Make lepton vector from muons and electrons
  for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
    zAverage = zAverage + CleanMuons->At(j)->BestTrk()->Z0();
  }
  for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
    zAverage = zAverage + CleanElectrons->At(j)->BestTrk()->Z0();
  }

  // Computing Z average (our primary vertex)
  if(leptons->GetEntries() > 0) zAverage = zAverage / leptons->GetEntries();

  // Z -> tautau -> emuX selection
  if (CleanMuons->GetEntries() == 1 && CleanElectrons->GetEntries() == 1 &&
      leptons->GetEntries() == 2 && leptons->At(0)->Charge()+leptons->At(1)->Charge() == 0 && 
      leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10){
    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));

    LoadBranch(fCaloJetName0);
    vector<Jet*> sortedJets;
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	 CleanJets->At(i)->Pt() > fPtJetCut){
        nCentralJets++;
    	Jet* jet_f = new Jet(CleanJets->At(i)->Px(),
    			     CleanJets->At(i)->Py(),
    			     CleanJets->At(i)->Pz(),
    			     CleanJets->At(i)->E() );
        jet_f->SetMatchedMCFlavor(CleanJets->At(i)->MatchedMCFlavor());
        jet_f->SetCombinedSecondaryVertexBJetTagsDisc(CleanJets->At(i)->CombinedSecondaryVertexBJetTagsDisc());
        jet_f->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJets->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
        jet_f->SetJetProbabilityBJetTagsDisc(CleanJets->At(i)->JetProbabilityBJetTagsDisc());
        jet_f->SetJetBProbabilityBJetTagsDisc(CleanJets->At(i)->JetBProbabilityBJetTagsDisc());
        jet_f->SetTrackCountingHighEffBJetTagsDisc(CleanJets->At(i)->TrackCountingHighEffBJetTagsDisc());
        jet_f->SetTrackCountingHighPurBJetTagsDisc(CleanJets->At(i)->TrackCountingHighPurBJetTagsDisc());
    	sortedJets.push_back(jet_f);
      }
    }
    for(UInt_t i=0; i<sortedJets.size(); i++){
      for(UInt_t j=i+1; j<sortedJets.size(); j++){
        if(sortedJets[i]->Pt() < sortedJets[j]->Pt()) {
          //swap i and j
    	  Jet* tempjet = sortedJets[i];
    	  sortedJets[i] = sortedJets[j];
    	  sortedJets[j] = tempjet;    
        }
      }
    }

    DiTauSystem ditau(leptons->At(0), leptons->At(1), caloMet);

    double deltaPhiMetLepton[2] = {fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi())),
    				   fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(1)->Phi()))};
    double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
      deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
    minDeltaPhiMetLepton = minDeltaPhiMetLepton * 180./TMath::Pi();

    double mTW[2] = {TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[0]))),
        	     TMath::Sqrt(2.0*leptons->At(1)->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[1])))};

    double deltaPhiLeptons = fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), 
    						 leptons->At(1)->Phi()))* 180./TMath::Pi();

    double deltaPhiDileptonMet = 180.-fabs(MathUtils::DeltaPhi(caloMet->Phi(), 
    						          dilepton.Phi()))* 180./TMath::Pi();

    Bool_t cuts[3] = {TMath::Min(mTW[0],mTW[1]) < 20.0, dilepton.Mass() > 12.0 && dilepton.Mass() < 80.0,
                      TMath::Min(TMath::Abs(cos(deltaPhiMetLepton[0])),TMath::Abs(cos(deltaPhiMetLepton[1]))) > 0.5};

    if(           cuts[1] && cuts[2]){
      hDZemSel[0]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
    }
    if(cuts[0] &&            cuts[2]){
      hDZemSel[1]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    }
    if(cuts[0] && cuts[1]           ){
      hDZemSel[2]->Fill(TMath::Min(TMath::Abs(cos(deltaPhiMetLepton[0])),TMath::Abs(cos(deltaPhiMetLepton[1]))),NNLOWeight->GetVal());
    }
    if(cuts[0] && cuts[1] && cuts[2]){
      hDZemSel[3]->Fill(TMath::Min((double)nCentralJets,9.499),NNLOWeight->GetVal());
      hDZemSel[4]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZemSel[4]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
        hDZemSel[5]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
      } 
      hDZemSel[6]->Fill(minDeltaPhiMetLepton,NNLOWeight->GetVal());
      hDZemSel[7]->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
      hDZemSel[8]->Fill(deltaPhiDileptonMet,NNLOWeight->GetVal());
      hDZemSel[9]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      hDZemSel[10]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
      if(nCentralJets > 0){
        hDZemSel[11]->Fill(TMath::Min(TMath::Max(sortedJets[0]->CombinedSecondaryVertexBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
        hDZemSel[12]->Fill(TMath::Min(TMath::Max(sortedJets[0]->CombinedSecondaryVertexMVABJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
        hDZemSel[13]->Fill(TMath::Min(TMath::Max(sortedJets[0]->JetProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
        hDZemSel[14]->Fill(TMath::Min(TMath::Max(sortedJets[0]->JetBProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
        hDZemSel[15]->Fill(TMath::Min(TMath::Max(sortedJets[0]->TrackCountingHighEffBJetTagsDisc(),-19.99999),19.9999),NNLOWeight->GetVal());
        hDZemSel[16]->Fill(TMath::Min(TMath::Max(sortedJets[0]->TrackCountingHighPurBJetTagsDisc(),-19.99999),19.9999),NNLOWeight->GetVal());
      }
    }
  }

  // Minimun Pt, Nleptons == 1 & Ntau == 1 requirements
  if (leptons->GetEntries() == 1 && CleanPFTaus->GetEntries() == 1 &&
      leptons->At(0)->Pt() > 20 && CleanPFTaus->At(0)->Pt() > 20){

    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(CleanPFTaus->At(0));

    // Sort and count the number of central Jets for vetoing
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      const Jet *jet = CleanJets->At(i);        
      if(TMath::Abs(jet->Eta()) >= 2.5) continue;
      if(jet->Pt() <= fPtJetCut) continue;
      Bool_t isTauOverlap = kFALSE;
      for (UInt_t j=0; j<CleanPFTaus->GetEntries(); ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanPFTaus->At(j)->Mom(),jet->Mom());  
        if (deltaR < 0.3) {
          isTauOverlap = kTRUE;
          break;
        }
      }
      if (isTauOverlap) continue;
 
      nCentralJets++;
    }

    int pairType = -1;
    if     (leptons->At(0)->ObjType() == kMuon    ) pairType = 0;
    else if(leptons->At(0)->ObjType() == kElectron) pairType = 1;
    else {
      cout << "Hey, this is not possible, leptonType: "
    	   << leptons->At(0)->ObjType() << endl;
    }

    double deltaRleptonHLT[1] = {999.};
    // HLT study
    if(objs){
      Int_t ents=objs->GetEntries();
      for(Int_t i=0;i<ents;++i) {
         const TriggerObject* to = objs->At(i);

         TString trName = to->TrigName();
         if(trName.Contains("HLT_IsoMu24_eta2p1_v"))								       hDZttHLT[0+100*pairType]->Fill(0.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Ele27_WP80_v"))								       hDZttHLT[0+100*pairType]->Fill(1.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"))				       hDZttHLT[0+100*pairType]->Fill(2.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"))				       hDZttHLT[0+100*pairType]->Fill(3.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Mu17_Mu8_v"))  								       hDZttHLT[0+100*pairType]->Fill(4.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Mu17_TkMu8_v"))								       hDZttHLT[0+100*pairType]->Fill(5.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) hDZttHLT[0+100*pairType]->Fill(6.,NNLOWeight->GetVal());

        if     (leptons->At(0)->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || 
	                                               to->Type() == TriggerObject::TriggerCluster)){
          double DeltaRAux = MathUtils::DeltaR(leptons->At(0)->Phi(), leptons->At(0)->Eta(),
                			       to->Phi(),to->Eta());
    	  if(DeltaRAux < deltaRleptonHLT[0]) deltaRleptonHLT[0] = DeltaRAux;
        }
        else if(leptons->At(0)->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
	                                                   to->Type() == TriggerObject::TriggerCluster  ||
							   to->Type() == TriggerObject::TriggerPhoton)
							   ){
          double DeltaRAux = MathUtils::DeltaR(leptons->At(0)->Phi(), leptons->At(0)->Eta(),
                			       to->Phi(),to->Eta());
    	  if(DeltaRAux < deltaRleptonHLT[0]) deltaRleptonHLT[0] = DeltaRAux;
        }
      } // Loop over HLT objects
      if     (leptons->At(0)->ObjType() == kMuon)
        hDZttHLT[1+100*pairType]->Fill(TMath::Min(deltaRleptonHLT[0],0.999),NNLOWeight->GetVal());
      else if(leptons->At(0)->ObjType() == kElectron)
        hDZttHLT[2+100*pairType]->Fill(TMath::Min(deltaRleptonHLT[0],0.999),NNLOWeight->GetVal());
      if      (deltaRleptonHLT[0] < 0.1)
        hDZttHLT[3+100*pairType]->Fill(0.0,NNLOWeight->GetVal());
      else
        hDZttHLT[3+100*pairType]->Fill(1.0,NNLOWeight->GetVal());
    } // if !objs

    int nCharge = CleanPFTaus->At(0)->Charge() + leptons->At(0)->Charge();
    hDZttSel[ 0+100*pairType]->Fill(nCharge,NNLOWeight->GetVal());
    if(nCharge != 0) hDZttSel[ 1+100*pairType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    // Ncharge == 0
    if(nCharge == 0){
      DiTauSystem ditau(leptons->At(0), CleanPFTaus->At(0), caloMet);
      double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi()));
      double mTW = TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
    		               (1.0 - cos(deltaPhiMetLepton)));
      double deltaPhiTauLepton = fabs(MathUtils::DeltaPhi(CleanPFTaus->At(0)->Phi(), leptons->At(0)->Phi()));
      double deltaEtaTauLepton = TMath::Abs(CleanPFTaus->At(0)->Eta()-leptons->At(0)->Eta());
      Bool_t cuts[1] = {mTW < 40};
      hDZttSel[ 3+100*pairType]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
      if(cuts[0]){
      	hDZttSel[ 8+100*pairType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
        if(TMath::Abs(dilepton.Mass()-91.1876) > 15.0){
	  hDZttSel[ 2+100*pairType]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
          hDZttSel[ 2+100*pairType]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
          hDZttSel[ 4+100*pairType]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
      	  hDZttSel[ 5+100*pairType]->Fill((double)nCentralJets,NNLOWeight->GetVal());
      	  hDZttSel[ 6+100*pairType]->Fill(deltaPhiTauLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
          hDZttSel[ 7+100*pairType]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
      	  hDZttSel[ 9+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      	  hDZttSel[10+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      	  hDZttSel[11+100*pairType]->Fill(TMath::Min(CleanPFTaus->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      	  hDZttSel[12+100*pairType]->Fill(TMath::Min(caloMet->MetSig(),19.999),NNLOWeight->GetVal());
      	  hDZttSel[13+100*pairType]->Fill(TMath::Min(caloMet->SumEt(),799.999),NNLOWeight->GetVal());
      	  hDZttSel[14+100*pairType]->Fill(TMath::Min(deltaEtaTauLepton,4.999),NNLOWeight->GetVal());
      	  hDZttSel[15+100*pairType]->Fill(TMath::Abs(leptons->At(0)->Eta()),NNLOWeight->GetVal());
      	  hDZttSel[16+100*pairType]->Fill(TMath::Abs(CleanPFTaus->At(0)->Eta()),NNLOWeight->GetVal());
        }
      }
    } // Ncharge == 0
  } // Minimun requirements for tautau->h mu/e
}
//--------------------------------------------------------------------------------------------------
void ZttEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fCaloJetName0, fCaloJet0);

  char sb[200];
  sprintf(sb,"hDZemSel_%d", 0); hDZemSel[ 0] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZemSel_%d", 1); hDZemSel[ 1] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZemSel_%d", 2); hDZemSel[ 2] = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDZemSel_%d", 3); hDZemSel[ 3] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZemSel_%d", 4); hDZemSel[ 4] = new TH1D(sb,sb,300,-1.0,2.0);
  sprintf(sb,"hDZemSel_%d", 5); hDZemSel[ 5] = new TH1D(sb,sb,200,0.0,400.0);
  sprintf(sb,"hDZemSel_%d", 6); hDZemSel[ 6] = new TH1D(sb,sb,180,0.0,180.0);
  sprintf(sb,"hDZemSel_%d", 7); hDZemSel[ 7] = new TH1D(sb,sb,180,0.0,180.0);
  sprintf(sb,"hDZemSel_%d", 8); hDZemSel[ 8] = new TH1D(sb,sb,180,0.0,180.0);
  sprintf(sb,"hDZemSel_%d", 9); hDZemSel[ 9] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZemSel_%d",10); hDZemSel[10] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDZemSel_%d",11); hDZemSel[11] = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDZemSel_%d",12); hDZemSel[12] = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDZemSel_%d",13); hDZemSel[13] = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDZemSel_%d",14); hDZemSel[14] = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDZemSel_%d",15); hDZemSel[15] = new TH1D(sb,sb,100,-20.0,20.0);
  sprintf(sb,"hDZemSel_%d",16); hDZemSel[16] = new TH1D(sb,sb,100,-20.0,20.0);

  for(int i=0; i<17; i++){
    AddOutput(hDZemSel[i]);
  }

  sprintf(sb,"hDZttMCSel_%d", 0); hDZttMCSel[ 0] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZttMCSel_%d", 1); hDZttMCSel[ 1] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZttMCSel_%d", 2); hDZttMCSel[ 2] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZttMCSel_%d", 3); hDZttMCSel[ 3] = new TH1D(sb,sb,4,-0.5,3.5);
  sprintf(sb,"hDZttMCSel_%d", 4); hDZttMCSel[ 4] = new TH1D(sb,sb,200,-2.0,2.0);
  sprintf(sb,"hDZttMCSel_%d", 5); hDZttMCSel[ 5] = new TH1D(sb,sb,90,0.0,180.0);
  sprintf(sb,"hDZttMCSel_%d", 6); hDZttMCSel[ 6] = new TH1D(sb,sb,300,-1.0,2.0);
  sprintf(sb,"hDZttMCSel_%d", 7); hDZttMCSel[ 7] = new TH1D(sb,sb,200,0.0,400.0);
  sprintf(sb,"hDZttMCSel_%d", 8); hDZttMCSel[ 8] = new TH1D(sb,sb,200,-2.0,2.0);
  sprintf(sb,"hDZttMCSel_%d", 9); hDZttMCSel[ 9] = new TH1D(sb,sb,90,0.0,180.0);
  sprintf(sb,"hDZttMCSel_%d",10); hDZttMCSel[10] = new TH1D(sb,sb,300,-1.0,2.0);
  sprintf(sb,"hDZttMCSel_%d",11); hDZttMCSel[11] = new TH1D(sb,sb,200,0.0,400.0);
  sprintf(sb,"hDZttMCSel_%d",12); hDZttMCSel[12] = new TH1D(sb,sb,200,-2.0,2.0);
  sprintf(sb,"hDZttMCSel_%d",13); hDZttMCSel[13] = new TH1D(sb,sb,90,0.0,180.0);
  sprintf(sb,"hDZttMCSel_%d",14); hDZttMCSel[14] = new TH1D(sb,sb,300,-1.0,2.0);
  sprintf(sb,"hDZttMCSel_%d",15); hDZttMCSel[15] = new TH1D(sb,sb,200,0.0,400.0);

  for(int i=0; i<16; i++){
    AddOutput(hDZttMCSel[i]);
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZttSel_%d",ind+ 0); hDZttSel[ind+ 0] = new TH1D(sb,sb,5,-2.5,2.5);
    sprintf(sb,"hDZttSel_%d",ind+ 1); hDZttSel[ind+ 1] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZttSel_%d",ind+ 2); hDZttSel[ind+ 2] = new TH1D(sb,sb,300,-1.0,2.0);
    sprintf(sb,"hDZttSel_%d",ind+ 3); hDZttSel[ind+ 3] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZttSel_%d",ind+ 4); hDZttSel[ind+ 4] = new TH1D(sb,sb,180,0.,180.);
    sprintf(sb,"hDZttSel_%d",ind+ 5); hDZttSel[ind+ 5] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZttSel_%d",ind+ 6); hDZttSel[ind+ 6] = new TH1D(sb,sb,180,0.,180.);
    sprintf(sb,"hDZttSel_%d",ind+ 7); hDZttSel[ind+ 7] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDZttSel_%d",ind+ 8); hDZttSel[ind+ 8] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDZttSel_%d",ind+ 9); hDZttSel[ind+ 9] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDZttSel_%d",ind+10); hDZttSel[ind+10] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZttSel_%d",ind+11); hDZttSel[ind+11] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZttSel_%d",ind+12); hDZttSel[ind+12] = new TH1D(sb,sb,200,0.,20.);
    sprintf(sb,"hDZttSel_%d",ind+13); hDZttSel[ind+13] = new TH1D(sb,sb,200,0.,800.);
    sprintf(sb,"hDZttSel_%d",ind+14); hDZttSel[ind+14] = new TH1D(sb,sb,100,0.,5.);
    sprintf(sb,"hDZttSel_%d",ind+15); hDZttSel[ind+15] = new TH1D(sb,sb,300,0.,3.);
    sprintf(sb,"hDZttSel_%d",ind+16); hDZttSel[ind+16] = new TH1D(sb,sb,300,0.,3.);
  }

  for(int i=0; i<17; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZttSel[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZttHLT_%d",ind+0);  hDZttHLT[ind+0]  = new TH1D(sb,sb,9,-0.5,8.5);
    sprintf(sb,"hDZttHLT_%d",ind+1);  hDZttHLT[ind+1]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZttHLT_%d",ind+2);  hDZttHLT[ind+2]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZttHLT_%d",ind+3);  hDZttHLT[ind+3]  = new TH1D(sb,sb,2,-0.5,1.5);
  }

  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZttHLT[i+j*100]);
    }
  }

}

//--------------------------------------------------------------------------------------------------
void ZttEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ZttEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
