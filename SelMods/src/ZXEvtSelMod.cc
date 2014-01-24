// $Id: ZXEvtSelMod.cc,v 1.11 2012/07/05 07:22:12 ceballos Exp $

#include "Analysis/SelMods/interface/ZXEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

using namespace mithep;
ClassImp(mithep::ZXEvtSelMod)

//-----------------------------------------------------------------------------
ZXEvtSelMod::ZXEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fCaloJetName0("AKt5Jets"),
  fCaloJet0(0),
  fUsePDFs(kFALSE),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fRho(0.0),
  fTrackName(Names::gkTrackBrn),
  fGsfTrackName("GsfTracks"),
  fTracks(0),
  fGsfTracks(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//-----------------------------------------------------------------------------
void ZXEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//-----------------------------------------------------------------------------
void ZXEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Get Generator Level information for matching
  ObjArray<MCParticle> *GenLeptons   = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCLeptonsName.Data()));

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);
  PFTauOArr  *CleanPFTaus    = GetObjThisEvt<PFTauOArr>(ModNames::gkCleanPFTausName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");
  FArrDouble *PDFArr = 0;
  if(fUsePDFs == kTRUE){
    PDFArr =  GetObjThisEvt<FArrDouble >("PDFWeights");
  //cout << PDFArr->GetEntries() << " " << PDFArr->At(0) << endl;
  }

  if(GenLeptons->GetEntries() > 0){
    hDZZGenSel[0]->Fill((double)GenLeptons->GetEntries(),NNLOWeight->GetVal());
    if(GenLeptons->GetEntries() == 4 &&
       GenLeptons->At(0)->AbsEta() < 5 && GenLeptons->At(1)->AbsEta() < 5 &&
       GenLeptons->At(2)->AbsEta() < 5 && GenLeptons->At(3)->AbsEta() < 5){
      double etaMax[2] = {-1, -1};
      double ptMin[2] = {9999, 9999};
      for(unsigned int i=0; i<GenLeptons->GetEntries(); i++){
        const MCParticle *geni = GenLeptons->At(i);
	if(geni->Pt() < ptMin[0]) ptMin[0] = geni->Pt();
        if(geni->AbsEta() > etaMax[0]) etaMax[0] = geni->AbsEta();
        for(unsigned int j=i+1; j<GenLeptons->GetEntries(); j++){
          const MCParticle *genj = GenLeptons->At(j);
	  if(geni->AbsPdgId() == genj->AbsPdgId() && geni->PdgId() != genj->PdgId()){
            CompositeParticle dilepton;
            dilepton.AddDaughter(geni);
            dilepton.AddDaughter(genj);
            hDZZGenSel[1]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
	    if(TMath::Abs(dilepton.Mass()-91.1876)<15){
              if(geni->Pt() < ptMin[1]) ptMin[1] = geni->Pt();
              if(geni->AbsEta() > etaMax[1]) etaMax[1] = geni->AbsEta();
              if(genj->Pt() < ptMin[1]) ptMin[1] = genj->Pt();
              if(genj->AbsEta() > etaMax[1]) etaMax[1] = genj->AbsEta();	      
	    } // |mll-mzz|<15
	  } // same-flavor, opposite-charge
	} // loop j
      } // loop i
      hDZZGenSel[2]->Fill(TMath::Min(ptMin[0],99.999),NNLOWeight->GetVal());
      if(etaMax[1] >= 0) hDZZGenSel[3]->Fill(TMath::Min(ptMin[1],99.999),NNLOWeight->GetVal());
      if(ptMin[0] >  5) hDZZGenSel[4]->Fill(etaMax[0],NNLOWeight->GetVal());
      if(ptMin[0] > 10) hDZZGenSel[5]->Fill(etaMax[0],NNLOWeight->GetVal());
      if(etaMax[1] >= 0 && ptMin[1] >  5) hDZZGenSel[6]->Fill(etaMax[1],NNLOWeight->GetVal());
      if(etaMax[1] >= 0 && ptMin[1] > 10) hDZZGenSel[7]->Fill(etaMax[1],NNLOWeight->GetVal());
    } // 4 leptons
  } // ZZ Gen analysis

  hDZXSel[0]->Fill((double)leptons->GetEntries(),NNLOWeight->GetVal());

  // HZZ Forward
  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries() >= 1){
    LoadBranch(fPFCandidatesName);
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
    fRho = fPileupEnergyDensity->At(0)->Rho();
    HZZForward(leptons, fPFCandidates, NNLOWeight->GetVal());

    PhotonOArr *CleanPhotons = GetObjThisEvt<PhotonOArr>(ModNames::gkCleanPhotonsName);
    HZZForward(leptons, CleanPhotons, NNLOWeight->GetVal());

    ParticleOArr *leptonsFakeable = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
    ParticleOArr *leptonsOnlyFake = new ObjArray<Particle>;
    for(UInt_t i=0; i<leptonsFakeable->GetEntries(); i++){
      Bool_t isOnlyFake = kTRUE;
      for(UInt_t j=0; j<leptons->GetEntries(); j++) {
	if(leptons->At(j) == leptonsFakeable->At(i)) {
          isOnlyFake = kFALSE;
	  break;
	}
      }
      if(isOnlyFake == kTRUE) leptonsOnlyFake->Add(leptonsFakeable->At(i));
    }
    leptonsOnlyFake->Sort();
    HZZForward(leptons, leptonsOnlyFake, NNLOWeight->GetVal());
    delete leptonsOnlyFake;
  }

  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries() == 2){
    LoadBranch(fTrackName);
    LoadBranch(fGsfTrackName);
    hDWGGenSel[0]->Fill((double)GenLeptons->GetEntries(),NNLOWeight->GetVal());
    int nLepFid = 0;
    double ptMin[1] = {9999};
    if(GenLeptons->GetEntries() >= 2){
      for(unsigned int i=0; i<GenLeptons->GetEntries(); i++){
        const MCParticle *geni = GenLeptons->At(i);
	if(geni->AbsEta() >= 2.6) continue;
	nLepFid++;
	if(geni->Pt() < ptMin[0]) ptMin[0] = geni->Pt();
        for(unsigned int j=i+1; j<GenLeptons->GetEntries(); j++){
          const MCParticle *genj = GenLeptons->At(j);
	  if(geni->AbsPdgId() == genj->AbsPdgId() && geni->PdgId() != genj->PdgId()){
            CompositeParticle dilepton;
            dilepton.AddDaughter(geni);
            dilepton.AddDaughter(genj);
            hDWGGenSel[1]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
	    if(dilepton.Mass() > 37 && dilepton.Mass() < 39) return;
	  } // same-flavor, opposite-charge
	} // loop j
      } // loop i
      hDWGGenSel[2]->Fill((double)nLepFid,NNLOWeight->GetVal());
      if(nLepFid == 2)
        hDWGGenSel[3]->Fill(TMath::Min(ptMin[0],49.999),NNLOWeight->GetVal()); 
      if(nLepFid >= 3)
        hDWGGenSel[4]->Fill(TMath::Min(ptMin[0],49.999),NNLOWeight->GetVal());
      if(nLepFid >=3){
	for (UInt_t i=0; i<CleanMuons->GetEntries(); i++) {
	  double dRMin = 9999.; int indexDrMin = -1;
	  for (UInt_t j=0; j<fTracks->GetEntries(); j++) {
            const Track *track = fTracks->At(j);
	    if(CleanMuons->At(i)->Charge() == track->Charge()) continue;
	    double dR = MathUtils::DeltaR(CleanMuons->At(i)->Mom(), track->Mom());
	    if(dR < dRMin) {dRMin = dR; indexDrMin = j;}
	  }
          GenericParticle *p = new GenericParticle(fTracks->At(indexDrMin)->Px(), fTracks->At(indexDrMin)->Py(), fTracks->At(indexDrMin)->Pz(), fTracks->At(indexDrMin)->P(), fTracks->At(indexDrMin)->Charge());
          CompositeParticle dilepton;
          dilepton.AddDaughter(CleanMuons->At(i));
          dilepton.AddDaughter(p);
	  hDWGGenSel[5]->Fill(TMath::Min(dilepton.Mass(),9.999),NNLOWeight->GetVal());
	  delete p;
	}
	for (UInt_t i=0; i<CleanElectrons->GetEntries(); i++) {
	  double dRMin = 9999.; int indexDrMin = -1;
	  for (UInt_t j=0; j<fTracks->GetEntries(); j++) {
            const Track *track = fTracks->At(j);
	    if(CleanElectrons->At(i)->Charge() == track->Charge()) continue;
	    double dR = MathUtils::DeltaR(CleanElectrons->At(i)->Mom(), track->Mom());
	    if(dR < dRMin) {dRMin = dR; indexDrMin = j;}
	  }
          GenericParticle *p = new GenericParticle(fTracks->At(indexDrMin)->Px(), fTracks->At(indexDrMin)->Py(), fTracks->At(indexDrMin)->Pz(), fTracks->At(indexDrMin)->P(), fTracks->At(indexDrMin)->Charge());
          CompositeParticle dilepton;
          dilepton.AddDaughter(CleanElectrons->At(i));
          dilepton.AddDaughter(p);
	  hDWGGenSel[6]->Fill(TMath::Min(dilepton.Mass(),9.999),NNLOWeight->GetVal());
	  delete p;
	}
      }
    } // Gen leptons >= 2
  } // WG study

  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries()+CleanPFTaus->GetEntries() < 3) return;
  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries()                           < 2) return;

  // WZ analysis
  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries() == 3){
    int nWZType = -1;
    vector<ChargedParticle*> leptonsWZ;
    if     (CleanMuons->GetEntries() == 2){
      leptonsWZ.push_back(CleanMuons->At(0));
      leptonsWZ.push_back(CleanMuons->At(1));
      leptonsWZ.push_back(CleanElectrons->At(0));
      nWZType = 0;
      if(leptonsWZ[0]->Charge() == leptonsWZ[1]->Charge()) nWZType = 4;
    }
    else if(CleanElectrons->GetEntries() == 2){
      leptonsWZ.push_back(CleanElectrons->At(0));
      leptonsWZ.push_back(CleanElectrons->At(1));
      leptonsWZ.push_back(CleanMuons->At(0));
      nWZType = 1;
      if(leptonsWZ[0]->Charge() == leptonsWZ[1]->Charge()) nWZType = 5;
    }
    else if(CleanMuons->GetEntries() == 3){
      if(TMath::Abs(CleanMuons->At(0)->Charge()+CleanMuons->At(1)->Charge()+
                    CleanMuons->At(2)->Charge())==1){
        nWZType = 2;
	if     (CleanMuons->At(0)->Charge() != CleanMuons->At(1)->Charge() &&
	        CleanMuons->At(0)->Charge() != CleanMuons->At(2)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanMuons->At(0));
          dileptonA.AddDaughter(CleanMuons->At(1));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanMuons->At(0));
          dileptonB.AddDaughter(CleanMuons->At(2));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanMuons->At(0));
            leptonsWZ.push_back(CleanMuons->At(1));
            leptonsWZ.push_back(CleanMuons->At(2));	    
	  }
	  else {
            leptonsWZ.push_back(CleanMuons->At(0));
            leptonsWZ.push_back(CleanMuons->At(2));
            leptonsWZ.push_back(CleanMuons->At(1));	    
	  }
	}
	else if(CleanMuons->At(1)->Charge() != CleanMuons->At(0)->Charge() &&
	        CleanMuons->At(1)->Charge() != CleanMuons->At(2)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanMuons->At(1));
          dileptonA.AddDaughter(CleanMuons->At(0));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanMuons->At(1));
          dileptonB.AddDaughter(CleanMuons->At(2));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanMuons->At(0));
            leptonsWZ.push_back(CleanMuons->At(1));
            leptonsWZ.push_back(CleanMuons->At(2));	    
	  }
	  else {
            leptonsWZ.push_back(CleanMuons->At(1));
            leptonsWZ.push_back(CleanMuons->At(2));
            leptonsWZ.push_back(CleanMuons->At(0));	    
	  }
	}
	else if(CleanMuons->At(2)->Charge() != CleanMuons->At(0)->Charge() &&
	        CleanMuons->At(2)->Charge() != CleanMuons->At(1)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanMuons->At(2));
          dileptonA.AddDaughter(CleanMuons->At(0));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanMuons->At(2));
          dileptonB.AddDaughter(CleanMuons->At(1));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanMuons->At(0));
            leptonsWZ.push_back(CleanMuons->At(2));
            leptonsWZ.push_back(CleanMuons->At(1));	    
	  }
	  else {
            leptonsWZ.push_back(CleanMuons->At(1));
            leptonsWZ.push_back(CleanMuons->At(2));
            leptonsWZ.push_back(CleanMuons->At(0));	    
	  }
	}
	else {
	  printf("Impossible in nWZType = 2\n");
	}
      }
      else {
        nWZType = 6;
      }
    }
    else if(CleanElectrons->GetEntries() == 3){
      if(TMath::Abs(CleanElectrons->At(0)->Charge()+CleanElectrons->At(1)->Charge()+
                    CleanElectrons->At(2)->Charge())==1){
        nWZType = 3;
	if     (CleanElectrons->At(0)->Charge() != CleanElectrons->At(1)->Charge() &&
	        CleanElectrons->At(0)->Charge() != CleanElectrons->At(2)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanElectrons->At(0));
          dileptonA.AddDaughter(CleanElectrons->At(1));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanElectrons->At(0));
          dileptonB.AddDaughter(CleanElectrons->At(2));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanElectrons->At(0));
            leptonsWZ.push_back(CleanElectrons->At(1));
            leptonsWZ.push_back(CleanElectrons->At(2));	    
	  }
	  else {
            leptonsWZ.push_back(CleanElectrons->At(0));
            leptonsWZ.push_back(CleanElectrons->At(2));
            leptonsWZ.push_back(CleanElectrons->At(1));	    
	  }
	}
	else if(CleanElectrons->At(1)->Charge() != CleanElectrons->At(0)->Charge() &&
	        CleanElectrons->At(1)->Charge() != CleanElectrons->At(2)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanElectrons->At(1));
          dileptonA.AddDaughter(CleanElectrons->At(0));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanElectrons->At(1));
          dileptonB.AddDaughter(CleanElectrons->At(2));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanElectrons->At(0));
            leptonsWZ.push_back(CleanElectrons->At(1));
            leptonsWZ.push_back(CleanElectrons->At(2));	    
	  }
	  else {
            leptonsWZ.push_back(CleanElectrons->At(1));
            leptonsWZ.push_back(CleanElectrons->At(2));
            leptonsWZ.push_back(CleanElectrons->At(0));	    
	  }
	}
	else if(CleanElectrons->At(2)->Charge() != CleanElectrons->At(0)->Charge() &&
	        CleanElectrons->At(2)->Charge() != CleanElectrons->At(1)->Charge()){
          CompositeParticle dileptonA;
          dileptonA.AddDaughter(CleanElectrons->At(2));
          dileptonA.AddDaughter(CleanElectrons->At(0));
          CompositeParticle dileptonB;
          dileptonB.AddDaughter(CleanElectrons->At(2));
          dileptonB.AddDaughter(CleanElectrons->At(1));
	  if(TMath::Abs(dileptonA.Mass()-91.1876) < TMath::Abs(dileptonB.Mass()-91.1876)){
            leptonsWZ.push_back(CleanElectrons->At(0));
            leptonsWZ.push_back(CleanElectrons->At(2));
            leptonsWZ.push_back(CleanElectrons->At(1));	    
	  }
	  else {
            leptonsWZ.push_back(CleanElectrons->At(1));
            leptonsWZ.push_back(CleanElectrons->At(2));
            leptonsWZ.push_back(CleanElectrons->At(0));	    
	  }
	}
	else {
	  printf("Impossible in nWZType = 3\n");
	}
      }
      else {
        nWZType = 7;
      }
    }

    // Start analysis after selecting WZ -> l0 l1 l2
    hDWZSel[90]->Fill((double) nWZType,NNLOWeight->GetVal());

    // Charge of the leptons should be opposite
    if (nWZType < 4){
      if(fUsePDFs == kTRUE){
        for(UInt_t i=0; i<PDFArr->GetEntries(); i++){
     	  hDWZPDF[i]->Fill(TMath::Min(caloMet->Pt(),199.999),PDFArr->At(i));
        }
      }

      CompositeParticle dilepton;
      dilepton.AddDaughter(leptonsWZ[0]);
      dilepton.AddDaughter(leptonsWZ[1]);
      CompositeParticle *trilepton = new CompositeParticle();
      trilepton->AddDaughter(leptonsWZ[0]);
      trilepton->AddDaughter(leptonsWZ[1]);
      trilepton->AddDaughter(leptonsWZ[2]);

      if(dilepton.Charge() != 0) 
         printf("Impossible, nchaZ != 0, nWZType ===> %d\n",nWZType);

      hDWZSel[ 0+100*nWZType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
      // Z mass requirement
      if(TMath::Abs(dilepton.Mass()-91.1876) < 20){
        hDWZSel[ 1+100*nWZType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
        // MET requirement
	if(caloMet->Pt() > 20){
          // Pt(W->ln) requirement & mTW
	  double deltaPhiln = fabs(MathUtils::DeltaPhi(leptonsWZ[2]->Phi(), caloMet->Phi()));
	  double mTW = TMath::Sqrt(2.0*leptonsWZ[2]->Pt()*caloMet->Pt()*
          			  (1.0 - cos(deltaPhiln)));
          hDWZSel[ 2+100*nWZType]->Fill(TMath::Min(leptonsWZ[2]->Pt(),199.999),NNLOWeight->GetVal());
	  hDWZSel[ 8+100*nWZType]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
	  if(mTW > 15.0){
	    // Sort and count the number of central Jets for vetoing
            double theDijetMass[2] = {999999999.0, 999999999.0};
	    UInt_t indexJet[2] = {0, 0};
	    vector<Jet*> sortedJets;
	    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
	      if(CleanJets->At(i)->AbsEta() < fEtaJetCut &&
	         CleanJets->At(i)->Pt() > fPtJetCut){
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
   	        sortedJets.push_back(jet_f);
	      }
	    }
            double nMaxMediumBTagJet = 0.0;
	    for(UInt_t i=0; i<sortedJets.size(); i++){
              if(nMaxMediumBTagJet < sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc())
                 nMaxMediumBTagJet = sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc();
	      for(UInt_t j=i+1; j<sortedJets.size(); j++){
      	  	CompositeParticle *dijet = new CompositeParticle();
      	  	dijet->AddDaughter(sortedJets[i]);
      	  	dijet->AddDaughter(sortedJets[j]);
      	  	if(theDijetMass[0] > 200){
      	  	  theDijetMass[0] = dijet->Mass();
		  indexJet[0] = i; indexJet[1] = j;
      	  	}
      	  	delete dijet;
		if(sortedJets[i]->Pt() < sortedJets[j]->Pt()) {
		  //swap i and j
        	  Jet* tempjet = sortedJets[i];
        	  sortedJets[i] = sortedJets[j];
        	  sortedJets[j] = tempjet;	
		}
	      }
	    }
            if(sortedJets.size() <= 1) {theDijetMass[0] = 0.001; theDijetMass[1] = 0.001;}
	    if(sortedJets.size()  > 1) {
      	      CompositeParticle *dijet = new CompositeParticle();
      	      dijet->AddDaughter(sortedJets[0]);
      	      dijet->AddDaughter(sortedJets[1]);
	      theDijetMass[1] = dijet->Mass();
	    }

	    double deltaRJetl = 999.;
	    double deltaRJetMet = 999.;
	    if(sortedJets.size() > 0){
	     for(int i=0; i<3; i++){
        	if(MathUtils::DeltaR(leptonsWZ[i]->Phi(), leptonsWZ[i]->Eta(),
	                	     sortedJets[0]->Phi(), sortedJets[0]->Eta()) < deltaRJetl)
		  deltaRJetl = MathUtils::DeltaR(leptonsWZ[i]->Phi(), leptonsWZ[i]->Eta(),
	                                	 sortedJets[0]->Phi(), sortedJets[0]->Eta()); 
	      }
	      deltaRJetMet = fabs(MathUtils::DeltaPhi(caloMet->Phi(), sortedJets[0]->Phi()));
	    }
	    double deltaPhiln = fabs(MathUtils::DeltaPhi(leptonsWZ[2]->Phi(), caloMet->Phi()));
	    double deltaRWl[2] = {MathUtils::DeltaR(leptonsWZ[0]->Phi(), leptonsWZ[0]->Eta(),
	                                	    leptonsWZ[2]->Phi(), leptonsWZ[2]->Eta()),
				  MathUtils::DeltaR(leptonsWZ[1]->Phi(), leptonsWZ[1]->Eta(),
	                                	    leptonsWZ[2]->Phi(), leptonsWZ[2]->Eta())};		    ;
            double deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(caloMet->Phi(), 
                                                       dilepton.Phi()))* 180./TMath::Pi();
            double deltaPhiTrileptonMet = fabs(MathUtils::DeltaPhi(caloMet->Phi(), 
                                                       trilepton->Phi()))* 180./TMath::Pi();
	    hDWZSel[ 3+100*nWZType]->Fill(TMath::Min(leptonsWZ[0]->Pt(),199.999),NNLOWeight->GetVal());
	    hDWZSel[ 4+100*nWZType]->Fill(TMath::Min(leptonsWZ[1]->Pt(),199.999),NNLOWeight->GetVal());
	    hDWZSel[ 5+100*nWZType]->Fill((double)sortedJets.size(),NNLOWeight->GetVal());
	    hDWZSel[ 6+100*nWZType]->Fill(caloMet->MetSig(),NNLOWeight->GetVal());
	    hDWZSel[ 7+100*nWZType]->Fill(caloMet->SumEt(),NNLOWeight->GetVal());
	    hDWZSel[ 9+100*nWZType]->Fill(deltaPhiln * 180./TMath::Pi(),NNLOWeight->GetVal());
	    hDWZSel[10+100*nWZType]->Fill(TMath::Min(deltaRWl[0],deltaRWl[1]),NNLOWeight->GetVal());
	    hDWZSel[11+100*nWZType]->Fill(TMath::Max(deltaRWl[0],deltaRWl[1]),NNLOWeight->GetVal());
	    hDWZSel[12+100*nWZType]->Fill(TMath::Min(caloMet->Pt()*deltaPhiln/4.,199.999),NNLOWeight->GetVal());
	    if(sortedJets.size() > 0){
	      hDWZSel[13+100*nWZType]->Fill(TMath::Min(sortedJets[0]->Pt(),199.99),NNLOWeight->GetVal());
	      hDWZSel[14+100*nWZType]->Fill(deltaRJetl,NNLOWeight->GetVal());
	      hDWZSel[15+100*nWZType]->Fill(deltaRJetMet * 180./TMath::Pi(),NNLOWeight->GetVal());
	      
	      if(sortedJets.size() > 1){
	        hDZHSel[0]->Fill(TMath::Min(TMath::Max(nMaxMediumBTagJet,0.00001),0.9999),NNLOWeight->GetVal());
		if(nMaxMediumBTagJet < 0.8){
	          hDZHSel[1]->Fill(TMath::Min(theDijetMass[1],199.99),NNLOWeight->GetVal());
                  if(theDijetMass[0] < 200){
	            hDZHSel[2]->Fill(TMath::Min(theDijetMass[0],199.99),NNLOWeight->GetVal());
	            hDZHSel[3]->Fill(TMath::Min(sortedJets[indexJet[0]]->Pt(),199.99),NNLOWeight->GetVal());
	            hDZHSel[4]->Fill(TMath::Min(sortedJets[indexJet[1]]->Pt(),199.99),NNLOWeight->GetVal());
	            hDZHSel[5]->Fill(TMath::Min((double)nWZType,3.4999),NNLOWeight->GetVal());
    		    dilepton.AddDaughter(leptons->At(0));
    		    dilepton.AddDaughter(leptons->At(1));
    		    TVector3 v0(sortedJets[indexJet[0]]->Px(),sortedJets[indexJet[0]]->Py(),sortedJets[indexJet[0]]->Pz());
    		    TVector3 v1(sortedJets[indexJet[1]]->Px(),sortedJets[indexJet[1]]->Py(),sortedJets[indexJet[1]]->Pz());
	            hDZHSel[6] ->Fill(TMath::Min(MathUtils::DeltaR(sortedJets[indexJet[0]]->Mom(), sortedJets[indexJet[1]]->Mom()),5.99),NNLOWeight->GetVal());
	            hDZHSel[7] ->Fill(fabs(MathUtils::DeltaPhi(sortedJets[indexJet[0]]->Phi(), sortedJets[indexJet[1]]->Phi())) * 180./TMath::Pi(),NNLOWeight->GetVal());
	            hDZHSel[8]->Fill(v0.Angle(v1) * 180./TMath::Pi(),NNLOWeight->GetVal());
		    if(sortedJets[indexJet[0]]->Pt() > 30 &&
		       fabs(MathUtils::DeltaPhi(sortedJets[indexJet[0]]->Phi(), sortedJets[indexJet[1]]->Phi())) * 180./TMath::Pi() < 178.0){
		      hDZHSel[9]->Fill(TMath::Min(theDijetMass[0],199.99),NNLOWeight->GetVal());
                      if(theDijetMass[0] > 50 && theDijetMass[0] < 110){
		        CompositeParticle *mh = new CompositeParticle();
                        mh->AddDaughter(leptonsWZ[2]);
                        mh->AddDaughter(caloMet);
                        mh->AddDaughter(sortedJets[indexJet[0]]);
                        mh->AddDaughter(sortedJets[indexJet[1]]);
		        hDZHSel[10]->Fill(TMath::Min(mh->Mass(),399.99),NNLOWeight->GetVal());
		        hDZHSel[11]->Fill(TMath::Min(mh->Mt(),399.99),NNLOWeight->GetVal());
		        delete mh;
		      }
		    }
		  }
		}
	      }
	    }
	    hDWZSel[16+100*nWZType]->Fill(180.-deltaPhiDileptonMet,NNLOWeight->GetVal());
	    hDWZSel[17+100*nWZType]->Fill(180.-deltaPhiTrileptonMet,NNLOWeight->GetVal());
	    delete trilepton;
	    for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
          } // mtw requirement
	} // MET requirement
      } // Z mass requirement
    } // Z Charge == 0
  } // WZ Selection

  // ZZ->4l analysis
  if(CleanMuons->GetEntries()+CleanElectrons->GetEntries() >= 4){
    int nZZType = -3;
    vector<ChargedParticle*> leptonsZZ;
    if     (CleanMuons->GetEntries() >= 2 && CleanElectrons->GetEntries() >= 2){
      nZZType = 0;
      leptonsZZ.push_back(CleanMuons->At(0));
      leptonsZZ.push_back(CleanMuons->At(1));
      leptonsZZ.push_back(CleanElectrons->At(0));
      leptonsZZ.push_back(CleanElectrons->At(1));
    }
    else if(CleanMuons->GetEntries() >= 4){
      CompositeParticle dileptonA;
      CompositeParticle dileptonB;
      if     (CleanMuons->At(0)->Charge() != CleanMuons->At(1)->Charge() &&
              CleanMuons->At(0)->Charge() != CleanMuons->At(2)->Charge()){
	nZZType = 1;
        dileptonA.AddDaughter(CleanMuons->At(0));
        dileptonA.AddDaughter(CleanMuons->At(1));
        dileptonB.AddDaughter(CleanMuons->At(0));
        dileptonB.AddDaughter(CleanMuons->At(2));
      }
      else if(CleanMuons->At(0)->Charge() != CleanMuons->At(1)->Charge() &&
              CleanMuons->At(0)->Charge() != CleanMuons->At(3)->Charge()){
	nZZType = 2;
        dileptonA.AddDaughter(CleanMuons->At(0));
        dileptonA.AddDaughter(CleanMuons->At(1));
        dileptonB.AddDaughter(CleanMuons->At(0));
        dileptonB.AddDaughter(CleanMuons->At(3));
      } 
      else if(CleanMuons->At(0)->Charge() != CleanMuons->At(2)->Charge() &&
              CleanMuons->At(0)->Charge() != CleanMuons->At(3)->Charge()){
	nZZType = 3;
        dileptonA.AddDaughter(CleanMuons->At(0));
        dileptonA.AddDaughter(CleanMuons->At(2));
        dileptonB.AddDaughter(CleanMuons->At(0));
        dileptonB.AddDaughter(CleanMuons->At(3));
      }
      else {
        nZZType = -1;
      }
      if     (nZZType == 1 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(3));
      }
      else if(nZZType == 1) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(3));
        nZZType = 4;
      }
      else if(nZZType == 2 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(3));
      }
      else if(nZZType == 2) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(3));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        nZZType = 5;
      }
      else if(nZZType == 3 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(3));
      }
      else if(nZZType == 3) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(3));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        nZZType = 6;
      }
    }
    else if(CleanElectrons->GetEntries() >= 4){
      CompositeParticle dileptonA;
      CompositeParticle dileptonB;
      if     (CleanElectrons->At(0)->Charge() != CleanElectrons->At(1)->Charge() &&
              CleanElectrons->At(0)->Charge() != CleanElectrons->At(2)->Charge()){
	nZZType = 7;
        dileptonA.AddDaughter(CleanElectrons->At(0));
        dileptonA.AddDaughter(CleanElectrons->At(1));
        dileptonB.AddDaughter(CleanElectrons->At(0));
        dileptonB.AddDaughter(CleanElectrons->At(2));
      }
      else if(CleanElectrons->At(0)->Charge() != CleanElectrons->At(1)->Charge() &&
              CleanElectrons->At(0)->Charge() != CleanElectrons->At(3)->Charge()){
	nZZType = 8;
        dileptonA.AddDaughter(CleanElectrons->At(0));
        dileptonA.AddDaughter(CleanElectrons->At(1));
        dileptonB.AddDaughter(CleanElectrons->At(0));
        dileptonB.AddDaughter(CleanElectrons->At(3));
      } 
      else if(CleanElectrons->At(0)->Charge() != CleanElectrons->At(2)->Charge() &&
              CleanElectrons->At(0)->Charge() != CleanElectrons->At(3)->Charge()){
	nZZType = 9;
        dileptonA.AddDaughter(CleanElectrons->At(0));
        dileptonA.AddDaughter(CleanElectrons->At(2));
        dileptonB.AddDaughter(CleanElectrons->At(0));
        dileptonB.AddDaughter(CleanElectrons->At(3));
      }
      else {
        nZZType = -2;
      }
      if     (nZZType == 7 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(3));
      }
      else if(nZZType == 7) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(3));
        nZZType = 10;
      }
      else if(nZZType == 8 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(3));
      }
      else if(nZZType == 8) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(3));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        nZZType = 11;
      }
      else if(nZZType == 9 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(3));
      }
      else if(nZZType == 9) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(3));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        nZZType = 12;
      }
    }
    else if(CleanMuons->GetEntries() >= 3 && CleanElectrons->GetEntries() >= 1){
      CompositeParticle dileptonA;
      CompositeParticle dileptonB;
      if     (CleanMuons->At(0)->Charge() != CleanMuons->At(1)->Charge() &&
              CleanMuons->At(0)->Charge() != CleanMuons->At(2)->Charge()){
        dileptonA.AddDaughter(CleanMuons->At(0));
        dileptonA.AddDaughter(CleanMuons->At(1));
        dileptonB.AddDaughter(CleanMuons->At(0));
        dileptonB.AddDaughter(CleanMuons->At(2));
        nZZType = 13;
      }
      else if(CleanMuons->At(1)->Charge() != CleanMuons->At(2)->Charge() &&
              CleanMuons->At(1)->Charge() != CleanMuons->At(0)->Charge()){
        dileptonA.AddDaughter(CleanMuons->At(1));
        dileptonA.AddDaughter(CleanMuons->At(2));
        dileptonB.AddDaughter(CleanMuons->At(1));
        dileptonB.AddDaughter(CleanMuons->At(0));
        nZZType = 14;
      }
      else {
        nZZType = -4;
      }

      if     (nZZType == 13 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(0));
      }
      else if(nZZType == 13) {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(0));
        nZZType = 14;
      }
      else if(nZZType == 14 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(0));
      }
      else if(nZZType == 14) {
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(0));
        nZZType = 14;
      }
    }
    else if(CleanElectrons->GetEntries() >= 3 && CleanMuons->GetEntries() >= 1){
      CompositeParticle dileptonA;
      CompositeParticle dileptonB;

      if     (CleanElectrons->At(0)->Charge() != CleanElectrons->At(1)->Charge() &&
              CleanElectrons->At(0)->Charge() != CleanElectrons->At(2)->Charge()){
        dileptonA.AddDaughter(CleanElectrons->At(0));
        dileptonA.AddDaughter(CleanElectrons->At(1));
        dileptonB.AddDaughter(CleanElectrons->At(0));
        dileptonB.AddDaughter(CleanElectrons->At(2));
        nZZType = 15;
      }
      else if(CleanElectrons->At(1)->Charge() != CleanElectrons->At(2)->Charge() &&
              CleanElectrons->At(1)->Charge() != CleanElectrons->At(0)->Charge()){
        dileptonA.AddDaughter(CleanElectrons->At(1));
        dileptonA.AddDaughter(CleanElectrons->At(2));
        dileptonB.AddDaughter(CleanElectrons->At(1));
        dileptonB.AddDaughter(CleanElectrons->At(0));
        nZZType = 16;
      }
      else {
        nZZType = -5;
      }

      if     (nZZType == 15 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanMuons->At(0));
      }
      else if(nZZType == 15) {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanMuons->At(0));
        nZZType = 16;
      }
      else if(nZZType == 16 &&
	      TMath::Abs(dileptonA.Mass()-91.1876) <
	      TMath::Abs(dileptonB.Mass()-91.1876)  ) {
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanMuons->At(0));
      }
      else if(nZZType == 16) {
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanMuons->At(0));
        nZZType = 16;
      }
    }

    hDZZSel[90]->Fill((double)nZZType,NNLOWeight->GetVal());
    if(nZZType >= 0 && nZZType <= 12){
      if(fUsePDFs == kTRUE){
        for(UInt_t i=0; i<PDFArr->GetEntries(); i++){
          hDZZPDF[i]->Fill(TMath::Min(caloMet->Pt(),199.999),PDFArr->At(i));
        }
      }

      CompositeParticle *dileptonZ1 = new CompositeParticle();
      dileptonZ1->AddDaughter(leptonsZZ[0]);
      dileptonZ1->AddDaughter(leptonsZZ[1]);
      CompositeParticle *dileptonZ2 = new CompositeParticle();
      dileptonZ2->AddDaughter(leptonsZZ[2]);
      dileptonZ2->AddDaughter(leptonsZZ[3]);
      hDZZSel[91]->Fill((double)(dileptonZ1->Charge()+dileptonZ2->Charge()),NNLOWeight->GetVal());

      CompositeParticle *particleH = new CompositeParticle();
      particleH->AddDaughter(leptonsZZ[0]);
      particleH->AddDaughter(leptonsZZ[1]);
      particleH->AddDaughter(leptonsZZ[2]);
      particleH->AddDaughter(leptonsZZ[3]);

      int simpleZZType = 0;
      if     (nZZType >= 1 && nZZType <=  6) simpleZZType = 1;
      else if(nZZType >= 7 && nZZType <= 12) simpleZZType = 2;

      hDZZSel[ 0+100*simpleZZType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 1+100*simpleZZType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 2+100*simpleZZType]->Fill(TMath::Min(leptons->At(2)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 3+100*simpleZZType]->Fill(TMath::Min(leptons->At(3)->Pt(),199.999),NNLOWeight->GetVal());
      if(TMath::Abs(dileptonZ1->Mass()-91.1876) < TMath::Abs(dileptonZ2->Mass()-91.1876)){
        hDZZSel[ 4+100*simpleZZType]->Fill(TMath::Min(dileptonZ1->Mass(),199.99),NNLOWeight->GetVal());
        hDZZSel[ 5+100*simpleZZType]->Fill(TMath::Min(dileptonZ2->Mass(),199.99),NNLOWeight->GetVal());
      }
      else {
        hDZZSel[ 4+100*simpleZZType]->Fill(TMath::Min(dileptonZ2->Mass(),199.99),NNLOWeight->GetVal());
        hDZZSel[ 5+100*simpleZZType]->Fill(TMath::Min(dileptonZ1->Mass(),199.99),NNLOWeight->GetVal());
      }
      if(dileptonZ1->Mass() > 12 && dileptonZ2->Mass() > 12){
        hDZZSel[ 6+100*simpleZZType]->Fill(TMath::Min(particleH->Mass(),999.99),NNLOWeight->GetVal());
      }
      delete dileptonZ1;
      delete dileptonZ2;
      delete particleH;
    } // nZZType >= 0 && <= 12
    else if(nZZType >= 13 && nZZType <= 16){
      CompositeParticle *dileptonZ1 = new CompositeParticle();
      dileptonZ1->AddDaughter(leptonsZZ[0]);
      dileptonZ1->AddDaughter(leptonsZZ[1]);
      CompositeParticle *dileptonZ2 = new CompositeParticle();
      dileptonZ2->AddDaughter(leptonsZZ[2]);
      dileptonZ2->AddDaughter(leptonsZZ[3]);
      hDZZSel[92]->Fill((double)(dileptonZ1->Charge()+dileptonZ2->Charge()),NNLOWeight->GetVal());

      CompositeParticle *particleH = new CompositeParticle();
      particleH->AddDaughter(leptonsZZ[0]);
      particleH->AddDaughter(leptonsZZ[1]);
      particleH->AddDaughter(leptonsZZ[2]);
      particleH->AddDaughter(leptonsZZ[3]);

      int simpleZZType = 3;
      if (nZZType >= 15 && nZZType <=  16) simpleZZType = 4;

      hDZZSel[ 0+100*simpleZZType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 1+100*simpleZZType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 2+100*simpleZZType]->Fill(TMath::Min(leptons->At(2)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 3+100*simpleZZType]->Fill(TMath::Min(leptons->At(3)->Pt(),199.999),NNLOWeight->GetVal());
      hDZZSel[ 4+100*simpleZZType]->Fill(TMath::Min(dileptonZ1->Mass(),199.99),NNLOWeight->GetVal());
      hDZZSel[ 5+100*simpleZZType]->Fill(TMath::Min(dileptonZ2->Mass(),199.99),NNLOWeight->GetVal());

      DiTauSystem ditau(leptonsZZ[2], leptonsZZ[3], caloMet);

      hDZZtt0Sel[ 0+(simpleZZType-3)*100]->Fill(TMath::Max(TMath::Min(ditau.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZZtt0Sel[ 0+(simpleZZType-3)*100]->Fill(TMath::Max(TMath::Min(ditau.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
        hDZZtt0Sel[ 1+(simpleZZType-3)*100]->Fill(TMath::Min(ditau.RecoMass(),399.99),NNLOWeight->GetVal());
        
	Jet *pt0 = new Jet(leptonsZZ[2]->Px()/ditau.XTau1(),
	                   leptonsZZ[2]->Py()/ditau.XTau1(),
	                   leptonsZZ[2]->Pz()/ditau.XTau1(),
			   leptonsZZ[2]->E()/ditau.XTau1());
	Jet *pt1 = new Jet(leptonsZZ[3]->Px()/ditau.XTau2(),
	                   leptonsZZ[3]->Py()/ditau.XTau2(),
	                   leptonsZZ[3]->Pz()/ditau.XTau2(),
			   leptonsZZ[3]->E()/ditau.XTau2());
        CompositeParticle *particleH2 = new CompositeParticle();
        particleH2->AddDaughter(leptonsZZ[0]);
        particleH2->AddDaughter(leptonsZZ[1]);
        particleH2->AddDaughter(pt0);
        particleH2->AddDaughter(pt1);
        hDZZtt0Sel[ 2+(simpleZZType-3)*100]->Fill(TMath::Min(particleH2->Mass(),999.99),NNLOWeight->GetVal());
        delete pt0;
        delete pt1;
        delete particleH2;
      }

      if(leptons->At(2)->Pt() > 15 &&
         dileptonZ1->Mass() > 12 && dileptonZ2->Mass() > 12){
        hDZZSel[ 6+100*simpleZZType]->Fill(TMath::Min(particleH->Mass(),999.99),NNLOWeight->GetVal());
      }
      delete dileptonZ1;
      delete dileptonZ2;
      delete particleH;
    } // nZZType >= 13 && <= 16
    leptonsZZ.clear();
  } // ZZ->4l analysis

  // ZZ->lltt analysis
  else if(CleanMuons->GetEntries()+CleanElectrons->GetEntries()+
     CleanPFTaus->GetEntries() == 4){
    int nZZType = -2;
    vector<Particle*> leptonsZZ;
    GenericParticle *tau0 = new GenericParticle();
    GenericParticle *tau1 = new GenericParticle();

    tau0->SetMom(CleanPFTaus->At(0)->Px(),
                 CleanPFTaus->At(0)->Py(),
                 CleanPFTaus->At(0)->Pz(),
                 CleanPFTaus->At(0)->E());
    tau0->SetCharge(CleanPFTaus->At(0)->Charge());

    if(CleanPFTaus->GetEntries() >= 2){
      tau1->SetMom(CleanPFTaus->At(1)->Px(),
        	   CleanPFTaus->At(1)->Py(),
        	   CleanPFTaus->At(1)->Pz(),
        	   CleanPFTaus->At(1)->E());
      tau1->SetCharge(CleanPFTaus->At(1)->Charge());
    }

    if(CleanMuons->GetEntries() == 1 && CleanElectrons->GetEntries() == 1 &&
       CleanPFTaus->GetEntries() == 2){
      nZZType = 0;
      if     (CleanMuons->At(0)->Charge()*CleanElectrons->At(0)->Charge() < 0 &&
              tau0->Charge()*tau1->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau0);
        leptonsZZ.push_back(tau1);
      }
      else if(CleanMuons->At(0)->Charge()*tau0->Charge() < 0 &&
              CleanElectrons->At(0)->Charge()*tau1->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(tau0);
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau1);
      }
      else if(CleanElectrons->At(0)->Charge()*tau0->Charge() < 0 &&
              CleanMuons->At(0)->Charge()*tau1->Charge() < 0){
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau0);
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(tau1);
      }
      else if(CleanMuons->At(0)->Charge()*tau1->Charge() < 0 &&
              CleanElectrons->At(0)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(tau1);
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau0);
      }
      else if(CleanElectrons->At(0)->Charge()*tau1->Charge() < 0 &&
              CleanMuons->At(0)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau1);
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(tau0);
      }
      else {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau0);
        leptonsZZ.push_back(tau1);
      }
    } // emhh
    else if(CleanMuons->GetEntries() == 2 &&
            CleanPFTaus->GetEntries() == 2){
      nZZType = 1;
      leptonsZZ.push_back(CleanMuons->At(0));
      leptonsZZ.push_back(CleanMuons->At(1));
      leptonsZZ.push_back(tau0);
      leptonsZZ.push_back(tau1);
    } // mmhh
    else if(CleanElectrons->GetEntries() == 2 &&
            CleanPFTaus->GetEntries() == 2){
      nZZType = 2;
      leptonsZZ.push_back(CleanElectrons->At(0));
      leptonsZZ.push_back(CleanElectrons->At(1));
      leptonsZZ.push_back(tau0);
      leptonsZZ.push_back(tau1);
    } // eehh
    else if(CleanMuons->GetEntries() == 2 && CleanElectrons->GetEntries() == 1 &&
            CleanPFTaus->GetEntries() == 1){
      nZZType = 3;
      leptonsZZ.push_back(CleanMuons->At(0));
      leptonsZZ.push_back(CleanMuons->At(1));
      leptonsZZ.push_back(CleanElectrons->At(0));
      leptonsZZ.push_back(tau0);
    } // mmeh
    else if(CleanMuons->GetEntries() == 1 && CleanElectrons->GetEntries() == 2 &&
            CleanPFTaus->GetEntries() == 1){
      nZZType = 4;
      leptonsZZ.push_back(CleanElectrons->At(0));
      leptonsZZ.push_back(CleanElectrons->At(1));
      leptonsZZ.push_back(CleanMuons->At(0));
      leptonsZZ.push_back(tau0);
    } // eemh
    else if(CleanMuons->GetEntries() == 3 &&
            CleanPFTaus->GetEntries() == 1){
      nZZType = 5;
      if     (CleanMuons->At(0)->Charge()*CleanMuons->At(1)->Charge() < 0 &&
              CleanMuons->At(2)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(tau0);
      }
      else if(CleanMuons->At(0)->Charge()*CleanMuons->At(2)->Charge() < 0 &&
              CleanMuons->At(1)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(tau0);
      }
      else if(CleanMuons->At(1)->Charge()*CleanMuons->At(2)->Charge() < 0 &&
              CleanMuons->At(0)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(tau0);
      }
      else {
        leptonsZZ.push_back(CleanMuons->At(0));
        leptonsZZ.push_back(CleanMuons->At(1));
        leptonsZZ.push_back(CleanMuons->At(2));
        leptonsZZ.push_back(tau0);	
      }
    } // mmmh
    else if(CleanElectrons->GetEntries() == 3 &&
            CleanPFTaus->GetEntries() == 1){
      nZZType = 6;
      if     (CleanElectrons->At(0)->Charge()*CleanElectrons->At(1)->Charge() < 0 &&
              CleanElectrons->At(2)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(tau0);
      }
      else if(CleanElectrons->At(0)->Charge()*CleanElectrons->At(2)->Charge() < 0 &&
              CleanElectrons->At(1)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(tau0);
      }
      else if(CleanElectrons->At(1)->Charge()*CleanElectrons->At(2)->Charge() < 0 &&
              CleanElectrons->At(0)->Charge()*tau0->Charge() < 0){
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(tau0);
      }
      else {
        leptonsZZ.push_back(CleanElectrons->At(0));
        leptonsZZ.push_back(CleanElectrons->At(1));
        leptonsZZ.push_back(CleanElectrons->At(2));
        leptonsZZ.push_back(tau0);
      }
    } // eeeh
    else {
      nZZType = -1;
      printf("Not possible combination: m(%d), e(%d), h(%d)\n", CleanElectrons->GetEntries(),
             CleanMuons->GetEntries(), CleanPFTaus->GetEntries());
      assert(0);
    }
    // types: 0(emhh), 1(mmhh), 2(eehh), 3(mmeh), 4(eemh), 5(mmmh), 6(eeeh)
    CompositeParticle *dileptonZ1 = new CompositeParticle();
    dileptonZ1->AddDaughter(leptonsZZ[0]);
    dileptonZ1->AddDaughter(leptonsZZ[1]);
    CompositeParticle *dileptonZ2 = new CompositeParticle();
    dileptonZ2->AddDaughter(leptonsZZ[2]);
    dileptonZ2->AddDaughter(leptonsZZ[3]);
    hDZXSel[1]->Fill((double)nZZType,NNLOWeight->GetVal());
    hDZZtt1Sel[ 0+100*nZZType]->Fill(dileptonZ1->Charge(),NNLOWeight->GetVal());
    hDZZtt1Sel[ 1+100*nZZType]->Fill(dileptonZ2->Charge(),NNLOWeight->GetVal());
    hDZZtt1Sel[ 2+100*nZZType]->Fill(dileptonZ1->Charge()+dileptonZ2->Charge(),NNLOWeight->GetVal());
    if(dileptonZ1->Charge() == 0 && dileptonZ2->Charge() == 0 &&
       dileptonZ1->Charge()+dileptonZ2->Charge() == 0){
      hDZZtt1Sel[ 3+100*nZZType]->Fill(TMath::Min(dileptonZ1->Mass(),199.99),NNLOWeight->GetVal());
      hDZZtt1Sel[ 4+100*nZZType]->Fill(TMath::Min(dileptonZ2->Mass(),199.99),NNLOWeight->GetVal());
      DiTauSystem ditau0(leptonsZZ[0], leptonsZZ[1], caloMet);
      DiTauSystem ditau1(leptonsZZ[2], leptonsZZ[3], caloMet);

      hDZZtt1Sel[ 5+100*nZZType]->Fill(TMath::Max(TMath::Min(ditau0.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZZtt1Sel[ 5+100*nZZType]->Fill(TMath::Max(TMath::Min(ditau0.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      hDZZtt1Sel[ 6+100*nZZType]->Fill(TMath::Max(TMath::Min(ditau1.XTau1(),1.99),-0.99),NNLOWeight->GetVal());
      hDZZtt1Sel[ 6+100*nZZType]->Fill(TMath::Max(TMath::Min(ditau1.XTau2(),1.99),-0.99),NNLOWeight->GetVal());
      if(dileptonZ1->Mass() > 12 && dileptonZ1->Mass() < 110 &&
         dileptonZ2->Mass() > 12 && dileptonZ2->Mass() < 110){
        hDZZtt1Sel[ 7+100*nZZType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      }
      if(ditau1.XTau1() > 0 && ditau1.XTau2() > 0){
        hDZZtt1Sel[ 8+100*nZZType]->Fill(TMath::Min(ditau1.RecoMass(),399.99),NNLOWeight->GetVal());

	Jet *pt2 = new Jet(leptonsZZ[2]->Px()/ditau1.XTau1(),
	                   leptonsZZ[2]->Py()/ditau1.XTau1(),
	                   leptonsZZ[2]->Pz()/ditau1.XTau1(),
			   leptonsZZ[2]->E()/ditau1.XTau1());
	Jet *pt3 = new Jet(leptonsZZ[3]->Px()/ditau1.XTau2(),
	                   leptonsZZ[3]->Py()/ditau1.XTau2(),
	                   leptonsZZ[3]->Pz()/ditau1.XTau2(),
			   leptonsZZ[3]->E()/ditau1.XTau2());
        CompositeParticle *particleH2 = new CompositeParticle();
        particleH2->AddDaughter(leptonsZZ[0]);
        particleH2->AddDaughter(leptonsZZ[1]);
        particleH2->AddDaughter(pt2);
        particleH2->AddDaughter(pt3);
        hDZZtt1Sel[ 9+100*nZZType]->Fill(TMath::Min(particleH2->Mass(),999.99),NNLOWeight->GetVal());
        delete pt2;
        delete pt3;
        delete particleH2;
      }
      if(ditau0.XTau1() > 0 && ditau0.XTau2() > 0 &&
         ditau1.XTau1() > 0 && ditau1.XTau2() > 0){
        hDZZtt1Sel[10+100*nZZType]->Fill(TMath::Min(ditau1.RecoMass(),399.99),NNLOWeight->GetVal());        
	Jet *pt0 = new Jet(leptonsZZ[0]->Px()/ditau1.XTau1(),
	                   leptonsZZ[0]->Py()/ditau1.XTau1(),
	                   leptonsZZ[0]->Pz()/ditau1.XTau1(),
			   leptonsZZ[0]->E()/ditau1.XTau1());
	Jet *pt1 = new Jet(leptonsZZ[1]->Px()/ditau1.XTau2(),
	                   leptonsZZ[1]->Py()/ditau1.XTau2(),
	                   leptonsZZ[1]->Pz()/ditau1.XTau2(),
			   leptonsZZ[1]->E()/ditau1.XTau2());
	Jet *pt2 = new Jet(leptonsZZ[2]->Px()/ditau0.XTau1(),
	                   leptonsZZ[2]->Py()/ditau0.XTau1(),
	                   leptonsZZ[2]->Pz()/ditau0.XTau1(),
			   leptonsZZ[2]->E()/ditau0.XTau1());
	Jet *pt3 = new Jet(leptonsZZ[3]->Px()/ditau0.XTau2(),
	                   leptonsZZ[3]->Py()/ditau0.XTau2(),
	                   leptonsZZ[3]->Pz()/ditau0.XTau2(),
			   leptonsZZ[3]->E()/ditau0.XTau2());
        CompositeParticle *particleH2 = new CompositeParticle();
        particleH2->AddDaughter(pt0);
        particleH2->AddDaughter(pt1);
        particleH2->AddDaughter(pt2);
        particleH2->AddDaughter(pt3);
        hDZZtt1Sel[11+100*nZZType]->Fill(TMath::Min(particleH2->Mass(),999.99),NNLOWeight->GetVal());
        delete pt0;
        delete pt1;
        delete pt2;
        delete pt3;
        delete particleH2;
      }     
    }
    delete dileptonZ1;
    delete dileptonZ2;
    delete tau0;
    delete tau1;
  }
}
//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fCaloJetName0,                  fCaloJet0);
  ReqBranch(fPFCandidatesName,              fPFCandidates);
  ReqEventObject(fPileupEnergyDensityName,  fPileupEnergyDensity, kTRUE);
  ReqBranch(fTrackName,       fTracks);
  ReqBranch(fGsfTrackName,    fGsfTracks);

  char sb[200];
  sprintf(sb,"hDZZGenSel_0"); hDZZGenSel[ 0] = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDZZGenSel_1"); hDZZGenSel[ 1] = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDZZGenSel_2"); hDZZGenSel[ 2] = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDZZGenSel_3"); hDZZGenSel[ 3] = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDZZGenSel_4"); hDZZGenSel[ 4] = new TH1D(sb,sb,50,0.0,5.0); 
  sprintf(sb,"hDZZGenSel_5"); hDZZGenSel[ 5] = new TH1D(sb,sb,50,0.0,5.0); 
  sprintf(sb,"hDZZGenSel_6"); hDZZGenSel[ 6] = new TH1D(sb,sb,50,0.0,5.0); 
  sprintf(sb,"hDZZGenSel_7"); hDZZGenSel[ 7] = new TH1D(sb,sb,50,0.0,5.0); 
  for(int i=0; i<8; i++){
    AddOutput(hDZZGenSel[i]);
  }

  sprintf(sb,"hDWGGenSel_0"); hDWGGenSel[ 0] = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDWGGenSel_1"); hDWGGenSel[ 1] = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDWGGenSel_2"); hDWGGenSel[ 2] = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDWGGenSel_3"); hDWGGenSel[ 3] = new TH1D(sb,sb,50,0.0,50.0); 
  sprintf(sb,"hDWGGenSel_4"); hDWGGenSel[ 4] = new TH1D(sb,sb,50,0.0,50.0); 
  sprintf(sb,"hDWGGenSel_5"); hDWGGenSel[ 5] = new TH1D(sb,sb,100,0.0,10.0); 
  sprintf(sb,"hDWGGenSel_6"); hDWGGenSel[ 6] = new TH1D(sb,sb,100,0.0,10.0); 
  for(int i=0; i<7; i++){
    AddOutput(hDWGGenSel[i]);
  }

  sprintf(sb,"hDZXSel_0"); hDZXSel[0] = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDZXSel_1"); hDZXSel[1] = new TH1D(sb,sb,9,-2.5,6.5); 
  AddOutput(hDZXSel[0]);
  AddOutput(hDZXSel[1]);

  // WZ histograms
  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWZSel_%d",ind+0);  hDWZSel[ind+0]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+1);  hDWZSel[ind+1]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+2);  hDWZSel[ind+2]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+3);  hDWZSel[ind+3]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+4);  hDWZSel[ind+4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+5);  hDWZSel[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDWZSel_%d",ind+6);  hDWZSel[ind+6]  = new TH1D(sb,sb,100,0.0,20.); 
    sprintf(sb,"hDWZSel_%d",ind+7);  hDWZSel[ind+7]  = new TH1D(sb,sb,200,0.0,800.);
    sprintf(sb,"hDWZSel_%d",ind+8);  hDWZSel[ind+8]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDWZSel_%d",ind+9);  hDWZSel[ind+9]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDWZSel_%d",ind+10); hDWZSel[ind+10] = new TH1D(sb,sb,100,0.0,5.);
    sprintf(sb,"hDWZSel_%d",ind+11); hDWZSel[ind+11] = new TH1D(sb,sb,100,0.0,5.);
    sprintf(sb,"hDWZSel_%d",ind+12); hDWZSel[ind+12] = new TH1D(sb,sb,100,0.0,200.); 
    sprintf(sb,"hDWZSel_%d",ind+13); hDWZSel[ind+13] = new TH1D(sb,sb,100,0.0,200);
    sprintf(sb,"hDWZSel_%d",ind+14); hDWZSel[ind+14] = new TH1D(sb,sb,100,0.0,5.); 
    sprintf(sb,"hDWZSel_%d",ind+15); hDWZSel[ind+15] = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDWZSel_%d",ind+16); hDWZSel[ind+16] = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDWZSel_%d",ind+17); hDWZSel[ind+17] = new TH1D(sb,sb,90,0.0,180.); 
  }

  for(int i=0; i<18; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDWZSel[i+j*100]);
    }
  }

  sprintf(sb,"hDWZSel_90"); hDWZSel[90] = new TH1D(sb,sb,8,-0.5,7.5); 
  AddOutput(hDWZSel[90]);

  // ZZ histograms
  for(int j=0; j<5; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZZSel_%d",ind+0);  hDZZSel[ind+0]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+1);  hDZZSel[ind+1]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+2);  hDZZSel[ind+2]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+3);  hDZZSel[ind+3]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+4);  hDZZSel[ind+4]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+5);  hDZZSel[ind+5]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDZZSel_%d",ind+6);  hDZZSel[ind+6]  = new TH1D(sb,sb,500,0.0,1000.);
  }

  for(int i=0; i<7; i++){
    for(int j=0; j<5; j++){
      AddOutput(hDZZSel[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZZtt0Sel_%d",ind+0);  hDZZtt0Sel[ind+0]  = new TH1D(sb,sb,150,-1.0,2.0);
    sprintf(sb,"hDZZtt0Sel_%d",ind+1);  hDZZtt0Sel[ind+1]  = new TH1D(sb,sb,200,0.0,400.0);
    sprintf(sb,"hDZZtt0Sel_%d",ind+2);  hDZZtt0Sel[ind+2]  = new TH1D(sb,sb,500,0.0,1000.);
  }

  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZZtt0Sel[i+j*100]);
    }
  }

  sprintf(sb,"hDZZSel_90"); hDZZSel[90] = new TH1D(sb,sb,22,-5.5,16.5); 
  sprintf(sb,"hDZZSel_91"); hDZZSel[91] = new TH1D(sb,sb,9,-4.5,4.5); 
  sprintf(sb,"hDZZSel_92"); hDZZSel[92] = new TH1D(sb,sb,9,-4.5,4.5); 
  AddOutput(hDZZSel[90]);
  AddOutput(hDZZSel[91]);
  AddOutput(hDZZSel[92]);

  for(int j=0; j<7; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 0); hDZZtt1Sel[ind+ 0] = new TH1D(sb,sb,5,-2.5,2.5);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 1); hDZZtt1Sel[ind+ 1] = new TH1D(sb,sb,5,-2.5,2.5);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 2); hDZZtt1Sel[ind+ 2] = new TH1D(sb,sb,9,-4.5,4.5);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 3); hDZZtt1Sel[ind+ 3] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 4); hDZZtt1Sel[ind+ 4] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 5); hDZZtt1Sel[ind+ 5] = new TH1D(sb,sb,300,-1,2.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 6); hDZZtt1Sel[ind+ 6] = new TH1D(sb,sb,300,-1,2.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 7); hDZZtt1Sel[ind+ 7] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 8); hDZZtt1Sel[ind+ 8] = new TH1D(sb,sb,200,0.0,400.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+ 9); hDZZtt1Sel[ind+ 9] = new TH1D(sb,sb,200,0.0,1000.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+10); hDZZtt1Sel[ind+10] = new TH1D(sb,sb,200,0.0,400.0);
    sprintf(sb,"hDZZtt1Sel_%d",ind+11); hDZZtt1Sel[ind+11] = new TH1D(sb,sb,200,0.0,1000.0);
  }

  for(int i=0; i<12; i++){
    for(int j=0; j<7; j++){
      AddOutput(hDZZtt1Sel[i+j*100]);
    }
  }

 for(int j=0; j<50; j++){
    sprintf(sb,"hDWZPDF_%d",j); hDWZPDF[j] = new TH1D(sb,sb,100,0,200);
    sprintf(sb,"hDZZPDF_%d",j); hDZZPDF[j] = new TH1D(sb,sb,100,0,200);
 }

 for(int j=0; j<50; j++){
   AddOutput(hDWZPDF[j]);
   AddOutput(hDZZPDF[j]);
 }

 sprintf(sb,"hDZHSel_%d", 0); hDZHSel[ 0] = new TH1D(sb,sb,100,0.0,1.);
 sprintf(sb,"hDZHSel_%d", 1); hDZHSel[ 1] = new TH1D(sb,sb,200,0,200);
 sprintf(sb,"hDZHSel_%d", 2); hDZHSel[ 2] = new TH1D(sb,sb,200,0,200);
 sprintf(sb,"hDZHSel_%d", 3); hDZHSel[ 3] = new TH1D(sb,sb,200,0,200);
 sprintf(sb,"hDZHSel_%d", 4); hDZHSel[ 4] = new TH1D(sb,sb,200,0,200);
 sprintf(sb,"hDZHSel_%d", 5); hDZHSel[ 5] = new TH1D(sb,sb,4,-0.5,3.5);
 sprintf(sb,"hDZHSel_%d", 6); hDZHSel[ 6] = new TH1D(sb,sb,60,0.0,6.0);
 sprintf(sb,"hDZHSel_%d", 7); hDZHSel[ 7] = new TH1D(sb,sb,180,0.0,180.0);
 sprintf(sb,"hDZHSel_%d", 8); hDZHSel[ 8] = new TH1D(sb,sb,180,0.0,180.0);
 sprintf(sb,"hDZHSel_%d", 9); hDZHSel[ 9] = new TH1D(sb,sb,200,0,200);
 sprintf(sb,"hDZHSel_%d",10); hDZHSel[10] = new TH1D(sb,sb,200,0,400);
 sprintf(sb,"hDZHSel_%d",11); hDZHSel[11] = new TH1D(sb,sb,200,0,400);

 for(int j=0; j<12; j++){
   AddOutput(hDZHSel[j]);
 }

 sprintf(sb,"hDHZZForwardSel_%d", 0); hDHZZForwardSel[ 0] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 1); hDHZZForwardSel[ 1] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 2); hDHZZForwardSel[ 2] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 3); hDHZZForwardSel[ 3] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 4); hDHZZForwardSel[ 4] = new TH1D(sb,sb,100,0.0,1.);
 sprintf(sb,"hDHZZForwardSel_%d", 5); hDHZZForwardSel[ 5] = new TH1D(sb,sb,100,0.0,1.);
 sprintf(sb,"hDHZZForwardSel_%d", 6); hDHZZForwardSel[ 6] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 7); hDHZZForwardSel[ 7] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 8); hDHZZForwardSel[ 8] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d", 9); hDHZZForwardSel[ 9] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d",10); hDHZZForwardSel[10] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d",11); hDHZZForwardSel[11] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d",12); hDHZZForwardSel[12] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardSel_%d",13); hDHZZForwardSel[13] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardSel_%d",14); hDHZZForwardSel[14] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardSel_%d",15); hDHZZForwardSel[15] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardSel_%d",16); hDHZZForwardSel[16] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardSel_%d",17); hDHZZForwardSel[17] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d",18); hDHZZForwardSel[18] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardSel_%d",19); hDHZZForwardSel[19] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardSel_%d",20); hDHZZForwardSel[20] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardSel_%d",21); hDHZZForwardSel[21] = new TH1D(sb,sb,5,-0.5,4.5);
 for(int j=0; j<22; j++){
   AddOutput(hDHZZForwardSel[j]);
 }

 sprintf(sb,"hDHZZForwardPSel_%d", 0); hDHZZForwardPSel[ 0] = new TH1D(sb,sb,80,0.0,0.04);
 sprintf(sb,"hDHZZForwardPSel_%d", 1); hDHZZForwardPSel[ 1] = new TH1D(sb,sb,80,0.0,0.04);
 sprintf(sb,"hDHZZForwardPSel_%d", 2); hDHZZForwardPSel[ 2] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 3); hDHZZForwardPSel[ 3] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 4); hDHZZForwardPSel[ 4] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 5); hDHZZForwardPSel[ 5] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 6); hDHZZForwardPSel[ 6] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 7); hDHZZForwardPSel[ 7] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 8); hDHZZForwardPSel[ 8] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d", 9); hDHZZForwardPSel[ 9] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d",10); hDHZZForwardPSel[10] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d",11); hDHZZForwardPSel[11] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d",12); hDHZZForwardPSel[12] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardPSel_%d",13); hDHZZForwardPSel[13] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardPSel_%d",14); hDHZZForwardPSel[14] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardPSel_%d",15); hDHZZForwardPSel[15] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardPSel_%d",16); hDHZZForwardPSel[16] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardPSel_%d",17); hDHZZForwardPSel[17] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d",18); hDHZZForwardPSel[18] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardPSel_%d",19); hDHZZForwardPSel[19] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardPSel_%d",20); hDHZZForwardPSel[20] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardPSel_%d",21); hDHZZForwardPSel[21] = new TH1D(sb,sb,5,-0.5,4.5);
 for(int j=0; j<22; j++){
   AddOutput(hDHZZForwardPSel[j]);
 }

 sprintf(sb,"hDHZZForwardMSel_%d", 0); hDHZZForwardMSel[ 0] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 1); hDHZZForwardMSel[ 1] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 2); hDHZZForwardMSel[ 2] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 3); hDHZZForwardMSel[ 3] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 4); hDHZZForwardMSel[ 4] = new TH1D(sb,sb,100,0.0,1.);
 sprintf(sb,"hDHZZForwardMSel_%d", 5); hDHZZForwardMSel[ 5] = new TH1D(sb,sb,100,0.0,1.);
 sprintf(sb,"hDHZZForwardMSel_%d", 6); hDHZZForwardMSel[ 6] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 7); hDHZZForwardMSel[ 7] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 8); hDHZZForwardMSel[ 8] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d", 9); hDHZZForwardMSel[ 9] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d",10); hDHZZForwardMSel[10] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d",11); hDHZZForwardMSel[11] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d",12); hDHZZForwardMSel[12] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardMSel_%d",13); hDHZZForwardMSel[13] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardMSel_%d",14); hDHZZForwardMSel[14] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardMSel_%d",15); hDHZZForwardMSel[15] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardMSel_%d",16); hDHZZForwardMSel[16] = new TH1D(sb,sb,600,0.0,600.);
 sprintf(sb,"hDHZZForwardMSel_%d",17); hDHZZForwardMSel[17] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d",18); hDHZZForwardMSel[18] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardMSel_%d",19); hDHZZForwardMSel[19] = new TH1D(sb,sb,200,0.0,200.);
 sprintf(sb,"hDHZZForwardMSel_%d",20); hDHZZForwardMSel[20] = new TH1D(sb,sb,50,0.0,5.0);
 sprintf(sb,"hDHZZForwardMSel_%d",21); hDHZZForwardMSel[21] = new TH1D(sb,sb,5,-0.5,4.5);
 for(int j=0; j<22; j++){
   AddOutput(hDHZZForwardMSel[j]);
 }

}
//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::HZZForward(const ParticleOArr *leptons, const PFCandidateCol *fPFCandidates, double weight) {
  
  if(leptons->GetEntries() <= 0) return;

  UInt_t nPFCandZ = 0;
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); i++) {
    const PFCandidate *pf = fPFCandidates->At(i);

    if(pf->HasTrk() || pf->Pt() <= 10 || pf->AbsEta() <= 3.0 || pf->AbsEta() >= 5.0) continue;

    Double_t ptSum = 0.0;
    for (UInt_t j=0; j<fPFCandidates->GetEntries(); j++) {
      const PFCandidate *pfcone = fPFCandidates->At(j);
      if(pf == pfcone) continue;
      Double_t dr = MathUtils::DeltaR(pfcone->Mom(), pf->Mom());
      if ( dr <  0.30 && dr >= 0.01 ) ptSum += pfcone->Pt();
    }
    ptSum = TMath::Max(ptSum - fRho * TMath::Pi() * 0.3 * 0.3, 0.0);
    if(ptSum/pf->Pt() >= 1.0) continue;

    double MinMassZ = 999.;
    int nMinMassZ = -1;
    for (UInt_t nl = 0; nl<leptons->GetEntries(); nl++) {
      CompositeParticle dilepton;
      dilepton.AddDaughter(leptons->At(nl));
      dilepton.AddDaughter(pf);
      if(leptons->At(nl)->ObjType() == kMuon    ) hDHZZForwardSel[0]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(leptons->At(nl)->ObjType() == kElectron) hDHZZForwardSel[1]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(ptSum/pf->Pt() < 0.10) {
        if(leptons->At(nl)->ObjType() == kMuon    ) hDHZZForwardSel[2]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
        if(leptons->At(nl)->ObjType() == kElectron) hDHZZForwardSel[3]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      }
      if(leptons->At(nl)->ObjType() == kMuon     && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardSel[4]->Fill(TMath::Min(ptSum/pf->Pt(),0.999),weight);
      if(leptons->At(nl)->ObjType() == kElectron && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardSel[5]->Fill(TMath::Min(ptSum/pf->Pt(),0.999),weight);
      
      if(leptons->At(nl)->ObjType() == kElectron && ptSum/pf->Pt() < 0.10 && leptons->GetEntries() == 3 &&
         TMath::Abs(dilepton.Mass()-91.1876) < MinMassZ){
	MinMassZ  = TMath::Abs(dilepton.Mass()-91.1876);
	nMinMassZ = nl;
      }
    } // loop over leptons
    if(nMinMassZ >= 0){
      CompositeParticle dilepton;;
      dilepton.AddDaughter(leptons->At(nMinMassZ));
      dilepton.AddDaughter(pf);
      hDHZZForwardSel[6]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      bool isGoodPair = TMath::Abs(dilepton.Mass()-91.1876)<15.0;
      if(isGoodPair == true){
        UInt_t nlep[2] = {0, 0};
        if     (nMinMassZ == 0) {nlep[0] = 1; nlep[1] = 2;}
        else if(nMinMassZ == 1) {nlep[0] = 0; nlep[1] = 2;}
        else if(nMinMassZ == 2) {nlep[0] = 1; nlep[1] = 0;}
        CompositeParticle dilepton2;;
        dilepton2.AddDaughter(leptons->At(nlep[0]));
        dilepton2.AddDaughter(leptons->At(nlep[1]));
        CompositeParticle fourlepton;;
        fourlepton.AddDaughter(leptons->At(0));
        fourlepton.AddDaughter(leptons->At(1));
        fourlepton.AddDaughter(leptons->At(2));
        fourlepton.AddDaughter(pf);
	int type = -1;
	if     (leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() == 0) type = 0;
	else if(leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() != 0) type = 1;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() == 0) type = 2;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() != 0) type = 3;
	else                                                                                                                             type = 4;
	if     (type == 0) hDHZZForwardSel[ 7]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 1) hDHZZForwardSel[ 8]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 2) hDHZZForwardSel[ 9]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 3) hDHZZForwardSel[10]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else               hDHZZForwardSel[11]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	if(dilepton2.Mass() < 91.1876+20.0){
	  if     (type == 0) hDHZZForwardSel[12]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 1) hDHZZForwardSel[13]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 2) hDHZZForwardSel[14]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 3) hDHZZForwardSel[15]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else               hDHZZForwardSel[16]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
          if(type == 0 || type == 2){
            hDHZZForwardSel[17]->Fill(TMath::Min(pf->Pt(),199.999),weight);
            hDHZZForwardSel[18]->Fill(TMath::Min(pf->AbsEta(),4.999),weight);
	  }
	  else {
            hDHZZForwardSel[19]->Fill(TMath::Min(pf->Pt(),199.999),weight);
            hDHZZForwardSel[20]->Fill(TMath::Min(pf->AbsEta(),4.999),weight);
	  }
	  nPFCandZ++;
        }
      }
    }
  } // loop over PFCandidates
  hDHZZForwardSel[21]->Fill(TMath::Min((double)nPFCandZ,4.499),weight);
}

//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::HZZForward(const ParticleOArr *leptons, const PhotonCol *fPhotons, double weight) {

  if(leptons->GetEntries() <= 0) return;

  UInt_t nPHCandZ = 0;
  for (UInt_t i=0; i<fPhotons->GetEntries(); i++) {
    const Photon *ph = fPhotons->At(i);

    if(ph->Pt() <= 10 || ph->AbsEta() > 3.0) continue;

    double MinMassZ = 999.;
    int nMinMassZ = -1;
    for (UInt_t nl = 0; nl<leptons->GetEntries(); nl++) {
      CompositeParticle dilepton;
      dilepton.AddDaughter(leptons->At(nl));
      dilepton.AddDaughter(ph);
      if(leptons->At(nl)->ObjType() == kMuon     && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardPSel[0]->Fill(TMath::Min(ph->CoviEtaiEta(),0.03999),weight);
      if(leptons->At(nl)->ObjType() == kElectron && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardPSel[1]->Fill(TMath::Min(ph->CoviEtaiEta(),0.03999),weight);
      if(leptons->At(nl)->ObjType() == kMuon    ) hDHZZForwardPSel[2]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(leptons->At(nl)->ObjType() == kElectron) hDHZZForwardPSel[3]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(leptons->At(nl)->ObjType() == kMuon     && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardPSel[4]->Fill(TMath::Min(ph->Pt(),199.999),weight);
      if(leptons->At(nl)->ObjType() == kElectron && TMath::Abs(dilepton.Mass()-91.1876)<15.0) hDHZZForwardPSel[5]->Fill(TMath::Min(ph->Pt(),199.999),weight);

      if(leptons->At(nl)->ObjType() == kElectron && leptons->GetEntries() == 3 &&
         TMath::Abs(dilepton.Mass()-91.1876) < MinMassZ){
	MinMassZ  = TMath::Abs(dilepton.Mass()-91.1876);
	nMinMassZ = nl;
      }
    } // loop over leptons
    if(nMinMassZ >= 0){
      CompositeParticle dilepton;;
      dilepton.AddDaughter(leptons->At(nMinMassZ));
      dilepton.AddDaughter(ph);
      hDHZZForwardPSel[6]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      bool isGoodPair = TMath::Abs(dilepton.Mass()-91.1876)<15.0;
      if(isGoodPair == true){
        UInt_t nlep[2] = {0, 0};
        if     (nMinMassZ == 0) {nlep[0] = 1; nlep[1] = 2;}
        else if(nMinMassZ == 1) {nlep[0] = 0; nlep[1] = 2;}
        else if(nMinMassZ == 2) {nlep[0] = 1; nlep[1] = 0;}
        CompositeParticle dilepton2;;
        dilepton2.AddDaughter(leptons->At(nlep[0]));
        dilepton2.AddDaughter(leptons->At(nlep[1]));
        CompositeParticle fourlepton;;
        fourlepton.AddDaughter(leptons->At(0));
        fourlepton.AddDaughter(leptons->At(1));
        fourlepton.AddDaughter(leptons->At(2));
        fourlepton.AddDaughter(ph);
	int type = -1;
	if     (leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() == 0) type = 0;
	else if(leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() != 0) type = 1;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() == 0) type = 2;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() != 0) type = 3;
	else                                                                                                                             type = 4;
	if     (type == 0) hDHZZForwardPSel[ 7]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 1) hDHZZForwardPSel[ 8]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 2) hDHZZForwardPSel[ 9]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 3) hDHZZForwardPSel[10]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else               hDHZZForwardPSel[11]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	if(dilepton2.Mass() < 91.1876+20.0){
	  if     (type == 0) hDHZZForwardPSel[12]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 1) hDHZZForwardPSel[13]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 2) hDHZZForwardPSel[14]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 3) hDHZZForwardPSel[15]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else               hDHZZForwardPSel[16]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
          if(type == 0 || type == 2){
            hDHZZForwardPSel[17]->Fill(TMath::Min(ph->Pt(),199.999),weight);
            hDHZZForwardPSel[18]->Fill(TMath::Min(ph->AbsEta(),4.999),weight);
	  }
	  else {
            hDHZZForwardPSel[19]->Fill(TMath::Min(ph->Pt(),199.999),weight);
            hDHZZForwardPSel[20]->Fill(TMath::Min(ph->AbsEta(),4.999),weight);
	  }
	  nPHCandZ++;
        }
      }
    }
  } // loop over PFCandidates
  hDHZZForwardPSel[21]->Fill(TMath::Min((double)nPHCandZ,4.499),weight);
}

//--------------------------------------------------------------------------------------------------
void ZXEvtSelMod::HZZForward(const ParticleOArr *leptons, const ParticleOArr *fakes, double weight) {

  if(leptons->GetEntries() <= 0) return;

  UInt_t nMUCandZ = 0;
  for (UInt_t i=0; i<fakes->GetEntries(); i++) {
    const Particle *fake = fakes->At(i);

    if(fake->Pt() <= 10 || fake->AbsEta() > 2.5) continue;
    if(fake->ObjType() != kMuon) continue; 

    Double_t ptSum = 0.0;
    for (UInt_t j=0; j<fPFCandidates->GetEntries(); j++) {
      const PFCandidate *pfcone = fPFCandidates->At(j);
      Double_t dr = MathUtils::DeltaR(pfcone->Mom(), fake->Mom());
      if ( dr <  0.30 && dr >= 0.03 ) ptSum += pfcone->Pt();
    }
    ptSum = TMath::Max(ptSum - fRho * TMath::Pi() * 0.3 * 0.3, 0.0);
 
    double MinMassZ = 999.;
    int nMinMassZ = -1;
    for (UInt_t nl = 0; nl<leptons->GetEntries(); nl++) {
      CompositeParticle dilepton;
      dilepton.AddDaughter(leptons->At(nl));
      dilepton.AddDaughter(fake);
      if(dilepton.Charge() != 0) continue;
      if(leptons->At(nl)->ObjType() == kMuon    ) hDHZZForwardMSel[0]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(leptons->At(nl)->ObjType() == kElectron) hDHZZForwardMSel[1]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      if(ptSum/fake->Pt() < 0.10) {
        if(leptons->At(nl)->ObjType() == kMuon    ) hDHZZForwardMSel[2]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
        if(leptons->At(nl)->ObjType() == kElectron) hDHZZForwardMSel[3]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      }

      if(leptons->At(nl)->ObjType() == kMuon && leptons->GetEntries() == 3 &&
         TMath::Abs(dilepton.Mass()-91.1876) < MinMassZ){
	MinMassZ  = TMath::Abs(dilepton.Mass()-91.1876);
	nMinMassZ = nl;
      }
    } // loop over leptons
    if(nMinMassZ >= 0){
      CompositeParticle dilepton;;
      dilepton.AddDaughter(leptons->At(nMinMassZ));
      dilepton.AddDaughter(fake);
      hDHZZForwardMSel[6]->Fill(TMath::Min(dilepton.Mass(),199.999),weight);
      bool isGoodPair = TMath::Abs(dilepton.Mass()-91.1876)<15.0;
      if(isGoodPair == true){
        UInt_t nlep[2] = {0, 0};
        if     (nMinMassZ == 0) {nlep[0] = 1; nlep[1] = 2;}
        else if(nMinMassZ == 1) {nlep[0] = 0; nlep[1] = 2;}
        else if(nMinMassZ == 2) {nlep[0] = 1; nlep[1] = 0;}
        CompositeParticle dilepton2;;
        dilepton2.AddDaughter(leptons->At(nlep[0]));
        dilepton2.AddDaughter(leptons->At(nlep[1]));
        CompositeParticle fourlepton;;
        fourlepton.AddDaughter(leptons->At(0));
        fourlepton.AddDaughter(leptons->At(1));
        fourlepton.AddDaughter(leptons->At(2));
        fourlepton.AddDaughter(fake);
	int type = -1;
	if     (leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() == 0) type = 0;
	else if(leptons->At(nlep[0])->ObjType() == kMuon     && leptons->At(nlep[1])->ObjType() == kMuon     && dilepton2.Charge() != 0) type = 1;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() == 0) type = 2;
	else if(leptons->At(nlep[0])->ObjType() == kElectron && leptons->At(nlep[1])->ObjType() == kElectron && dilepton2.Charge() != 0) type = 3;
	else                                                                                                                             type = 4;
	if     (type == 0) hDHZZForwardMSel[ 7]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 1) hDHZZForwardMSel[ 8]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 2) hDHZZForwardMSel[ 9]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else if(type == 3) hDHZZForwardMSel[10]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	else               hDHZZForwardMSel[11]->Fill(TMath::Min(dilepton2.Mass(),199.999),weight);
	if(dilepton2.Mass() < 91.1876+20.0){
	  if     (type == 0) hDHZZForwardMSel[12]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 1) hDHZZForwardMSel[13]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 2) hDHZZForwardMSel[14]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else if(type == 3) hDHZZForwardMSel[15]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
	  else               hDHZZForwardMSel[16]->Fill(TMath::Min(fourlepton.Mass(),599.999),weight);
          if(type == 0 || type == 2){
            hDHZZForwardMSel[17]->Fill(TMath::Min(fake->Pt(),199.999),weight);
            hDHZZForwardMSel[18]->Fill(TMath::Min(fake->AbsEta(),4.999),weight);
	    hDHZZForwardMSel[4]->Fill(TMath::Min(ptSum/fake->Pt(),0.999),weight);
	  }
	  else {
            hDHZZForwardMSel[19]->Fill(TMath::Min(fake->Pt(),199.999),weight);
            hDHZZForwardMSel[20]->Fill(TMath::Min(fake->AbsEta(),4.999),weight);
	    hDHZZForwardMSel[5]->Fill(TMath::Min(ptSum/fake->Pt(),0.999),weight);
	  }
	  nMUCandZ++;
        }
      }
    }
  } // loop over PFCandidates
  hDHZZForwardMSel[21]->Fill(TMath::Min((double)nMUCandZ,4.499),weight);
}
