// $Id: IsoStudy.cc,v 1.6 2012/01/26 13:17:56 ceballos Exp $

#include "Analysis/SelMods/interface/IsoStudy.h"
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
ClassImp(mithep::IsoStudy)

//--------------------------------------------------------------------------------------------------
IsoStudy::IsoStudy(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsFastSim(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("randomName1"),
  fCleanJetsName("randomName2"),
  fVertexName(ModNames::gkGoodVertexesName),
  fVertices(0),
  fMuonName(Names::gkMuonBrn),
  fMuons(0),
  fElectronName(Names::gkElectronBrn),
  fElectrons(0),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void IsoStudy::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void IsoStudy::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //MuonOArr  *CleanMuonsFakeable        = GetObjThisEvt<MuonOArr>("CleanMuonsFakeable");
  //ElectronOArr *CleanElectronsFakeable = GetObjThisEvt<ElectronOArr>("CleanElectronsFakeable");
  ParticleOArr *leptonsFakeable        = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  ParticleOArr *leptons                = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);

  // No need to go further if no fakeable lepton is found in the event
  // No need to go further if Pt <= 20
  if(leptonsFakeable->GetEntries() >= 2 && leptonsFakeable->At(0)->Pt() > 20 &&
     leptonsFakeable->GetEntries() - leptons->GetEntries() > 0){

    TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

    MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
    ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
    JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
    MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
    const Met *caloMet           = CleanMet->At(0);
    fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

    LoadBranch(fMuonName);
    LoadBranch(fElectronName);
    LoadBranch(fEvtHdrName);

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

    // No need to go further if no good reco lepton is found in the event
    if(leptonsOnlyFake->GetEntries() > 0 && leptonsOnlyFake->GetEntries() <= 10){
      hDisoPresel[0]->Fill(TMath::Min((double)leptons->GetEntries(),4.499) +
                	 5*TMath::Min((double)leptonsOnlyFake->GetEntries(),4.499),NNLOWeight->GetVal());

      Double_t metCorrection[3] = {0.0, 0.0, 0.0};
      Double_t theIso[10] = {999., 999., 999., 999., 999., 999., 999., 999., 999., 999.};
      Double_t theD0[10]  = {999., 999., 999., 999., 999., 999., 999., 999., 999., 999.};
      Int_t idLepVar[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      for(UInt_t i=0; i<leptonsOnlyFake->GetEntries(); i++){
	if     (leptonsOnlyFake->At(i)->ObjType() == kMuon){
          for (UInt_t m=0; m<fMuons->GetEntries(); ++m) {
            const Muon *mu = fMuons->At(m);
	    if(mu->Pt() == leptonsOnlyFake->At(i)->Pt()) {
	      theIso[i] = 1.0 * mu->IsoR03SumPt() + 
                          1.0 * mu->IsoR03EmEt()  +
		          1.0 * mu->IsoR03HadEt();

              const Track *mt = mu->BestTrk();
	      if(mt){
      		for(UInt_t i0 = 0; i0<fVertices->GetEntries(); i0++) {
		  if(fVertices->At(i0)->NTracks() > 0){
      	            Double_t pD0 = mt->D0Corrected(*fVertices->At(i0));
      	            theD0[i] = TMath::Abs(pD0);
      		    break;
		  }
		}
              }

              Double_t eta = TMath::Abs(mu->Eta());
	      if(mu->IsTrackerMuon()) eta = TMath::Abs(mu->TrackerTrk()->Eta());

              if(mu->HasGlobalTrk() == kTRUE                              ) idLepVar[i]+=1;
	      if(mu->BestTrk()->NHits()                               > 10) idLepVar[i]+=2;
	      if(mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight)) idLepVar[i]+=4;
	      if(mu->NSegments()                                       > 0) idLepVar[i]+=8;

	      if(eta > 2.4 || mu->IsTrackerMuon() == kFALSE || mu->HasGlobalTrk() == kFALSE ||
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight)  == kFALSE || mu->Pt() < 10 ||
		 mu->BestTrk()->NHits() <= 10 || mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() >= 10 ||
		 theD0[i] >= 0.2) {
		metCorrection[0] = metCorrection[0] + mu->Px() - mu->EmEnergy()*mu->Px()/mu->P() - mu->HadEnergy()*mu->Px()/mu->P();
		metCorrection[1] = metCorrection[1] + mu->Py() - mu->EmEnergy()*mu->Py()/mu->P() - mu->HadEnergy()*mu->Py()/mu->P();
	      }
	      break;
	    }
	  }
	}
	else if(leptonsOnlyFake->At(i)->ObjType() == kElectron){
          for (UInt_t e=0; e<fElectrons->GetEntries(); ++e) {
            const Electron *el = fElectrons->At(e);
	    if(el->Pt() == leptonsOnlyFake->At(i)->Pt()) {
	      theIso[i] = el->TrackIsolationDr03() + 
                          el->EcalRecHitIsoDr03()  +
			  el->HcalTowerSumEtDr03() - 1.0;

    	      for(UInt_t i0 = 0; i0<fVertices->GetEntries(); i0++) {
    		if(fVertices->At(i0)->NTracks() > 0){
		  Double_t pD0 = el->GsfTrk()->D0Corrected(*fVertices->At(i0));
    		  theD0[i] = TMath::Abs(pD0);
    	          break;
		}
	      }

              if(el->IsEB()){
        	if(el->HadronicOverEm()                             < 0.050) idLepVar[i]+=1;
        	if(el->CoviEtaiEta()                                < 0.010) idLepVar[i]+=2;
        	if(TMath::Abs(el->DeltaPhiSuperClusterTrackAtVtx()) < 0.020) idLepVar[i]+=4;
        	if(TMath::Abs(el->DeltaEtaSuperClusterTrackAtVtx()) < 0.006) idLepVar[i]+=8;
              }
	      else {
        	if(el->HadronicOverEm()                             < 0.025) idLepVar[i]+=1;
        	if(el->CoviEtaiEta()                                < 0.030) idLepVar[i]+=2;
        	if(TMath::Abs(el->DeltaPhiSuperClusterTrackAtVtx()) < 0.020) idLepVar[i]+=4;
        	if(TMath::Abs(el->DeltaEtaSuperClusterTrackAtVtx()) < 0.006) idLepVar[i]+=8;
	      }

	      break;
	    }
	  }        
	}
      }
      metCorrection[2] = sqrt((caloMet->Px()-metCorrection[0])*(caloMet->Px()-metCorrection[0])+
                              (caloMet->Py()-metCorrection[1])*(caloMet->Py()-metCorrection[1]));

      if(metCorrection[2] > 30.0){
	hDisoPresel[1]->Fill(TMath::Min((double)leptons->GetEntries(),4.499) +
                	   5*TMath::Min((double)leptonsOnlyFake->GetEntries(),4.499),NNLOWeight->GetVal());
      }

      // Sort and count the number of central Jets for vetoing
      int nCentralJets = 0;
      for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
	const Jet *jet = CleanJets->At(i);        
	if(TMath::Abs(jet->Eta()) < fEtaJetCut &&
	   jet->Pt() > fPtJetCut){
	  double deltaRMin = 999.0;
          for(UInt_t il=0; il<leptonsOnlyFake->GetEntries(); il++){
	    double dR = MathUtils::DeltaR(leptonsOnlyFake->At(il)->Mom(), jet->Mom());
	    if(dR < deltaRMin) deltaRMin = dR;
	  }
	  if(deltaRMin >= 0.5) nCentralJets++;	  
	  hDisoPresel[2]->Fill(deltaRMin,NNLOWeight->GetVal());
	}
      }
      hDisoPresel[3]->Fill(TMath::Min((double)nCentralJets,9.4999),NNLOWeight->GetVal());
      int lType = -1;
      if        (leptons->GetEntries() == 1 && leptonsOnlyFake->GetEntries() >= 1){
	if (leptons->At(0)->ObjType() == kMuon && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 0;
	else if(leptons->At(0)->ObjType() == kElectron && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 1;
	else if(leptons->At(0)->ObjType() == kElectron && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 2;
	else if(leptons->At(0)->ObjType() == kMuon && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 3;
	else {
	  cout << "Hey, this is not possible, leptonTypes1: "
    	       << leptons->At(0)->ObjType() << " - " 
               << leptonsOnlyFake->At(0)->ObjType() << endl;
	}
      }
      else if(leptons->GetEntries() == 0 && leptonsOnlyFake->GetEntries() >= 2){
	if     (leptonsOnlyFake->At(0)->ObjType() == kMuon     && leptonsOnlyFake->At(1)->ObjType() == kMuon)
	  lType = 4;
	else if(leptonsOnlyFake->At(0)->ObjType() == kElectron && leptonsOnlyFake->At(1)->ObjType() == kElectron)
	  lType = 5;
	else if((leptonsOnlyFake->At(0)->ObjType() == kMuon     && leptonsOnlyFake->At(1)->ObjType() == kElectron) ||
        	(leptonsOnlyFake->At(0)->ObjType() == kElectron && leptonsOnlyFake->At(1)->ObjType() == kMuon))
	  lType = 6;
	else {
	  cout << "Hey, this is not possible, leptonTypes2: "
    	       << leptonsOnlyFake->At(0)->ObjType() << " - " 
               << leptonsOnlyFake->At(1)->ObjType() << endl;
	}
      }
      else if(leptons->GetEntries() == 2 && leptonsOnlyFake->GetEntries() >= 1){
	if (CleanMuons->GetEntries() == 2 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 7;
	else if (CleanMuons->GetEntries() == 2 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 8;
	else if (CleanElectrons->GetEntries() == 2 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 9;
	else if (CleanElectrons->GetEntries() == 2 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 10;
	else if (CleanMuons->GetEntries() == 1 && CleanElectrons->GetEntries() == 1 && 
        	 leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 11;
	else if (CleanMuons->GetEntries() == 1 && CleanElectrons->GetEntries() == 1 && 
        	 leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 12;
	else {
	  cout << "Hey, this is not possible, leptonTypes3" << endl;
	}
      }
      else if(leptons->GetEntries() == 3 && leptonsOnlyFake->GetEntries() >= 1){
	if (CleanMuons->GetEntries() == 3 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 13;
	else if (CleanMuons->GetEntries() == 3 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 14;
	else if (CleanMuons->GetEntries() == 2 && CleanElectrons->GetEntries() == 1 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 15;
	else if (CleanMuons->GetEntries() == 2 && CleanElectrons->GetEntries() == 1 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 16;
	else if (CleanElectrons->GetEntries() == 3 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 17;
	else if (CleanElectrons->GetEntries() == 3 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 18;
	else if (CleanElectrons->GetEntries() == 2 && CleanMuons->GetEntries() == 1 && leptonsOnlyFake->At(0)->ObjType() == kMuon)
	  lType = 19;
	else if (CleanElectrons->GetEntries() == 2 && CleanMuons->GetEntries() == 1 && leptonsOnlyFake->At(0)->ObjType() == kElectron)
	  lType = 20;
	else {
	  cout << "Hey, this is not possible, leptonTypes4" << endl;
	}
      }
      hDisoPresel[4]->Fill(lType,NNLOWeight->GetVal());

      if     (lType >= 0 && lType <= 3 && leptons->At(0)->Pt() > 20 && leptonsOnlyFake->At(0)->Pt() > 20){
	CompositeParticle dilepton;
	dilepton.AddDaughter(leptons->At(0));
	dilepton.AddDaughter(leptonsOnlyFake->At(0));
	if(dilepton.Mass() > 12.0){
	  double deltaPhiMetLepton[2] = {fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi())),
    	  				 fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptonsOnlyFake->At(0)->Phi()))};

	  double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
	     deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
	  hDiso2lSel[ 0+10*lType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso2lSel[ 1+10*lType]->Fill(TMath::Min(metCorrection[2],199.999),NNLOWeight->GetVal());
	  hDiso2lSel[ 2+10*lType]->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
	  if(theIso[0]/leptonsOnlyFake->At(0)->Pt() < 0.15)
            hDiso2lSel[ 3+10*lType]->Fill((double)idLepVar[0],NNLOWeight->GetVal());
	  if(idLepVar[0] == 15)
            hDiso2lSel[ 4+10*lType]->Fill(TMath::Min(theIso[0]/leptonsOnlyFake->At(0)->Pt(),1.999),NNLOWeight->GetVal());
	  hDiso2lSel[ 5+10*lType]->Fill(TMath::Min(theD0[0],0.0499),NNLOWeight->GetVal());
	  if(metCorrection[2] > 30 && dilepton.Mass() > 12 && dilepton.Charge() == 0){ 
	    hDiso2lSel[ 6+10*lType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
	    hDiso2lSel[ 7+10*lType]->Fill(TMath::Min(leptonsOnlyFake->At(0)->Pt(),199.999),NNLOWeight->GetVal());
	    hDiso2lSel[ 8+10*lType]->Fill(nCentralJets,NNLOWeight->GetVal());
	    hDiso2lSel[ 9+10*lType]->Fill(minDeltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	  }
	}
      }
      else if(lType >= 4 && lType <= 6 && leptonsOnlyFake->At(1)->Pt() > 20){
	CompositeParticle dilepton;
	dilepton.AddDaughter(leptonsOnlyFake->At(0));
	dilepton.AddDaughter(leptonsOnlyFake->At(1));
	double deltaPhiMetLepton[2] = {fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptonsOnlyFake->At(0)->Phi())),
    				       fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptonsOnlyFake->At(1)->Phi()))};

	double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
	   deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
	hDiso2lSel[ 0+10*lType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
	hDiso2lSel[ 1+10*lType]->Fill(TMath::Min(metCorrection[2],199.999),NNLOWeight->GetVal());
	hDiso2lSel[ 2+10*lType]->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
	if(theIso[1]/leptonsOnlyFake->At(1)->Pt() < 0.15)
          hDiso2lSel[ 3+10*lType]->Fill((double)idLepVar[1],NNLOWeight->GetVal());
	if(idLepVar[1] == 15)
          hDiso2lSel[ 4+10*lType]->Fill(TMath::Min(theIso[1]/leptonsOnlyFake->At(1)->Pt(),1.999),NNLOWeight->GetVal());
	hDiso2lSel[ 5+10*lType]->Fill(TMath::Min(theD0[1],0.0499),NNLOWeight->GetVal());
	if(metCorrection[2] > 30 && dilepton.Mass() > 12 && dilepton.Charge() == 0){ 
	  hDiso2lSel[ 6+10*lType]->Fill(TMath::Min(leptonsOnlyFake->At(0)->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso2lSel[ 7+10*lType]->Fill(TMath::Min(leptonsOnlyFake->At(1)->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso2lSel[ 8+10*lType]->Fill(nCentralJets,NNLOWeight->GetVal());
	  hDiso2lSel[ 9+10*lType]->Fill(minDeltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	}
      }
      else if(lType >= 7 && lType <= 12){
	CompositeParticle trilepton;
	trilepton.AddDaughter(leptons->At(0));
	trilepton.AddDaughter(leptons->At(1));
	trilepton.AddDaughter(leptonsOnlyFake->At(0));
        /*
	printf("CANDIDATE3l: %d %d %d %d | %f %f %f %d %f | %f %f %f %d %f | %f %f %f %d %f - %f\n",
    	    lType,fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
            leptons->At(0)->Pt(),leptons->At(0)->Eta(),leptons->At(0)->Phi(),leptons->At(0)->ObjType(),leptons->At(0)->Charge(),
            leptons->At(1)->Pt(),leptons->At(1)->Eta(),leptons->At(1)->Phi(),leptons->At(1)->ObjType(),leptons->At(1)->Charge(),
            leptonsOnlyFake->At(0)->Pt(),leptonsOnlyFake->At(0)->Eta(),leptonsOnlyFake->At(0)->Phi(),leptonsOnlyFake->At(0)->ObjType(),leptonsOnlyFake->At(0)->Charge(),
            metCorrection[2]);
	*/
	double deltaPhiMetLepton[3] = {fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi())),
    				       fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(1)->Phi())),
    				       fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptonsOnlyFake->At(0)->Phi()))};

	double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
	   deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
	if(minDeltaPhiMetLepton > deltaPhiMetLepton[2]) minDeltaPhiMetLepton = deltaPhiMetLepton[2];
	TVector3 v0(leptons->At(0)->Px()        ,leptons->At(0)->Py()        ,leptons->At(0)->Pz());
	TVector3 v1(leptons->At(1)->Px()        ,leptons->At(1)->Py()        ,leptons->At(1)->Pz());
	TVector3 v2(leptonsOnlyFake->At(0)->Px(),leptonsOnlyFake->At(0)->Py(),leptonsOnlyFake->At(0)->Pz());
	hDiso3lSel[ 0+10*(lType-7)]->Fill(TMath::Max(TMath::Max(v0.Angle(v1),v0.Angle(v2)),v1.Angle(v2)) * 180./TMath::Pi(),NNLOWeight->GetVal());
	hDiso3lSel[ 1+10*(lType-7)]->Fill(TMath::Min(metCorrection[2],199.999),NNLOWeight->GetVal());
	hDiso3lSel[ 2+10*(lType-7)]->Fill(trilepton.Charge(),NNLOWeight->GetVal());
	hDiso3lSel[ 3+10*(lType-7)]->Fill((double)idLepVar[0],NNLOWeight->GetVal());
	if(metCorrection[2] > 30 && abs(trilepton.Charge()) == 1){
	  hDiso3lSel[ 4+10*(lType-7)]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso3lSel[ 5+10*(lType-7)]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso3lSel[ 6+10*(lType-7)]->Fill(TMath::Min(leptonsOnlyFake->At(0)->Pt(),199.999),NNLOWeight->GetVal());
	  hDiso3lSel[ 7+10*(lType-7)]->Fill(minDeltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());

          CompositeParticle dilepton0;
          dilepton0.AddDaughter(leptons->At(0));
          dilepton0.AddDaughter(leptons->At(1));
          CompositeParticle dilepton1;
          dilepton1.AddDaughter(leptons->At(0));
          dilepton1.AddDaughter(leptonsOnlyFake->At(0));
          CompositeParticle dilepton2;
          dilepton2.AddDaughter(leptons->At(1));
          dilepton2.AddDaughter(leptonsOnlyFake->At(0));

	  hDiso3lSel[ 8+10*(lType-7)]->Fill(TMath::Min(dilepton0.Mass(),199.999),NNLOWeight->GetVal());
	  hDiso3lSel[ 9+10*(lType-7)]->Fill(TMath::Min(dilepton1.Mass(),199.999),NNLOWeight->GetVal());
	  hDiso3lSel[ 9+10*(lType-7)]->Fill(TMath::Min(dilepton2.Mass(),199.999),NNLOWeight->GetVal());
	}
      }
    } // leptonsOnlyFake >= 0
    delete leptonsOnlyFake;
  } // At least one fakeable object with pt>20
}
//--------------------------------------------------------------------------------------------------
void IsoStudy::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMuonName,     fMuons);
  ReqBranch(fElectronName, fElectrons);
  ReqBranch(fEvtHdrName,   fEventHeader);

  char sb[200];
  sprintf(sb,"hDisoPresel_%d",0); hDisoPresel[0] = new TH1D(sb,sb,50,-0.5,49.5); 
  sprintf(sb,"hDisoPresel_%d",1); hDisoPresel[1] = new TH1D(sb,sb,50,-0.5,49.5); 
  sprintf(sb,"hDisoPresel_%d",2); hDisoPresel[2] = new TH1D(sb,sb,100,-0.5,9.5);
  sprintf(sb,"hDisoPresel_%d",3); hDisoPresel[3] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDisoPresel_%d",4); hDisoPresel[4] = new TH1D(sb,sb,21,-0.5,20.5);

  for(int i=0; i<5; i++){
    AddOutput(hDisoPresel[i]);
  }

  for(int j=0; j<7; j++){
    int ind = 10 * j;
    sprintf(sb,"hDiso2lSel_%d",ind+ 0); hDiso2lSel[ind+ 0] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 1); hDiso2lSel[ind+ 1] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 2); hDiso2lSel[ind+ 2] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 3); hDiso2lSel[ind+ 3] = new TH1D(sb,sb,16,-0.5,15.5);
    sprintf(sb,"hDiso2lSel_%d",ind+ 4); hDiso2lSel[ind+ 4] = new TH1D(sb,sb,200,0.0,2.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 5); hDiso2lSel[ind+ 5] = new TH1D(sb,sb,50,0.0,0.05);
    sprintf(sb,"hDiso2lSel_%d",ind+ 6); hDiso2lSel[ind+ 6] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 7); hDiso2lSel[ind+ 7] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso2lSel_%d",ind+ 8); hDiso2lSel[ind+ 8] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDiso2lSel_%d",ind+ 9); hDiso2lSel[ind+ 9] = new TH1D(sb,sb,90,0.0,180.);
  }

  for(int j=0; j<7; j++){
    for(int i=0; i<10; i++){
      AddOutput(hDiso2lSel[i+j*10]);
    }
  }

  for(int j=0; j<6; j++){
    int ind = 10 * j;
    sprintf(sb,"hDiso3lSel_%d",ind+ 0); hDiso3lSel[ind+ 0] = new TH1D(sb,sb,360,0.0,180.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 1); hDiso3lSel[ind+ 1] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 2); hDiso3lSel[ind+ 2] = new TH1D(sb,sb,7,-3.5,3.5);
    sprintf(sb,"hDiso3lSel_%d",ind+ 3); hDiso3lSel[ind+ 3] = new TH1D(sb,sb,16,-0.5,15.5);
    sprintf(sb,"hDiso3lSel_%d",ind+ 4); hDiso3lSel[ind+ 4] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 5); hDiso3lSel[ind+ 5] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 6); hDiso3lSel[ind+ 6] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 7); hDiso3lSel[ind+ 7] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 8); hDiso3lSel[ind+ 8] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDiso3lSel_%d",ind+ 9); hDiso3lSel[ind+ 9] = new TH1D(sb,sb,200,0.0,200.);
  }

  for(int j=0; j<6; j++){
    for(int i=0; i<10; i++){
      AddOutput(hDiso3lSel[i+j*10]);
    }
  }

}

//--------------------------------------------------------------------------------------------------
void IsoStudy::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void IsoStudy::Terminate()
{
  // Run finishing code on the client computer
}
