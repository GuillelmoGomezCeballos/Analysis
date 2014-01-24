// $Id: WlnFakeSelMod.cc,v 1.5 2012/04/18 14:59:41 ceballos Exp $

#include "Analysis/SelMods/interface/WlnFakeSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/StableDataCol.h"
#include "MitAna/DataTree/interface/MetCol.h"

using namespace mithep;
ClassImp(mithep::WlnFakeSelMod)

//--------------------------------------------------------------------------------------------------
WlnFakeSelMod::WlnFakeSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fTrackName(Names::gkTrackBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fTracks(0),
  fApplyConversionFilter(kFALSE),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WlnFakeSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WlnFakeSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "WlnFakeSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  //Get Generator Level information for matching
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  //JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet	       = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet	       = CleanMet->At(0);
  MCParticleOArr *GenRadPhotons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCRadPhotonsName);
  MCParticleOArr *GenISRPhotons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCISRPhotonsName);
  MCParticleOArr *GenLeptons    = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
  MCParticleOArr *GenAllLeptons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCAllLeptonsName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  // At least one identified lepton in the event
  hDWlnFakeSel[50]->Fill(leptons->GetEntries(),NNLOWeight->GetVal());
  if(leptons->GetEntries() != 1)  return;
  if(leptons->At(0)->Pt() <= 20)  return;
  if(caloMet->Pt() <= 30)  return;

  int theFirstGoodLepton = 0;
  ParticleOArr *leptonsFakeable        = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  int NleptonsOnlyFake = 0;
  for(UInt_t i=0; i<leptonsFakeable->GetEntries(); i++){
    Bool_t isOnlyFake = kTRUE;
    for(UInt_t j=0; j<leptons->GetEntries(); j++) {
      if(leptons->At(j) == leptonsFakeable->At(i)) {
  	isOnlyFake = kFALSE;
        break;
      }
    }
    if(isOnlyFake == kTRUE) NleptonsOnlyFake++;
  }
  if(NleptonsOnlyFake <= 0) return;
  
  double theQ = leptons->At(0)->Charge();

  LoadBranch(fTrackName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  
  int NFakeableNoIso = 0;
  int NFakeableIso = 0;
  int index = 0;
  for (UInt_t i=0; i<fTracks->GetEntries(); i++) {
    const Track *track = fTracks->At(i);
    
    hDWlnFakeSel[51]->Fill(TMath::Min(MathUtils::DeltaR(track->Phi(), track->Eta(), 
                           leptons->At(theFirstGoodLepton)->Phi(), leptons->At(theFirstGoodLepton)->Eta()),5.999),NNLOWeight->GetVal());
    if(MathUtils::DeltaR(track->Phi(), track->Eta(),
       leptons->At(theFirstGoodLepton)->Phi(), leptons->At(theFirstGoodLepton)->Eta()) < 0.1)
       continue;

    if(track->Pt() <= 10) continue;
    if(fabs(track->Eta()) >= 2.5) continue;

    fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
    double d0_real = 99999;
    for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
      if(fVertices->At(i0)->NTracks() > 0){
        double pD0 = track->D0Corrected(*fVertices->At(i0));
        d0_real = TMath::Abs(pD0);
        break;
      }
    }
    if(d0_real >= 0.025) continue;

    // Apply Conversion Filter
    Bool_t isGoodConversion = kFALSE;
    if (fApplyConversionFilter) {
      LoadBranch(fConversionBranchName);
      for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {
        
	Bool_t ConversionMatchFound = kFALSE;
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
	               (fConversions->At(ifc)->Daughter(d))->Trk();
          if (track == trk) {
            ConversionMatchFound = kTRUE;
            break;
          }
        }

        // if match between the e-track and one of the conversion legs
        if (ConversionMatchFound == kTRUE){
          isGoodConversion =  (fConversions->At(ifc)->Prob() > 0.0005) &&
                              (fConversions->At(ifc)->Lxy() > 0) &&
                              (fConversions->At(ifc)->Lz() > 0);

          if (isGoodConversion == kTRUE) {
	    for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
              const Track *trk = dynamic_cast<const ChargedParticle*>
	                   (fConversions->At(ifc)->Daughter(d))->Trk();
            	      
              if (trk) {
            	if (!(trk->NHits() > 8 && trk->Prob() > 0.005) && d == 0)
            	  isGoodConversion = kFALSE;
              
            	const StableData *sd = dynamic_cast<const StableData*>
		                  (fConversions->At(ifc)->DaughterDat(d));
            	if (sd->NWrongHits() != 0)
            	  isGoodConversion = kFALSE;
              
              } else {
            	isGoodConversion = kFALSE;
              }
            }
	  }
        }

        if (isGoodConversion == kTRUE) break;

      } // loop over all conversions 
      
    }
    if (isGoodConversion == kTRUE) continue;

    bool isGenLepton = kFALSE;
    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
      MCParticle *gen = GenLeptons->At(j);
      hDWlnFakeSel[52]->Fill(TMath::Min(MathUtils::DeltaR(gen->Phi(), gen->Eta(), track->Phi(), track->Eta()),5.999),NNLOWeight->GetVal());
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), track->Phi(), track->Eta()) < 0.1)
        isGenLepton = kTRUE;
    }
    if (isGenLepton == kTRUE) continue;
 
    UInt_t trackGenType = 0;
    // look for real W/Z -> l leptons
    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
      MCParticle *gen = GenLeptons->At(j);
      if(MathUtils::DeltaR(gen->Mom(), track->Mom()) < 0.05) {
        trackGenType = 1;
        break;
      }
    }
    // look for other leptons
    if(trackGenType == 0){
      for (UInt_t j=0; j<GenAllLeptons->GetEntries(); ++j) {
    	MCParticle *gen = GenAllLeptons->At(j);
    	if(MathUtils::DeltaR(gen->Mom(), track->Mom()) < 0.05) {
          trackGenType = 2;
          break;
        }
      }
    }
    // look for ISR photons
    if(trackGenType == 0){
      for (UInt_t j=0; j<GenISRPhotons->GetEntries(); ++j) {
    	MCParticle *gen = GenISRPhotons->At(j);
    	if(MathUtils::DeltaR(gen->Mom(), track->Mom()) < 0.05) {
          trackGenType = 3;
          break;
        }
      }
    }
    // look for Rad photons
    if(trackGenType == 0){
      for (UInt_t j=0; j<GenRadPhotons->GetEntries(); ++j) {
    	MCParticle *gen = GenRadPhotons->At(j);
    	if(MathUtils::DeltaR(gen->Mom(), track->Mom()) < 0.05) {
          trackGenType = 4;
          break;
        }
      }
    }
    hDWlnFakeGenType[0]->Fill(trackGenType,NNLOWeight->GetVal());

    NFakeableNoIso++;

    index = 0;
    if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[0 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[1 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[2 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[3 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());

    if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[4 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[5 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[6 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[7 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());

    Double_t trackIsolation = IsolationTools::TrackIsolation(track, 0.3, 0.015, 1.0, 0.5, fTracks);
    if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[80+index]->Fill(TMath::Min(trackIsolation,99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[81+index]->Fill(TMath::Min(trackIsolation,99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[82+index]->Fill(TMath::Min(trackIsolation,99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[83+index]->Fill(TMath::Min(trackIsolation,99.999),NNLOWeight->GetVal());

    bool isOneMuon[3] = {kFALSE, kFALSE, kFALSE};
    double isoMuon[3] = {-1.0, -1.0, -1.0};
    for (UInt_t j=0; j<fMuons->GetEntries(); ++j) {
      const Track *muonTrack = fMuons->At(j)->Trk();
      if(!muonTrack) continue;
      if(MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         muonTrack->Charge() * track->Charge() > 0){
	isOneMuon[0] = kTRUE;
        isoMuon[0] = fMuons->At(j)->IsoR03SumPt() + fMuons->At(j)->IsoR03EmEt() + fMuons->At(j)->IsoR03HadEt();
      }

      if(!fMuons->At(j)->GlobalTrk()) continue;
      if(MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         muonTrack->Charge() * track->Charge() > 0){
	isOneMuon[1] = kTRUE;
        isoMuon[1] = fMuons->At(j)->IsoR03SumPt() + fMuons->At(j)->IsoR03EmEt() + fMuons->At(j)->IsoR03HadEt();
      }

      if(! myMuonTools.IsGood(fMuons->At(j), MuonTools::kTMOneStationLoose) ||
	 ! myMuonTools.IsGood(fMuons->At(j), MuonTools::kTM2DCompatibilityLoose)) continue;
      if(MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         muonTrack->Charge() * track->Charge() > 0){
	isOneMuon[2] = kTRUE;
        isoMuon[2] = fMuons->At(j)->IsoR03SumPt() + fMuons->At(j)->IsoR03EmEt() + fMuons->At(j)->IsoR03HadEt();
      }
    }
   
    if(isOneMuon[0] == kTRUE){
      hDWlnFakeGenType[1]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 100;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[80+index]->Fill(TMath::Min(isoMuon[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[81+index]->Fill(TMath::Min(isoMuon[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[82+index]->Fill(TMath::Min(isoMuon[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[83+index]->Fill(TMath::Min(isoMuon[0],99.999),NNLOWeight->GetVal());
    }
    if(isOneMuon[1] == kTRUE){
      hDWlnFakeGenType[2]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 200;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[80+index]->Fill(TMath::Min(isoMuon[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[81+index]->Fill(TMath::Min(isoMuon[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[82+index]->Fill(TMath::Min(isoMuon[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[83+index]->Fill(TMath::Min(isoMuon[1],99.999),NNLOWeight->GetVal());
    }
    if(isOneMuon[2] == kTRUE){
      hDWlnFakeGenType[3]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 300;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[80+index]->Fill(TMath::Min(isoMuon[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[81+index]->Fill(TMath::Min(isoMuon[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[82+index]->Fill(TMath::Min(isoMuon[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[83+index]->Fill(TMath::Min(isoMuon[2],99.999),NNLOWeight->GetVal());
    }

    bool isOneElectron[3] = {kFALSE, kFALSE, kFALSE};
    double isoElectron[3] = {-1.0, -1.0, -1.0};
    for (UInt_t j=0; j<fElectrons->GetEntries(); ++j) {
      const Track *electronTrack = electronTrack = fElectrons->At(j)->BestTrk();
      if(!electronTrack) continue;
      if(MathUtils::DeltaR(electronTrack->Phi(), electronTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         electronTrack->Charge() * track->Charge() > 0){
	isOneElectron[0] = kTRUE;
        isoElectron[0] = fElectrons->At(j)->TrackIsolationDr03() + fElectrons->At(j)->EcalRecHitIsoDr04();
      }

      if(fElectrons->At(j)->PassLooseID()){
        if(MathUtils::DeltaR(electronTrack->Phi(), electronTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
           electronTrack->Charge() * track->Charge() > 0){
	  isOneElectron[1] = kTRUE;
          isoElectron[1] = fElectrons->At(j)->TrackIsolationDr03() + fElectrons->At(j)->EcalRecHitIsoDr04();
        }
      }

      if(fElectrons->At(j)->PassTightID()){
        if(MathUtils::DeltaR(electronTrack->Phi(), electronTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
           electronTrack->Charge() * track->Charge() > 0){
	  isOneElectron[2] = kTRUE;
          isoElectron[2] = fElectrons->At(j)->TrackIsolationDr03() + fElectrons->At(j)->EcalRecHitIsoDr04();
        }
      }
    }

    if(isOneElectron[0] == kTRUE){
      hDWlnFakeGenType[4]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 100;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[90+index]->Fill(TMath::Min(isoElectron[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[91+index]->Fill(TMath::Min(isoElectron[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[92+index]->Fill(TMath::Min(isoElectron[0],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[93+index]->Fill(TMath::Min(isoElectron[0],99.999),NNLOWeight->GetVal());
    }
    if(isOneElectron[1] == kTRUE){
      hDWlnFakeGenType[5]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 200;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[90+index]->Fill(TMath::Min(isoElectron[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[91+index]->Fill(TMath::Min(isoElectron[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[92+index]->Fill(TMath::Min(isoElectron[1],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[93+index]->Fill(TMath::Min(isoElectron[1],99.999),NNLOWeight->GetVal());
    }
    if(isOneElectron[2] == kTRUE){
      hDWlnFakeGenType[6]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 300;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[90+index]->Fill(TMath::Min(isoElectron[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[91+index]->Fill(TMath::Min(isoElectron[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[92+index]->Fill(TMath::Min(isoElectron[2],99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[93+index]->Fill(TMath::Min(isoElectron[2],99.999),NNLOWeight->GetVal());
    }
    if(trackIsolation >= 10) continue;
    hDWlnFakeGenType[7]->Fill(trackGenType,NNLOWeight->GetVal());
    NFakeableIso++;

    index = 100;
    if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[0 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[1 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[2 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[3 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());

    if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[4 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[5 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[6 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[7 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());

    bool isMuon = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {  
      const Track *muonTrack = CleanMuons->At(j)->Trk();
      if(MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         muonTrack->Charge() * track->Charge() > 0) isMuon = kTRUE;
      else if(MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), track->Phi(), track->Eta()) < 0.1)
         printf("goodMu: %f %f %f, trkMu: %f %f %f\n",muonTrack->Eta(), muonTrack->Phi(), muonTrack->Pt(),
	                                              track->Eta(), track->Phi(), track->Pt());        
    }
    if(isMuon == kTRUE){
      hDWlnFakeGenType[8]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 200;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[0 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[1 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[2 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[3 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());

      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[4 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[5 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[6 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[7 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    }

    bool isElectron = kFALSE;
    for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {  
      const Track *electronTrack = electronTrack = CleanElectrons->At(j)->BestTrk();
      if(MathUtils::DeltaR(electronTrack->Phi(), electronTrack->Eta(), track->Phi(), track->Eta()) < 0.1 &&
         electronTrack->Charge() * track->Charge() > 0) isElectron = kTRUE;
      else if(MathUtils::DeltaR(electronTrack->Phi(), electronTrack->Eta(), track->Phi(), track->Eta()) < 0.1)
         printf("goodEl: %f %f %f, trkEl: %f %f %f\n",electronTrack->Eta(), electronTrack->Phi(), electronTrack->Pt(),
	                                              track->Eta(), track->Phi(), track->Pt());
    }
    if(isElectron == kTRUE){
      hDWlnFakeGenType[9]->Fill(trackGenType,NNLOWeight->GetVal());
      index = 300;
      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[0 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[1 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[2 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[3 +index]->Fill(TMath::Min(track->Pt(),99.999),NNLOWeight->GetVal());

      if     (theQ * track->Charge() > 0  && theQ > 0) hDWlnFakeSel[4 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() > 0  && theQ < 0) hDWlnFakeSel[5 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ > 0) hDWlnFakeSel[6 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
      else if(theQ * track->Charge() < 0  && theQ < 0) hDWlnFakeSel[7 +index]->Fill(TMath::Min(fabs(track->Eta()),2.499),NNLOWeight->GetVal());
    }
  }
  hDWlnFakeSel[53]->Fill(NFakeableNoIso,NNLOWeight->GetVal());
  hDWlnFakeSel[54]->Fill(NFakeableIso,NNLOWeight->GetVal());
}
//--------------------------------------------------------------------------------------------------
void WlnFakeSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fTrackName,            fTracks);
  ReqBranch(fMuonName,             fMuons);
  ReqBranch(fElectronName,         fElectrons);
  ReqBranch(fConversionBranchName, fConversions);

  char sb[200];
  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWlnFakeSel_%d",ind+0);  hDWlnFakeSel[ind+0]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+1);  hDWlnFakeSel[ind+1]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+2);  hDWlnFakeSel[ind+2]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+3);  hDWlnFakeSel[ind+3]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+4);  hDWlnFakeSel[ind+4] = new TH1D(sb,sb,100,0.,2.5);
    sprintf(sb,"hDWlnFakeSel_%d",ind+5);  hDWlnFakeSel[ind+5] = new TH1D(sb,sb,100,0.,2.5);
    sprintf(sb,"hDWlnFakeSel_%d",ind+6);  hDWlnFakeSel[ind+6] = new TH1D(sb,sb,100,0.,2.5);
    sprintf(sb,"hDWlnFakeSel_%d",ind+7);  hDWlnFakeSel[ind+7] = new TH1D(sb,sb,100,0.,2.5);
  }

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    for(int i=0; i<8; i++){
      AddOutput(hDWlnFakeSel[ind+i]);
    }
  }

  sprintf(sb,"hDWlnFakeSel_50");  hDWlnFakeSel[50]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWlnFakeSel_51");  hDWlnFakeSel[51]  = new TH1D(sb,sb,600,0.0,6.0);
  sprintf(sb,"hDWlnFakeSel_52");  hDWlnFakeSel[52]  = new TH1D(sb,sb,600,0.0,6.0);
  sprintf(sb,"hDWlnFakeSel_53");  hDWlnFakeSel[53]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWlnFakeSel_54");  hDWlnFakeSel[54]  = new TH1D(sb,sb,10,-0.5,9.5);
  for(int j=50; j<=54; j++){
    AddOutput(hDWlnFakeSel[j]);
  }

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWlnFakeSel_%d",ind+80);  hDWlnFakeSel[ind+80]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+81);  hDWlnFakeSel[ind+81]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+82);  hDWlnFakeSel[ind+82]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+83);  hDWlnFakeSel[ind+83]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+90);  hDWlnFakeSel[ind+90]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+91);  hDWlnFakeSel[ind+91]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+92);  hDWlnFakeSel[ind+92]  = new TH1D(sb,sb,100,0.,100.);
    sprintf(sb,"hDWlnFakeSel_%d",ind+93);  hDWlnFakeSel[ind+93]  = new TH1D(sb,sb,100,0.,100.);
  }
  for(int j=0; j<4; j++){
    int ind = 100 * j;
    for(int i=80; i<84; i++){
      AddOutput(hDWlnFakeSel[ind+i]);
    }
    for(int i=90; i<94; i++){
      AddOutput(hDWlnFakeSel[ind+i]);
    }
  }

  for(int j=0; j<10; j++){
    sprintf(sb,"hDWlnFakeGenType_%d",j); hDWlnFakeGenType[j] = new TH1D(sb,sb,5,-0.5,4.5);
    AddOutput(hDWlnFakeGenType[j]);
  }  
}

//--------------------------------------------------------------------------------------------------
void WlnFakeSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WlnFakeSelMod::Terminate()
{
  // Run finishing code on the client computer
}
