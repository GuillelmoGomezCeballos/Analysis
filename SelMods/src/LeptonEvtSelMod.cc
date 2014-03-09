// $Id: LeptonEvtSelMod.cc,v 1.27 2012/04/18 14:59:41 ceballos Exp $

#include "Analysis/SelMods/interface/LeptonEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/MetCol.h"

using namespace mithep;
ClassImp(mithep::LeptonEvtSelMod)

//--------------------------------------------------------------------------------------------------
LeptonEvtSelMod::LeptonEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsData(kFALSE),
  fMetName(Names::gkCaloMetBrn),
  fTrackName(Names::gkTrackBrn),
  fTrack0Name("ConversionInOutTracks"),
  fTrack1Name("ConversionOutInTracks"),
  fTrack2Name("GsfTracks"),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fAllVertexName("random"),
  fVertexName(ModNames::gkGoodVertexesName),
  fConversionName(Names::gkMvfConversionBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fBeamSpot(0),
  fMuonTools(0),
  fMuonIDMVA(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void LeptonEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void LeptonEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Get Generator Level information for matching
  MCParticleOArr *GenLeptons       = GetObjThisEvt<MCParticleOArr>(fMCLeptonsName);
  ElectronOArr *CleanElectrons     = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons            = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");
  MCParticleOArr *GenTaus          = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCTausName);
  ParticleOArr *leptonsFakeable    = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  MetOArr *CleanMet	           = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet	           = CleanMet->At(0);
  Bool_t RemoveLeptons = kFALSE;
  if(GenLeptons && GenLeptons->GetEntries() == 2){
    CompositeParticle dilepton;
    dilepton.AddDaughter(GenLeptons->At(0));
    dilepton.AddDaughter(GenLeptons->At(1));
    if(dilepton.Mass() <= 12) RemoveLeptons = kTRUE;
  }
  if(RemoveLeptons == kTRUE) return;
  
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  if(fVertices->GetEntries() <= 0) return;
  int nVertices = (int)(fVertices->GetEntries()/4.); if(nVertices > 4) nVertices = 4;

  // Lepton efficiencies study
  for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
    MCParticle *gen = GenLeptons->At(j);
    if(gen->AbsEta() > 2.4) continue;
    if(gen->Pt() < 10.0) continue;

    if(gen->Is(MCParticle::kMu) == kTRUE){ // Muons
      hDLepSel[221]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[231]->Fill(gen->AbsEta(),NNLOWeight->GetVal());

      bool isRecoMuon = kFALSE;
      for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
	const Muon *mu = fMuons->At(i);
	if(!(mu->HasGlobalTrk() && mu->IsTrackerMuon())) continue;
	if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), mu->Phi(), mu->Eta()) < 0.1){
          isRecoMuon = kTRUE;
          break;
	}
      }
      if(isRecoMuon == kTRUE){
	hDLepSel[222]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
	hDLepSel[232]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	bool isGoodMuon = kFALSE;
	for (UInt_t i=0; i<CleanMuons->GetEntries(); i++) {
          const Muon *mu = CleanMuons->At(i);
          if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), mu->Phi(), mu->Eta()) < 0.1){
            isGoodMuon = kTRUE;
            break;
          }
	}
	if(isGoodMuon == kTRUE){
          hDLepSel[223]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
          hDLepSel[233]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	}
      }
    } // Muons loop
    else if(gen->Is(MCParticle::kEl) == kTRUE){ // Electrons
      hDLepSel[224]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[234]->Fill(gen->AbsEta(),NNLOWeight->GetVal());

      bool isRecoElectron = kFALSE;
      Bool_t passElId     = kFALSE;
      for (UInt_t i=0; i<fElectrons->GetEntries(); i++) {
	const Electron *el = fElectrons->At(i);
	if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), el->Phi(), el->Eta()) < 0.1){
          isRecoElectron = kTRUE;
          ElectronIDMod *electronID = new ElectronIDMod();
    	  electronID->SetIDType("VBTFWorkingPoint80Id");
    	  electronID->SetIsoType("TrackJuraSliding");
    	  electronID->Setup();
          passElId = electronID->PassIDCut(el, ElectronTools::kVBTFWorkingPoint80Id, fVertices->At(0));
	  delete electronID;
          break;
	}
      }
      if(isRecoElectron == kTRUE){
	hDLepSel[225]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
	hDLepSel[235]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	if(passElId){
          hDLepSel[226]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
          hDLepSel[236]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	  bool isGoodElectron = kFALSE;
	  for (UInt_t i=0; i<CleanElectrons->GetEntries(); i++) {
            Electron *cleanel = CleanElectrons->At(i);
            if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), cleanel->Phi(), cleanel->Eta()) < 0.1){
              isGoodElectron = kTRUE;
              break;
            }
	  }
	  if(isGoodElectron == kTRUE){
            hDLepSel[227]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
            hDLepSel[237]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	  }
        } // pass Id
      } // pass RECO
    } // Electrons loop
  } // Lepton loop

  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fTrackName);
  LoadBranch(fTrack0Name);
  LoadBranch(fTrack1Name);
  LoadBranch(fTrack2Name);
  LoadBranch(fPileupEnergyDensityName);
  LoadBranch(fBeamSpotName);
  //if(fMuons->GetEntries() + fElectrons->GetEntries() < 2) {return;}
  LoadBranch(fAllVertexName);
  LoadBranch(fPFCandidatesName);
  double isoAux;
  const PileupEnergyDensity *rho =  fPileupEnergyDensity->At(0);

  // Muon loop
  for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
    const Muon *mu = fMuons->At(i);

    if(mu->BestTrk() == 0) continue;
    if(mu->Pt() <= 10) continue;
    if(mu->AbsEta() >= 2.5) continue;

    bool isGenTau = kFALSE;
    if(GenTaus){
      for (UInt_t j=0; j<GenTaus->GetEntries(); j++) {
        MCParticle *genTau = GenTaus->At(j);
        if(MathUtils::DeltaR(genTau->Mom(), mu->Mom()) < 0.1){
          isGenTau = kTRUE;
	  break;
        }
      }
    }
    if(isGenTau == kTRUE) continue;

    bool isGenLepton = kFALSE;
    int nMu = -1;
    if(GenLeptons){
      for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
        MCParticle *gen = GenLeptons->At(j);
        if(!gen->Is(MCParticle::kMu)) continue;
        if(gen->Charge() != mu->Charge()) continue;
        if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), mu->Phi(), mu->Eta()) < 0.1){
          isGenLepton = kTRUE;
	  nMu = j;
	  break;
        }
      }
    }
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[200]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[210]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }
    bool passDataCut = kTRUE;
    if(fIsData == kTRUE){
      Bool_t isZ = kFALSE;
      for(UInt_t j=0; j<leptonsFakeable->GetEntries(); j++){
	if(mu->ObjType() != leptonsFakeable->At(j)->ObjType()) continue;
	CompositeParticle dilepton;
	dilepton.AddDaughter(mu);
	dilepton.AddDaughter(leptonsFakeable->At(j));
	if(TMath::Abs(dilepton.Mass() - 91.1876) < 15.0) {
          isZ = kTRUE;
          break;
	}
      }
      if(isZ == kTRUE){
        isGenLepton = kTRUE;
	passDataCut = kTRUE;
      }
      else {
        isGenLepton = kFALSE;
        double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), mu->Phi()));
        double mTW = TMath::Sqrt(2.0*mu->Pt()*caloMet->Pt()*
    				(1.0 - cos(deltaPhiMetLepton)));
	if(caloMet->Pt() >= 20.0 || mTW >= 15.0) passDataCut = kFALSE;
	else                                     passDataCut = kTRUE;
      }
    }
    if(passDataCut == kFALSE) continue;
    int type = 0;
    if(mu->IsGlobalMuon())     type = 1;
    if(mu->IsTrackerMuon())    type = type + 2;
    if(mu->IsStandaloneMuon()) type = type + 4;
    hDLepSel[0+(int)(isGenLepton==kTRUE)]->Fill((double) type,NNLOWeight->GetVal());
    if(mu->IsGlobalMuon()     && !mu->HasGlobalTrk())     cout << "IsGlobalMuon and not HasGlobalTrk!" << endl;
    if(mu->IsTrackerMuon()    && !mu->HasTrackerTrk())    cout << "IsTrackerMuon and not HasTrackerTrk!" << endl;
    if(mu->IsStandaloneMuon() && !mu->HasStandaloneTrk()) cout << "IsStandaloneMuon and not HasStandaloneTrk!" << endl;

    if(mu->IsGlobalMuon()){
      hDLepSel[2+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof(),19.999),NNLOWeight->GetVal());
      hDLepSel[4+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->NValidHits(),9.499),NNLOWeight->GetVal());
      hDLepSel[6+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->NSegments(),9.499),NNLOWeight->GetVal());    
    }
    if(mu->IsTrackerMuon()){
      hDLepSel[8+(int)(isGenLepton==kTRUE)]->Fill((double)mu->Quality().Quality(MuonQuality::TMLastStationTight),NNLOWeight->GetVal());
    }
    bool loosepass = (mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
	             (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
	            (mu->IsTrackerMuon() &&
	             mu->Quality().Quality(MuonQuality::TMLastStationTight));
    if(loosepass == kFALSE) continue;
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[201]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[211]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }
    
    hDLepSel[10+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->NTrkLayersHit(),19.499),NNLOWeight->GetVal());
    hDLepSel[12+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->IsPFMuon(),1.499),NNLOWeight->GetVal());
    hDLepSel[14+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->BestTrk()->NPixelHits(),9.499),NNLOWeight->GetVal());
    hDLepSel[16+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(mu->BestTrk()->PtErr()/mu->BestTrk()->Pt(),0.999),NNLOWeight->GetVal());
    hDLepSel[18+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min((double)mu->TrkKink(),19.999),NNLOWeight->GetVal());

    bool idpass = mu->BestTrk() != 0 &&
	          mu->NTrkLayersHit() > 5 &&
		  mu->IsPFMuon() == kTRUE &&
		  mu->BestTrk()->NPixelHits() > 0 &&
		  mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 1.0 &&
		  mu->TrkKink() < 20.0;
    if(!idpass) continue;
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[202]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[212]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }

    Double_t d0_real = TMath::Abs(mu->BestTrk()->D0Corrected(*fVertices->At(0)));
    Double_t dz_real = TMath::Abs(mu->BestTrk()->DzCorrected(*fVertices->At(0)));
    hDLepSel[20+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(d0_real,0.0999),NNLOWeight->GetVal());
    hDLepSel[22+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(dz_real,0.9999),NNLOWeight->GetVal());

    if(dz_real >= 0.100) continue;
    if(d0_real >= 0.100) continue;
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[203]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[213]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }

    hDLepSel[24+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->D0PV()),0.0999),NNLOWeight->GetVal());
    hDLepSel[26+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->D0PVBS()),0.0999),NNLOWeight->GetVal());
    hDLepSel[28+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->D0PVSignificance()),9.9999),NNLOWeight->GetVal());
    hDLepSel[30+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->Ip3dPVSignificance()),9.9999),NNLOWeight->GetVal());
    hDLepSel[32+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->D0PVBSSignificance()),9.9999),NNLOWeight->GetVal());
    hDLepSel[34+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(TMath::Abs(mu->Ip3dPVBSSignificance()),9.9999),NNLOWeight->GetVal());

    if(d0_real >= 0.020) continue;
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[204]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[214]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }

    ElectronOArr *tempElectrons = new  ElectronOArr;
    MuonOArr	 *tempMuons	= new  MuonOArr;
    double MVAValue = fMuonIDMVA->MVAValue(mu, fVertices->At(0), fMuonTools, fPFCandidates,
    				           fPileupEnergyDensity, MuonTools::kMuEAFall11MC, tempElectrons, tempMuons,kFALSE);
    delete tempElectrons;
    delete tempMuons;

    hDLepSel[36+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min(mu->IsoR03SumPt()/(mu->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    hDLepSel[38+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min((mu->IsoR03SumPt() + mu->IsoR03EmEt())/(mu->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    hDLepSel[40+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min((mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt())/(mu->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    hDLepSel[42+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min((mu->IsoR03SumPt() + TMath::Max(mu->IsoR03EmEt() + mu->IsoR03HadEt() - rho->Rho() * TMath::Pi() * 0.3 * 0.3, 0.0))/(mu->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    hDLepSel[44+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min((IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0), 0.1, 1.0, 0.3, 0.0))/(mu->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    hDLepSel[46+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min(MVAValue,0.999),0.000),NNLOWeight->GetVal());

    hDLepSel2D[0+(int)(isGenLepton==kTRUE)]->Fill(TMath::Min(mu->BestTrk()->PtErr()/mu->BestTrk()->Pt(),0.399),TMath::Min(mu->Pt(),399.999));
    if (mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() >= 0.3) continue;
    hDLepSel2D[2+(int)(isGenLepton==kTRUE)]->Fill(TMath::Max(TMath::Min(MVAValue,0.999),0.000),TMath::Min(mu->Pt(),399.999));
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[205]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[215]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }

    Double_t MVACut = -1.0;
    Double_t eta = mu->AbsEta();
    if     (mu->Pt() <  20 && eta <  1.479) MVACut = 0.86;
    else if(mu->Pt() <  20 && eta >= 1.479) MVACut = 0.82;
    else if(mu->Pt() >= 20 && eta <  1.479) MVACut = 0.82;
    else if(mu->Pt() >= 20 && eta >= 1.479) MVACut = 0.86;

    double dRMin = 999.;
    for (UInt_t nm=0; nm<fMuons->GetEntries(); nm++) {
      const Muon *mu2 = fMuons->At(nm);
      if(mu->BestTrk() == 0) continue;
      if(mu == mu2) continue;
      if(MathUtils::DeltaR(mu->Mom(), mu2->Mom()) < dRMin) dRMin = MathUtils::DeltaR(mu->Mom(), mu2->Mom());
    }

    hDLepSel2D[4+(int)(isGenLepton==kTRUE)]->Fill(dRMin,TMath::Min(mu->Pt(),399.999));
    if (MVAValue <= MVACut) continue;
    if(isGenLepton == kTRUE && fIsData == kFALSE) {
      hDLepSel[206]->Fill(TMath::Min(GenLeptons->At(nMu)->Pt(),199.999),NNLOWeight->GetVal());
      hDLepSel[216]->Fill(GenLeptons->At(nMu)->AbsEta(),NNLOWeight->GetVal());
    }
    hDLepSel2D[6+(int)(isGenLepton==kTRUE)]->Fill(dRMin,TMath::Min(mu->Pt(),399.999));

  } // Muon loop

  // Electron loop
  for (UInt_t i=0; i<fElectrons->GetEntries(); i++) {
    const Electron *el = fElectrons->At(i);
    if(el->SCluster() == 0) continue;
    if(el->Pt() <= 10) continue;
    if(el->AbsEta() >= 2.5) continue;

    bool isGenTau = kFALSE;
    for (UInt_t j=0; j<GenTaus->GetEntries(); j++) {
      MCParticle *genTau = GenTaus->At(j);
      if(MathUtils::DeltaR(genTau->Mom(), el->Mom()) < 0.1){
        isGenTau = kTRUE;
	break;
      }
    }
    if(isGenTau == kTRUE) continue;

    bool isGenLepton = kFALSE;
    for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
      MCParticle *gen = GenLeptons->At(j);
      if(!gen->Is(MCParticle::kEl)) continue;
      if(gen->Charge() != el->Charge()) continue;
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), el->Phi(), el->Eta()) < 0.1){
        isGenLepton = kTRUE;
	break;
      }
    }

    bool passDataCut = kTRUE;
    if(fIsData == kTRUE){
      Bool_t isZ = kFALSE;
      for(UInt_t j=0; j<leptonsFakeable->GetEntries(); j++){
	if(el->ObjType() != leptonsFakeable->At(j)->ObjType()) continue;
	CompositeParticle dilepton;
	dilepton.AddDaughter(el);
	dilepton.AddDaughter(leptonsFakeable->At(j));
	if(TMath::Abs(dilepton.Mass() - 91.1876) < 15.0) {
          isZ = kTRUE;
          break;
	}
      }
      if(isZ == kTRUE){
        isGenLepton = kTRUE;
	passDataCut = kTRUE;
      }
      else {
        isGenLepton = kFALSE;
        double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), el->Phi()));
        double mTW = TMath::Sqrt(2.0*el->Pt()*caloMet->Pt()*
    				(1.0 - cos(deltaPhiMetLepton)));
	if(caloMet->Pt() >= 20.0 || mTW >= 15.0) passDataCut = kFALSE;
	else                                     passDataCut = kTRUE;
      }
    }
    if(passDataCut == kFALSE) continue;

    Double_t dz_real = TMath::Abs(el->BestTrk()->DzCorrected(*fVertices->At(0)));
    if(dz_real >= 0.100) return;

    if(el->TrackIsolationDr03() < 5.0){
      if(isGenLepton == kTRUE) hDLepSel[150]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      else                     hDLepSel[151]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[152]->Fill(TMath::Max(TMath::Min(el->IDLikelihood(),0.99999),0.),NNLOWeight->GetVal());
      else                     hDLepSel[153]->Fill(TMath::Max(TMath::Min(el->IDLikelihood(),0.99999),0.),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE && el->PassTightID()) hDLepSel[154]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      else if(el->PassTightID())                    hDLepSel[155]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE && el->PassLooseID()) hDLepSel[156]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      else if(el->PassLooseID())                    hDLepSel[157]->Fill(TMath::Min(el->Pt(),99.999),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[158]->Fill(TMath::Min(el->TrackIsolationDr03(),19.999),NNLOWeight->GetVal());
      else                     hDLepSel[159]->Fill(TMath::Min(el->TrackIsolationDr03(),19.999),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[160]->Fill(TMath::Min(el->CaloIsolation(),19.999),NNLOWeight->GetVal());
      else                     hDLepSel[161]->Fill(TMath::Min(el->CaloIsolation(),19.999),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[162]->Fill(TMath::Max(TMath::Min(el->EcalRecHitIsoDr03(),19.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[163]->Fill(TMath::Max(TMath::Min(el->EcalRecHitIsoDr03(),19.999),0.000),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[164]->Fill(TMath::Max(TMath::Min(el->HcalTowerSumEtDr03(),19.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[165]->Fill(TMath::Max(TMath::Min(el->HcalTowerSumEtDr03(),19.999),0.000),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[166]->Fill(TMath::Max(TMath::Min(el->HcalTowerSumEtDr03()+el->EcalRecHitIsoDr03(),19.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[167]->Fill(TMath::Max(TMath::Min(el->HcalTowerSumEtDr03()+el->EcalRecHitIsoDr03(),19.999),0.000),NNLOWeight->GetVal());
    }

    if     (isGenLepton == kTRUE && el->IsEcalDriven()) hDLepSel[171]->Fill(TMath::Min(el->Pt(),199.99),NNLOWeight->GetVal());
    else if(isGenLepton == kTRUE                      ) hDLepSel[172]->Fill(TMath::Min(el->Pt(),199.99),NNLOWeight->GetVal());
    else if(                        el->IsEcalDriven()) hDLepSel[173]->Fill(TMath::Min(el->Pt(),199.99),NNLOWeight->GetVal());
    else                                                hDLepSel[174]->Fill(TMath::Min(el->Pt(),199.99),NNLOWeight->GetVal());
    if     (isGenLepton == kTRUE && el->IsEcalDriven()) hDLepSel[175]->Fill(el->AbsEta(),NNLOWeight->GetVal());
    else if(isGenLepton == kTRUE                      ) hDLepSel[176]->Fill(el->AbsEta(),NNLOWeight->GetVal());
    else if(                        el->IsEcalDriven()) hDLepSel[177]->Fill(el->AbsEta(),NNLOWeight->GetVal());
    else                                                hDLepSel[178]->Fill(el->AbsEta(),NNLOWeight->GetVal());

    // Begin Electron id study
    LoadBranch(fConversionName);
    ElectronIDMod *electronID = new ElectronIDMod();
    electronID->SetIDType("CustomTight");
    electronID->SetIsoType("TrackJuraSliding");
    electronID->Setup();
    Bool_t passElId     = electronID->PassIDCut(el, ElectronTools::kCustomIdTight, fVertices->At(0));
    electronID->SetCombIsoCut(0.10);
    Bool_t passIso      = electronID->PassIsolationCut(el, ElectronTools::kTrackJuraSliding, fTracks, fVertices->At(0), rho->Rho());
    if(el->Pt() < 20) electronID->SetCombIsoCut(0.10); // same cut as high pt leptons now
    Bool_t passIsoLowPt = electronID->PassIsolationCut(el, ElectronTools::kTrackJuraSliding, fTracks, fVertices->At(0), rho->Rho());
    Bool_t passD0       = ElectronTools::PassD0Cut(el, *&fVertices, kTRUE);

    Bool_t passConvVetoNoFit = TMath::Abs(el->ConvPartnerDCotTheta()) >= 0.02 ||  TMath::Abs(el->ConvPartnerDist()) >= 0.02;

    Bool_t passConvVeto_v0 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-6,   2.0,   kFALSE, kTRUE );	//(This was what you guys have already been testing, so we know it's too loose, but might as well leave it as a benchmark)
    Bool_t passConvVeto_v1 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-6,   2.0,   kFALSE, kFALSE);	//attempts matching to all conversions without regard to double-counting)
    Bool_t passConvVeto_v2 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-6,   2.0,   kTRUE,  kFALSE);	//also attempt matching through ckf track in addition to gsf track
    Bool_t passConvVeto_v3 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-6,  -1.0,   kTRUE,  kFALSE);	//relax lxy cut
    Bool_t passConvVeto_v4 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-6,  -999.9, kTRUE,  kFALSE);	//remove lxy cut
    Bool_t passConvVeto_v5 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 1,   1e-10, -999.9, kTRUE,  kFALSE);	//relax probability cut
    Bool_t passConvVeto_v6 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 2,   1e-10, -999.9, kTRUE,  kFALSE);	//relax hits before vtx cut
    Bool_t passConvVeto_v7 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 999, 1e-10, -999.9, kTRUE,  kFALSE);	//remove hits before vtx cut (this will probably be too tight)
    Bool_t passConvVeto_v8 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kFALSE, kTRUE );	// v0 with NextraHits=0
    Bool_t passConvVeto_v9 = ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kFALSE, kFALSE);	// v1 with NextraHits=0
    Bool_t passConvVeto_v10= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kTRUE,  kFALSE);	// v2 with NextraHits=0 <---- default
    Bool_t passConvVeto_v11= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,  -1.0,   kTRUE,  kFALSE);	// v3 with NextraHits=0
    Bool_t passConvVeto_v12= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,  -999.9, kTRUE,  kFALSE);	// v4 with NextraHits=0
    Bool_t passConvVeto_v13= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kTRUE,  kFALSE, 0.2); // v11 with pttrkmin>0.2
    Bool_t passConvVeto_v14= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kTRUE,  kFALSE, 0.5); // v11 with pttrkmin>0.5
    Bool_t passConvVeto_v15= ElectronTools::PassConversionFilter(el, fConversions, fBeamSpot->At(0), 0,   1e-6,   2.0,   kTRUE,  kFALSE, 1.0); // v11 with pttrkmin>1.0

/*
    Double_t eOverP	  = el->ESuperClusterOverP();
    Double_t fBrem	  = el->FBrem();
    //Double_t eSeedOverPin = el->ESeedClusterOverPIn(); 
    //Double_t hOverE	  = el->HadronicOverEm();
    //Double_t sigmaee	  = el->CoviEtaiEta();
    //Double_t deltaPhiIn   = TMath::Abs(el->DeltaPhiSuperClusterTrackAtVtx());
    //Double_t deltaEtaIn   = TMath::Abs(el->DeltaEtaSuperClusterTrackAtVtx());
    Int_t cat = 2;
    if ((el->IsEB() && fBrem<0.06) || (el->IsEE() && fBrem<0.1)) 
      cat=1;
    else if (eOverP < 1.2 && eOverP > 0.8) 
      cat=0;

    Int_t eb = 1;
    if (el->IsEB()) 
      eb = 0;

    cat = cat + 3*eb;
    int catNew = ElectronTools::Classify(el);

    if(fBrem >= 1) fBrem = 0.999;
*/
    double betaForCiC = IsolationTools::BetaE(fTracks, el, fVertices->At(0), 0.0, 0.2, 0.4, 0.02); 
    Int_t result0 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 0, betaForCiC);
    Int_t result1 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 1, betaForCiC);
    Int_t result2 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 2, betaForCiC);
    Int_t result3 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 3, betaForCiC);
    Int_t result4 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 4, betaForCiC);
    Int_t result5 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 5, betaForCiC);
    Int_t result6 = ElectronTools::PassTightId(el, *&fVertices, fConversions, 6, betaForCiC);
    Double_t scEt = el->Pt(); // el->SCluster()->Et(); 

    ElectronIDMod *electronIDVBTF80 = new ElectronIDMod();
    electronIDVBTF80->SetIDType("VBTFWorkingPoint80Id");
    electronIDVBTF80->SetIsoType("TrackJuraSliding");
    electronIDVBTF80->Setup();
    Bool_t passElIdVBTF80 = electronIDVBTF80->PassIDCut(el, ElectronTools::kVBTFWorkingPoint80Id, fVertices->At(0));
    delete electronIDVBTF80;
    ElectronIDMod *electronIDVBTFLowPt = new ElectronIDMod();
    electronIDVBTFLowPt->SetIDType("VBTFWorkingPointLowPtId");
    electronIDVBTFLowPt->SetIsoType("TrackJuraSliding");
    electronIDVBTFLowPt->Setup();
    Bool_t passElIdVBTFLowPt = electronIDVBTFLowPt->PassIDCut(el, ElectronTools::kVBTFWorkingPoint80Id, fVertices->At(0));
    delete electronIDVBTFLowPt;

    Bool_t  idLikcut = (el->Pt() >  15 && el->IDLikelihood() > 0.75) ||
                       (el->Pt() <= 15 && ElectronTools::PassCustomID(el, ElectronTools::kVBTFWorkingPoint80Id));

    int theHistoEle = 0;
    if     (isGenLepton == kTRUE && el->IsEB()) theHistoEle =   0;
    else if(isGenLepton == kTRUE	      ) theHistoEle =  50;
    else if(			    el->IsEB()) theHistoEle = 100;
    else					theHistoEle = 150;

    hDCutEleSel[theHistoEle+ 0]->Fill(result0,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 1]->Fill(result1,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 2]->Fill(result2,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 3]->Fill(result3,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 4]->Fill(result4,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 5]->Fill(result5,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 6]->Fill(result6,NNLOWeight->GetVal());
    hDCutEleSel[theHistoEle+ 7]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(passElIdVBTF80 && passIso && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0)
                                       hDCutEleSel[theHistoEle+ 8]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(idLikcut && passIso && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0)
                                       hDCutEleSel[theHistoEle+ 9]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(passElId && passIso && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0)
                                       hDCutEleSel[theHistoEle+10]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(passElIdVBTF80 && passIso && passD0 && passConvVetoNoFit)
                                       hDCutEleSel[theHistoEle+11]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(passElIdVBTFLowPt && passIsoLowPt && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0)
                                       hDCutEleSel[theHistoEle+12]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());

    if(passElIdVBTFLowPt && passIso && passD0 && el->CorrectedNExpectedHitsInner() <= 0){
      hDCutEleSel[theHistoEle+13]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());
      if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[ 0]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      else if(isGenLepton == kTRUE		) {hDEleConvSel[ 1]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      else if(  		      el->IsEB()) {hDEleConvSel[ 2]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      else					  {hDEleConvSel[ 3]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      if(passConvVetoNoFit){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[ 4]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[ 5]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[ 6]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[ 7]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v0){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[ 8]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[ 9]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[10]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[11]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v1){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[12]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[13]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[14]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[15]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v2){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[16]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[17]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[18]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[19]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v3){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[20]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[21]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[22]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[23]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v4){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[24]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[25]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[26]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[27]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v5){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[28]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[29]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[30]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[31]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v6){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[32]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[33]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[34]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[35]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v7){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[36]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[37]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[38]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[39]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v8){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[40]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[41]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[42]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[43]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v9){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[44]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[45]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[46]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[47]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v10){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[48]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[49]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[50]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[51]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v11){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[52]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[53]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[54]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[55]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v12){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[56]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[57]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[58]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[59]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v13){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[60]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[61]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[62]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[63]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v14){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[64]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[65]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[66]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[67]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
      if(passConvVeto_v15){
        if     (isGenLepton == kTRUE && el->IsEB()) {hDEleConvSel[68]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(isGenLepton == kTRUE		  ) {hDEleConvSel[69]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else if(			el->IsEB()) {hDEleConvSel[70]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
        else					    {hDEleConvSel[71]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());}
      }
    }

    if(passElIdVBTF80 && passIso           && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0){
      hDCutEleSel[theHistoEle+14]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());
      Double_t d0_real = TMath::Abs(el->BestTrk()->D0Corrected(*fVertices->At(0)));
      if(isGenLepton == kTRUE) {hDD0LepSel[ 0]->Fill(TMath::Min(d0_real,0.0999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[10]->Fill(TMath::Min(d0_real,0.0999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 1]->Fill(TMath::Min(TMath::Abs(el->D0PV()),0.0999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[11]->Fill(TMath::Min(TMath::Abs(el->D0PV()),0.0999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 2]->Fill(TMath::Min(TMath::Abs(el->D0PVCkf()),0.0999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[12]->Fill(TMath::Min(TMath::Abs(el->D0PVCkf()),0.0999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 3]->Fill(TMath::Min(TMath::Abs(el->D0PVSignificance()),9.9999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[13]->Fill(TMath::Min(TMath::Abs(el->D0PVSignificance()),9.9999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 4]->Fill(TMath::Min(TMath::Abs(el->Ip3dPVSignificance()),9.9999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[14]->Fill(TMath::Min(TMath::Abs(el->Ip3dPVSignificance()),9.9999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 5]->Fill(TMath::Min(TMath::Abs(el->D0PVCkfSignificance()),9.9999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[15]->Fill(TMath::Min(TMath::Abs(el->D0PVCkfSignificance()),9.9999),NNLOWeight->GetVal());}
      if(isGenLepton == kTRUE) {hDD0LepSel[ 6]->Fill(TMath::Min(TMath::Abs(el->Ip3dPVCkfSignificance()),9.9999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[16]->Fill(TMath::Min(TMath::Abs(el->Ip3dPVCkfSignificance()),9.9999),NNLOWeight->GetVal());}
      double d0_real_1rst = 0.0; if(fAllVertices->At(0)) d0_real_1rst = el->BestTrk()->D0Corrected(*fAllVertices->At(0));
      if(isGenLepton == kTRUE) {hDD0LepSel[ 7]->Fill(TMath::Min(TMath::Abs(d0_real_1rst),0.0999),NNLOWeight->GetVal());}
      else		       {hDD0LepSel[17]->Fill(TMath::Min(TMath::Abs(d0_real_1rst),0.0999),NNLOWeight->GetVal());}
    }

    if(                  passIso && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0){
      hDCutEleSel[theHistoEle+15]->Fill(TMath::Min(scEt,199.99),NNLOWeight->GetVal());
    }

    delete electronID;
    bool ElNoIsoCuts = passElIdVBTF80 && passD0 && passConvVetoNoFit && el->CorrectedNExpectedHitsInner() <= 0;
    if(!ElNoIsoCuts) continue;

    isoAux = el->TrackIsolationDr03();
    if(isGenLepton == kTRUE) hDLepSel[100]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[101]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    isoAux = el->TrackIsolationDr03() + el->EcalRecHitIsoDr03();
    if(el->SCluster()->AbsEta() < 1.479) isoAux = el->TrackIsolationDr03() + TMath::Max(el->EcalRecHitIsoDr03() - 1.0, 0.0);
    if(isGenLepton == kTRUE) hDLepSel[102]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[103]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    double isostd =                               el->TrackIsolationDr03() +            el->EcalRecHitIsoDr03()             + el->HcalTowerSumEtDr03();
    if(el->SCluster()->AbsEta() < 1.479) isostd = el->TrackIsolationDr03() + TMath::Max(el->EcalRecHitIsoDr03() - 1.0, 0.0) + el->HcalTowerSumEtDr03();
    if(isGenLepton == kTRUE) hDLepSel[104]->Fill(TMath::Max(TMath::Min(isostd/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[105]->Fill(TMath::Max(TMath::Min(isostd/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    double isorho =                               el->TrackIsolationDr03() + TMath::Max(           el->EcalRecHitIsoDr03()             + el->HcalTowerSumEtDr03() - rho->Rho() * TMath::Pi() * 0.3 * 0.3, 0.0);
    if(el->SCluster()->AbsEta() < 1.479) isorho = el->TrackIsolationDr03() + TMath::Max(TMath::Max(el->EcalRecHitIsoDr03() - 1.0, 0.0) + el->HcalTowerSumEtDr03() - rho->Rho() * TMath::Pi() * 0.3 * 0.3, 0.0);
    if(isGenLepton == kTRUE) hDLepSel[106]->Fill(TMath::Max(TMath::Min(isorho/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[107]->Fill(TMath::Max(TMath::Min(isorho/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    double isopf1 = IsolationTools::PFElectronIsolation(el, fPFCandidates, fVertices->At(0), 0.1, 1.0, 0.3, 0.0);
    if(isGenLepton == kTRUE) hDLepSel[108]->Fill(TMath::Max(TMath::Min(isopf1/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[109]->Fill(TMath::Max(TMath::Min(isopf1/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    double isopf2 = IsolationTools::PFElectronIsolation(el, fPFCandidates, fVertices->At(0), 0.1, 1.0, 0.4, 0.0);
    if(isGenLepton == kTRUE) hDLepSel[110]->Fill(TMath::Max(TMath::Min(isopf2/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else                     hDLepSel[111]->Fill(TMath::Max(TMath::Min(isopf2/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    int theHisto = -1;
    if     (isGenLepton == kTRUE  && el->AbsEta() <  1.479) theHisto =  0;
    else if(isGenLepton == kTRUE  && el->AbsEta() >= 1.479) theHisto =  5;
    else if(isGenLepton == kFALSE && el->AbsEta() <  1.479) theHisto = 10;
    else if(isGenLepton == kFALSE && el->AbsEta() >= 1.479) theHisto = 15;

    if     (el->Pt() < 20)
    hDIsoELepSel0[theHisto+nVertices+0]->Fill(TMath::Max(TMath::Min(isostd/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else if(el->Pt() < 35)
    hDIsoELepSel1[theHisto+nVertices+0]->Fill(TMath::Max(TMath::Min(isostd/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    if     (el->Pt() < 20)
    hDIsoELepSel0[theHisto+nVertices+20]->Fill(TMath::Max(TMath::Min(isorho/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else if(el->Pt() < 35)
    hDIsoELepSel1[theHisto+nVertices+20]->Fill(TMath::Max(TMath::Min(isorho/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    if     (el->Pt() < 20)
    hDIsoELepSel0[theHisto+nVertices+40]->Fill(TMath::Max(TMath::Min(isopf1/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else if(el->Pt() < 35)
    hDIsoELepSel1[theHisto+nVertices+40]->Fill(TMath::Max(TMath::Min(isopf1/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    if     (el->Pt() < 20)
    hDIsoELepSel0[theHisto+nVertices+60]->Fill(TMath::Max(TMath::Min(isopf2/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else if(el->Pt() < 35)
    hDIsoELepSel1[theHisto+nVertices+60]->Fill(TMath::Max(TMath::Min(isopf2/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());

    isoAux = isostd;
    if(el->Pt() < 20){
      if(isGenLepton == kTRUE) hDLepSel[112]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[113]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      if(isGenLepton == kTRUE) hDLepSel[114]->Fill(TMath::Max(TMath::Min(isoAux/(20.00000-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[115]->Fill(TMath::Max(TMath::Min(isoAux/(20.00000-0.0),0.999),0.000),NNLOWeight->GetVal());
    } else if(el->Pt() < 25){
      if(isGenLepton == kTRUE) hDLepSel[116]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[117]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    } else if(el->Pt() < 30){
      if(isGenLepton == kTRUE) hDLepSel[118]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[119]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    } else if(el->Pt() < 40){
      if(isGenLepton == kTRUE) hDLepSel[120]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[121]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    } else if(el->Pt() < 60){
      if(isGenLepton == kTRUE) hDLepSel[122]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[123]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    } else {
      if(isGenLepton == kTRUE) hDLepSel[124]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
      else                     hDLepSel[125]->Fill(TMath::Max(TMath::Min(isoAux/(el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    }
    if(isGenLepton == kTRUE) hDLepSel[126]->Fill(TMath::Max(TMath::Min(isoAux / (el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    else		     hDLepSel[127]->Fill(TMath::Max(TMath::Min(isoAux / (el->Pt()-0.0),0.999),0.000),NNLOWeight->GetVal());
    if(isGenLepton == kTRUE) hDLepSel[128]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());
    else		     hDLepSel[129]->Fill(TMath::Min(el->Pt(),99.99),NNLOWeight->GetVal());

    if(isoAux > 0 && isoAux < 20.0){
      hDLepSel2D[8+(int)(isGenLepton==kTRUE)]->Fill(isoAux,el->Pt());
    }

    if(isoAux / (el->Pt()-0.0) < 1.0){
      int theHisto = 130;
      if(isGenLepton == kFALSE) theHisto = 131;
      hDLepSel[theHisto]->Fill(0.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.05) hDLepSel[theHisto]->Fill(1.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.10) hDLepSel[theHisto]->Fill(2.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.15) hDLepSel[theHisto]->Fill(3.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.20) hDLepSel[theHisto]->Fill(4.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.25) hDLepSel[theHisto]->Fill(5.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.30) hDLepSel[theHisto]->Fill(6.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.35) hDLepSel[theHisto]->Fill(7.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.40) hDLepSel[theHisto]->Fill(8.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.45) hDLepSel[theHisto]->Fill(9.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()-10.0)*0.50) hDLepSel[theHisto]->Fill(10.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.05) hDLepSel[theHisto]->Fill(11.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.10) hDLepSel[theHisto]->Fill(12.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.15) hDLepSel[theHisto]->Fill(13.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.20) hDLepSel[theHisto]->Fill(14.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.25) hDLepSel[theHisto]->Fill(15.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.30) hDLepSel[theHisto]->Fill(16.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.35) hDLepSel[theHisto]->Fill(17.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.40) hDLepSel[theHisto]->Fill(18.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.45) hDLepSel[theHisto]->Fill(19.0,NNLOWeight->GetVal());
      if(isoAux < (el->Pt()- 0.0)*0.50) hDLepSel[theHisto]->Fill(20.0,NNLOWeight->GetVal());
      if(isoAux <  1) hDLepSel[theHisto]->Fill(21.0,NNLOWeight->GetVal());
      if(isoAux <  2) hDLepSel[theHisto]->Fill(22.0,NNLOWeight->GetVal());
      if(isoAux <  3) hDLepSel[theHisto]->Fill(23.0,NNLOWeight->GetVal());
      if(isoAux <  4) hDLepSel[theHisto]->Fill(24.0,NNLOWeight->GetVal());
      if(isoAux <  5) hDLepSel[theHisto]->Fill(25.0,NNLOWeight->GetVal());
      if(isoAux <  6) hDLepSel[theHisto]->Fill(26.0,NNLOWeight->GetVal());
      if(isoAux <  7) hDLepSel[theHisto]->Fill(27.0,NNLOWeight->GetVal());
      if(isoAux <  8) hDLepSel[theHisto]->Fill(28.0,NNLOWeight->GetVal());
      if(isoAux <  9) hDLepSel[theHisto]->Fill(29.0,NNLOWeight->GetVal());
      if(isoAux < 10) hDLepSel[theHisto]->Fill(30.0,NNLOWeight->GetVal());
    } // iso vs. Pt study
  } // Electron loop

}
//--------------------------------------------------------------------------------------------------
void LeptonEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fAllVertexName,          fAllVertices);
  ReqBranch(fTrackName,              fTracks);
  ReqBranch(fTrack0Name,             fTracks0);
  ReqBranch(fTrack1Name,             fTracks1);
  ReqBranch(fTrack2Name,             fTracks2);
  ReqBranch(fMuonName,               fMuons);
  ReqBranch(fElectronName,           fElectrons);
  ReqBranch(fConversionName,         fConversions);
  ReqBranch(fPFCandidatesName,       fPFCandidates);
  ReqBranch(fPileupEnergyDensityName,fPileupEnergyDensity);
  ReqBranch(fBeamSpotName,           fBeamSpot);

  fMuonTools = new MuonTools();
  fMuonIDMVA = new MuonIDMVA();
  std::vector<std::string> muonidiso_weightfiles;
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_lowpt.weights.xml"))));
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_highpt.weights.xml"))));
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_lowpt.weights.xml"))));
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_highpt.weights.xml"))));
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_tracker.weights.xml"))));
  muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_global.weights.xml"))));
  fMuonIDMVA->Initialize("MuonIso_BDTG_IsoRings",
  		    MuonIDMVA::kIsoRingsV0,
  		    kTRUE,
  		    muonidiso_weightfiles,RhoUtilities::CMS_RHO_RHOKT6PFJETS);

  char sb[200];

  // Isolation
  for(int i=0; i<80; i++){
    sprintf(sb,"hDIsoELepSel0_%d", i);  hDIsoELepSel0[i]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDIsoELepSel1_%d", i);  hDIsoELepSel1[i]  = new TH1D(sb,sb,100,0.0,1.0);
    AddOutput(hDIsoELepSel0[i]);
    AddOutput(hDIsoELepSel1[i]);
  }

  sprintf(sb,"hDLepSel_%d", 0);  hDLepSel[ 0]  = new TH1D(sb,sb,8,-0.5,7.5);
  sprintf(sb,"hDLepSel_%d", 1);  hDLepSel[ 1]  = new TH1D(sb,sb,8,-0.5,7.5);
  sprintf(sb,"hDLepSel_%d", 2);  hDLepSel[ 2]  = new TH1D(sb,sb,100,0.0,20.0);
  sprintf(sb,"hDLepSel_%d", 3);  hDLepSel[ 3]  = new TH1D(sb,sb,100,0.0,20.0);
  sprintf(sb,"hDLepSel_%d", 4);  hDLepSel[ 4]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d", 5);  hDLepSel[ 5]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d", 6);  hDLepSel[ 6]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d", 7);  hDLepSel[ 7]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d", 8);  hDLepSel[ 8]  = new TH1D(sb,sb,2,-0.5,1.5);
  sprintf(sb,"hDLepSel_%d", 9);  hDLepSel[ 9]  = new TH1D(sb,sb,2,-0.5,1.5);
  sprintf(sb,"hDLepSel_%d",10);  hDLepSel[10]  = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDLepSel_%d",11);  hDLepSel[11]  = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDLepSel_%d",12);  hDLepSel[12]  = new TH1D(sb,sb,2,-0.5,1.5);
  sprintf(sb,"hDLepSel_%d",13);  hDLepSel[13]  = new TH1D(sb,sb,2,-0.5,1.5);
  sprintf(sb,"hDLepSel_%d",14);  hDLepSel[14]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d",15);  hDLepSel[15]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDLepSel_%d",16);  hDLepSel[16]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",17);  hDLepSel[17]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",18);  hDLepSel[18]  = new TH1D(sb,sb,40,0.0,20.0);
  sprintf(sb,"hDLepSel_%d",19);  hDLepSel[19]  = new TH1D(sb,sb,40,0.0,20.0);
  sprintf(sb,"hDLepSel_%d",20);  hDLepSel[20]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",21);  hDLepSel[21]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",22);  hDLepSel[22]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",23);  hDLepSel[23]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",24);  hDLepSel[24]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",25);  hDLepSel[25]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",26);  hDLepSel[26]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",27);  hDLepSel[27]  = new TH1D(sb,sb,100,0.0,0.1);
  sprintf(sb,"hDLepSel_%d",28);  hDLepSel[28]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",29);  hDLepSel[29]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",30);  hDLepSel[30]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",31);  hDLepSel[31]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",32);  hDLepSel[32]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",33);  hDLepSel[33]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",34);  hDLepSel[34]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",35);  hDLepSel[35]  = new TH1D(sb,sb,100,0.0,10.0);
  sprintf(sb,"hDLepSel_%d",36);  hDLepSel[36]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",37);  hDLepSel[37]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",38);  hDLepSel[38]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",39);  hDLepSel[39]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",40);  hDLepSel[40]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",41);  hDLepSel[41]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",42);  hDLepSel[42]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",43);  hDLepSel[43]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",44);  hDLepSel[44]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",45);  hDLepSel[45]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",46);  hDLepSel[46]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",47);  hDLepSel[47]  = new TH1D(sb,sb,100,0.0,1.0);

  for(int i=0; i<48; i++){
    AddOutput(hDLepSel[i]);
  }

  // Electrons
  sprintf(sb,"hDLepSel_%d", 0+100);  hDLepSel[ 0+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 1+100);  hDLepSel[ 1+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 2+100);  hDLepSel[ 2+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 3+100);  hDLepSel[ 3+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 4+100);  hDLepSel[ 4+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 5+100);  hDLepSel[ 5+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 6+100);  hDLepSel[ 6+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 7+100);  hDLepSel[ 7+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 8+100);  hDLepSel[ 8+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d", 9+100);  hDLepSel[ 9+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",10+100);  hDLepSel[10+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",11+100);  hDLepSel[11+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",12+100);  hDLepSel[12+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",13+100);  hDLepSel[13+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",14+100);  hDLepSel[14+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",15+100);  hDLepSel[15+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",16+100);  hDLepSel[16+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",17+100);  hDLepSel[17+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",18+100);  hDLepSel[18+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",19+100);  hDLepSel[19+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",20+100);  hDLepSel[20+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",21+100);  hDLepSel[21+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",22+100);  hDLepSel[22+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",23+100);  hDLepSel[23+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",24+100);  hDLepSel[24+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",25+100);  hDLepSel[25+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",26+100);  hDLepSel[26+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",27+100);  hDLepSel[27+100]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDLepSel_%d",28+100);  hDLepSel[28+100]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",29+100);  hDLepSel[29+100]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",30+100);  hDLepSel[30+100]  = new TH1D(sb,sb,41,-0.5,40.5);
  sprintf(sb,"hDLepSel_%d",31+100);  hDLepSel[31+100]  = new TH1D(sb,sb,41,-0.5,40.5);

  for(int i=0; i<32; i++){
    AddOutput(hDLepSel[i+100]);
  }

  sprintf(sb,"hDLepSel_%d",150);  hDLepSel[150]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",151);  hDLepSel[151]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",152);  hDLepSel[152]  = new TH1D(sb,sb,200,0.,1.0);
  sprintf(sb,"hDLepSel_%d",153);  hDLepSel[153]  = new TH1D(sb,sb,200,0.,1.0);
  sprintf(sb,"hDLepSel_%d",154);  hDLepSel[154]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",155);  hDLepSel[155]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",156);  hDLepSel[156]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",157);  hDLepSel[157]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDLepSel_%d",158);  hDLepSel[158]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",159);  hDLepSel[159]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",160);  hDLepSel[160]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",161);  hDLepSel[161]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",162);  hDLepSel[162]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",163);  hDLepSel[163]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",164);  hDLepSel[164]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",165);  hDLepSel[165]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",166);  hDLepSel[166]  = new TH1D(sb,sb,200,0.,20.0);
  sprintf(sb,"hDLepSel_%d",167);  hDLepSel[167]  = new TH1D(sb,sb,200,0.,20.0);

  for(int i=150; i<168; i++){
    AddOutput(hDLepSel[i]);
  }

  // Ecal/tracker driven
  sprintf(sb,"hDLepSel_%d",171);  hDLepSel[171]  = new TH1D(sb,sb,100,0.,100.0);
  sprintf(sb,"hDLepSel_%d",172);  hDLepSel[172]  = new TH1D(sb,sb,100,0.,100.0);
  sprintf(sb,"hDLepSel_%d",173);  hDLepSel[173]  = new TH1D(sb,sb,100,0.,100.0);
  sprintf(sb,"hDLepSel_%d",174);  hDLepSel[174]  = new TH1D(sb,sb,100,0.,100.0);
  sprintf(sb,"hDLepSel_%d",175);  hDLepSel[175]  = new TH1D(sb,sb,100,0.0,2.5);
  sprintf(sb,"hDLepSel_%d",176);  hDLepSel[176]  = new TH1D(sb,sb,100,0.0,2.5);
  sprintf(sb,"hDLepSel_%d",177);  hDLepSel[177]  = new TH1D(sb,sb,100,0.0,2.5);
  sprintf(sb,"hDLepSel_%d",178);  hDLepSel[178]  = new TH1D(sb,sb,100,0.0,2.5);

  for(int i=171; i<=178; i++){
    AddOutput(hDLepSel[i]);
  }

  // Lepton Id
  sprintf(sb,"hDLepSel_%d",200);  hDLepSel[200]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",201);  hDLepSel[201]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",202);  hDLepSel[202]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",203);  hDLepSel[203]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",204);  hDLepSel[204]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",205);  hDLepSel[205]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",206);  hDLepSel[206]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",210);  hDLepSel[210]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",211);  hDLepSel[211]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",212);  hDLepSel[212]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",213);  hDLepSel[213]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",214);  hDLepSel[214]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",215);  hDLepSel[215]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",216);  hDLepSel[216]  = new TH1D(sb,sb,96,0,2.4);

  for(int i=0; i<=6; i++){
    AddOutput(hDLepSel[i+200]);
    AddOutput(hDLepSel[i+210]);
  }

  sprintf(sb,"hDLepSel_%d",221);  hDLepSel[221]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",222);  hDLepSel[222]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",223);  hDLepSel[223]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",224);  hDLepSel[224]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",225);  hDLepSel[225]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",226);  hDLepSel[226]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",227);  hDLepSel[227]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDLepSel_%d",231);  hDLepSel[231]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",232);  hDLepSel[232]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",233);  hDLepSel[233]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",234);  hDLepSel[234]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",235);  hDLepSel[235]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",236);  hDLepSel[236]  = new TH1D(sb,sb,96,0,2.4);
  sprintf(sb,"hDLepSel_%d",237);  hDLepSel[237]  = new TH1D(sb,sb,96,0,2.4);

  for(int i=1; i<=7; i++){
    AddOutput(hDLepSel[i+220]);
    AddOutput(hDLepSel[i+230]);
  }

  // 2D
  sprintf(sb,"hDLepSel2D_%d",0);   hDLepSel2D[0] = new TH2D(sb,sb,40,0.0,0.4,40,0,400);
  sprintf(sb,"hDLepSel2D_%d",1);   hDLepSel2D[1] = new TH2D(sb,sb,40,0.0,0.4,40,0,400); 
  sprintf(sb,"hDLepSel2D_%d",2);   hDLepSel2D[2] = new TH2D(sb,sb,50,0.0,1.0,40,0,400);
  sprintf(sb,"hDLepSel2D_%d",3);   hDLepSel2D[3] = new TH2D(sb,sb,50,0.0,1.0,40,0,400); 
  sprintf(sb,"hDLepSel2D_%d",4);   hDLepSel2D[4] = new TH2D(sb,sb,80,0.0,4.0,40,0,400);
  sprintf(sb,"hDLepSel2D_%d",5);   hDLepSel2D[5] = new TH2D(sb,sb,80,0.0,4.0,40,0,400); 
  sprintf(sb,"hDLepSel2D_%d",6);   hDLepSel2D[6] = new TH2D(sb,sb,80,0.0,4.0,40,0,400);
  sprintf(sb,"hDLepSel2D_%d",7);   hDLepSel2D[7] = new TH2D(sb,sb,80,0.0,4.0,40,0,400); 
  sprintf(sb,"hDLepSel2D_%d",8);   hDLepSel2D[8] = new TH2D(sb,sb,50,0.0,20.0,50,0,100); 
  sprintf(sb,"hDLepSel2D_%d",9);   hDLepSel2D[9] = new TH2D(sb,sb,50,0.0,20.0,50,0,100); 

  for(int i=0; i<10; i++){
    AddOutput(hDLepSel2D[i]);
  }

  // Electron id Study
  for(int i=0; i<=6; i++){
    for(int j=0; j<4; j++){
      sprintf(sb,"hDCutEleSel_%d",i+50*j); hDCutEleSel[i+50*j] = new TH1D(sb,sb,16,-0.5,15.5);
      AddOutput(hDCutEleSel[i+50*j]);
    }
  }

  for(int i=7; i<16; i++){
    for(int j=0; j<4; j++){
      sprintf(sb,"hDCutEleSel_%d",i+50*j); hDCutEleSel[i+50*j] = new TH1D(sb,sb,200,0.0,200.0);
      AddOutput(hDCutEleSel[i+50*j]);
    }
  }

  for(int j=0; j<72; j++){
    sprintf(sb,"hDEleConvSel_%d", j); hDEleConvSel[j] = new TH1D(sb,sb,100,0.0,100);
    AddOutput(hDEleConvSel[j]);
  }

  // d0 studies
  for(int j=0; j<2; j++){
    sprintf(sb,"hDD0LepSel_%d", 0+10*j); hDD0LepSel[ 0+10*j] = new TH1D(sb,sb,100,0.0,0.1);
    sprintf(sb,"hDD0LepSel_%d", 1+10*j); hDD0LepSel[ 1+10*j] = new TH1D(sb,sb,100,0.0,0.1);
    sprintf(sb,"hDD0LepSel_%d", 2+10*j); hDD0LepSel[ 2+10*j] = new TH1D(sb,sb,100,0.0,0.1);
    sprintf(sb,"hDD0LepSel_%d", 3+10*j); hDD0LepSel[ 3+10*j] = new TH1D(sb,sb,100,0.,10);
    sprintf(sb,"hDD0LepSel_%d", 4+10*j); hDD0LepSel[ 4+10*j] = new TH1D(sb,sb,100,0.,10);
    sprintf(sb,"hDD0LepSel_%d", 5+10*j); hDD0LepSel[ 5+10*j] = new TH1D(sb,sb,100,0.,10);
    sprintf(sb,"hDD0LepSel_%d", 6+10*j); hDD0LepSel[ 6+10*j] = new TH1D(sb,sb,100,0.,10);
    sprintf(sb,"hDD0LepSel_%d", 7+10*j); hDD0LepSel[ 7+10*j] = new TH1D(sb,sb,80,0.,40);
  }

  for(int j=0; j<2; j++){
    AddOutput(hDD0LepSel[ 0+10*j]);
    AddOutput(hDD0LepSel[ 1+10*j]);
    AddOutput(hDD0LepSel[ 2+10*j]);
    AddOutput(hDD0LepSel[ 3+10*j]);
    AddOutput(hDD0LepSel[ 4+10*j]);
    AddOutput(hDD0LepSel[ 5+10*j]);
    AddOutput(hDD0LepSel[ 6+10*j]);
    AddOutput(hDD0LepSel[ 7+10*j]);
  }
}

//--------------------------------------------------------------------------------------------------
void LeptonEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
  delete fMuonTools;
  delete fMuonIDMVA;
}

//--------------------------------------------------------------------------------------------------
void LeptonEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
