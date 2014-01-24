// $Id: WWEvtSelMod.cc,v 1.50 2013/01/09 10:20:13 ceballos Exp $

#include "Analysis/SelMods/interface/WWEvtSelMod.h"
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/TrackJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"

using namespace mithep;
ClassImp(mithep::WWEvtSelMod)

//--------------------------------------------------------------------------------------------------
WWEvtSelMod::WWEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsFastSim(kFALSE),
  fIsData(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("randomMet"),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fTrackName(Names::gkTrackBrn),
  fAllVertexName("random"),
  fVertexName(ModNames::gkGoodVertexesName),
  fCleanJetsNoPtCutName("random"),
  fMuons(0),
  fElectrons(0),
  fPFMetName("PFMet"),
  fPFMetStd(0),
  fTCMetName("TCMet"),
  fTCMetStd(0),
  fPileupInfoName(Names::gkPileupInfoBrn),
  fPileupInfo(0),
  fCaloJetName0("AKt5Jets"),
  fCaloJetName1("AKt7Jets"),
  fCaloJetName2("Kt4Jets"),
  fCaloJetName3("Kt6Jets"),
  fTrackJetName0("AKt5TrackJets"),
  fPFJetName0("AKt5PFJets"),
  fPFJetName1("AKt7PFJets"),
  fPFJetName2("Kt4PFJets"),
  fPFJetName3("Kt6PFJets"),
  fCaloJet0(0),
  fCaloJet1(0),
  fCaloJet2(0),
  fCaloJet3(0),
  fTrackJet0(0),
  fPFJet0(0),
  fPFJet1(0),
  fPFJet2(0),
  fPFJet3(0),
  fJetScaleSyst(0.0),
  fUsePDFs(kFALSE),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fBeamSpot(0),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fConversions(0),
  fMCPartName(Names::gkMCPartBrn),
  fParticles(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fIntRadius(0.0),
  fIsOldSelection(kFALSE),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");
  FArrDouble *PDFArr = 0;
  if(fUsePDFs == kTRUE){
    PDFArr =  GetObjThisEvt<FArrDouble >("PDFWeights");
  //cout << PDFArr->GetEntries() << " " << PDFArr->At(0) << endl;
  }

  //Get Generator Level information for matching
  //ObjArray<MCParticle> *GenNeutrinos = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCNeutrinosName.Data()));

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons  = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons         = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJetsNoPtCut     = GetObjThisEvt<JetOArr>(fCleanJetsNoPtCutName);
  ParticleOArr *leptons         = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet             = GetObjThisEvt<MetOArr>(fMetName);
  const Met *cleanMet           = CleanMet->At(0);
  MCParticleOArr *GenRadPhotons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCRadPhotonsName);
  MCParticleOArr *GenISRPhotons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCISRPhotonsName);
  MCParticleOArr *GenLeptons    = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
  MCParticleOArr *GenAllLeptons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCAllLeptonsName);
  MCParticleOArr *GenBosons     = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCBosonsName);

  LoadBranch(fPFCandidatesName);
  LoadBranch(fBeamSpotName);
  LoadBranch(fConversionBranchName);

  MetOArr *GenMet              = GetObjThisEvt<MetOArr>(ModNames::gkMCMETName);

  LoadBranch(fPFMetName);
  const PFMet *PFMetStd        = fPFMetStd->At(0);

  LoadBranch(fTCMetName);
  const Met *TCMetStd          = fTCMetStd->At(0);

  //PFCandidateOArr *pfPileUp   = GetObjThisEvt<PFCandidateOArr>("PFPileUp");
  //PFCandidateOArr *pfNoPileUp = GetObjThisEvt<PFCandidateOArr>("PFNoPileUp");
  //cout << "PF: " << pfPileUp->GetEntries() << " " << pfNoPileUp->GetEntries() << " " << fPFCandidates->GetEntries() << endl;

  LoadBranch(fEvtHdrName);
  LoadBranch(fAllVertexName);
  LoadBranch(fCaloJetName0);
  LoadBranch(fCaloJetName1);
  LoadBranch(fCaloJetName2);
  LoadBranch(fCaloJetName3);
  LoadBranch(fTrackJetName0);
  LoadBranch(fPFJetName0);
  LoadBranch(fPFJetName1);
  LoadBranch(fPFJetName2);
  LoadBranch(fPFJetName3);
  LoadBranch(fPileupEnergyDensityName);
  if(fIsData == kFALSE){
    if (GenLeptons->GetEntries() == 2 &&
        GenLeptons->At(0)->AbsEta() < 2.5 && GenLeptons->At(0)->Pt() > 10.0 &&
        GenLeptons->At(1)->AbsEta() < 2.5 && GenLeptons->At(1)->Pt() > 10.0) {
      CompositeParticle dilepton;
      dilepton.AddDaughter(GenLeptons->At(0));
      dilepton.AddDaughter(GenLeptons->At(1));
      hDwwPresel[25]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    }
    //const PileupEnergyDensity *rho =  fPileupEnergyDensity->At(0);
    //printf(" %f %d %d %d\n",rho->Rho(),fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
    LoadEventObject(fMCPartName, fParticles);
    LoadBranch(fPileupInfoName);
    for(unsigned int npu=0; npu<fPileupInfo->GetEntries(); npu++){
      if     (fPileupInfo->At(npu)->GetBunchCrossing() == -1) hDwwPresel[22]->Fill(TMath::Min((double)fPileupInfo->At(npu)->GetPU_NumInteractions(),29.4999),NNLOWeight->GetVal());
      else if(fPileupInfo->At(npu)->GetBunchCrossing() ==  0) hDwwPresel[23]->Fill(TMath::Min((double)fPileupInfo->At(npu)->GetPU_NumInteractions(),29.4999),NNLOWeight->GetVal());
      else if(fPileupInfo->At(npu)->GetBunchCrossing() == +1) hDwwPresel[24]->Fill(TMath::Min((double)fPileupInfo->At(npu)->GetPU_NumInteractions(),29.4999),NNLOWeight->GetVal());
      //else cout << "Warning, weird PU! " << fPileupInfo->At(npu)->GetBunchCrossing() << " "<< fPileupInfo->At(npu)->GetPU_NumInteractions() <<endl;
    }
  }
  
  // Z analysis at GEN/RECO level
  int nZBoson = -1;
  for(UInt_t i=0; i<GenBosons->GetEntries(); i++) {
    if(GenBosons->At(i)->Is(MCParticle::kZ)) {
      nZBoson = i;
      break;
    }
  }
  if(nZBoson != -1) {
    hDGenZwwSel[0]->Fill(TMath::Min(GenBosons->At(nZBoson)->Mass(),199.999),NNLOWeight->GetVal());
    int passLeptons = 0;
    for(UInt_t i=0; i<GenAllLeptons->GetEntries(); i++) {
      if(GenAllLeptons->At(i)->Pt() > 10.0 && GenAllLeptons->At(i)->AbsEta() < 2.5) passLeptons++;
    }
    if(passLeptons >= 2) hDGenZwwSel[5]->Fill(TMath::Min(GenBosons->At(nZBoson)->Mass(),199.999),NNLOWeight->GetVal());
    int zType = -1;
    if     (GenBosons->At(nZBoson)->Mass() > 10 && GenBosons->At(nZBoson)->Mass() <= 20) zType = 0;
    else if(GenBosons->At(nZBoson)->Mass() > 20)                                         zType = 1;
    if(zType >= 0){
      hDGenZwwSel[5*zType+ 1]->Fill(TMath::Min((double)passLeptons,4.499),NNLOWeight->GetVal());
      if(passLeptons >= 2){
        hDGenZwwSel[5*zType+ 2]->Fill(TMath::Min(PFMetStd->Pt(),99.999),NNLOWeight->GetVal());
	if(PFMetStd->Pt() > 15) {
	  double jetPt = 0.0;
	  if(CleanJetsNoPtCut->GetEntries() > 0) jetPt = CleanJetsNoPtCut->At(0)->Pt();
	  hDGenZwwSel[5*zType+ 3]->Fill(TMath::Min(jetPt,99.999),NNLOWeight->GetVal());
	  if(jetPt > 20) {
	    hDGenZwwSel[5*zType+ 4]->Fill(TMath::Min(CleanJetsNoPtCut->At(0)->AbsEta(),4.999),NNLOWeight->GetVal());
	  }
	}
      }
    }
  }
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  if(fVertices->GetEntries() == 0 || fVertices->GetEntries() == 0) return;
  if(fAllVertices->At(0) != fVertices->At(0)) return;

  double zAverage = 0.0;
  double zDiffMax = 0.0;
  std::vector<double> leptonsDz;

  LoadBranch(fMuonName);
  ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
  for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
    const Muon *mu = fMuons->At(i);
    if(mu->BestTrk() == 0) continue;
    if(!MuonTools::PassSoftMuonCut(mu, fVertices, 0.2)) continue;
    
    bool isCleanMuon = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      if(fMuons->At(i) == CleanMuons->At(j) &&
  	 CleanMuons->At(j)->Pt() > 10) isCleanMuon = kTRUE;
    }
    if(isCleanMuon == kFALSE) DirtyMuons->Add(mu);
  }
  // Make lepton vector from muons and electrons
  for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
    double pDz = CleanMuons->At(j)->BestTrk()->DzCorrected(*fVertices->At(0));
    leptonsDz.push_back(pDz);
  }

  LoadBranch(fElectronName);
  for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
    double pDz = CleanElectrons->At(j)->GsfTrk()->DzCorrected(*fVertices->At(0));
    leptonsDz.push_back(pDz);
  }

  for(UInt_t t=0; t<leptonsDz.size(); t++) {
    zAverage = zAverage + leptonsDz[t];
    for(UInt_t i=t+1; i<leptonsDz.size(); i++) {
      if(TMath::Abs(leptonsDz[t]-leptonsDz[i]) > zDiffMax) zDiffMax = TMath::Abs(leptonsDz[t]-leptonsDz[i]);
    }
  }
  leptonsDz.clear();
  // Sort and count the number of central Jets for vetoing
  int nCentralJets = 0;
  double etjetmax  = 0.0;
  double etajetmax = 4.999;
  vector<Jet*> sortedJets;
  vector<Jet*> sortedJetsAll;
  vector<Jet*> sortedJetsLowPt;
  if(CleanJetsNoPtCut->GetEntries() > 0) {
    etjetmax  = CleanJetsNoPtCut->At(0)->Pt();
    etajetmax = CleanJetsNoPtCut->At(0)->Eta();
  }

  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
    if(CleanJetsNoPtCut->At(i)->RawMom().Pt() <= 7) continue;
    Jet* jet_a = new Jet(CleanJetsNoPtCut->At(i)->Px()*(1.0+fJetScaleSyst),
   			 CleanJetsNoPtCut->At(i)->Py()*(1.0+fJetScaleSyst),
   			 CleanJetsNoPtCut->At(i)->Pz()*(1.0+fJetScaleSyst),
   			 CleanJetsNoPtCut->At(i)->E() *(1.0+fJetScaleSyst));
    int nCloseStdJet = -1;
    double deltaRMin = 999.;
    for(UInt_t nj=0; nj<fCaloJet0->GetEntries(); nj++){
      const CaloJet *jet = fCaloJet0->At(nj);
      Double_t deltaR = MathUtils::DeltaR(jet_a->Mom(),jet->Mom());
      if(deltaR < deltaRMin) {
   	nCloseStdJet = nj;
   	deltaRMin = deltaR;
      }
    }
    jet_a->SetMatchedMCFlavor(CleanJetsNoPtCut->At(i)->MatchedMCFlavor());
    jet_a->SetCombinedSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexBJetTagsDisc());
    jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
    jet_a->SetJetProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetProbabilityBJetTagsDisc());
    jet_a->SetJetBProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetBProbabilityBJetTagsDisc());
    jet_a->SetTrackCountingHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc());
    jet_a->SetTrackCountingHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighPurBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighEffBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighPurBJetTagsDisc());
    sortedJetsAll.push_back(jet_a);

    if(nCloseStdJet >= 0){
      hDwwBTagCheck[0]->Fill(TMath::Max(TMath::Min(CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc(),17.999),-1.999),NNLOWeight->GetVal());
      hDwwBTagCheck[1]->Fill(TMath::Max(TMath::Min(fCaloJet0->At(nCloseStdJet)->TrackCountingHighEffBJetTagsDisc(),17.999),-1.999),NNLOWeight->GetVal());
    }
  }

  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
    if(TMath::Abs(CleanJetsNoPtCut->At(i)->Eta()) < fEtaJetCut &&
       CleanJetsNoPtCut->At(i)->Pt()*(1.0+fJetScaleSyst) > fPtJetCut){
      Jet* jet_b = new Jet(CleanJetsNoPtCut->At(i)->Px()*(1.0+fJetScaleSyst),
     			   CleanJetsNoPtCut->At(i)->Py()*(1.0+fJetScaleSyst),
   			   CleanJetsNoPtCut->At(i)->Pz()*(1.0+fJetScaleSyst),
   			   CleanJetsNoPtCut->At(i)->E() *(1.0+fJetScaleSyst));
      nCentralJets++;
      sortedJets.push_back(jet_b);
    }
  }
  for(UInt_t i=0; i<sortedJetsAll.size(); i++){
    bool overlap = kFALSE;
    for(UInt_t j=0; j<sortedJets.size(); j++){
      if(sortedJetsAll[i]->Pt() == sortedJets[j]->Pt() ||
        (sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc() == sortedJets[j]->CombinedSecondaryVertexBJetTagsDisc() &&
	 sortedJetsAll[i]->JetBProbabilityBJetTagsDisc()	 == sortedJets[j]->JetBProbabilityBJetTagsDisc() &&
	 sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc()	 == sortedJets[j]->TrackCountingHighPurBJetTagsDisc())
        ) {
        sortedJets[j]->SetMatchedMCFlavor(sortedJetsAll[i]->MatchedMCFlavor());
        sortedJets[j]->SetCombinedSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetCombinedSecondaryVertexMVABJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc());
        sortedJets[j]->SetJetProbabilityBJetTagsDisc(sortedJetsAll[i]->JetProbabilityBJetTagsDisc());
        sortedJets[j]->SetJetBProbabilityBJetTagsDisc(sortedJetsAll[i]->JetBProbabilityBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighEffBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighPurBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighEffBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighEffBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighPurBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighPurBJetTagsDisc());        
  	overlap = kTRUE;
        break;
      }
    }
    if(overlap == kFALSE){
      sortedJetsLowPt.push_back(sortedJetsAll[i]);
    }
  }

  MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, fVertices->At(0), 0.1, 8.0, 5.0, fIntRadius);
  double pMET[3] = {metTools.GetProjectedMet(leptons,PFMetStd),
  		    metTools.GetProjectedTrackMet(leptons),
  		    metTools.GetProjectedMet(leptons)};
  //Met* newMET = new Met();
  //newMET->SetMex(metTools.GetMinimumTrackMet(PFMetStd).Px());
  //newMET->SetMey(metTools.GetMinimumTrackMet(PFMetStd).Py());
  //newMET->SetSumEt(PFMetStd->SumEt());
  Met* newMET = new Met(PFMetStd->Px(), PFMetStd->Py());
  newMET->SetSumEt(PFMetStd->SumEt());

  if(fPrintDebug == kTRUE && fMuons    ->GetEntries() > 50) printf("WARNING, more than 50     muon cands\n");
  if(fPrintDebug == kTRUE && fElectrons->GetEntries() > 50) printf("WARNING, more than 50 electron cands\n");
  
  Int_t Nmu[8]  = {0, 0, 0, 0, 0, 0, 0, 0};
  Int_t Qmu[8][50];
  Double_t Ptmu[8][50];
  for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
    const Muon *mu = fMuons->At(i);
    if(mu->BestTrk() == 0) continue;

    Qmu[0][Nmu[0]] = mu->Charge(); Ptmu[0][Nmu[0]] = mu->Pt(); Nmu[0]++;

    if(mu->BestTrk()->Pt() <= 10.0) continue;
    if(TMath::Abs(mu->BestTrk()->Eta()) >= 2.4) continue;
    if(!(mu->HasGlobalTrk() || mu->IsTrackerMuon())) continue;
    if(fVertices->GetEntries() == 0) continue;
    Qmu[1][Nmu[1]] = mu->Charge(); Ptmu[1][Nmu[1]] = mu->Pt(); Nmu[1]++;

    const Track *mt = mu->BestTrk();
    Double_t d0_real = 99999.;
    Double_t dz_real = 99999.;
    if (mt) {
      d0_real = TMath::Abs(mt->D0Corrected(*fVertices->At(0)));
      dz_real = TMath::Abs(mt->DzCorrected(*fVertices->At(0)));
    }

    Double_t totalIso =  IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0), 0.1, 1.0, 0.3, 0.0, fIntRadius);
    if(fPrintDebug == kTRUE){
    double ptErrOverPt = mu->BestTrk()->PtErr()/mu->BestTrk()->Pt();
    //if(mu->HasGlobalTrk()) ptErrOverPt = mu->BestTrk()->PtErr()/mu->GlobalTrk()->Pt();
    printf("mu(%d): %f %f %f %f - %f %d %d %d %f - %d %d %d - %d %d %f - %f %f - %f | ",Nmu[1],mu->Pt(),mu->Eta(),mu->Phi(),mu->Charge(),
                                        mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof(),mu->NSegments(),mu->NMatches(),mu->NValidHits(),mu->TrkKink(),
					mu->HasGlobalTrk(),mu->IsTrackerMuon(),mu->Quality().Quality(MuonQuality::TMLastStationTight),
					mu->NTrkLayersHit(),mu->BestTrk()->NPixelHits(),ptErrOverPt,
					d0_real,dz_real,totalIso/mu->BestTrk()->Pt());
    }

    Qmu[2][Nmu[2]] = mu->Charge(); Ptmu[2][Nmu[2]] = mu->Pt(); Nmu[2]++;

    double fD0Cut = 0.02;
    if (mu->Pt() <= 20.0) fD0Cut = 0.01;
    Bool_t passD0cut = MuonTools::PassD0Cut(mu, fVertices, fD0Cut, 0);
    if (!passD0cut) continue;
    Qmu[3][Nmu[3]] = mu->Charge(); Ptmu[3][Nmu[3]] = mu->Pt(); Nmu[3]++;

    Bool_t passDZcut = MuonTools::PassDZCut(mu, fVertices, 0.1, 0);
    if (!passDZcut) continue;
    Qmu[4][Nmu[4]] = mu->Charge(); Ptmu[4][Nmu[4]] = mu->Pt(); Nmu[4]++;

    Bool_t passId = ((mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
	             (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
	             (mu->IsTrackerMuon() &&
	              mu->Quality().Quality(MuonQuality::TMLastStationTight))) &&
		      mu->BestTrk() != 0 &&
                      mu->NTrkLayersHit() > 5 &&
		      mu->IsPFMuon() == kTRUE &&
                      mu->BestTrk()->NPixelHits() > 0 &&
                      mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
		      mu->TrkKink() < 20.0;
    if(fIsOldSelection == true){
      passId = ((mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
     		(mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
     		(mu->IsTrackerMuon() &&
     		 mu->Quality().Quality(MuonQuality::TMLastStationTight))) &&
    	 	 mu->BestTrk() != 0 &&
     		 mu->BestTrk()->NHits() > 10 &&
    	 	 mu->BestTrk()->NPixelHits() > 0 &&
    	 	 mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1 &&
    	 	 mu->TrkKink() < 20.0;
    }
    if(!passId) continue;
    Qmu[5][Nmu[5]] = mu->Charge(); Ptmu[5][Nmu[5]] = mu->Pt(); Nmu[5]++;

    // trick to change the signal region cut
    Double_t pfIsoCutValue = 9999;
    if (mu->AbsEta() < 1.479) {
      if (mu->Pt() > 20) {
    	pfIsoCutValue = 0.13;
      } else {
    	pfIsoCutValue = 0.06;
      }
    } else {
      if (mu->Pt() > 20) {
        pfIsoCutValue = 0.09;
      } else {
        pfIsoCutValue = 0.05;
      }
    }

    Bool_t passIso = kFALSE;
    for (UInt_t nmu=0; nmu<CleanMuons->GetEntries(); nmu++) {   
      if(mu == CleanMuons->At(nmu)) { passIso = kTRUE; break;}
    }
    //if(fIsOldSelection == true){
    //  if((totalIso <  mu->Pt()*pfIsoCutValue && !passIso) ||
    //     (totalIso >= mu->Pt()*pfIsoCutValue &&  passIso)) {printf("DANGER: %f %f %d %d\n",totalIso/mu->Pt(),pfIsoCutValue,totalIso<mu->Pt()*pfIsoCutValue,passIso); assert(0);}
    //}

    if(!passIso) continue;
    Qmu[6][Nmu[6]] = mu->Charge(); Ptmu[6][Nmu[6]] = mu->Pt(); Nmu[6]++;

   if(PFMetStd->Pt() > 20) {Qmu[7][Nmu[7]] = mu->Charge(); Ptmu[7][Nmu[7]] = mu->Pt(); Nmu[7]++;};
  }
  Int_t Nel[8]  = {0, 0, 0, 0, 0, 0, 0, 0};
  Int_t Qel[8][50];
  Double_t Ptel[8][50];
  for (UInt_t i=0; i<fElectrons->GetEntries(); i++) {
    const Electron *el = fElectrons->At(i);
    if(el->SCluster() == 0) continue;

    Qel[0][Nel[0]] = el->Charge(); Ptel[0][Nel[0]] = el->Pt(); Nel[0]++;

    if(el->Pt() <= 10.0)   continue;
    if(el->AbsEta() >= 2.5)   continue;
    if(fVertices->GetEntries() == 0) continue;
    Qel[1][Nel[1]] = el->Charge(); Ptel[1][Nel[1]] = el->Pt(); Nel[1]++;

    Double_t d0_real = TMath::Abs(el->GsfTrk()->D0Corrected(*fVertices->At(0)));
    Double_t dz_real = TMath::Abs(el->GsfTrk()->DzCorrected(*fVertices->At(0)));

    Bool_t  passConvVeto = ElectronTools::PassConversionFilter(el, fConversions, 
                                                               fBeamSpot->At(0), 0, 1e-6, 2.0, kTRUE, kFALSE);      
    ElectronOArr *tempIsoElectrons = new  ElectronOArr;
    MuonOArr	 *tempIsoMuons     = new  MuonOArr;
    PFCandidateCol *fPFNoPileUpCands = GetObjThisEvt<PFCandidateCol>("PFNoPileUp");
    Double_t IsoOverPt = IsolationTools::PFElectronIsolation2012(el, fVertices->At(0), fPFNoPileUpCands, 
     fPileupEnergyDensity, ElectronTools::kEleEANoCorr, tempIsoElectrons, tempIsoMuons, 0.4, kFALSE);
    delete tempIsoElectrons;
    delete tempIsoMuons;
    Double_t eta = el->SCluster()->AbsEta();
    Double_t IsoCut = -1;
    if (el->Pt() <  20 && eta <  0.800  	      ) IsoCut = 0.150;
    if (el->Pt() <  20 && eta >= 0.800 && eta < 1.479 ) IsoCut = 0.150;
    if (el->Pt() <  20 && eta >= 1.479  	      ) IsoCut = 0.150;
    if (el->Pt() >= 20 && eta <  0.800  	      ) IsoCut = 0.150;
    if (el->Pt() >= 20 && eta >= 0.800 && eta < 1.479 ) IsoCut = 0.150;
    if (el->Pt() >= 20 && eta >= 1.479  	      ) IsoCut = 0.150;

    Double_t totalIso = IsolationTools::PFElectronIsolation(el, fPFCandidates, fVertices->At(0), 0.1, 1.0, 0.4, fIntRadius);
    Double_t pfIsoCutValue = 9999;
    if (el->SCluster()->AbsEta() < 1.479) {
      if (el->Pt() > 20) {
    	pfIsoCutValue = 0.13;
      } else {
    	pfIsoCutValue = 0.13;
      }
    } else {
      if (el->Pt() > 20) {
    	pfIsoCutValue = 0.09;
      } else {
    	pfIsoCutValue = 0.09;
      }
    }
    bool passIso = IsoOverPt < IsoCut;
    if(fIsOldSelection == true) {
      passIso = totalIso < (el->Pt()*pfIsoCutValue);
    }

    Double_t hOverE       = el->HadronicOverEm();
    Double_t sigmaee      = el->CoviEtaiEta();
    Double_t deltaPhiIn   = TMath::Abs(el->DeltaPhiSuperClusterTrackAtVtx());
    Double_t deltaEtaIn   = TMath::Abs(el->DeltaEtaSuperClusterTrackAtVtx());
    Double_t eOverP       = el->ESuperClusterOverP();
    Double_t fBrem        = el->FBrem();
    if(fPrintDebug == kTRUE){
    printf("el(%d): %f %f %f %f - %f %f %f - %d - %f %f %f %f %f %f FO: %d | ",Nel[1],el->Pt(),el->Eta(),el->Phi(),el->Charge(),
                                        d0_real,dz_real,IsoOverPt,
                                        passConvVeto,
					hOverE,sigmaee,deltaPhiIn,deltaEtaIn,eOverP,fBrem,ElectronTools::PassCustomID(el, ElectronTools::kVBTFWorkingPointFakeableId));
    }
 
    if (!passIso) continue;
    Qel[2][Nel[2]] = el->Charge(); Ptel[2][Nel[2]] = el->Pt(); Nel[2]++;

    if (passConvVeto == kFALSE) continue;
    if (el->CorrectedNExpectedHitsInner() != 0) continue;
    Qel[3][Nel[3]] = el->Charge(); Ptel[3][Nel[3]] = el->Pt(); Nel[3]++;

    Bool_t passD0cut = ElectronTools::PassD0Cut(el, fVertices, 0.02, 0);
    if (!passD0cut) continue;
    Qel[4][Nel[4]] = el->Charge(); Ptel[4][Nel[4]] = el->Pt(); Nel[4]++;

    Bool_t passDZcut = ElectronTools::PassDZCut(el, fVertices, 0.1, 0);
    if (!passDZcut) continue;
    Qel[5][Nel[5]] = el->Charge(); Ptel[5][Nel[5]] = el->Pt(); Nel[5]++;

    Bool_t passId = kFALSE;
    for (UInt_t nel=0; nel<CleanElectrons->GetEntries(); nel++) {   
      if(el == CleanElectrons->At(nel)) { passId = kTRUE; break;}
    }
    if(!passId) continue;
    Qel[6][Nel[6]] = el->Charge(); Ptel[6][Nel[6]] = el->Pt(); Nel[6]++;

    if(PFMetStd->Pt() > 20) {Qel[7][Nel[7]] = el->Charge(); Ptel[7][Nel[7]] = el->Pt(); Nel[7]++;};
  }
  if(fPrintDebug == kTRUE){
  if(Nmu[0] > 0 || Nel[0] > 0) printf(" %f %d %d %d\n",PFMetStd->Pt(),fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
  }
  Int_t isQ[2];
  Int_t compType[3] = {0, 0, 0};
  for(int i=0; i<8; i++){
    if(Nmu[i]>= 2) {
      isQ[0] = 0; isQ[1] = 0; double isPt20 = kFALSE;
      for(int k=0; k<Nmu[i]; k++){
        if(Qmu[i][k] > 0) isQ[0] = +1;
        if(Qmu[i][k] < 0) isQ[1] = -1;
	if(Ptmu[i][k] > 20) isPt20 = kTRUE;
      }
      if(isQ[0]*isQ[1] < 0 && isPt20 == kTRUE) {hDwwXS[0]->Fill((double)i,NNLOWeight->GetVal()); compType[0] = i+1;}
    }
    if(Nel[i]>= 2) {
      isQ[0] = 0; isQ[1] = 0; double isPt20 = kFALSE;
      for(int k=0; k<Nel[i]; k++){
        if(Qel[i][k] > 0) isQ[0] = +1;
        if(Qel[i][k] < 0) isQ[1] = -1;
	if(Ptel[i][k] > 20) isPt20 = kTRUE;
      }
      if(isQ[0]*isQ[1] < 0 && isPt20 == kTRUE) {hDwwXS[1]->Fill((double)i,NNLOWeight->GetVal()); compType[1] = i+1;}
    }
    if(Nmu[i]>= 1 && Nel[i] >=1) {
      isQ[0] = 0; isQ[1] = 0;
      for(int k0=0; k0<Nmu[i]; k0++){
        for(int k1=0; k1<Nel[i]; k1++){
          if(Qmu[i][k0]        *   Qel[i][k1] < 0 &&
	   (Ptmu[i][k0] > 20.0 || Ptel[i][k1] > 20.0) &&
	    Ptmu[i][k0] <= Ptel[i][k1]) isQ[0] = 1;
          if(Qmu[i][k0]        *   Qel[i][k1] < 0 &&
	   (Ptmu[i][k0] > 20.0 || Ptel[i][k1] > 20.0) &&
	    Ptmu[i][k0] >  Ptel[i][k1]) isQ[1] = 1;
        }
      }
      if(isQ[0] == 1) {hDwwXS[2]->Fill((double)i,NNLOWeight->GetVal()); compType[2] = i+1;}
      if(isQ[1] == 1) {hDwwXS[3]->Fill((double)i,NNLOWeight->GetVal()); compType[2] = i+1;}
    }
    if(Nmu[i]+Nel[i] >= 2){
      isQ[0] = 0;
      if(Nmu[i] > 0){
        for(int k0=0; k0<Nmu[i]; k0++){
          for(int k1=0; k1<Nel[i]; k1++){
            if(Qmu[i][k0]        *   Qel[i][k1] < 0 &&
	     (Ptmu[i][k0] > 20.0 || Ptel[i][k1] > 20.0)) isQ[0] = 1;
          }
          for(int k1=0; k1<Nmu[i]; k1++){
            if(Qmu[i][k0]        *   Qmu[i][k1] < 0 &&
	     (Ptmu[i][k0] > 20.0 || Ptmu[i][k1] > 20.0)) isQ[0] = 1;
          }
	}
      }
      if(Nel[i] > 0 && isQ[0] == 0){
        for(int k0=0; k0<Nel[i]; k0++){
          for(int k1=0; k1<Nel[i]; k1++){
            if(Qel[i][k0]        *   Qel[i][k1] < 0 &&
	     (Ptel[i][k0] > 20.0 || Ptel[i][k1] > 20.0)) isQ[0] = 1;
          }
          for(int k1=0; k1<Nmu[i]; k1++){
            if(Qel[i][k0]        *   Qel[i][k1] < 0 &&
	     (Ptel[i][k0] > 20.0 || Ptel[i][k1] > 20.0)) isQ[0] = 1;
          }
	}
      }
      if(isQ[0] == 1) hDwwXS[4]->Fill((double)i,NNLOWeight->GetVal());
    }
  }
  if(fPrintDebug == kTRUE){
  if(compType[0] > 0 || compType[1] > 0 || compType[2] > 0) printf("TYPE: %d %d %d %d %d %d\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
  								   compType[0],compType[1],compType[2]);
  }
  // Begin skim level plots
  hDwwSkim[0]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
  hDwwSkim[1]->Fill(TMath::Min((double)nCentralJets,9.499),NNLOWeight->GetVal());
  // End skim level plots

  if(leptons->GetEntries() >= 2){
    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));
    TVector3 v0(leptons->At(0)->Px()	    ,leptons->At(0)->Py()	 ,leptons->At(0)->Pz());
    TVector3 v1(leptons->At(1)->Px()	    ,leptons->At(1)->Py()	 ,leptons->At(1)->Pz());

    hDwwSkim[2]->Fill(v0.Angle(v1) * 180./TMath::Pi(),NNLOWeight->GetVal());
    /*
    printf("CANDIDATE: %d %d %d | %f %f %f %d %f | %f %f %f %d %f - %d %f\n",
    	    fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
            leptons->At(0)->Pt(),leptons->At(0)->Eta(),leptons->At(0)->Phi(),leptons->At(0)->ObjType(),leptons->At(0)->Charge(),
            leptons->At(1)->Pt(),leptons->At(1)->Eta(),leptons->At(1)->Phi(),leptons->At(1)->ObjType(),leptons->At(1)->Charge(),
            leptons->GetEntries(),dilepton.Mass());
    */
  }

  hDwwPresel[0]->Fill(TMath::Min((double)leptons->GetEntries(),9.499),NNLOWeight->GetVal());
  if(leptons->GetEntries() >= 1)
    hDwwPresel[1]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
  if(leptons->GetEntries() >= 2 && leptons->At(0)->Pt() > 20)
    hDwwPresel[2]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
  // Minimun Pt, Nleptons==2 requirements
  if (leptons->GetEntries() >= 2 &&
      leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10){

    if(fUsePDFs == kTRUE){
      for(UInt_t i=0; i<PDFArr->GetEntries(); i++){
        hDWWPDF[i]->Fill(TMath::Min(newMET->Pt(),199.999),PDFArr->At(i)*NNLOWeight->GetVal());
      }
    }

    CompositeParticle dilepton;;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));

    int pairType = -1;
    if (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon )
      pairType = 0;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron)
      pairType = 1;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon)
      pairType = 2;
    else if(leptons->At(0)->ObjType() == kMuon     && leptons->At(1)->ObjType() == kElectron)
      pairType = 3;
    else {
      cout << "Hey, this is not possible, leptonTypes: "
    	   << leptons->At(0)->ObjType() << " - " 
           << leptons->At(1)->ObjType() << endl;
    }

    // Check consistency with ntuple
    bool isMixedState = kFALSE;
    if(leptons->GetEntries() == 2){
      if((leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon    )||
    	 (leptons->At(0)->ObjType() == kMuon	 && leptons->At(1)->ObjType() == kElectron)){
    	isMixedState = kTRUE;
      }
    }
    if ((leptons->GetEntries() == 2 && PFMetStd->Pt() > 20 &&
    	 leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10) ||
    	leptons->GetEntries() == 3 ||
    	isMixedState == kTRUE){ 
      hDwwPresel[13]->Fill(0.0,NNLOWeight->GetVal());
      if(leptons->At(1)->Pt() > 20 && leptons->GetEntries() == 2){
    	hDwwPresel[13]->Fill(1.0,NNLOWeight->GetVal());
    	if(dilepton.Charge() == 0){
    	  hDwwPresel[13]->Fill(2.0,NNLOWeight->GetVal());
    	  if(zDiffMax < 0.2){
    	    hDwwPresel[13]->Fill(3.0,NNLOWeight->GetVal());
    	    if(dilepton.Mass() > 12){
    	      hDwwPresel[13]->Fill(4.0,NNLOWeight->GetVal());
    	      if(pairType >= 2 || TMath::Abs(dilepton.Mass()-91.1876) > 15){
    		hDwwPresel[13]->Fill(5.0,NNLOWeight->GetVal());
    	      }
    	    }
    	  }
    	}
      }
    }
    UInt_t leptonGenType[2] = {0, 0};
    // look for real W/Z -> l leptons
    for(UInt_t i=0; i<2; i++) {
      for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
        MCParticle *gen = GenLeptons->At(j);
        if(MathUtils::DeltaR(gen->Mom(), leptons->At(i)->Mom()) < 0.05) {
	  leptonGenType[i] = 1;
	  break;
	}
      }
      // look for other leptons
      if(leptonGenType[i] == 0){
      	for (UInt_t j=0; j<GenAllLeptons->GetEntries(); j++) {
      	  MCParticle *gen = GenAllLeptons->At(j);
      	  if(MathUtils::DeltaR(gen->Mom(), leptons->At(i)->Mom()) < 0.05) {
	    leptonGenType[i] = 2;
	    break;
	  }
      	}
      }
      // look for ISR photons
      if(leptonGenType[i] == 0){
      	for (UInt_t j=0; j<GenISRPhotons->GetEntries(); j++) {
      	  MCParticle *gen = GenISRPhotons->At(j);
      	  if(MathUtils::DeltaR(gen->Mom(), leptons->At(i)->Mom()) < 0.05) {
	    leptonGenType[i] = 3;
	    break;
	  }
      	}
      }
      // look for Rad photons
      if(leptonGenType[i] == 0){
      	for (UInt_t j=0; j<GenRadPhotons->GetEntries(); j++) {
      	  MCParticle *gen = GenRadPhotons->At(j);
      	  if(MathUtils::DeltaR(gen->Mom(), leptons->At(i)->Mom()) < 0.05) {
	    leptonGenType[i] = 4;
	    break;
	  }
      	}
      }
    } // loop 0,1

    hDwwFake[0+100*pairType]->Fill(leptonGenType[0],NNLOWeight->GetVal());
    hDwwFake[1+100*pairType]->Fill(leptonGenType[1],NNLOWeight->GetVal());
    hDwwFake[2+100*pairType]->Fill(leptonGenType[0]+5*leptonGenType[1],NNLOWeight->GetVal());

    hDwwPresel[3]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    hDwwPresel[4]->Fill((double)nCentralJets,NNLOWeight->GetVal());
    hDwwPresel[5]->Fill((double)dilepton.Charge(),NNLOWeight->GetVal());

    if(dilepton.Mass() < -5.0){
      printf("**************************************\n");
      printf("%f %d %f\n",dilepton.Mass(),nCentralJets,dilepton.Charge());
      printf("%f %f | %f %f | %f %f | %f %f | %d %d\n",leptons->At(0)->Pt(),leptons->At(1)->Pt(),
                                                       leptons->At(0)->Eta(),leptons->At(1)->Eta(),
                                                       leptons->At(0)->Phi(),leptons->At(1)->Phi(),
                                                       leptons->At(0)->Charge(),leptons->At(1)->Charge(),
                                                       leptons->At(0)->ObjType(),leptons->At(1)->ObjType());
      for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
        double isoAux = 1.0 * CleanMuons->At(j)->IsoR03SumPt() + 
                        1.0 * CleanMuons->At(j)->IsoR03EmEt()  +
  	         	1.0 * CleanMuons->At(j)->IsoR03HadEt();
        printf("m(%d): %f - %f %f %f ==> %f\n",j,CleanMuons->At(j)->Pt(),
               CleanMuons->At(j)->IsoR03SumPt(),CleanMuons->At(j)->IsoR03EmEt(),CleanMuons->At(j)->IsoR03HadEt(),
               isoAux);
      }
      for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
        double isoAux = 1.0 * CleanElectrons->At(j)->TrackIsolationDr03() + 
                        1.0 * CleanElectrons->At(j)->EcalRecHitIsoDr04();
        printf("e(%d): %f - %f %f ==> %f\n",j,CleanElectrons->At(j)->Pt(),
               CleanElectrons->At(j)->TrackIsolationDr03(),CleanElectrons->At(j)->EcalRecHitIsoDr04(),
               isoAux);
      }
      printf("**************************************\n");
    }

    if(dilepton.Mass() > 12 && dilepton.Charge() == 0 && leptons->GetEntries() == 2) {
      if(GenMet->GetEntries() > 0){
        hDwwMET[0]->Fill(TMath::Min(GenMet->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      }
      hDwwMET[1]->Fill(TMath::Min(cleanMet->Pt(),199.999),NNLOWeight->GetVal());
      hDwwMET[2]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
      hDwwMET[3]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
      hDwwMET[4]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
    
      int nVertices = (int)(fVertices->GetEntries()/3.); if(nVertices > 4) nVertices = 4;
      //int nPU = (int)(fPileupInfo->At(0)->GetPU_NumInteractions()/3.); if(nPU > 4) nPU = 4;
      hDwwMET[5 +nVertices]->Fill(TMath::Max(TMath::Min(PFMetStd->Px(),99.999),-99.99),NNLOWeight->GetVal());
      hDwwMET[10+nVertices]->Fill(TMath::Max(TMath::Min(PFMetStd->Py(),99.999),-99.99),NNLOWeight->GetVal());
      hDwwMET[15+nVertices]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
      hDwwMET[20+nVertices]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
      if(nCentralJets == 0){
        hDwwMET[25+nVertices]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[30+nVertices]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[35+nVertices]->Fill(TMath::Min(metTools.GetCorrectedMet().Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[40+nVertices]->Fill(TMath::Min(metTools.GetMinimumMet(PFMetStd).Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[45+nVertices]->Fill(TMath::Min(metTools.GetCorrectedTrackMet().Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[50+nVertices]->Fill(TMath::Min(metTools.GetMinimumTrackMet(PFMetStd).Pt(),199.999),NNLOWeight->GetVal());
        hDwwMET[55+nVertices]->Fill(TMath::Min(pMET[0],199.999),NNLOWeight->GetVal());
        hDwwMET[60+nVertices]->Fill(TMath::Min(pMET[1],199.999),NNLOWeight->GetVal());
        hDwwMET[65+nVertices]->Fill(TMath::Min(pMET[2],199.999),NNLOWeight->GetVal());
        hDwwMET[70+nVertices]->Fill(TMath::Min(TMath::Min(pMET[0],pMET[1]),199.999),NNLOWeight->GetVal());
        hDwwMET[75+nVertices]->Fill(TMath::Min(TMath::Min(pMET[0],pMET[2]),199.999),NNLOWeight->GetVal());
        hDwwMET[80+nVertices]->Fill(TMath::Min(TMath::Min(pMET[1],pMET[2]),199.999),NNLOWeight->GetVal());
      }

      // Begin jet studies
      UInt_t   nCentralJetsColl[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t maxPtJetsColl[10]    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      for(UInt_t nj=0; nj<fCaloJet0->GetEntries(); nj++){
	const CaloJet *jet = fCaloJet0->At(nj);

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
          UInt_t n = CleanElectrons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
            if (deltaR < 0.3) {
              isElectronOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
          UInt_t n = CleanMuons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
            if (deltaR < 0.3) {
              isMuonOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut/2.0) nCentralJetsColl[0]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[0]) maxPtJetsColl[0] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fCaloJet1->GetEntries(); nj++){
	const CaloJet *jet = fCaloJet1->At(nj);        

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
          UInt_t n = CleanElectrons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
            if (deltaR < 0.3) {
              isElectronOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
          UInt_t n = CleanMuons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
            if (deltaR < 0.3) {
              isMuonOverlap = kTRUE;
              break;		   
            }	 
          }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut/2.0) nCentralJetsColl[1]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[1]) maxPtJetsColl[1] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fCaloJet2->GetEntries(); nj++){
	const CaloJet *jet = fCaloJet2->At(nj);        

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
          UInt_t n = CleanElectrons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
            if (deltaR < 0.3) {
              isElectronOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
          UInt_t n = CleanMuons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
            if (deltaR < 0.3) {
              isMuonOverlap = kTRUE;
              break;		   
            }	 
          }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut/2.0) nCentralJetsColl[2]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[2]) maxPtJetsColl[2] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fCaloJet3->GetEntries(); nj++){
	const CaloJet *jet = fCaloJet3->At(nj);        

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
          UInt_t n = CleanElectrons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
            if (deltaR < 0.3) {
              isElectronOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
          UInt_t n = CleanMuons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
            if (deltaR < 0.3) {
              isMuonOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut/2.0) nCentralJetsColl[3]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[3]) maxPtJetsColl[3] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fTrackJet0->GetEntries(); nj++){
	const TrackJet *jet = fTrackJet0->At(nj);        

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
          UInt_t n = CleanElectrons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
            if (deltaR < 0.3) {
              isElectronOverlap = kTRUE;
              break;		   
            }
          }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
          UInt_t n = CleanMuons->GetEntries();
          for (UInt_t j=0; j<n; j++) {
            Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
            if (deltaR < 0.3) {
              isMuonOverlap = kTRUE;
              break;
            }	 
          }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut*2./3.) nCentralJetsColl[5]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[5]) maxPtJetsColl[5] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      int nGoodJets[2] = {0, 0};
      for(UInt_t nj=0; nj<fPFJet0->GetEntries(); nj++){
	const PFJet *jet = fPFJet0->At(nj);	 

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
    	  UInt_t n = CleanElectrons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
    	    if (deltaR < 0.3) {
    	      isElectronOverlap = kTRUE;
    	      break;
    	    }
    	  }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
    	  UInt_t n = CleanMuons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
    	    if (deltaR < 0.3) {
    	      isMuonOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isMuonOverlap) continue;
	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut) nCentralJetsColl[6]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[6]) maxPtJetsColl[6] = jet->Pt()*(1.0+fJetScaleSyst);

	double totalE = jet->ChargedEmEnergy()+jet->ChargedHadronEnergy()+jet->NeutralEmEnergy()+jet->NeutralHadronEnergy();	  
        bool isGoodJet = jet->NeutralEmEnergy()/totalE < 0.99 && jet->NeutralHadronEnergy()/totalE < 0.99 &&
	        	 jet->NConstituents() > 1;
	if(jet->AbsEta() < 2.4) isGoodJet = isGoodJet && jet->ChargedEmEnergy()/totalE < 0.99 && jet->ChargedHadronEnergy()/totalE > 0 &&
	        					 jet->ChargedMultiplicity() > 0;
	if(jet->Pt() > fPtJetCut		      ) nGoodJets[0]++;
	if(jet->Pt() > fPtJetCut && isGoodJet == kTRUE) nGoodJets[1]++;
	if(1){
          if(jet->Pt() > fPtJetCut){
	    if(jet->AbsEta() >= 2.4){
	      hDwwJetVar[ 0]->Fill(TMath::Min(jet->NeutralEmEnergy()/totalE,0.999),NNLOWeight->GetVal()); // < 0.99
	      hDwwJetVar[ 1]->Fill(TMath::Min(jet->NeutralHadronEnergy()/totalE,0.999),NNLOWeight->GetVal()); // < 0.99
	      hDwwJetVar[ 2]->Fill(TMath::Min((double)jet->NConstituents(),49.499),NNLOWeight->GetVal()); // > 1
	    } else {
	      hDwwJetVar[ 3]->Fill(TMath::Min(jet->NeutralEmEnergy()/totalE,0.999),NNLOWeight->GetVal()); // < 0.99
	      hDwwJetVar[ 4]->Fill(TMath::Min(jet->NeutralHadronEnergy()/totalE,0.999),NNLOWeight->GetVal()); // < 0.99
	      hDwwJetVar[ 5]->Fill(TMath::Min((double)jet->NConstituents(),49.499),NNLOWeight->GetVal()); // > 1
	      hDwwJetVar[ 6]->Fill(TMath::Min(jet->ChargedEmEnergy()/totalE,0.999),NNLOWeight->GetVal()); // < 0.99
	      hDwwJetVar[ 7]->Fill(TMath::Min(jet->ChargedHadronEnergy()/totalE,0.999),NNLOWeight->GetVal()); // > 0
	      hDwwJetVar[ 8]->Fill(TMath::Min((double)jet->ChargedMultiplicity(),49.499),NNLOWeight->GetVal()); // > 0
	    }
	  }
	  if(jet->AbsEta() < 2.4) hDwwJetVar[ 9]->Fill(TMath::Min(jet->Pt(),199.999),NNLOWeight->GetVal());
	  else  		  hDwwJetVar[10]->Fill(TMath::Min(jet->Pt(),199.999),NNLOWeight->GetVal());
	  if(jet->AbsEta() < 2.4 && isGoodJet == kTRUE) hDwwJetVar[11]->Fill(TMath::Min(jet->Pt(),199.999),NNLOWeight->GetVal());
	  else if(                  isGoodJet == kTRUE) hDwwJetVar[12]->Fill(TMath::Min(jet->Pt(),199.999),NNLOWeight->GetVal());
        }
      }
      hDwwJetVar[13]->Fill(TMath::Min((double)nGoodJets[0],9.499),NNLOWeight->GetVal());
      hDwwJetVar[14]->Fill(TMath::Min((double)nGoodJets[1],9.499),NNLOWeight->GetVal());
      for(UInt_t nj=0; nj<fPFJet1->GetEntries(); nj++){
	const PFJet *jet = fPFJet1->At(nj);	 

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
    	  UInt_t n = CleanElectrons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
    	    if (deltaR < 0.3) {
    	      isElectronOverlap = kTRUE;
    	      break;		 
    	    }
    	  }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
    	  UInt_t n = CleanMuons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
    	    if (deltaR < 0.3) {
    	      isMuonOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut) nCentralJetsColl[7]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[7]) maxPtJetsColl[7] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fPFJet2->GetEntries(); nj++){
	const PFJet *jet = fPFJet2->At(nj);	 

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
    	  UInt_t n = CleanElectrons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
    	    if (deltaR < 0.3) {
    	      isElectronOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
    	  UInt_t n = CleanMuons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
    	    if (deltaR < 0.3) {
    	      isMuonOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isMuonOverlap) continue;

	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut) nCentralJetsColl[8]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[8]) maxPtJetsColl[8] = jet->Pt()*(1.0+fJetScaleSyst);
      }
      for(UInt_t nj=0; nj<fPFJet3->GetEntries(); nj++){
	const PFJet *jet = fPFJet3->At(nj);	 

	if(TMath::Abs(jet->Eta()) >= fEtaJetCut) continue;

	Bool_t isElectronOverlap = kFALSE;
	if (CleanElectrons) {
    	  UInt_t n = CleanElectrons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
                                              CleanElectrons->At(j)->SCluster()->Eta(), 
                                              jet->Phi(), jet->Eta());
    	    if (deltaR < 0.3) {
    	      isElectronOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isElectronOverlap) continue;
	Bool_t isMuonOverlap = kFALSE;
	if (CleanMuons) {
    	  UInt_t n = CleanMuons->GetEntries();
    	  for (UInt_t j=0; j<n; j++) {
    	    Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
    	    if (deltaR < 0.3) {
    	      isMuonOverlap = kTRUE;
    	      break;		 
    	    }    
    	  }
	}
	if (isMuonOverlap) continue;
	if(jet->Pt()*(1.0+fJetScaleSyst) > fPtJetCut) nCentralJetsColl[9]++;
	if(jet->Pt()*(1.0+fJetScaleSyst) > maxPtJetsColl[9]) maxPtJetsColl[9] = jet->Pt()*(1.0+fJetScaleSyst);       
      }
      for(UInt_t nj=0; nj<10; nj++) {
	hDwwJet[    nj+100*pairType]->Fill(TMath::Min((double)nCentralJetsColl[nj],9.499),NNLOWeight->GetVal());
	hDwwJet[ 10+nj+100*pairType]->Fill(TMath::Min(maxPtJetsColl[nj],199.999),NNLOWeight->GetVal());
      }
      if(newMET->Pt() > 40.0){
	for(UInt_t nj=0; nj<10; nj++) {
          hDwwJet[ 20+nj+100*pairType]->Fill(TMath::Min((double)nCentralJetsColl[nj],9.499),NNLOWeight->GetVal());
          hDwwJet[ 30+nj+100*pairType]->Fill(TMath::Min(maxPtJetsColl[nj],199.999),NNLOWeight->GetVal());
	}
      }
      // End jet studies
    }

    if(dilepton.Mass() > 12 && dilepton.Charge() == 0) {
      hDwwMET[85]->Fill(TMath::Max(TMath::Min(PFMetStd->Px()-newMET->Px(),99.999),-99.99),NNLOWeight->GetVal());
      hDwwMET[86]->Fill(TMath::Max(TMath::Min(PFMetStd->Py()-newMET->Py(),99.999),-99.99),NNLOWeight->GetVal());
      hDwwMET[87]->Fill(TMath::Max(TMath::Min(PFMetStd->Pt()-newMET->Pt(),99.999),-99.99),NNLOWeight->GetVal());
      hDwwMET[88]->Fill(fabs(MathUtils::DeltaPhi(newMET->Phi(), PFMetStd->Phi()))*180/TMath::Pi(),NNLOWeight->GetVal());
    }
    // Angle between MET and closest lepton
    double deltaPhiMetLepton[4] = {fabs(MathUtils::DeltaPhi(newMET->Phi(), leptons->At(0)->Phi())),
    				   fabs(MathUtils::DeltaPhi(newMET->Phi(), leptons->At(1)->Phi())),
				   TMath::Pi()-fabs(MathUtils::DeltaPhi(newMET->Phi(), leptons->At(0)->Phi())),
				   TMath::Pi()-fabs(MathUtils::DeltaPhi(newMET->Phi(), leptons->At(1)->Phi()))};

    double mTW[2] = {TMath::Sqrt(2.0*leptons->At(0)->Pt()*newMET->Pt()*
    				(1.0 - cos(deltaPhiMetLepton[0]))),
    		     TMath::Sqrt(2.0*leptons->At(1)->Pt()*newMET->Pt()*
    				(1.0 - cos(deltaPhiMetLepton[1])))};

    double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
      deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

    double METdeltaPhilEt = TMath::Min(pMET[0],pMET[1]);

    if((pairType >= 2 || TMath::Abs(dilepton.Mass()-91.1876) > 15) &&
      nCentralJets == 0 &&
      dilepton.Mass() > 12 && dilepton.Charge() == 0){
      hDwwDeltaPhiMetLeptonMet->Fill(minDeltaPhiMetLepton * 180./TMath::Pi(),TMath::Min(newMET->Pt(),199.999),NNLOWeight->GetVal());
      hDwwFake[3+100*pairType]->Fill(leptonGenType[0],NNLOWeight->GetVal());
      hDwwFake[4+100*pairType]->Fill(leptonGenType[1],NNLOWeight->GetVal());
      hDwwFake[5+100*pairType]->Fill(leptonGenType[0]+5*leptonGenType[1],NNLOWeight->GetVal());
    }

    // Njet studies
    if(dilepton.Mass() > 12 && dilepton.Charge() == 0 && leptons->GetEntries() == 2 &&
       DirtyMuons->GetEntries() == 0 && leptons->At(1)->Pt() >= 20){
      int nPU = 0;
      if(fIsData == kFALSE){
        nPU = (int)(fPileupInfo->At(0)->GetPU_NumInteractions()/5.); if(nPU > 3) nPU = 3;
      }
      hDwwJetSel[0+100*nPU]->Fill(TMath::Min(etjetmax,399.999),NNLOWeight->GetVal());
      hDwwJetSel[1+100*nPU]->Fill(etajetmax,NNLOWeight->GetVal());
      hDwwJetSel[2+100*nPU]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
      if(METdeltaPhilEt > 20 && (pairType >= 2 || METdeltaPhilEt > 37.0 + fVertices->GetEntries()/2.0)){
        double deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(newMET->Phi(), 
        						 dilepton.Phi()))* 180./TMath::Pi();
        hDwwJetSel[3+100*pairType]->Fill(180.-deltaPhiDileptonMet,NNLOWeight->GetVal());
        hDwwJetSel[4+100*pairType]->Fill(TMath::Min(newMET->Pt()/dilepton.Pt(),3.99),NNLOWeight->GetVal());
        hDwwSelAlphaEP0->Fill(TMath::Min(newMET->Pt()/dilepton.Pt(),3.99),180.-deltaPhiDileptonMet);
        hDwwJetSel[5+100*pairType]->Fill(TMath::Min(TMath::Abs(dilepton.Mass()-91.1876),199.999),NNLOWeight->GetVal());  
        if(pairType == 2 || TMath::Abs(dilepton.Mass()-91.1876) > 15){
          hDwwJetSel[6+100*pairType]->Fill(minDeltaPhiMetLepton*180.0/TMath::Pi(),NNLOWeight->GetVal());
          hDwwJetSel[7+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());		
          hDwwJetSel[8+100*pairType]->Fill(TMath::Min(etjetmax,399.999),NNLOWeight->GetVal());
          hDwwJetSel[9+100*pairType]->Fill(etajetmax,NNLOWeight->GetVal());
        }
      }
    }

    LoadBranch(fTrackName);
    ObjArray<Track> *CleanTracks = new ObjArray<Track>;
    if(dilepton.Mass() > 12 && fVertices->GetEntries()>= 1 && leptons->GetEntries() == 2){
      for(UInt_t i=0; i<fTracks->GetEntries(); i++){
    	double pD0 = 0.0;
    	for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
    	  if(fVertices->At(i0)->NTracks() > 0){
    	    pD0 = fTracks->At(i)->D0Corrected(*fVertices->At(i0));
    	    break;
    	  }
    	}
        hDwwSelD0Phi->Fill(fTracks->At(i)->Phi0() * 180./TMath::Pi(),pD0,NNLOWeight->GetVal());
      }
      hDwwPresel[ 6]->Fill(TMath::Min((double)CleanTracks->GetEntries(),199.4999),NNLOWeight->GetVal());
      hDwwPresel[ 7]->Fill((double)DirtyMuons->GetEntries(),NNLOWeight->GetVal());
      hDwwPresel[ 8]->Fill(TMath::Min(newMET->MetSig(),19.999),NNLOWeight->GetVal());
      hDwwPresel[ 9]->Fill(TMath::Min(newMET->SumEt(),999.999),NNLOWeight->GetVal());
      hDwwPresel[10]->Fill(TMath::Min(zDiffMax,0.999),NNLOWeight->GetVal());
      hDwwPresel[11]->Fill(TMath::Min((double)fVertices->GetEntries(),19.4999),NNLOWeight->GetVal());
      hDwwPresel[12]->Fill(TMath::Min((double)fAllVertices->GetEntries(),19.4999),NNLOWeight->GetVal());
      if     (fAllVertices->GetEntries() == 0)                       hDwwPresel[14]->Fill(0.0,NNLOWeight->GetVal());
      else if(fVertices->GetEntries() == 0)                          hDwwPresel[14]->Fill(1.0,NNLOWeight->GetVal());
      else if(fAllVertices->GetEntries() == fVertices->GetEntries()) hDwwPresel[14]->Fill(2.0,NNLOWeight->GetVal());
      else if(fAllVertices->At(0) == fVertices->At(0))               hDwwPresel[14]->Fill(3.0,NNLOWeight->GetVal());
      else                                                           hDwwPresel[14]->Fill(4.0,NNLOWeight->GetVal());

      Int_t closestVtx[2] = {0, 0};
      Double_t distVtx[2] = {999.9, 999.9};
      Double_t distVtx0[2] = {999.9, 999.9};
      for (UInt_t j=0; j<TMath::Min((double)leptons->GetEntries(), 2.0); j++) {
	for (UInt_t i=0; i<CleanMuons->GetEntries(); i++) {
          const Muon *mu =  CleanMuons->At(i);
	  if(leptons->At(j) == mu){
	    distVtx0[j] = mu->BestTrk()->DzCorrected(*fVertices->At(0));
	    for(UInt_t nv=0; nv<fVertices->GetEntries(); nv++){
	      double dz = TMath::Abs(mu->BestTrk()->DzCorrected(*fVertices->At(nv)));
	      if(dz < distVtx[j]) {
	        closestVtx[j] = nv;
		distVtx[j] = dz;
              }
	    }
	    break;
	  }
	}
	for (UInt_t i=0; i<CleanElectrons->GetEntries(); i++) {
          const Electron *el =  CleanElectrons->At(i);
          if(leptons->At(j) == el){
	    distVtx0[j] = el->BestTrk()->DzCorrected(*fVertices->At(0));
	    for(UInt_t nv=0; nv<fVertices->GetEntries(); nv++){
	      double dz = TMath::Abs(el->BestTrk()->DzCorrected(*fVertices->At(nv)));
	      if(dz < distVtx[j]) {
	        closestVtx[j] = nv;
		distVtx[j] = dz;
              }
	    }
	    break;
	  }
	}
      }
      if     (closestVtx[0] == 0 && closestVtx[1] == 0) hDwwPresel[15]->Fill(0.0,NNLOWeight->GetVal());
      else if(closestVtx[0]      == closestVtx[1]     ) hDwwPresel[15]->Fill(1.0,NNLOWeight->GetVal());
      else if(closestVtx[0] >= 1 && closestVtx[1] == 0) hDwwPresel[15]->Fill(2.0,NNLOWeight->GetVal());
      else if(closestVtx[0] == 0 && closestVtx[1] >= 1) hDwwPresel[15]->Fill(3.0,NNLOWeight->GetVal());
      else                                              hDwwPresel[15]->Fill(4.0,NNLOWeight->GetVal());
      if(closestVtx[0] == closestVtx[1]){
        hDwwPresel[16]->Fill(TMath::Min(TMath::Max(distVtx[0],distVtx[1]),0.999),NNLOWeight->GetVal());      
        hDwwPresel[17]->Fill(TMath::Min(TMath::Abs(distVtx0[0]-distVtx0[1]),0.999),NNLOWeight->GetVal());      
      }
      else {
        hDwwPresel[18]->Fill(TMath::Min(TMath::Max(distVtx[0],distVtx[1]),0.999),NNLOWeight->GetVal());      
        hDwwPresel[19]->Fill(TMath::Min(TMath::Abs(distVtx0[0]-distVtx0[1]),0.999),NNLOWeight->GetVal());      
      }
      hDwwPresel[20]->Fill(0.0,NNLOWeight->GetVal());
      if(closestVtx[0] == closestVtx[1]) hDwwPresel[20]->Fill(1.0,NNLOWeight->GetVal());
      if(TMath::Max(distVtx[0],distVtx[1]) < 0.2 && TMath::Abs(distVtx0[0]-distVtx0[1]) < 0.2) hDwwPresel[20]->Fill(2.0,NNLOWeight->GetVal());
      if(closestVtx[0] == closestVtx[1] && 
         TMath::Max(distVtx[0],distVtx[1]) < 0.2 && TMath::Abs(distVtx0[0]-distVtx0[1]) < 0.2) hDwwPresel[20]->Fill(3.0,NNLOWeight->GetVal());
      if(TMath::Max(distVtx[0],distVtx[1]) < 0.1 && TMath::Abs(distVtx0[0]-distVtx0[1]) < 0.1) hDwwPresel[20]->Fill(4.0,NNLOWeight->GetVal());
      if(closestVtx[0] == closestVtx[1] && 
         TMath::Max(distVtx[0],distVtx[1]) < 0.1 && TMath::Abs(distVtx0[0]-distVtx0[1]) < 0.1) hDwwPresel[20]->Fill(5.0,NNLOWeight->GetVal());
    } // End Zll studies

    // Preselection level
    Bool_t PreselPtCut = kTRUE;
    if(leptons->At(0)->Pt() <= 20) PreselPtCut = kFALSE;
    if(leptons->At(1)->Pt() <= 10) PreselPtCut = kFALSE;
    //if(leptons->At(1)->ObjType() == kElectron && leptons->At(1)->Pt() <= 15) PreselPtCut = kFALSE;
    
    hDwwPresel[21]->Fill(TMath::Min(zDiffMax,0.999),NNLOWeight->GetVal());
    if(PreselPtCut == kTRUE &&
       dilepton.Charge() == 0 && zDiffMax < 10000.0 && PFMetStd->Pt() > 20){

      if(fPrintDebug == kTRUE){
	printf("EVT: %d %d %d %d %d %f %f %f %f %d %f %d %f\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
	                                        pairType,leptons->GetEntries(),newMET->Pt(),
	                                        dilepton.Mass(),TMath::Abs(dilepton.Mass()-91.1876),
	                                        METdeltaPhilEt,nCentralJets,etjetmax,
						DirtyMuons->GetEntries(),zDiffMax);
      }
      hDwwSel[22+100*pairType]->Fill(TMath::Min(zDiffMax,4.999),NNLOWeight->GetVal());
      if(dilepton.Mass() > 12 && (pairType >= 2 || dilepton.Mass() > 20)){
        hDwwSel[23+100*pairType]->Fill(TMath::Max(TMath::Min(0.0,0.999),0.000),NNLOWeight->GetVal());
        hDwwSel[24+100*pairType]->Fill(TMath::Max(TMath::Min(0.0,0.999),0.000),NNLOWeight->GetVal());
        hDwwXS[pairType]->Fill(8.0,NNLOWeight->GetVal());
        hDwwXS[4]       ->Fill(8.0,NNLOWeight->GetVal());
        if(fPrintDebug == kTRUE) {
	  if     (pairType == 0) printf("TYPE: %d %d %d 8 0 0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	  else if(pairType == 1) printf("TYPE: %d %d %d 0 8 0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	  else if(pairType == 2) printf("TYPE: %d %d %d 0 0 8\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
        }
	if(pairType >= 2 || TMath::Abs(dilepton.Mass()-91.1876) > 15){
	  if(fIsData == kTRUE){
	      printf("EVTDATA: %d %d %d %d %d %f %f %f %f %d %f %d %f %f %f %f\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
					      pairType,leptons->GetEntries(),newMET->Pt(),
					      dilepton.Mass(),TMath::Abs(dilepton.Mass()-91.1876),
					      METdeltaPhilEt,nCentralJets,etjetmax,DirtyMuons->GetEntries(),
	    				      leptons->At(0)->Pt(),leptons->At(0)->Eta(),
	    				      leptons->At(1)->Pt(),leptons->At(1)->Eta());
	  }
          hDwwXS[pairType]->Fill(9.0,NNLOWeight->GetVal());
          hDwwXS[4]       ->Fill(9.0,NNLOWeight->GetVal());
          if(fPrintDebug == kTRUE) {
	    if     (pairType == 0) printf("TYPE: %d %d %d 9 0 0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	    else if(pairType == 1) printf("TYPE: %d %d %d 0 9 0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	    else if(pairType == 2) printf("TYPE: %d %d %d 0 0 9\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
          }
          // Delta phi between the 2 leptons in degrees
          double deltaPhiLeptons = fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), 
    	  					       leptons->At(1)->Phi()))* 180./TMath::Pi();

          double deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(newMET->Phi(), 
      						           dilepton.Phi()))* 180./TMath::Pi();

          if(nCentralJets == 0){
            hDwwSel[14+100*0]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
            hDwwSel[16+100*0]->Fill(TMath::Min(newMET->Pt(),199.999),NNLOWeight->GetVal());
          }
          if(METdeltaPhilEt > 20 &&
	    (pairType >= 2 || METdeltaPhilEt > 37.0 + fVertices->GetEntries()/2.0)){
	    Double_t Beta = 1.0;
	    Double_t ptjetmaxBeta = 0.0;
            for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
              for(UInt_t nj=0; nj<fPFJet0->GetEntries(); nj++){
                const PFJet *pfjet = fPFJet0->At(nj);
    	        Double_t deltaR = MathUtils::DeltaR(CleanJetsNoPtCut->At(i)->Mom(),pfjet->Mom());
		if(deltaR < 0.05){
	  	  Beta = JetTools::Beta(pfjet, fVertices->At(0), 0.2);
		  if(ptjetmaxBeta < CleanJetsNoPtCut->At(i)->Pt()*Beta) {
		    ptjetmaxBeta = CleanJetsNoPtCut->At(i)->Pt()*Beta;
		  }
		  break;
		}
	      } // PFjets loop
	    } // jets loop
            if(DirtyMuons->GetEntries() == 0 && leptons->GetEntries() == 2){
	      hDwwSel[4+0*pairType]->Fill(TMath::Min(ptjetmaxBeta,199.999),NNLOWeight->GetVal());
              ObjArray<Jet> *JetDummy = new ObjArray<Jet>;
              ParticleOArr *PFCandidateSet = new ParticleOArr;
	      for(int i = 0;i<int(fPFCandidates->GetEntries());i++){
		if (MathUtils::DeltaR(fPFCandidates->At(i)->Mom(),leptons->At(0)->Mom()) > 0.05 && 
		    MathUtils::DeltaR(fPFCandidates->At(i)->Mom(),leptons->At(1)->Mom()) > 0.05){
                  PFCandidateSet->Add(fPFCandidates->At(i));
	        }
              }
      	      double etjetmax  = 0.0;
      	      if(CleanJetsNoPtCut->GetEntries() > 0) etjetmax  = CleanJetsNoPtCut->At(0)->Pt();
              double varJ[6] = {JetTools::NJettiness(CleanTracks, JetDummy), JetTools::NJettiness(CleanTracks, JetDummy, kTRUE),
	                        JetTools::NJettiness(CleanJetsNoPtCut, JetDummy), JetTools::NJettiness(CleanJetsNoPtCut, JetDummy, kTRUE),
	                        JetTools::NJettiness(PFCandidateSet, JetDummy), JetTools::NJettiness(PFCandidateSet, JetDummy, kTRUE)};
	      delete JetDummy;
	      delete PFCandidateSet;
              hDwwSel[25+0*pairType]->Fill(TMath::Min(varJ[0],199.999),NNLOWeight->GetVal());
              hDwwSel[26+0*pairType]->Fill(TMath::Min(varJ[1],199.999),NNLOWeight->GetVal());
              hDwwSel[27+0*pairType]->Fill(TMath::Min(varJ[2],199.999),NNLOWeight->GetVal());
              hDwwSel[28+0*pairType]->Fill(TMath::Min(varJ[3],4.999),NNLOWeight->GetVal());
              hDwwSel[29+0*pairType]->Fill(TMath::Min(varJ[4],199.999),NNLOWeight->GetVal());
              hDwwSel[30+0*pairType]->Fill(TMath::Min(varJ[5],199.999),NNLOWeight->GetVal());
              hDwwSel[31+0*pairType]->Fill(TMath::Min(etjetmax,199.999),NNLOWeight->GetVal());
              hDwwSel[32+0*pairType]->Fill(TMath::Min((double)nCentralJets,9.499),NNLOWeight->GetVal());
	    }

            hDwwSel[19+100*pairType]->Fill(180.-deltaPhiDileptonMet,NNLOWeight->GetVal());
            hDwwSel[20+100*pairType]->Fill(TMath::Min(newMET->Pt()/dilepton.Pt(),3.99),NNLOWeight->GetVal());
            hDwwXS[pairType]->Fill(10.0,NNLOWeight->GetVal());
            hDwwXS[4]       ->Fill(10.0,NNLOWeight->GetVal());
            if(fPrintDebug == kTRUE) {
	      if     (pairType == 0) printf("TYPE: %d %d %d 10  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	      else if(pairType == 1) printf("TYPE: %d %d %d  0 10  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	      else if(pairType >= 2) printf("TYPE: %d %d %d  0  0 10\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
            }
            hDwwSelAlphaEP1->Fill(TMath::Min(newMET->Pt()/dilepton.Pt(),3.99),180.-deltaPhiDileptonMet);
            // b-tagging study
	    double maxBtag[7] = {-99999., -99999., -99999., -99999., -99999., -99999., -99999.};
	    double imaxBtag[7] = {-1, -1, -1, -1, -1, -1, -1};
            for(UInt_t i=0; i<sortedJetsLowPt.size(); i++){
      	      double dZAverageJetPt = 0.0;
      	      double sumJetPt = 0.0;
      	      double jetPt = 0.0;
              for(UInt_t iPF=0; iPF<fPFJet0->GetEntries(); iPF++){							
      	        const PFJet *jet = fPFJet0->At(iPF);									
      	        if(MathUtils::DeltaR(jet->Mom(),sortedJetsLowPt[i]->Mom()) < 0.01){
      	          jetPt = jet->Pt();
      	          for (UInt_t npf=0; npf<jet->NPFCands();npf++) {
      	            const PFCandidate *pf = jet->PFCand(npf);
      	            if(pf->BestTrk()) {
      	              dZAverageJetPt = dZAverageJetPt + pf->Pt()*pf->Pt()*pf->BestTrk()->DzCorrected(*fVertices->At(0));
      	              sumJetPt = sumJetPt + pf->Pt()*pf->Pt();
      	            }
      	          }
      	          if(sumJetPt > 0) dZAverageJetPt = TMath::Abs(dZAverageJetPt)/sumJetPt;
      	          break;
      	        }
      	      } // loop over PF jets
      	      if(dZAverageJetPt < 2.0 && jetPt > 10){
	      	if(sortedJetsLowPt[i]->CombinedSecondaryVertexBJetTagsDisc() > maxBtag[0]){
	      	  maxBtag[0]  = sortedJetsLowPt[i]->CombinedSecondaryVertexBJetTagsDisc();
	      	  imaxBtag[0] = i;
	      	}
	      	if(sortedJetsLowPt[i]->CombinedSecondaryVertexMVABJetTagsDisc() > maxBtag[1]){
	      	  maxBtag[1]  = sortedJetsLowPt[i]->CombinedSecondaryVertexMVABJetTagsDisc();
	      	  imaxBtag[1] = i;
	      	}
	      	if(sortedJetsLowPt[i]->JetProbabilityBJetTagsDisc() > maxBtag[2]){
	      	  maxBtag[2]  = sortedJetsLowPt[i]->JetProbabilityBJetTagsDisc();
	      	  imaxBtag[2] = i;
	      	}
	      	if(sortedJetsLowPt[i]->JetBProbabilityBJetTagsDisc() > maxBtag[3]){
	      	  maxBtag[3]  = sortedJetsLowPt[i]->JetBProbabilityBJetTagsDisc();
	      	  imaxBtag[3] = i;
	      	}
	      	if(sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc() > maxBtag[4]){
	      	  maxBtag[4]  = sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc();
	      	  imaxBtag[4] = i;
	      	}
	      	if(sortedJetsLowPt[i]->TrackCountingHighPurBJetTagsDisc() > maxBtag[5]){
	      	  maxBtag[5]  = sortedJetsLowPt[i]->TrackCountingHighPurBJetTagsDisc();
	      	  imaxBtag[5] = i;
	      	}
	      	if(sortedJetsLowPt[i]->SimpleSecondaryVertexBJetTagsDisc() > maxBtag[6]){
	      	  maxBtag[6]  = sortedJetsLowPt[i]->SimpleSecondaryVertexBJetTagsDisc();
	      	  imaxBtag[6] = i;
	      	}
	      }
	    }
	    if(leptons->GetEntries() == 2){
	      if(imaxBtag[0] != -1){
		hDwwBTag[ 0+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[0],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 1+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[1],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 2+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[2],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 3+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[3],0.00001),3.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 4+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[4],-4.99999),14.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 5+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[5],-4.99999),14.9999),NNLOWeight->GetVal());
	        hDwwBTag[ 6+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[6],0.00001),4.9999),NNLOWeight->GetVal());
	      }
	      else {
	        hDwwBTag[ 0+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[ 1+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[ 2+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[ 3+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[ 4+30*TMath::Min(nCentralJets,2)]->Fill(-4.9999,NNLOWeight->GetVal());
	        hDwwBTag[ 5+30*TMath::Min(nCentralJets,2)]->Fill(-4.9999,NNLOWeight->GetVal());
	        hDwwBTag[ 6+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	      }
	      ObjArray<Muon> *DirtyMuonsNoJet = new ObjArray<Muon>;
	      for(UInt_t i=0; i<DirtyMuons->GetEntries(); i++){
	        double DeltaRMin = 999.;
		for(UInt_t j=0; j<sortedJets.size(); j++){
	          if(MathUtils::DeltaR(DirtyMuons->At(i)->Mom(),sortedJets[j]->Mom()) < DeltaRMin)
		    DeltaRMin = MathUtils::DeltaR(DirtyMuons->At(i)->Mom(),sortedJets[j]->Mom());
		}
	        hDwwBTag[7+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(DeltaRMin,5.999),NNLOWeight->GetVal());
		if(DeltaRMin > 0.3) DirtyMuonsNoJet->Add(DirtyMuons->At(i));
	      }
	      hDwwBTag[8+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min((double)DirtyMuons->GetEntries(),4.4999),NNLOWeight->GetVal());
	      hDwwBTag[9+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min((double)DirtyMuonsNoJet->GetEntries(),4.4999),NNLOWeight->GetVal());
	      if(DirtyMuonsNoJet->GetEntries() > 0){
	        hDwwBTag[10+30*TMath::Min(nCentralJets,2)]->Fill(0.9999,NNLOWeight->GetVal());
	        hDwwBTag[11+30*TMath::Min(nCentralJets,2)]->Fill(0.9999,NNLOWeight->GetVal());
	        hDwwBTag[12+30*TMath::Min(nCentralJets,2)]->Fill(0.9999,NNLOWeight->GetVal());
	        hDwwBTag[13+30*TMath::Min(nCentralJets,2)]->Fill(3.9999,NNLOWeight->GetVal());
	        hDwwBTag[14+30*TMath::Min(nCentralJets,2)]->Fill(14.9999,NNLOWeight->GetVal());
	        hDwwBTag[15+30*TMath::Min(nCentralJets,2)]->Fill(14.9999,NNLOWeight->GetVal());	      
	        hDwwBTag[16+30*TMath::Min(nCentralJets,2)]->Fill(4.9999,NNLOWeight->GetVal());
	      }
	      else if(imaxBtag[0] != -1){
	        hDwwBTag[10+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[0],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[11+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[1],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[12+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[2],0.00001),0.9999),NNLOWeight->GetVal());
	        hDwwBTag[13+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[3],0.00001),3.9999),NNLOWeight->GetVal());
	        hDwwBTag[14+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[4],-4.99999),14.9999),NNLOWeight->GetVal());
	        hDwwBTag[15+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[5],-4.99999),14.9999),NNLOWeight->GetVal());
	        hDwwBTag[16+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Max(maxBtag[6],0.00001),4.9999),NNLOWeight->GetVal());
	      }
	      else {
	        hDwwBTag[10+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[11+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[12+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[13+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	        hDwwBTag[14+30*TMath::Min(nCentralJets,2)]->Fill(-4.9999,NNLOWeight->GetVal());
	        hDwwBTag[15+30*TMath::Min(nCentralJets,2)]->Fill(-4.9999,NNLOWeight->GetVal());
	        hDwwBTag[16+30*TMath::Min(nCentralJets,2)]->Fill(0.0,NNLOWeight->GetVal());
	      }

	      hDwwBTag[20]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());
	      if(maxBtag[4] > 2.1){
	        hDwwBTag[17+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(sortedJetsLowPt[imaxBtag[4]]->Pt(),29.9999),NNLOWeight->GetVal());
		for(UInt_t i=0; i<fPFJet0->GetEntries(); i++){
		  const PFJet *jet = fPFJet0->At(i); 
		  if(MathUtils::DeltaR(jet->Mom(),sortedJetsLowPt[imaxBtag[4]]->Mom()) < 0.01){
                    double dZAverageJet = 0.0; double dZAverageJetPt = 0.0;
		    int nChaJet = 0; double sumPt = 0.0;
		    for (UInt_t npf=0; npf<jet->NPFCands();npf++) {   
                      const PFCandidate *pf = jet->PFCand(npf);
		      if(pf->BestTrk()) {
		        dZAverageJet   = dZAverageJet   +                   pf->BestTrk()->DzCorrected(*fVertices->At(0));
			nChaJet++;
		        dZAverageJetPt = dZAverageJetPt + pf->Pt()*pf->Pt()*pf->BestTrk()->DzCorrected(*fVertices->At(0));
			sumPt = sumPt + pf->Pt()*pf->Pt();
		      }
		    }
		    if(nChaJet > 0) dZAverageJet   = dZAverageJet  /nChaJet;
		    if(sumPt   > 0) dZAverageJetPt = dZAverageJetPt/sumPt;
		    hDwwBTag[18+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Abs(dZAverageJet),9.999),NNLOWeight->GetVal());
		    hDwwBTag[19+30*TMath::Min(nCentralJets,2)]->Fill(TMath::Min(TMath::Abs(dZAverageJetPt),9.999),NNLOWeight->GetVal());
		    break;
		  }
		}
	        hDwwBTag[21]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());
		if(TMath::Abs(sortedJetsLowPt[imaxBtag[4]]->MatchedMCFlavor()) == 5){
		  hDwwBTag[22]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());
		}
	      }
	      if(DirtyMuonsNoJet->GetEntries() > 0){
	        hDwwBTag[23]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());	        
	      }
	      if(maxBtag[4] > 2.1 || DirtyMuonsNoJet->GetEntries() > 0){
	        hDwwBTag[24]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());	        
	      }
	      if(maxBtag[4] > 2.1 || DirtyMuonsNoJet->GetEntries() > 0 || maxBtag[6] > 0.0){
	        hDwwBTag[25]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());	        
	      }
	      if(DirtyMuons->GetEntries() > 0){
	        hDwwBTag[26]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());	        
	      }
	      if(maxBtag[4] <= 2.1 && DirtyMuonsNoJet->GetEntries() == 0){
	        hDwwBTag[27]->Fill(TMath::Min((double)nCentralJets,5.499),NNLOWeight->GetVal());
		if(nCentralJets == 0){
		  hDwwBTag[28]->Fill((double)pairType,NNLOWeight->GetVal());
		}
	        if(nCentralJets == 1){
                  double deltaPhiLLJet  = fabs(MathUtils::DeltaPhi(dilepton.Phi(), sortedJets[0]->Phi()));
		  hDwwBTag[50]->Fill(deltaPhiLLJet*180./TMath::Pi(),NNLOWeight->GetVal());
	          if(deltaPhiLLJet*180./TMath::Pi() < 160.0){
		    hDwwBTag[51]->Fill(TMath::Min(TMath::Max(sortedJets[0]->CombinedSecondaryVertexBJetTagsDisc()   ,0.00001),0.9999),NNLOWeight->GetVal());
	            hDwwBTag[52]->Fill(TMath::Min(TMath::Max(sortedJets[0]->CombinedSecondaryVertexMVABJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
	            hDwwBTag[53]->Fill(TMath::Min(TMath::Max(sortedJets[0]->JetProbabilityBJetTagsDisc()            ,0.00001),0.9999),NNLOWeight->GetVal());
	            hDwwBTag[54]->Fill(TMath::Min(TMath::Max(sortedJets[0]->JetBProbabilityBJetTagsDisc()           ,0.00001),3.9999),NNLOWeight->GetVal());
	            hDwwBTag[55]->Fill(TMath::Min(TMath::Max(sortedJets[0]->TrackCountingHighEffBJetTagsDisc()      ,-4.99999),14.9999),NNLOWeight->GetVal());
	            hDwwBTag[56]->Fill(TMath::Min(TMath::Max(sortedJets[0]->TrackCountingHighPurBJetTagsDisc()      ,-4.99999),14.9999),NNLOWeight->GetVal());
	            hDwwBTag[57]->Fill(TMath::Min(TMath::Max(sortedJets[0]->SimpleSecondaryVertexBJetTagsDisc()     ,0.00001),4.9999),NNLOWeight->GetVal());
	          }
		}
	      }
              delete DirtyMuonsNoJet;
	      double NBjets[2] = {0, 0};
	      for(UInt_t i=0; i<sortedJetsAll.size(); i++){
                if(TMath::Abs(sortedJetsAll[i]->MatchedMCFlavor()) == 5                                    ) NBjets[0]++;
                if(TMath::Abs(sortedJetsAll[i]->MatchedMCFlavor()) == 5 && sortedJetsAll[i]->AbsEta() < 2.5) NBjets[1]++;
	      }
	      hDwwBTag[58]->Fill(TMath::Min(nCentralJets,4)*10+TMath::Min(NBjets[0],7.499),NNLOWeight->GetVal());
	      hDwwBTag[59]->Fill(TMath::Min(nCentralJets,4)*10+TMath::Min(NBjets[1],7.499),NNLOWeight->GetVal());
	    } // nleptons == 2

	    if(nCentralJets == 0){
	    //if(ptjetmaxBeta <= fPtJetCut){
              hDwwXS[pairType]->Fill(11.0,NNLOWeight->GetVal());
              hDwwXS[4]       ->Fill(11.0,NNLOWeight->GetVal());
              if(fPrintDebug == kTRUE) {
	        if     (pairType == 0) printf("TYPE: %d %d %d 11  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        else if(pairType == 1) printf("TYPE: %d %d %d  0 11  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        else if(pairType == 2) printf("TYPE: %d %d %d  0  0 11\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        if     (pairType == 0) printf("TYPE: %d %d %d 11  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        printf("BEFORESOFTMUON: %d %d %d | M: ",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
                for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
                  const Muon *mu = fMuons->At(i);
		  printf("m(%d): %f %f %f - ",i,mu->Pt(),mu->Eta(),mu->Phi());
		}
	        printf(" | E: ");
                for (UInt_t i=0; i<fElectrons->GetEntries(); i++) {
                  const Electron *el = fElectrons->At(i);
		  printf("e(%d): %f %f %f - ",i,el->Pt(),el->Eta(),el->Phi());
		}
	        printf(" | J: ");
                for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
                const Jet *jet = CleanJetsNoPtCut->At(i);
		  printf("j(%d): %f %f %f %f %f - ",i,jet->Pt(),jet->Eta(),jet->Phi(),jet->RawMom().Pt(),jet->TrackCountingHighEffBJetTagsDisc());
		}
	        printf("\n");
              }
              double deltaPhiLLJet = 0.0;
              if(sortedJetsAll.size() > 0 && sortedJetsAll[0]->Pt() > 15.0 && (pairType == 0 || pairType == 1)){
	        deltaPhiLLJet = fabs(MathUtils::DeltaPhi(dilepton.Phi(), sortedJetsAll[0]->Phi()))*180.0/TMath::Pi();
              }
	      if(deltaPhiLLJet < 165.0){
                hDwwXS[pairType]->Fill(12.0,NNLOWeight->GetVal());
                hDwwXS[4]	->Fill(12.0,NNLOWeight->GetVal());
		if(DirtyMuons->GetEntries() == 0){
                  hDwwXS[pairType]->Fill(13.0,NNLOWeight->GetVal());
                  hDwwXS[4]       ->Fill(13.0,NNLOWeight->GetVal());
                  if(fPrintDebug == kTRUE) {
	            if     (pairType == 0) printf("TYPE: %d %d %d 12  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	            else if(pairType == 1) printf("TYPE: %d %d %d  0 12  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	            else if(pairType == 2) printf("TYPE: %d %d %d  0  0 12\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
                  }
                  if(leptons->GetEntries() == 2){
                    hDwwXS[pairType]->Fill(14.0,NNLOWeight->GetVal());
                    hDwwXS[4]       ->Fill(14.0,NNLOWeight->GetVal());
                    if(fPrintDebug == kTRUE) {
	              if     (pairType == 0) printf("TYPE: %d %d %d 13  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	              else if(pairType == 1) printf("TYPE: %d %d %d  0 13  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	              else if(pairType == 2) printf("TYPE: %d %d %d  0  0 13\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
                    }
                    if(maxBtag[4] < 2.1){
	      	      for(UInt_t nj=0; nj<fPFJet0->GetEntries(); nj++){
			const PFJet *jet = fPFJet0->At(nj);      

			Bool_t isElectronOverlap = kFALSE;
			if (CleanElectrons) {
    		    	  UInt_t n = CleanElectrons->GetEntries();
    		    	  for (UInt_t j=0; j<n; j++) {
    		    	    Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(),
              	    					        CleanElectrons->At(j)->SCluster()->Eta(), 
              	    					        jet->Phi(), jet->Eta());
		    	    hDwwOverlap[0]->Fill(TMath::Min(deltaR,1.999),NNLOWeight->GetVal());
    		    	    if (deltaR < 0.3) {
		    	      hDwwOverlap[1]->Fill(TMath::Min((double)jet->ChargedMultiplicity(),9.499),NNLOWeight->GetVal());
			      if(jet->ChargedMultiplicity() >= 3){
		    		hDwwOverlap[2]->Fill(TMath::Min(CleanElectrons->At(j)->E()/(jet->ChargedHadronEnergy()+jet->ChargedEmEnergy()+jet->MuonEnergy()),0.999),NNLOWeight->GetVal());
		    		hDwwOverlap[3]->Fill(TMath::Min(jet->ChargedHadronEnergy()/CleanElectrons->At(j)->E(),0.999),NNLOWeight->GetVal());
			      }
    		    	      isElectronOverlap = kTRUE;
    		    	      break;
    		    	    }
    		    	  }
			}
			if (isElectronOverlap) continue;
			Bool_t isMuonOverlap = kFALSE;
			if (CleanMuons) {
    		    	  UInt_t n = CleanMuons->GetEntries();
    		    	  for (UInt_t j=0; j<n; j++) {
    		    	    Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());
		    	    hDwwOverlap[4]->Fill(TMath::Min(deltaR,1.999),NNLOWeight->GetVal());
    		    	    if (deltaR < 0.3) {
		    	      hDwwOverlap[5]->Fill(TMath::Min((double)jet->ChargedMultiplicity(),9.499),NNLOWeight->GetVal());
			      if(jet->ChargedMultiplicity() >= 3){
		    		hDwwOverlap[6]->Fill(TMath::Min(CleanMuons->At(j)->E()/(jet->ChargedHadronEnergy()+jet->ChargedEmEnergy()+jet->MuonEnergy()),0.999),NNLOWeight->GetVal());
		    		hDwwOverlap[7]->Fill(TMath::Min(jet->ChargedHadronEnergy()/CleanMuons->At(j)->E(),0.999),NNLOWeight->GetVal());
		    	      }
    		    	      isMuonOverlap = kTRUE;
    		    	      break;	       
    		    	    }
    		    	  }
			}
			if (isMuonOverlap) continue;
	      	      }
                      hDwwXS[pairType]->Fill(15.0,NNLOWeight->GetVal());
                      hDwwXS[4]       ->Fill(15.0,NNLOWeight->GetVal());
		      // H(160)->WW->2l2n selection
                      if(leptons->At(0)->Pt() > 30.0) {
                        hDwwXS[pairType]->Fill(16.0,NNLOWeight->GetVal());
                        hDwwXS[4]       ->Fill(16.0,NNLOWeight->GetVal());
                        if(leptons->At(1)->Pt() > 25.0) {
                          hDwwXS[pairType]->Fill(17.0,NNLOWeight->GetVal());
                          hDwwXS[4]       ->Fill(17.0,NNLOWeight->GetVal());
			  if(dilepton.Mass() < 50.0) {
                            hDwwXS[pairType]->Fill(18.0,NNLOWeight->GetVal());
                            hDwwXS[4]       ->Fill(18.0,NNLOWeight->GetVal());
			    double metFraction[2] = {1.0,0.0};
			    double mtHiggs = JetTools::MtHiggs(leptons,newMET,metFraction,7);
			    if(mtHiggs > 90.0 && mtHiggs < 160.0) {
                              hDwwXS[pairType]->Fill(19.0,NNLOWeight->GetVal());
                              hDwwXS[4]       ->Fill(19.0,NNLOWeight->GetVal());
			      if(deltaPhiLeptons < 60.0) {
                                hDwwXS[pairType]->Fill(20.0,NNLOWeight->GetVal());
                                hDwwXS[4]       ->Fill(20.0,NNLOWeight->GetVal());
			        if(dilepton.Pt() > 45) {
                                  hDwwXS[pairType]->Fill(21.0,NNLOWeight->GetVal());
                                  hDwwXS[4]       ->Fill(21.0,NNLOWeight->GetVal());
			        } // dilepton pt cut
			      } // deltaPhill cut
			    } // mt cut
			  } // massll cut
			} // ptmin cut
	              } // ptmax cut
                      if(fPrintDebug == kTRUE) {
	        	if     (pairType == 0) printf("TYPE: %d %d %d 14  0  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        	else if(pairType == 1) printf("TYPE: %d %d %d  0 14  0\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
	        	else if(pairType == 2) printf("TYPE: %d %d %d  0  0 14\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec());
                      }
                      /*
                      if(pairType==0) {
			printf("********************************************************************************\n");
			printf("EVT: %d %d %d %d %d %f %f %f %f %d %f %d %f\n",fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
		        					pairType,leptons->GetEntries(),newMET->Pt(),
		        					dilepton.Mass(),TMath::Abs(dilepton.Mass()-91.1876),
		        					METdeltaPhilEt,nCentralJets,etjetmax,
		        					DirtyMuons->GetEntries(),zDiffMax);
			printf("LEP: %f %f %f | %f %f %f\n",leptons->At(0)->Pt(),leptons->At(0)->Eta(),leptons->At(0)->Phi(),
		        				    leptons->At(1)->Pt(),leptons->At(1)->Eta(),leptons->At(1)->Phi());
    			for (UInt_t ip=0; ip<fParticles->GetEntries(); ip++) {
    		          const MCParticle *p = fParticles->At(ip);
                          p->Print("l");
			}
			printf("********************************************************************************\n");
		      }
		      */
        	      hDwwFake[6+100*pairType]->Fill(leptonGenType[0],NNLOWeight->GetVal());
        	      hDwwFake[7+100*pairType]->Fill(leptonGenType[1],NNLOWeight->GetVal());
        	      hDwwFake[8+100*pairType]->Fill(leptonGenType[0]+5*leptonGenType[1],NNLOWeight->GetVal());

		      hDwwSel[ 0+  0*pairType]->Fill(deltaPhiMetLepton[0] * 180./TMath::Pi(),NNLOWeight->GetVal());
		      hDwwSel[ 1+  0*pairType]->Fill(deltaPhiMetLepton[1] * 180./TMath::Pi(),NNLOWeight->GetVal());
		      hDwwSel[ 2+  0*pairType]->Fill(TMath::Min(deltaPhiMetLepton[0],deltaPhiMetLepton[1]) * 180./TMath::Pi(),NNLOWeight->GetVal());
		      hDwwSel[ 3+  0*pairType]->Fill(TMath::Max(deltaPhiMetLepton[0],deltaPhiMetLepton[1]) * 180./TMath::Pi(),NNLOWeight->GetVal());
		      hDwwSel[ 9+100*pairType]->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
        	      hDwwSel[10+100*pairType]->Fill(TMath::Min(dilepton.Mt(),399.999),NNLOWeight->GetVal());
        	      hDwwSel[11+100*pairType]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
        	      hDwwSel[12+100*pairType]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
        	      hDwwSel[13+100*pairType]->Fill(minDeltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
        	      hDwwSel[15+100*pairType]->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
        	      hDwwSel[17+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        	      hDwwSel[18+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
        	      DiTauSystem ditau(leptons->At(0), leptons->At(1), newMET);
        	      if(ditau.XTau1() > 0 && ditau.XTau2() > 0){
			if(ditau.RecoMass() <= 0) printf("RecoMass <= 0! == %f\n",ditau.RecoMass());
			hDwwSel[21+100*pairType]->Fill(TMath::Min(ditau.RecoMass(),199.99),NNLOWeight->GetVal());
		      }
        	      int metBins = 100;
		      double minMass = 9000.;
		      double nBin[2] = {-1, -1};
        	      for(int i=0; i<=metBins; i++){
        		for(int j=0; j<=metBins; j++){
			  double metAuxPx[2] = {newMET->Px() * i/metBins*1.0,
	                        		newMET->Px() * (metBins-i)/metBins*1.0};
			  double metAuxPy[2] = {newMET->Py() * j/metBins*1.0,
	                        		newMET->Py() * (metBins-j)/metBins*1.0};
			  double phiAux[2] = {atan2(metAuxPy[0],metAuxPx[0]),
	                        	      atan2(metAuxPy[1],metAuxPx[1])};
			  double ptAux[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]),
	                		     sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1])};
        		  double deltaPhiMetLepton[2] = {fabs(MathUtils::DeltaPhi(phiAux[0], leptons->At(0)->Phi())),
    				        		 fabs(MathUtils::DeltaPhi(phiAux[1], leptons->At(1)->Phi()))};
        		  double mTW[2] = {TMath::Sqrt(2.0*leptons->At(0)->Pt()*ptAux[0]*
    	  	  		        	      (1.0 - cos(deltaPhiMetLepton[0]))),
        	        		   TMath::Sqrt(2.0*leptons->At(1)->Pt()*ptAux[1]*
    				        	      (1.0 - cos(deltaPhiMetLepton[1])))};
			  if(minMass > TMath::Max(mTW[0],mTW[1])){
			    minMass = TMath::Max(mTW[0],mTW[1]);
			    nBin[0] = i;
			    nBin[1] = j;
			  }
        		}
		      }
		      double metAuxPx[2] = {newMET->Px() * nBin[0]/metBins*1.0,
					    newMET->Px() * (metBins-nBin[0])/metBins*1.0};
		      double metAuxPy[2] = {newMET->Py() * nBin[1]/metBins*1.0,
					    newMET->Py() * (metBins-nBin[1])/metBins*1.0};

        	      double T[2] = {80.40*80.40/2.0+leptons->At(0)->Px()*metAuxPx[0]+leptons->At(0)->Py()*metAuxPy[0],
				     80.40*80.40/2.0+leptons->At(1)->Px()*metAuxPx[1]+leptons->At(1)->Py()*metAuxPy[1]};
		      double S[2] = {metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0],
				     metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]};
		      double B[2] = {T[0]*T[0]*leptons->At(0)->Pz()*leptons->At(0)->Pz()-(leptons->At(0)->P()*leptons->At(0)->P()*S[0]-T[0]*T[0])*(leptons->At(0)->P()*leptons->At(0)->P()-leptons->At(0)->Pz()*leptons->At(0)->Pz()),
				     T[1]*T[1]*leptons->At(1)->Pz()*leptons->At(1)->Pz()-(leptons->At(1)->P()*leptons->At(1)->P()*S[1]-T[1]*T[1])*(leptons->At(1)->P()*leptons->At(1)->P()-leptons->At(1)->Pz()*leptons->At(1)->Pz())};
		      for(int i=0; i<2; i++) {if(B[i] > 0) B[i] = sqrt(B[i]); else B[i] = 0.0;}
		      double PzNeuP[2] = {(T[0]*leptons->At(0)->Pz()+B[0])/(leptons->At(0)->P()*leptons->At(0)->P()-leptons->At(0)->Pz()*leptons->At(0)->Pz()),
	                		  (T[1]*leptons->At(1)->Pz()+B[1])/(leptons->At(1)->P()*leptons->At(1)->P()-leptons->At(1)->Pz()*leptons->At(1)->Pz())};
		      double PzNeuN[2] = {(T[0]*leptons->At(0)->Pz()-B[0])/(leptons->At(0)->P()*leptons->At(0)->P()-leptons->At(0)->Pz()*leptons->At(0)->Pz()),
	                		  (T[1]*leptons->At(1)->Pz()-B[1])/(leptons->At(1)->P()*leptons->At(1)->P()-leptons->At(1)->Pz()*leptons->At(1)->Pz())};

		      double np[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]+PzNeuP[0]*PzNeuP[0]),
				      sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]+PzNeuP[1]*PzNeuP[1])};
		      double nn[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]+PzNeuN[0]*PzNeuN[0]),
				      sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]+PzNeuN[1]*PzNeuN[1])};
		      double massH[5];
        	      massH[0]= (np[0]+leptons->At(0)->P()+np[1]+leptons->At(1)->P())*(np[0]+leptons->At(0)->P()+np[1]+leptons->At(1)->P())
			     -(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())*(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())
			     -(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())*(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())
			     -(PzNeuP[0]  +leptons->At(0)->Pz()+PzNeuP[1]  +leptons->At(1)->Pz())*(PzNeuP[0]  +leptons->At(0)->Py()+PzNeuP[1]  +leptons->At(1)->Pz());
        	      massH[1]= (np[0]+leptons->At(0)->P()+nn[1]+leptons->At(1)->P())*(np[0]+leptons->At(0)->P()+nn[1]+leptons->At(1)->P())
			     -(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())*(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())
			     -(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())*(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())
			     -(PzNeuP[0]  +leptons->At(0)->Pz()+PzNeuN[1]  +leptons->At(1)->Pz())*(PzNeuP[0]  +leptons->At(0)->Py()+PzNeuN[1]  +leptons->At(1)->Pz());
        	      massH[2]= (nn[0]+leptons->At(0)->P()+np[1]+leptons->At(1)->P())*(nn[0]+leptons->At(0)->P()+np[1]+leptons->At(1)->P())
			     -(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())*(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())
			     -(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())*(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())
			     -(PzNeuN[0]  +leptons->At(0)->Pz()+PzNeuP[1]  +leptons->At(1)->Pz())*(PzNeuN[0]  +leptons->At(0)->Py()+PzNeuP[1]  +leptons->At(1)->Pz());
        	      massH[3]= (nn[0]+leptons->At(0)->P()+nn[1]+leptons->At(1)->P())*(nn[0]+leptons->At(0)->P()+nn[1]+leptons->At(1)->P())
			     -(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())*(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())
			     -(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())*(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())
			     -(PzNeuN[0]  +leptons->At(0)->Pz()+PzNeuN[1]  +leptons->At(1)->Pz())*(PzNeuN[0]  +leptons->At(0)->Py()+PzNeuN[1]  +leptons->At(1)->Pz());

        	      double metFraction[2] = {nBin[0]/metBins*1.0,nBin[1]/metBins*1.0};

		      // New version
		      double metAuxPz[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0])/leptons->At(0)->Pt()*leptons->At(0)->Pz(),
		                            sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1])/leptons->At(1)->Pt()*leptons->At(1)->Pz()};
                      double ENeu[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]+metAuxPz[0]*metAuxPz[0]),
		                	sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]+metAuxPz[1]*metAuxPz[1])};		  
        	      massH[4]= (ENeu[0]+leptons->At(0)->P()+ENeu[1]+leptons->At(1)->P())*(ENeu[0]+leptons->At(0)->P()+ENeu[1]+leptons->At(1)->P())
			     -(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())*(metAuxPx[0]+leptons->At(0)->Px()+metAuxPx[1]+leptons->At(1)->Px())
			     -(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())*(metAuxPy[0]+leptons->At(0)->Py()+metAuxPy[1]+leptons->At(1)->Py())
			     -(metAuxPz[0]+leptons->At(0)->Pz()+metAuxPz[1]+leptons->At(1)->Pz())*(metAuxPz[0]+leptons->At(0)->Py()+metAuxPz[1]+leptons->At(1)->Pz());

		      for(int i=0; i<5; i++) {if(massH[i] > 0) massH[i] = sqrt(massH[i]); else massH[i] = 0.0;}

        	      hDwwMassSel[ 0]->Fill(TMath::Min(minMass,199.999),NNLOWeight->GetVal());
        	      hDwwMassSel[ 1]->Fill(TMath::Min(massH[0],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[ 2]->Fill(TMath::Min(massH[1],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[ 3]->Fill(TMath::Min(massH[2],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[ 4]->Fill(TMath::Min(massH[3],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[ 5]->Fill(TMath::Min(massH[4],499.999),NNLOWeight->GetVal());
		      if(minMass > 50.0){
			hDwwMassSel[ 6]->Fill(TMath::Min(massH[0],499.999),NNLOWeight->GetVal());
			hDwwMassSel[ 7]->Fill(TMath::Min(massH[1],499.999),NNLOWeight->GetVal());
			hDwwMassSel[ 8]->Fill(TMath::Min(massH[2],499.999),NNLOWeight->GetVal());
			hDwwMassSel[ 9]->Fill(TMath::Min(massH[3],499.999),NNLOWeight->GetVal());
			hDwwMassSel[10]->Fill(TMath::Min(massH[4],499.999),NNLOWeight->GetVal());
		      }
		      double mthiggs[8] = {JetTools::MtHiggs(leptons,newMET,metFraction,0),
		                	   JetTools::MtHiggs(leptons,newMET,metFraction,1),
		                	   JetTools::MtHiggs(leptons,newMET,metFraction,2),
					   JetTools::MtHiggs(leptons,newMET,metFraction,3),
					   JetTools::MtHiggs(leptons,newMET,metFraction,4),
					   JetTools::MtHiggs(leptons,newMET,metFraction,5),
					   JetTools::MtHiggs(leptons,newMET,metFraction,6),
					   JetTools::MtHiggs(leptons,newMET,metFraction,7)};

        	      hDwwMassSel[11]->Fill(TMath::Min(mthiggs[0],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[12]->Fill(TMath::Min(mthiggs[1],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[13]->Fill(TMath::Min(mthiggs[2],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[14]->Fill(TMath::Min(mthiggs[3],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[15]->Fill(TMath::Min(mthiggs[4],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[16]->Fill(TMath::Min(mthiggs[5],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[17]->Fill(TMath::Min(mthiggs[6],499.999),NNLOWeight->GetVal());
        	      hDwwMassSel[18]->Fill(TMath::Min(mthiggs[7],499.999),NNLOWeight->GetVal());
		      if(dilepton.Mass() < 100.0 && deltaPhiLeptons < 120.0){
        		hDwwMassSel[ 0+20]->Fill(TMath::Min(minMass,199.999),NNLOWeight->GetVal());
        		hDwwMassSel[ 1+20]->Fill(TMath::Min(massH[0],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[ 2+20]->Fill(TMath::Min(massH[1],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[ 3+20]->Fill(TMath::Min(massH[2],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[ 4+20]->Fill(TMath::Min(massH[3],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[ 5+20]->Fill(TMath::Min(massH[4],499.999),NNLOWeight->GetVal());
			if(minMass > 50.0){
		          hDwwMassSel[ 6+20]->Fill(TMath::Min(massH[0],499.999),NNLOWeight->GetVal());
		          hDwwMassSel[ 7+20]->Fill(TMath::Min(massH[1],499.999),NNLOWeight->GetVal());
		          hDwwMassSel[ 8+20]->Fill(TMath::Min(massH[2],499.999),NNLOWeight->GetVal());
		          hDwwMassSel[ 9+20]->Fill(TMath::Min(massH[3],499.999),NNLOWeight->GetVal());
		          hDwwMassSel[10+20]->Fill(TMath::Min(massH[4],499.999),NNLOWeight->GetVal());
			}
        		hDwwMassSel[11+20]->Fill(TMath::Min(mthiggs[0],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[12+20]->Fill(TMath::Min(mthiggs[1],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[13+20]->Fill(TMath::Min(mthiggs[2],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[14+20]->Fill(TMath::Min(mthiggs[3],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[15+20]->Fill(TMath::Min(mthiggs[4],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[16+20]->Fill(TMath::Min(mthiggs[5],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[17+20]->Fill(TMath::Min(mthiggs[6],499.999),NNLOWeight->GetVal());
        		hDwwMassSel[18+20]->Fill(TMath::Min(mthiggs[7],499.999),NNLOWeight->GetVal());
		      }

                      double varX[5] = {JetTools::M_r(leptons), JetTools::Beta_r(leptons), JetTools::M_r_t(leptons, newMET),
		                	JetTools::Razor(leptons,newMET), JetTools::CosineOmega(leptons->At(0),leptons->At(1))};
                      hDwwSel[5+0*pairType]->Fill(TMath::Max(TMath::Min(varX[1],  1.999),-1.999),NNLOWeight->GetVal());
                      hDwwSel[6+0*pairType]->Fill(TMath::Max(TMath::Min(varX[0],499.999), 0.000),NNLOWeight->GetVal());
		      if(TMath::Abs(varX[1]) <= 1){
                	hDwwSel[7+0*pairType]->Fill(TMath::Max(TMath::Min(varX[0],499.999), 0.000),NNLOWeight->GetVal());
                	hDwwSel[8+0*pairType]->Fill(TMath::Max(TMath::Min(varX[2],499.999), 0.000),NNLOWeight->GetVal());
                	hDwwSel[33+0*pairType]->Fill(TMath::Max(TMath::Min(varX[3],  1.999), 0.000),NNLOWeight->GetVal());
                	hDwwSel[34+0*pairType]->Fill(TMath::Max(TMath::Min(varX[4],  0.999),-0.999),NNLOWeight->GetVal());
		      }
                    } // btagmax < 2.1
                  } // NLeptons == 2
		} // DirtyMuons == 0
	      } // deltaPhiLLJet cut
	    } // Njets == 0
	  } // projected MET > 20 (35) for emu (ee/em)
	} // |mz-mll|>15 for eemm
      } // mll > 12
    } // MET > 20 && ptl2 > 20 && q == 0 && |Z-Zvert|<1
    delete CleanTracks;
  } // Minimun Pt, Nleptons>=2 requirements
  delete newMET;
  delete DirtyMuons;
  for(UInt_t i=0; i<sortedJets.size();      i++) delete sortedJets[i];
  for(UInt_t i=0; i<sortedJetsAll.size();   i++) delete sortedJetsAll[i];
}
//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fAllVertexName, fAllVertices);
  ReqBranch(fEvtHdrName, fEventHeader);
  ReqBranch(fCaloJetName0, fCaloJet0);
  ReqBranch(fCaloJetName1, fCaloJet1);
  ReqBranch(fCaloJetName2, fCaloJet2);
  ReqBranch(fCaloJetName3, fCaloJet3);
  ReqBranch(fTrackJetName0,fTrackJet0);
  ReqBranch(fPFJetName0, fPFJet0);
  ReqBranch(fPFJetName1, fPFJet1);
  ReqBranch(fPFJetName2, fPFJet2);
  ReqBranch(fPFJetName3, fPFJet3);
  ReqBranch(fPFMetName,   fPFMetStd);
  ReqBranch(fTCMetName,   fTCMetStd);
  ReqBranch(fMuonName,    fMuons);
  ReqBranch(fElectronName,fElectrons);
  ReqBranch(fTrackName,   fTracks);
  ReqBranch(fPFCandidatesName, fPFCandidates);
  ReqBranch(fConversionBranchName, fConversions);
  ReqBranch(fBeamSpotName, fBeamSpot);
  ReqBranch(fPileupEnergyDensityName,fPileupEnergyDensity);

  if(fIsData == kFALSE){
    ReqEventObject(fMCPartName, fParticles, kTRUE);
    ReqBranch(fPileupInfoName,fPileupInfo);
  }

  char sb[200];
  sprintf(sb,"hDwwBTagCheck_%d",0);  hDwwBTagCheck[0]  = new TH1D(sb,sb,100,-2.0,18.0);
  sprintf(sb,"hDwwBTagCheck_%d",1);  hDwwBTagCheck[1]  = new TH1D(sb,sb,100,-2.0,18.0);

  for(int i=0; i<2; i++){
    AddOutput(hDwwBTagCheck[i]);
  }

  sprintf(sb,"hDwwXS_%d",0);  hDwwXS[0]  = new TH1D(sb,sb,22,-0.5,21.5);
  sprintf(sb,"hDwwXS_%d",1);  hDwwXS[1]  = new TH1D(sb,sb,22,-0.5,21.5);
  sprintf(sb,"hDwwXS_%d",2);  hDwwXS[2]  = new TH1D(sb,sb,22,-0.5,21.5);
  sprintf(sb,"hDwwXS_%d",3);  hDwwXS[3]  = new TH1D(sb,sb,22,-0.5,21.5);
  sprintf(sb,"hDwwXS_%d",4);  hDwwXS[4]  = new TH1D(sb,sb,22,-0.5,21.5);

  for(int i=0; i<5; i++){
    AddOutput(hDwwXS[i]);
  }

  sprintf(sb,"hDwwSkim_%d",0);  hDwwSkim[0]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwSkim_%d",1);  hDwwSkim[1]  = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDwwSkim_%d",2);  hDwwSkim[2]  = new TH1D(sb,sb,360,0.0,180.0); 

  for(int i=0; i<3; i++){
    AddOutput(hDwwSkim[i]);
  }

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDwwSel_%d",ind+0);  hDwwSel[ind+0]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+1);  hDwwSel[ind+1]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+2);  hDwwSel[ind+2]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+3);  hDwwSel[ind+3]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+4);  hDwwSel[ind+4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+5);  hDwwSel[ind+5] = new TH1D(sb,sb,100,-2.0,2.0);
    sprintf(sb,"hDwwSel_%d",ind+6);  hDwwSel[ind+6] = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+7);  hDwwSel[ind+7] = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+8);  hDwwSel[ind+8] = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+9);  hDwwSel[ind+9]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+10); hDwwSel[ind+10] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDwwSel_%d",ind+11); hDwwSel[ind+11] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+12); hDwwSel[ind+12] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+13); hDwwSel[ind+13] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwSel_%d",ind+14); hDwwSel[ind+14] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+15); hDwwSel[ind+15] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDwwSel_%d",ind+16); hDwwSel[ind+16] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+17); hDwwSel[ind+17] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+18); hDwwSel[ind+18] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+19); hDwwSel[ind+19] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwSel_%d",ind+20); hDwwSel[ind+20] = new TH1D(sb,sb,80,0.0,4.0);
    sprintf(sb,"hDwwSel_%d",ind+21); hDwwSel[ind+21] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDwwSel_%d",ind+22); hDwwSel[ind+22] = new TH1D(sb,sb,100,0.0,5.0);
    sprintf(sb,"hDwwSel_%d",ind+23); hDwwSel[ind+23] = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDwwSel_%d",ind+24); hDwwSel[ind+24] = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDwwSel_%d",ind+25); hDwwSel[ind+25] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+26); hDwwSel[ind+26] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+27); hDwwSel[ind+27] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+28); hDwwSel[ind+28] = new TH1D(sb,sb,100,0.0,5.);
    sprintf(sb,"hDwwSel_%d",ind+29); hDwwSel[ind+29] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+30); hDwwSel[ind+30] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+31); hDwwSel[ind+31] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+32); hDwwSel[ind+32] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwSel_%d",ind+33); hDwwSel[ind+33] = new TH1D(sb,sb,100,0.0,2.0);
    sprintf(sb,"hDwwSel_%d",ind+34); hDwwSel[ind+34] = new TH1D(sb,sb,100,-1.0,1.0);
  }

  for(int i=0; i<35; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDwwSel[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 20 * j;
    sprintf(sb,"hDwwMassSel_%d",ind+ 0); hDwwMassSel[ind+ 0]  = new TH1D(sb,sb,500,0.0,200.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 1); hDwwMassSel[ind+ 1]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 2); hDwwMassSel[ind+ 2]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 3); hDwwMassSel[ind+ 3]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 4); hDwwMassSel[ind+ 4]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 5); hDwwMassSel[ind+ 5]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 6); hDwwMassSel[ind+ 6]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 7); hDwwMassSel[ind+ 7]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 8); hDwwMassSel[ind+ 8]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+ 9); hDwwMassSel[ind+ 9]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+10); hDwwMassSel[ind+10]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+11); hDwwMassSel[ind+11]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+12); hDwwMassSel[ind+12]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+13); hDwwMassSel[ind+13]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+14); hDwwMassSel[ind+14]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+15); hDwwMassSel[ind+15]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+16); hDwwMassSel[ind+16]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+17); hDwwMassSel[ind+17]  = new TH1D(sb,sb,500,0.0,500.);
    sprintf(sb,"hDwwMassSel_%d",ind+18); hDwwMassSel[ind+18]  = new TH1D(sb,sb,500,0.0,500.);
 }

  for(int i=0; i<19; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDwwMassSel[i+j*20]);
    }
  }

  sprintf(sb,"hDwwSelAlphaEP0");
       hDwwSelAlphaEP0 = new TH2D(sb,sb,80,0.0,4.0,72,0.0,180.0);  
  AddOutput(hDwwSelAlphaEP0);
  sprintf(sb,"hDwwSelAlphaEP1");
       hDwwSelAlphaEP1 = new TH2D(sb,sb,80,0.0,4.0,72,0.0,180.0);  
  AddOutput(hDwwSelAlphaEP1);

  sprintf(sb,"hDwwDeltaPhiMetLeptonMet");
       hDwwDeltaPhiMetLeptonMet = new TH2D(sb,sb,36,0.0,180.0,40,0.0,200.0);  
  AddOutput(hDwwDeltaPhiMetLeptonMet);

  sprintf(sb,"hDwwJetVar_%d",0);  hDwwJetVar[0]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",1);  hDwwJetVar[1]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",2);  hDwwJetVar[2]  = new TH1D(sb,sb,50,-0.5,49.5);
  sprintf(sb,"hDwwJetVar_%d",3);  hDwwJetVar[3]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",4);  hDwwJetVar[4]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",5);  hDwwJetVar[5]  = new TH1D(sb,sb,50,-0.5,49.5);
  sprintf(sb,"hDwwJetVar_%d",6);  hDwwJetVar[6]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",7);  hDwwJetVar[7]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwJetVar_%d",8);  hDwwJetVar[8]  = new TH1D(sb,sb,50,-0.5,49.5);
  sprintf(sb,"hDwwJetVar_%d",9);  hDwwJetVar[9]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwJetVar_%d",10); hDwwJetVar[10] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwJetVar_%d",11); hDwwJetVar[11] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwJetVar_%d",12); hDwwJetVar[12] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwJetVar_%d",13); hDwwJetVar[13] = new TH1D(sb,sb,8,-0.5,7.5);
  sprintf(sb,"hDwwJetVar_%d",14); hDwwJetVar[14] = new TH1D(sb,sb,8,-0.5,7.5);
  for(int j=0; j<15; j++){
    AddOutput(hDwwJetVar[j]);
  }


  sprintf(sb,"hDwwMET_%d",0);  hDwwMET[0]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",1);  hDwwMET[1]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",2);  hDwwMET[2]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",3);  hDwwMET[3]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",4);  hDwwMET[4]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",5);  hDwwMET[5]  = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",6);  hDwwMET[6]  = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",7);  hDwwMET[7]  = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",8);  hDwwMET[8]  = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",9);  hDwwMET[9]  = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",10); hDwwMET[10] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",11); hDwwMET[11] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",12); hDwwMET[12] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",13); hDwwMET[13] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",14); hDwwMET[14] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",15); hDwwMET[15] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",16); hDwwMET[16] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",17); hDwwMET[17] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",18); hDwwMET[18] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",19); hDwwMET[19] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",20); hDwwMET[20] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",21); hDwwMET[21] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",22); hDwwMET[22] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",23); hDwwMET[23] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",24); hDwwMET[24] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",25); hDwwMET[25] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",26); hDwwMET[26] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",27); hDwwMET[27] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",28); hDwwMET[28] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",29); hDwwMET[29] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",30); hDwwMET[30] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",31); hDwwMET[31] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",32); hDwwMET[32] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",33); hDwwMET[33] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",34); hDwwMET[34] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",35); hDwwMET[35] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",36); hDwwMET[36] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",37); hDwwMET[37] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",38); hDwwMET[38] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",39); hDwwMET[39] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",40); hDwwMET[40] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",41); hDwwMET[41] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",42); hDwwMET[42] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",43); hDwwMET[43] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",44); hDwwMET[44] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",45); hDwwMET[45] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",46); hDwwMET[46] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",47); hDwwMET[47] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",48); hDwwMET[48] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",49); hDwwMET[49] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",50); hDwwMET[50] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",51); hDwwMET[51] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",52); hDwwMET[52] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",53); hDwwMET[53] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",54); hDwwMET[54] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",55); hDwwMET[55] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",56); hDwwMET[56] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",57); hDwwMET[57] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",58); hDwwMET[58] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",59); hDwwMET[59] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",60); hDwwMET[60] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",61); hDwwMET[61] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",62); hDwwMET[62] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",63); hDwwMET[63] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",64); hDwwMET[64] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",65); hDwwMET[65] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",66); hDwwMET[66] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",67); hDwwMET[67] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",68); hDwwMET[68] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",69); hDwwMET[69] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",70); hDwwMET[70] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",71); hDwwMET[71] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",72); hDwwMET[72] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",73); hDwwMET[73] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",74); hDwwMET[74] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",75); hDwwMET[75] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",76); hDwwMET[76] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",77); hDwwMET[77] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",78); hDwwMET[78] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",79); hDwwMET[79] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",80); hDwwMET[80] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",81); hDwwMET[81] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",82); hDwwMET[82] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",83); hDwwMET[83] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",84); hDwwMET[84] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwMET_%d",85); hDwwMET[85] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",86); hDwwMET[86] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",87); hDwwMET[87] = new TH1D(sb,sb,200,-100.0,100.0);
  sprintf(sb,"hDwwMET_%d",88); hDwwMET[88] = new TH1D(sb,sb,180,0.0,180.0);

  for(int i=0; i<89; i++){
    AddOutput(hDwwMET[i]);
  }

  sprintf(sb,"hDGenZwwSel_%d",0);   hDGenZwwSel[ 0]  = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDGenZwwSel_%d",1);   hDGenZwwSel[ 1]  = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDGenZwwSel_%d",2);   hDGenZwwSel[ 2]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDGenZwwSel_%d",3);   hDGenZwwSel[ 3]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDGenZwwSel_%d",4);   hDGenZwwSel[ 4]  = new TH1D(sb,sb,100,0.0,5.0); 
  sprintf(sb,"hDGenZwwSel_%d",5);   hDGenZwwSel[ 5]  = new TH1D(sb,sb,200,0.0,200.0); 
  sprintf(sb,"hDGenZwwSel_%d",6);   hDGenZwwSel[ 6]  = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDGenZwwSel_%d",7);   hDGenZwwSel[ 7]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDGenZwwSel_%d",8);   hDGenZwwSel[ 8]  = new TH1D(sb,sb,100,0.0,100.0); 
  sprintf(sb,"hDGenZwwSel_%d",9);   hDGenZwwSel[ 9]  = new TH1D(sb,sb,100,0.0,5.0); 

  for(int i=0; i<10; i++){
    AddOutput(hDGenZwwSel[i]);
  }

  sprintf(sb,"hDwwPresel_%d",0);   hDwwPresel[ 0]  = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDwwPresel_%d",1);   hDwwPresel[ 1]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",2);   hDwwPresel[ 2]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",3);   hDwwPresel[ 3]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",4);   hDwwPresel[ 4]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDwwPresel_%d",5);   hDwwPresel[ 5]  = new TH1D(sb,sb,9,-4.5,4.5); 
  sprintf(sb,"hDwwPresel_%d",6);   hDwwPresel[ 6]  = new TH1D(sb,sb,100,-0.5,199.5); 
  sprintf(sb,"hDwwPresel_%d",7);   hDwwPresel[ 7]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDwwPresel_%d",8);   hDwwPresel[ 8]  = new TH1D(sb,sb,200,0.0,20.0);
  sprintf(sb,"hDwwPresel_%d",9);   hDwwPresel[ 9]  = new TH1D(sb,sb,200,0.0,1000.0);
  sprintf(sb,"hDwwPresel_%d",10);  hDwwPresel[10]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",11);  hDwwPresel[11]  = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDwwPresel_%d",12);  hDwwPresel[12]  = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDwwPresel_%d",13);  hDwwPresel[13]  = new TH1D(sb,sb,6,-0.5,5.5); 
  sprintf(sb,"hDwwPresel_%d",14);  hDwwPresel[14]  = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDwwPresel_%d",15);  hDwwPresel[15]  = new TH1D(sb,sb,5,-0.5,4.5); 
  sprintf(sb,"hDwwPresel_%d",16);  hDwwPresel[16]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",17);  hDwwPresel[17]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",18);  hDwwPresel[18]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",19);  hDwwPresel[19]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",20);  hDwwPresel[20]  = new TH1D(sb,sb,6,-0.5,5.5); 
  sprintf(sb,"hDwwPresel_%d",21);  hDwwPresel[21]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwPresel_%d",22);  hDwwPresel[22]  = new TH1D(sb,sb,30,-0.5,29.5); 
  sprintf(sb,"hDwwPresel_%d",23);  hDwwPresel[23]  = new TH1D(sb,sb,30,-0.5,29.5); 
  sprintf(sb,"hDwwPresel_%d",24);  hDwwPresel[24]  = new TH1D(sb,sb,30,-0.5,29.5); 
  sprintf(sb,"hDwwPresel_%d",25);  hDwwPresel[25]  = new TH1D(sb,sb,200,0.0,200.0);

  for(int i=0; i<26; i++){
    AddOutput(hDwwPresel[i]);
  }

  sprintf(sb,"hDwwSelD0Phi");
       hDwwSelD0Phi = new TH2D(sb,sb,45,-180.0,180.0,100,-0.2,0.2);  
  AddOutput(hDwwSelD0Phi);

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDwwJet_%d",ind+0);  hDwwJet[ind+0]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+1);  hDwwJet[ind+1]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+2);  hDwwJet[ind+2]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+3);  hDwwJet[ind+3]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+4);  hDwwJet[ind+4]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+5);  hDwwJet[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+6);  hDwwJet[ind+6]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+7);  hDwwJet[ind+7]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+8);  hDwwJet[ind+8]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+9);  hDwwJet[ind+9]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+10); hDwwJet[ind+10] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+11); hDwwJet[ind+11] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+12); hDwwJet[ind+12] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+13); hDwwJet[ind+13] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+14); hDwwJet[ind+14] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+15); hDwwJet[ind+15] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+16); hDwwJet[ind+16] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+17); hDwwJet[ind+17] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+18); hDwwJet[ind+18] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+19); hDwwJet[ind+19] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+20); hDwwJet[ind+20] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+21); hDwwJet[ind+21] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+22); hDwwJet[ind+22] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+23); hDwwJet[ind+23] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+24); hDwwJet[ind+24] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+25); hDwwJet[ind+25] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+26); hDwwJet[ind+26] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+27); hDwwJet[ind+27] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+28); hDwwJet[ind+28] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+29); hDwwJet[ind+29] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDwwJet_%d",ind+30); hDwwJet[ind+30] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+31); hDwwJet[ind+31] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+32); hDwwJet[ind+32] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+33); hDwwJet[ind+33] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+34); hDwwJet[ind+34] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+35); hDwwJet[ind+35] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+36); hDwwJet[ind+36] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+37); hDwwJet[ind+37] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+38); hDwwJet[ind+38] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJet_%d",ind+39); hDwwJet[ind+39] = new TH1D(sb,sb,200,0.0,200.);
  }
  for(int i=0; i<40; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDwwJet[i+j*100]);
    }
  }

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDwwJetSel_%d",ind+0);  hDwwJetSel[ind+0] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDwwJetSel_%d",ind+1);  hDwwJetSel[ind+1] = new TH1D(sb,sb,100,-5.0,5.0);
    sprintf(sb,"hDwwJetSel_%d",ind+2);  hDwwJetSel[ind+2] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJetSel_%d",ind+3);  hDwwJetSel[ind+3] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwJetSel_%d",ind+4);  hDwwJetSel[ind+4] = new TH1D(sb,sb,100,0.0,5.0);
    sprintf(sb,"hDwwJetSel_%d",ind+5);  hDwwJetSel[ind+5] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJetSel_%d",ind+6);  hDwwJetSel[ind+6] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwJetSel_%d",ind+7);  hDwwJetSel[ind+7] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDwwJetSel_%d",ind+8);  hDwwJetSel[ind+8] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDwwJetSel_%d",ind+9);  hDwwJetSel[ind+9] = new TH1D(sb,sb,100,-5.0,5.0);
  }
  for(int i=0; i<10; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDwwJetSel[i+j*100]);
    }
  }

  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDwwFake_%d",ind+ 0); hDwwFake[ind+ 0] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 1); hDwwFake[ind+ 1] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 2); hDwwFake[ind+ 2] = new TH1D(sb,sb,25,-0.5,24.5);
    sprintf(sb,"hDwwFake_%d",ind+ 3); hDwwFake[ind+ 3] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 4); hDwwFake[ind+ 4] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 5); hDwwFake[ind+ 5] = new TH1D(sb,sb,25,-0.5,24.5);
    sprintf(sb,"hDwwFake_%d",ind+ 6); hDwwFake[ind+ 6] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 7); hDwwFake[ind+ 7] = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwFake_%d",ind+ 8); hDwwFake[ind+ 8] = new TH1D(sb,sb,25,-0.5,24.5);
  }
  for(int i=0; i<9; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDwwFake[i+j*100]);
    }
  }

  for(int j=0; j<50; j++){
    sprintf(sb,"hDWWPDF_%d",j); hDWWPDF[j] = new TH1D(sb,sb,100,0,200);
  }

  for(int j=0; j<50; j++){
    AddOutput(hDWWPDF[j]);
  }

  for(int j=0; j<3; j++){
    int ind = 30 * j;
    sprintf(sb,"hDwwBTag_%d",ind+ 0);  hDwwBTag[ind+ 0]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+ 1);  hDwwBTag[ind+ 1]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+ 2);  hDwwBTag[ind+ 2]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+ 3);  hDwwBTag[ind+ 3]  = new TH1D(sb,sb,200,0.0,4.);
    sprintf(sb,"hDwwBTag_%d",ind+ 4);  hDwwBTag[ind+ 4]  = new TH1D(sb,sb,200,-5.0,15.0);
    sprintf(sb,"hDwwBTag_%d",ind+ 5);  hDwwBTag[ind+ 5]  = new TH1D(sb,sb,200,-5.0,15.0);
    sprintf(sb,"hDwwBTag_%d",ind+ 6);  hDwwBTag[ind+ 6]  = new TH1D(sb,sb,200,0.0,5.);
    sprintf(sb,"hDwwBTag_%d",ind+ 7);  hDwwBTag[ind+ 7]  = new TH1D(sb,sb,120,0.0,6.0);
    sprintf(sb,"hDwwBTag_%d",ind+ 8);  hDwwBTag[ind+ 8]  = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwBTag_%d",ind+ 9);  hDwwBTag[ind+ 9]  = new TH1D(sb,sb,5,-0.5,4.5);
    sprintf(sb,"hDwwBTag_%d",ind+10);  hDwwBTag[ind+10]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+11);  hDwwBTag[ind+11]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+12);  hDwwBTag[ind+12]  = new TH1D(sb,sb,200,0.0,1.);
    sprintf(sb,"hDwwBTag_%d",ind+13);  hDwwBTag[ind+13]  = new TH1D(sb,sb,200,0.0,4.);
    sprintf(sb,"hDwwBTag_%d",ind+14);  hDwwBTag[ind+14]  = new TH1D(sb,sb,200,-5.0,15.0);
    sprintf(sb,"hDwwBTag_%d",ind+15);  hDwwBTag[ind+15]  = new TH1D(sb,sb,200,-5.0,15.0);
    sprintf(sb,"hDwwBTag_%d",ind+16);  hDwwBTag[ind+16]  = new TH1D(sb,sb,200,0.0,5.);
    sprintf(sb,"hDwwBTag_%d",ind+17);  hDwwBTag[ind+17]  = new TH1D(sb,sb,30,0.0,30.0);
    sprintf(sb,"hDwwBTag_%d",ind+18);  hDwwBTag[ind+18]  = new TH1D(sb,sb,100,0,10);
    sprintf(sb,"hDwwBTag_%d",ind+19);  hDwwBTag[ind+19]  = new TH1D(sb,sb,100,0,10);
  }
  int ind = 0;
  sprintf(sb,"hDwwBTag_%d",ind+20);  hDwwBTag[ind+20]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+21);  hDwwBTag[ind+21]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+22);  hDwwBTag[ind+22]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+23);  hDwwBTag[ind+23]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+24);  hDwwBTag[ind+24]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+25);  hDwwBTag[ind+25]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+26);  hDwwBTag[ind+26]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+27);  hDwwBTag[ind+27]  = new TH1D(sb,sb,6,-0.5,5.5);
  sprintf(sb,"hDwwBTag_%d",ind+28);  hDwwBTag[ind+28]  = new TH1D(sb,sb,4,-0.5,3.5);
  sprintf(sb,"hDwwBTag_%d",ind+50);  hDwwBTag[ind+50]  = new TH1D(sb,sb,180,0.0,180.0);
  sprintf(sb,"hDwwBTag_%d",ind+51);  hDwwBTag[ind+51]  = new TH1D(sb,sb,200,0.0,1.);
  sprintf(sb,"hDwwBTag_%d",ind+52);  hDwwBTag[ind+52]  = new TH1D(sb,sb,200,0.0,1.);
  sprintf(sb,"hDwwBTag_%d",ind+53);  hDwwBTag[ind+53]  = new TH1D(sb,sb,200,0.0,1.);
  sprintf(sb,"hDwwBTag_%d",ind+54);  hDwwBTag[ind+54]  = new TH1D(sb,sb,200,0.0,4.);
  sprintf(sb,"hDwwBTag_%d",ind+55);  hDwwBTag[ind+55]  = new TH1D(sb,sb,200,-5.0,15.0);
  sprintf(sb,"hDwwBTag_%d",ind+56);  hDwwBTag[ind+56]  = new TH1D(sb,sb,200,-5.0,15.0);
  sprintf(sb,"hDwwBTag_%d",ind+57);  hDwwBTag[ind+57]  = new TH1D(sb,sb,200,0.0,5.);
  sprintf(sb,"hDwwBTag_%d",ind+58);  hDwwBTag[ind+58]  = new TH1D(sb,sb,50,-0.5,49.5);
  sprintf(sb,"hDwwBTag_%d",ind+59);  hDwwBTag[ind+59]  = new TH1D(sb,sb,50,-0.5,49.5);
  for(int i=0; i<=19; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDwwBTag[i+j*30]);
    }
  }
  for(int i=20; i<=28; i++){
    AddOutput(hDwwBTag[i]);
  }
  for(int i=50; i<=59; i++){
    AddOutput(hDwwBTag[i]);
  }

  sprintf(sb,"hDwwOverlap_%d",0);  hDwwOverlap[ 0]  = new TH1D(sb,sb,100,0.0,2.0);
  sprintf(sb,"hDwwOverlap_%d",1);  hDwwOverlap[ 1]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDwwOverlap_%d",2);  hDwwOverlap[ 2]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwOverlap_%d",3);  hDwwOverlap[ 3]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwOverlap_%d",4);  hDwwOverlap[ 4]  = new TH1D(sb,sb,100,0.0,2.0);
  sprintf(sb,"hDwwOverlap_%d",5);  hDwwOverlap[ 5]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDwwOverlap_%d",6);  hDwwOverlap[ 6]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDwwOverlap_%d",7);  hDwwOverlap[ 7]  = new TH1D(sb,sb,100,0.0,1.0);
  for(int i=0; i<8; i++){
    AddOutput(hDwwOverlap[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
