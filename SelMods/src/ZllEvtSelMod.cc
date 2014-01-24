// $Id: ZllEvtSelMod.cc,v 1.13 2012/04/24 11:26:27 ceballos Exp $

#include "Analysis/SelMods/interface/ZllEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpot.h"

using namespace mithep;
ClassImp(mithep::ZllEvtSelMod)

//--------------------------------------------------------------------------------------------------
ZllEvtSelMod::ZllEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("randomMet"),
  fTrackName(Names::gkTrackBrn),
  fGsfTrackColName("GsfTracks"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fPFMetName("PFMet"),
  fPFMetStd(0),
  fTCMetName("TCMet"),
  fTCMetStd(0),
  fBarrelSCName("BarrelSuperClusters"),
  fEndcapSCName("EndcapSuperClusters"),
  fBarrelSC(0),
  fEndcapSC(0),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fVertexName(ModNames::gkGoodVertexesName),
  fVertices(0),
  fBeamSpotName("BeamSpot"),
  fBeamSpot(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ZllEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ZllEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Obtain all the good objects from the event cleaning modul
  MCParticleOArr *GenLeptons   = GetObjThisEvt<MCParticleOArr>(fMCLeptonsName);
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);

  MuonOArr  *CleanMuonsNoD0         	= GetObjThisEvt<MuonOArr>("CleanMuonsNoD0");
  ElectronOArr *CleanElectronsNoD0  	= GetObjThisEvt<ElectronOArr>("CleanElectronsNoD0");
  //ParticleOArr *leptonsNoD0         	= GetObjThisEvt<ParticleOArr>("MergedLeptonsNoD0");
  MuonOArr  *CleanMuonsNoId             = GetObjThisEvt<MuonOArr>("CleanMuonsNoId");
  ElectronOArr *CleanElectronsNoId      = GetObjThisEvt<ElectronOArr>("CleanElectronsNoId");
  //ParticleOArr *leptonsNoId           = GetObjThisEvt<ParticleOArr>("MergedLeptonsNoId");
  MuonOArr  *CleanMuonsNoIso            = GetObjThisEvt<MuonOArr>("CleanMuonsNoIso");
  ElectronOArr *CleanElectronsNoIso     = GetObjThisEvt<ElectronOArr>("CleanElectronsNoIso");
  //ParticleOArr *leptonsNoIso            = GetObjThisEvt<ParticleOArr>("MergedLeptonsNoIso");
  ElectronOArr *CleanElectronsNoConvF   = GetObjThisEvt<ElectronOArr>("CleanElectronsNoConvF");

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  LoadBranch(fEvtHdrName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fTrackName);
  LoadBranch(fBeamSpotName);
  //Get Generator Level information for matching
  //ObjArray<MCParticle> *GenLeptons = 
  //  dynamic_cast<ObjArray<MCParticle>*> (FindObjThisEvt(fMCLeptonsName.Data()));

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    printf("ZllEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *cleanMet          = CleanMet->At(0);

  LoadBranch(fPFMetName);
  const PFMet *PFMetStd        = fPFMetStd->At(0);

  LoadBranch(fTCMetName);
  const Met *TCMetStd          = fTCMetStd->At(0);

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

  // Z->ee official analysis
  if(CleanElectrons->GetEntries() >= 2 && 
     CleanElectrons->At(0)->SCluster()->Et() > 20 &&
     CleanElectrons->At(1)->SCluster()->Et() > 20){
    CompositeParticle diE;
    diE.AddDaughter(CleanElectrons->At(0));
    diE.AddDaughter(CleanElectrons->At(1));
    if(diE.Mass() > 60 && diE.Mass() < 120 && diE.Charge() == 0){
      if     (( CleanElectrons->At(0)->IsEB() &&  CleanElectrons->At(1)->IsEB())) hDZllZee[0]->Fill(diE.Mass(),NNLOWeight->GetVal());
      else if(( CleanElectrons->At(0)->IsEB() && !CleanElectrons->At(1)->IsEB()) ||
              (!CleanElectrons->At(0)->IsEB() &&  CleanElectrons->At(1)->IsEB())) hDZllZee[1]->Fill(diE.Mass(),NNLOWeight->GetVal());
      else if((!CleanElectrons->At(0)->IsEB() && !CleanElectrons->At(1)->IsEB())) hDZllZee[2]->Fill(diE.Mass(),NNLOWeight->GetVal());
      else printf("What is the status?\n");
    }
  }

  // Minimun Pt, Nleptons>=2 requirements
  if (leptons->GetEntries() >= 2 &&
      leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10){

    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));

    // Sort and count the number of central Jets for vetoing
    double deltaPhiMetLepton = 200.;
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	 CleanJets->At(i)->Pt() > fPtJetCut){
        nCentralJets++;
        if(deltaPhiMetLepton > fabs(MathUtils::DeltaPhi(cleanMet->Phi(), CleanJets->At(i)->Phi())))
	   deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(cleanMet->Phi(), CleanJets->At(i)->Phi()));
      }
    }

    int pairType = -1;
    if (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon )
      pairType = 0;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron)
      pairType = 1;
    else if((leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon) || 
            (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kElectron))
      pairType = 2;
    else {
      cout << "Hey, this is not possible, leptonTypes: "
    	   << leptons->At(0)->ObjType() << " - " 
           << leptons->At(1)->ObjType() << endl;
    }

    bool isGenLepton[2] = {false, false};
    bool isGenElectron[2] = {false, false};
    bool isCorrectQ[2] = {false, false};
    for(int nl=0; nl<2; nl++){
      for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
    	MCParticle *gen = GenLeptons->At(j);
    	if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), leptons->At(nl)->Phi(), leptons->At(nl)->Eta()) < 0.1){
          isGenLepton[nl] = true;
    	  if(leptons->At(nl)->ObjType() == kElectron) isGenElectron[nl] = true;
	  if(gen->Charge() * leptons->At(nl)->Charge() > 0) isCorrectQ[nl] = true;
          break;
        }
      }
    }

    // HLT study
    if(objs){
      Int_t ents=objs->GetEntries();
      double deltaRleptonHLT[2] = {999., 999.};
      for(Int_t i=0;i<ents;++i) {
         const TriggerObject* to = objs->At(i);

         TString trName = to->TrigName();
         if(trName.Contains("HLT_IsoMu24_eta2p1_v")) hDZllHLT[0+100*pairType]->Fill(0.,NNLOWeight->GetVal());
         if(trName.Contains("HLT_Ele27_WP80_v"))     hDZllHLT[0+100*pairType]->Fill(1.,NNLOWeight->GetVal());

        for(int nl=0; nl<2; nl++){
          if     (leptons->At(nl)->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || 
	                                                  to->Type() == TriggerObject::TriggerCluster)){
            double DeltaRAux = MathUtils::DeltaR(leptons->At(nl)->Phi(), leptons->At(nl)->Eta(),
        					 to->Phi(),to->Eta());
    	    if(DeltaRAux < deltaRleptonHLT[nl]) deltaRleptonHLT[nl] = DeltaRAux;
          }
          else if(leptons->At(nl)->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
	                                                      to->Type() == TriggerObject::TriggerCluster  ||
							      to->Type() == TriggerObject::TriggerPhoton)
							      ){
            double DeltaRAux = MathUtils::DeltaR(leptons->At(nl)->Phi(), leptons->At(nl)->Eta(),
        					 to->Phi(),to->Eta());
    	    if(DeltaRAux < deltaRleptonHLT[nl]) deltaRleptonHLT[nl] = DeltaRAux;
          }
        }
        hDZllHLT[1+100*pairType]->Fill((double)to->Type(),NNLOWeight->GetVal());
      } // Loop over HLT objects
      for(int nl=0; nl<2; nl++){
        if     (leptons->At(nl)->ObjType() == kMuon)
          hDZllHLT[2+100*pairType]->Fill(TMath::Min(deltaRleptonHLT[nl],0.999),NNLOWeight->GetVal());
        else if(leptons->At(nl)->ObjType() == kElectron)
          hDZllHLT[3+100*pairType]->Fill(TMath::Min(deltaRleptonHLT[nl],0.999),NNLOWeight->GetVal());
      }
      if      (deltaRleptonHLT[0] < 0.1 && deltaRleptonHLT[1] < 0.1)
        hDZllHLT[4+100*pairType]->Fill(0.0,NNLOWeight->GetVal());
      else if(deltaRleptonHLT[0] < 0.1 || deltaRleptonHLT[1] < 0.1)
        hDZllHLT[4+100*pairType]->Fill(1.0,NNLOWeight->GetVal());
      else
        hDZllHLT[4+100*pairType]->Fill(2.0,NNLOWeight->GetVal());

      int triggerGenType = 0;
      if((isGenLepton[0] == true  && deltaRleptonHLT[0] < 0.1) ||
         (isGenLepton[1] == true  && deltaRleptonHLT[1] < 0.1)) triggerGenType = triggerGenType + 1;
      if((isGenLepton[0] == false && deltaRleptonHLT[0] < 0.1) ||
         (isGenLepton[1] == false && deltaRleptonHLT[1] < 0.1)) triggerGenType = triggerGenType + 2;
      if(isGenElectron[0] == false && isGenElectron[1] == false){
        hDZllHLT[5+100*pairType]->Fill(triggerGenType,NNLOWeight->GetVal());
      }
      else {
        hDZllHLT[6+100*pairType]->Fill(triggerGenType,NNLOWeight->GetVal());
      }

    } // if !objs

    // Begin wrong sign study
    if(leptons->GetEntries() == 2 && TMath::Abs(dilepton.Mass()-91.1876) < 30 &&
       leptons->At(0)->AbsEta() < 2.5 && leptons->At(1)->AbsEta() < 2.5 &&
      (pairType ==0 || pairType == 1)) {
      int nBin0 = (int) (leptons->At(0)->AbsEta() * 2.0);
      int nBin1 = (int) (leptons->At(1)->AbsEta() * 2.0);
      int nhSS = 0;      
      if     ((nBin0 == 0 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(0.0,NNLOWeight->GetVal());
      else if((nBin0 == 0 && nBin1 == 1)||
              (nBin0 == 1 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(1.0,NNLOWeight->GetVal());
      else if((nBin0 == 0 && nBin1 == 2)||
              (nBin0 == 2 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(2.0,NNLOWeight->GetVal());
      else if((nBin0 == 0 && nBin1 == 3)||
              (nBin0 == 3 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(3.0,NNLOWeight->GetVal());
      else if((nBin0 == 0 && nBin1 == 4)||
              (nBin0 == 4 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(4.0,NNLOWeight->GetVal());
      else if((nBin0 == 1 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(5.0,NNLOWeight->GetVal());
      else if((nBin0 == 1 && nBin1 == 2)||
              (nBin0 == 2 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(6.0,NNLOWeight->GetVal());
      else if((nBin0 == 1 && nBin1 == 3)||
              (nBin0 == 3 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(7.0,NNLOWeight->GetVal());
      else if((nBin0 == 1 && nBin1 == 4)||
              (nBin0 == 4 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(8.0,NNLOWeight->GetVal());
      else if((nBin0 == 2 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(9.0,NNLOWeight->GetVal());
      else if((nBin0 == 2 && nBin1 == 3)||
              (nBin0 == 3 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(10.0,NNLOWeight->GetVal());
      else if((nBin0 == 2 && nBin1 == 4)||
              (nBin0 == 4 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(11.0,NNLOWeight->GetVal());
      else if((nBin0 == 3 && nBin1 == 3)) hDZllSS[ nhSS+100*pairType]->Fill(12.0,NNLOWeight->GetVal());
      else if((nBin0 == 3 && nBin1 == 4)||
              (nBin0 == 4 && nBin1 == 3)) hDZllSS[ nhSS+100*pairType]->Fill(13.0,NNLOWeight->GetVal());
      else if((nBin0 == 4 && nBin1 == 4)) hDZllSS[ nhSS+100*pairType]->Fill(14.0,NNLOWeight->GetVal());
      hDZllSS2D->Fill(leptons->At(0)->AbsEta(),leptons->At(1)->AbsEta());
      if(dilepton.Charge() != 0) {
        nhSS = 1;
	if     ((nBin0 == 0 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(0.0,NNLOWeight->GetVal());
	else if((nBin0 == 0 && nBin1 == 1)||
        	(nBin0 == 1 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(1.0,NNLOWeight->GetVal());
	else if((nBin0 == 0 && nBin1 == 2)||
        	(nBin0 == 2 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(2.0,NNLOWeight->GetVal());
	else if((nBin0 == 0 && nBin1 == 3)||
        	(nBin0 == 3 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(3.0,NNLOWeight->GetVal());
	else if((nBin0 == 0 && nBin1 == 4)||
        	(nBin0 == 4 && nBin1 == 0)) hDZllSS[ nhSS+100*pairType]->Fill(4.0,NNLOWeight->GetVal());
	else if((nBin0 == 1 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(5.0,NNLOWeight->GetVal());
	else if((nBin0 == 1 && nBin1 == 2)||
        	(nBin0 == 2 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(6.0,NNLOWeight->GetVal());
	else if((nBin0 == 1 && nBin1 == 3)||
        	(nBin0 == 3 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(7.0,NNLOWeight->GetVal());
	else if((nBin0 == 1 && nBin1 == 4)||
        	(nBin0 == 4 && nBin1 == 1)) hDZllSS[ nhSS+100*pairType]->Fill(8.0,NNLOWeight->GetVal());
	else if((nBin0 == 2 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(9.0,NNLOWeight->GetVal());
	else if((nBin0 == 2 && nBin1 == 3)||
        	(nBin0 == 3 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(10.0,NNLOWeight->GetVal());
	else if((nBin0 == 2 && nBin1 == 4)||
        	(nBin0 == 4 && nBin1 == 2)) hDZllSS[ nhSS+100*pairType]->Fill(11.0,NNLOWeight->GetVal());
	else if((nBin0 == 3 && nBin1 == 3)) hDZllSS[ nhSS+100*pairType]->Fill(12.0,NNLOWeight->GetVal());
	else if((nBin0 == 3 && nBin1 == 4)||
        	(nBin0 == 4 && nBin1 == 3)) hDZllSS[ nhSS+100*pairType]->Fill(13.0,NNLOWeight->GetVal());
	else if((nBin0 == 4 && nBin1 == 4)) hDZllSS[ nhSS+100*pairType]->Fill(14.0,NNLOWeight->GetVal());
      }
      if(isCorrectQ[0] == false) {
        hDZllSS[ 2+100*pairType]->Fill(leptons->At(0)->AbsEta(),NNLOWeight->GetVal());
        hDZllSS[ 3+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      }
      if(isCorrectQ[1] == false) {
        hDZllSS[ 2+100*pairType]->Fill(leptons->At(1)->AbsEta(),NNLOWeight->GetVal());
        hDZllSS[ 3+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
      }
      hDZllSS[ 4+100*pairType]->Fill(leptons->At(0)->AbsEta(),NNLOWeight->GetVal());
      hDZllSS[ 5+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      hDZllSS[ 4+100*pairType]->Fill(leptons->At(1)->AbsEta(),NNLOWeight->GetVal());
      hDZllSS[ 5+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
    }
    // End wrong sign study

    if(pairType == 0 || pairType == 1){
      hDZllMET[ 0+100*pairType]->Fill((double)leptons->GetEntries(),NNLOWeight->GetVal());
      if(leptons->GetEntries() == 2 && dilepton.Charge() == 0)
            hDZllMET[ 1+100*pairType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
      // No more than 2 isolated good leptons
      if(leptons->GetEntries() == 2 && dilepton.Charge() == 0 &&
	 TMath::Abs(dilepton.Mass()-91.1876) < 30){

	// MET study for nCentralJets == 0
	if     (nCentralJets == 0){
	  double METCorrPhi = 0;
	  double METCorrPt  = 0;
 	  double deltaPhiCorMetLepton[2] = {fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(),METCorrPhi) * 180./TMath::Pi()),
	  				    fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(),METCorrPhi) * 180./TMath::Pi())};
 	  double deltaPtCorMetLepton[2] = {TMath::Abs(leptons->At(0)->Pt() - METCorrPt),
	  				   TMath::Abs(leptons->At(1)->Pt() - METCorrPt)};
         if(cleanMet->Pt() < 25.0){
            hDZllMET[ 2+100*pairType]->Fill(TMath::Max(TMath::Min(0.0,49.999),-49.999),NNLOWeight->GetVal());
            hDZllMET[ 4+100*pairType]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
            hDZllMET[ 5+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());        
	    hDZllMET[ 8+100*pairType]->Fill(fabs(MathUtils::DeltaPhi(dilepton.Phi(),METCorrPhi)) * 180./TMath::Pi(),NNLOWeight->GetVal());
	    hDZllMET[36+100*pairType]->Fill(TMath::Min(deltaPhiCorMetLepton[0],deltaPhiCorMetLepton[1]),NNLOWeight->GetVal());
	    hDZllMET[37+100*pairType]->Fill(TMath::Max(deltaPhiCorMetLepton[0],deltaPhiCorMetLepton[1]),NNLOWeight->GetVal());
	    if(deltaPhiCorMetLepton[0] < deltaPhiCorMetLepton[1]){
	      hDZllMET[40+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[0],49.999),NNLOWeight->GetVal());
	      hDZllMET[41+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[1],49.999),NNLOWeight->GetVal());
            }
	    else {
	      hDZllMET[40+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[1],49.999),NNLOWeight->GetVal());
	      hDZllMET[41+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[0],49.999),NNLOWeight->GetVal());
	    }
	  }
	  else {
            hDZllMET[ 3+100*pairType]->Fill(TMath::Max(TMath::Min(0.0,49.999),-49.999),NNLOWeight->GetVal());
            hDZllMET[ 6+100*pairType]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
            hDZllMET[ 7+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
	    hDZllMET[ 9+100*pairType]->Fill(fabs(MathUtils::DeltaPhi(dilepton.Phi(),METCorrPhi)) * 180./TMath::Pi(),NNLOWeight->GetVal());
	    hDZllMET[38+100*pairType]->Fill(TMath::Min(deltaPhiCorMetLepton[0],deltaPhiCorMetLepton[1]),NNLOWeight->GetVal());
	    hDZllMET[39+100*pairType]->Fill(TMath::Max(deltaPhiCorMetLepton[0],deltaPhiCorMetLepton[1]),NNLOWeight->GetVal());
	    if(deltaPhiCorMetLepton[0] < deltaPhiCorMetLepton[1]){
	      hDZllMET[42+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[0],49.999),NNLOWeight->GetVal());
	      hDZllMET[43+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[1],49.999),NNLOWeight->GetVal());
            }
	    else {
	      hDZllMET[42+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[1],49.999),NNLOWeight->GetVal());
	      hDZllMET[43+100*pairType]->Fill(TMath::Min(deltaPtCorMetLepton[0],49.999),NNLOWeight->GetVal());
	    }
            if(pairType == 0 && cleanMet->Pt() > 25.0){
	      int noMuonCorrL = -1;
	      if     (deltaPhiCorMetLepton[0] > 178 && deltaPtCorMetLepton[0] < 5) noMuonCorrL = 0;
	      else if(deltaPhiCorMetLepton[1] > 178 && deltaPtCorMetLepton[1] < 5) noMuonCorrL = 1;
	      if(noMuonCorrL >= 0){
		printf("LARGE MET: %6.3f %6.3f %6.3f | %2d %11d %6d | %2d | %6.3f %6.3f %6.3f / %6.3f %6.3f %6.3f ==> ",
	        	0.0,PFMetStd->Pt(),TCMetStd->Pt(),
	        	fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),
			fMuons->GetEntries(),leptons->At(0)->Pt(),leptons->At(0)->Eta(),leptons->At(0)->Phi(),
		                             leptons->At(1)->Pt(),leptons->At(1)->Eta(),leptons->At(1)->Phi());
        	for(UInt_t i=0; i<fMuons->GetEntries(); i++) printf(" | %6.3f %6.3f %6.3f",fMuons->At(i)->Pt(),fMuons->At(i)->Eta(),fMuons->At(i)->Phi());
		printf("\n");
	        int noMuonCorrM[2] = {-1, -1};
	        for(UInt_t i=0; i<fMuons->GetEntries(); i++) {
	          if(leptons->At(0) == fMuons->At(i)) {
	            noMuonCorrM[0] = i;
	          }
	          if(leptons->At(1) == fMuons->At(i)) {
	            noMuonCorrM[1] = i;
	          }
	        }
		for(UInt_t i=0; i<2; i++) {
		  const Muon *mu = fMuons->At(noMuonCorrM[i]);
    		  Double_t d0_real = 1e30;
    		  for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
                    if(fVertices->At(i0)->NTracks() > 0){
		      Double_t pD0 = mu->BestTrk()->D0Corrected(*fVertices->At(i0));
    		      d0_real = TMath::Abs(pD0);
    		      break;
		    }
		  }
		  Double_t pD0BS = mu->BestTrk()->D0Corrected(*fBeamSpot->At(0));
		  int type = 4;
	  	  if	 (mu->HasGlobalTrk()	) type = 0;
	  	  else if(mu->IsTrackerMuon()	) type = 1;
	  	  else if(mu->HasStandaloneTrk()) type = 2;
	  	  else if(mu->IsCaloMuon()	) type = 3;
	  	  else  			  type = 4;
	  	  hDZllMuon[ 0]->Fill(type,NNLOWeight->GetVal());
                  hDZllMuon[ 1]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
           	  hDZllMuon[ 2]->Fill(TMath::Min(mu->P(),99.999),NNLOWeight->GetVal());
            	  hDZllMuon[ 3]->Fill(mu->Eta(),NNLOWeight->GetVal());
            	  hDZllMuon[ 4]->Fill(mu->Phi() * 180./TMath::Pi(),NNLOWeight->GetVal());
            	  hDZllMuon[ 5]->Fill(TMath::Min(mu->EmEnergy(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[ 6]->Fill(TMath::Min(mu->HadEnergy(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[ 7]->Fill(TMath::Min(mu->IsoR03SumPt(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[ 8]->Fill(TMath::Min(mu->IsoR03EmEt(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[ 9]->Fill(TMath::Min(mu->IsoR03HadEt(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[10]->Fill(TMath::Min(mu->IsoR03HoEt(),19.999),NNLOWeight->GetVal());
            	  hDZllMuon[11]->Fill(TMath::Min((double)mu->IsoR03NTracks(),9.499),NNLOWeight->GetVal());
            	  hDZllMuon[12]->Fill(TMath::Min((double)mu->IsoR03NJets(),9.499),NNLOWeight->GetVal());
            	  hDZllMuon[13]->Fill(TMath::Min((double)mu->NChambers(),19.499),NNLOWeight->GetVal());
            	  hDZllMuon[14]->Fill(TMath::Min((double)mu->NSegments(),19.499),NNLOWeight->GetVal());
            	  hDZllMuon[15]->Fill(TMath::Min((double)mu->LastHit(),19.499),NNLOWeight->GetVal());
	    	  hDZllMuon[16]->Fill(TMath::Min((double)mu->BestTrk()->NHits(),29.499),NNLOWeight->GetVal());
            	  hDZllMuon[17]->Fill(TMath::Min(mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof(),19.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight))  hDZllMuon[18]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TMOneStationLoose))	 hDZllMuon[19]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TMOneStationTight))	 hDZllMuon[20]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TMOneStationLoose))	 hDZllMuon[21]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TMLastStationTight))	 hDZllMuon[22]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose)) hDZllMuon[23]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
	    	  if(mu->Quality().Quality(MuonQuality::TM2DCompatibilityTight)) hDZllMuon[24]->Fill(TMath::Min(mu->Pt(),99.999),NNLOWeight->GetVal());
            	  hDZllMuon[25]->Fill(TMath::Min(TMath::Abs(d0_real),0.0499),NNLOWeight->GetVal());
            	  hDZllMuon[26]->Fill(TMath::Min(TMath::Abs(pD0BS),0.0499),NNLOWeight->GetVal());
	        }
	      } // noMuonCorrL >= 0
	    } // pairType == 10 && cleanMet->Pt() > 25.0
	  }
          hDZllMET[10+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
          hDZllMET[21+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMET[27+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
	}
	else if(nCentralJets == 1){
          hDZllMET[11+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
          hDZllMET[22+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMET[28+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
	}
	else if(nCentralJets == 2){
          hDZllMET[12+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
          hDZllMET[23+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMET[29+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
	}
	else if(nCentralJets == 3){
          hDZllMET[13+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
          hDZllMET[24+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMET[30+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
	}
	else if(nCentralJets >= 4){
          hDZllMET[14+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
          hDZllMET[25+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMET[31+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());
	}
	hDZllMET[15+100*pairType]->Fill(TMath::Min(0.0,199.999),NNLOWeight->GetVal());
	hDZllMET[26+100*pairType]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
	hDZllMET[32+100*pairType]->Fill(TMath::Min(TCMetStd->Pt(),199.999),NNLOWeight->GetVal());

	hDZllMET[16+100*pairType]->Fill((double)nCentralJets,NNLOWeight->GetVal());
	hDZllMET[17+100*pairType]->Fill(leptons->At(0)->Pt(),NNLOWeight->GetVal());
	hDZllMET[18+100*pairType]->Fill(leptons->At(1)->Pt(),NNLOWeight->GetVal());
	hDZllMET[19+100*pairType]->Fill(cleanMet->MetSig(),NNLOWeight->GetVal());
	hDZllMET[20+100*pairType]->Fill(cleanMet->SumEt(),NNLOWeight->GetVal());

	for(UInt_t i=0; i<2; i++){
          if(deltaPhiMetLepton > fabs(MathUtils::DeltaPhi(cleanMet->Phi(), leptons->At(i)->Phi())))
	     deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(cleanMet->Phi(), leptons->At(i)->Phi()));
	}
	if(nCentralJets == 0){
          hDZllMET[33+100*pairType]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	}
	else if(nCentralJets == 1){
          hDZllMET[34+100*pairType]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	}
	else {
          hDZllMET[35+100*pairType]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	}
      } // Nleptons == 2 , Ncharge == 0 and TMath::Abs(mll-mz)<30 requirements
    } // Minimun Pt and Nleptons >= 2 requirements
  } // Either two muons or two electrons

  // Begin tag and probe study
  // Minimun Pt, Nleptons>=1 requirements
  if (leptons->GetEntries() >= 1 && leptons->At(0)->Pt() > 20){

    ObjArray<Track> *GoodTracks = new ObjArray<Track>;
    for(unsigned int i = 0; i < fTracks->GetEntries(); i++) {
      const mithep::Track* pTrack = fTracks->At(i);
      if(pTrack->Pt() <= 10) continue;
      if(TMath::Abs(pTrack->Eta()) >= 2.4) continue;
      GoodTracks->Add(pTrack);
    }
    LoadBranch(fGsfTrackColName);
    ObjArray<Track> *GoodGsfTracks = new ObjArray<Track>;
    for(unsigned int i = 0; i < fGsfTracks->GetEntries(); i++) {
      const mithep::Track* pTrack = fGsfTracks->At(i);
      if(pTrack->Pt() <= 10) continue;
      if(TMath::Abs(pTrack->Eta()) >= 2.4) continue;
      if(!pTrack->SCluster()) continue;
      GoodGsfTracks->Add(pTrack);
    }
    LoadBranch(fBarrelSCName);
    LoadBranch(fEndcapSCName);
    ObjArray<SuperCluster> *GoodSC = new ObjArray<SuperCluster>;
    for(unsigned int i = 0; i < fBarrelSC->GetEntries(); i++) {
      const SuperCluster* sc = fBarrelSC->At(i);
      double theP  = sqrt(sc->Point().X() * sc->Point().X() + sc->Point().Y() * sc->Point().Y() + sc->Point().Z() * sc->Point().Z());
      double thePt = sc->Energy() * sqrt(sc->Point().X() * sc->Point().X() + sc->Point().Y() * sc->Point().Y())/theP;
      if(thePt <= 10) continue;
      if(TMath::Abs(sc->Eta()) >= 2.4) continue;
      GoodSC->Add(sc);
    }
    for(unsigned int i = 0; i < fEndcapSC->GetEntries(); i++) {
      const SuperCluster* sc = fEndcapSC->At(i);
      double theP  = sqrt(sc->Point().X() * sc->Point().X() + sc->Point().Y() * sc->Point().Y() + sc->Point().Z() * sc->Point().Z());
      double thePt = sc->Energy() * sqrt(sc->Point().X() * sc->Point().X() + sc->Point().Y() * sc->Point().Y())/theP;
      if(thePt <= 10) continue;
      if(TMath::Abs(sc->Eta()) >= 2.4) continue;
      GoodSC->Add(sc);
    }

    // Muons
    int LeptonType = 0;
    for (UInt_t i=0; i<CleanMuons->GetEntries(); i++) {
      const Muon *mu = CleanMuons->At(i);
      if(mu->HasGlobalTrk()     == kFALSE ||
         mu->HasStandaloneTrk() == kFALSE ||
	 mu->HasTrackerTrk()    == kFALSE) continue;
      if(mu->Pt() <= 20) continue;
      bool isTriggerLepton = kFALSE;
      for(UInt_t nt=0; nt<objs->GetEntries(); nt++) {
        const TriggerObject* to = objs->At(nt);
        if((to->Type() == TriggerObject::TriggerMuon || 
            to->Type() == TriggerObject::TriggerCluster) &&
           MathUtils::DeltaR(mu->Phi(), mu->Eta(), to->Phi(),to->Eta()) < 0.1){
	  isTriggerLepton = kTRUE;
	  break;
	}
      }
      if(isTriggerLepton == kFALSE) continue;

      // Looking at GEN level
      for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
        MCParticle *gen = GenLeptons->At(j);
        if(gen->AbsEta() > 2.4) continue;
        if(gen->Pt() < 10.0) continue;
	if(gen->Is(MCParticle::kMu) == kFALSE) continue;
        if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), mu->Phi(), mu->Eta()) < 0.1) continue;
        hDZllTPGen[0]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
        hDZllTPGen[1]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
        for (UInt_t jm=0; jm<fMuons->GetEntries(); jm++) {
          const Muon *cleanmu = fMuons->At(jm);
	  if(cleanmu->HasGlobalTrk() == kFALSE) continue;
	  if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), cleanmu->Phi(), cleanmu->Eta()) < 0.1 &&
	     gen->Charge() * cleanmu->Charge() > 0){
            hDZllTPGen[2]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
            hDZllTPGen[3]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	    break;
	  }
	}
        for (UInt_t jm=0; jm<CleanMuons->GetEntries(); jm++) {
          const Muon *cleanmu = CleanMuons->At(jm);
	  if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), cleanmu->Phi(), cleanmu->Eta()) < 0.1 &&
	     gen->Charge() * cleanmu->Charge() > 0){
            hDZllTPGen[4]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
            hDZllTPGen[5]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
	    break;
	  }
	}
      }

      // Looking at RECO level
      int nInter[7] = {0, 0, 0, 0, 0, 0, 0};
      int indexInter[7] = {-1, -1, -1, -1, -1, -1, -1};
      double ptInter[7] = {-1, -1, -1, -1, -1, -1, -1};
      for (UInt_t j=0; j<GoodTracks->GetEntries(); j++) {
        const Track *track = GoodTracks->At(j);
        if(track->Pt() <= 10) continue;
        if(TMath::Abs(track->Eta()) >= 2.4) continue;
        if(MathUtils::DeltaR(track->Phi(), track->Eta(), mu->TrackerTrk()->Phi(), mu->TrackerTrk()->Eta()) < 0.1 ||
           track == mu->TrackerTrk()) continue;
        GenericParticle *p = new GenericParticle(track->Px(), track->Py(), track->Pz(), track->P(), track->Charge());
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(p);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[0]++;
	  if(track->Pt() > ptInter[0]){
	    ptInter[0]    = track->Pt();
	    indexInter[0] = j;
	  }
	}
	delete p;
      }
      for (UInt_t j=0; j<fMuons->GetEntries(); j++) {
        const Muon *muProbe = fMuons->At(j);
        if(muProbe->Pt() <= 10) continue;
        if(muProbe->AbsEta() >= 2.4) continue;
        if(muProbe == mu) continue;
	if(muProbe->HasStandaloneTrk() == kFALSE) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[1]++;
	  if(muProbe->Pt() > ptInter[1]){
	    ptInter[1]    = muProbe->Pt();
	    indexInter[1] = j;
	  }
	  if(muProbe->HasGlobalTrk() == kTRUE){
	    nInter[2]++;
	    if(muProbe->Pt() > ptInter[2]){
	      ptInter[2]    = muProbe->Pt();
	      indexInter[2] = j;
	    }
	  }
	}
      }
      for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
        const Muon *muProbe = CleanMuons->At(j);
        if(muProbe->Pt() <= 10) continue;
        if(muProbe->AbsEta() >= 2.4) continue;
        if(muProbe == mu) continue;

        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[3]++;
	  if(muProbe->Pt() > ptInter[3]){
	    ptInter[3]    = muProbe->Pt();
	    indexInter[3] = j;
	  }
	}
      }
      for (UInt_t j=0; j<CleanMuonsNoD0->GetEntries(); j++) {
        const Muon *muProbe = CleanMuonsNoD0->At(j);
        if(muProbe->Pt() <= 10) continue;
        if(muProbe->AbsEta() >= 2.4) continue;
        if(muProbe == mu) continue;

        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[4]++;
	  if(muProbe->Pt() > ptInter[4]){
	    ptInter[4]    = muProbe->Pt();
	    indexInter[4] = j;
	  }
	}
      }
      for (UInt_t j=0; j<CleanMuonsNoId->GetEntries(); j++) {
        const Muon *muProbe = CleanMuonsNoId->At(j);
        if(muProbe->Pt() <= 10) continue;
        if(muProbe->AbsEta() >= 2.4) continue;
        if(muProbe == mu) continue;

        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[5]++;
	  if(muProbe->Pt() > ptInter[5]){
	    ptInter[5]    = muProbe->Pt();
	    indexInter[5] = j;
	  }
	}
      }
      for (UInt_t j=0; j<CleanMuonsNoIso->GetEntries(); j++) {
        const Muon *muProbe = CleanMuonsNoIso->At(j);
        if(muProbe->Pt() <= 10) continue;
        if(muProbe->AbsEta() >= 2.4) continue;
        if(muProbe == mu) continue;

        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[6]++;
	  if(muProbe->Pt() > ptInter[6]){
	    ptInter[6]    = muProbe->Pt();
	    indexInter[6] = j;
	  }
	}
      }
      hDZllTPIni[ 0+LeptonType*100]->Fill(nInter[0],NNLOWeight->GetVal());
      hDZllTPIni[ 1+LeptonType*100]->Fill(nInter[1],NNLOWeight->GetVal());
      hDZllTPIni[ 2+LeptonType*100]->Fill(nInter[2],NNLOWeight->GetVal());
      hDZllTPIni[ 3+LeptonType*100]->Fill(nInter[3],NNLOWeight->GetVal());
      hDZllTPIni[ 4+LeptonType*100]->Fill(nInter[4],NNLOWeight->GetVal());
      hDZllTPIni[ 5+LeptonType*100]->Fill(nInter[5],NNLOWeight->GetVal());
      hDZllTPIni[ 6+LeptonType*100]->Fill(nInter[6],NNLOWeight->GetVal());
      if(nInter[0] > 0){ // standalone muon efficiency
        const Track *track = GoodTracks->At(indexInter[0]);
        GenericParticle *p = new GenericParticle(track->Px(), track->Py(), track->Pz(), track->P(), track->Charge());
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(p);
	int theHisto = 0;
	if(mu->Charge()+track->Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<fMuons->GetEntries(); j++) {
          const Muon *cleanmu = fMuons->At(j);
	  if(cleanmu->HasTrackerTrk() == kFALSE) continue;
	  if(track == cleanmu->TrackerTrk() && cleanmu->HasStandaloneTrk()){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[ 0+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[ 1+LeptonType*100+theHisto]->Fill(TMath::Min(track->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[ 2+LeptonType*100+theHisto]->Fill(TMath::Abs(track->Eta()),NNLOWeight->GetVal());    
        }
	else {
	  hDZllTP[ 6+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[ 7+LeptonType*100+theHisto]->Fill(TMath::Min(track->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[ 8+LeptonType*100+theHisto]->Fill(TMath::Abs(track->Eta()),NNLOWeight->GetVal());
	}
	delete p;
      }
      if(nInter[1] > 0){ // global muon efficiency
        const Muon *muProbe = fMuons->At(indexInter[1]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<fMuons->GetEntries(); j++) {
          const Muon *cleanmu = fMuons->At(j);
	  if(cleanmu == muProbe && cleanmu->HasTrackerTrk() && cleanmu->HasGlobalTrk()){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[12+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[13+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[14+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[18+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[19+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[20+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[2] > 0){ // full muon id efficiency
        const Muon *muProbe = fMuons->At(indexInter[2]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	if(dilepton.Charge() == 0 && TMath::Abs(dilepton.Mass()-91.1876) < 15.0){
	  hDZllMuSel[0]->Fill(TMath::Min((double)muProbe->BestTrk()->NHits(),29.499),NNLOWeight->GetVal());
          hDZllMuSel[1]->Fill(TMath::Min(muProbe->BestTrk()->Chi2()/muProbe->BestTrk()->Ndof(),19.999),NNLOWeight->GetVal());
	  if(muProbe->Quality().Quality(MuonQuality::GlobalMuonPromptTight))  hDZllMuSel[2]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMuSel[3]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllMuSel[4]->Fill(TMath::Min(muProbe->AbsEta(),2.4999),NNLOWeight->GetVal());
          Double_t totalIso = 1.0 * muProbe->IsoR03SumPt() + 
                              1.0 * muProbe->IsoR03EmEt() + 
                              1.0 * muProbe->IsoR03HadEt();
          hDZllMuSel[5]->Fill(TMath::Min(totalIso/muProbe->Pt(),0.999),NNLOWeight->GetVal());
          hDZllMuSel[8]->Fill(TMath::Min(muProbe->EmEnergy(),19.999),NNLOWeight->GetVal());
          hDZllMuSel[9]->Fill(TMath::Min(muProbe->HadEnergy(),19.999),NNLOWeight->GetVal());
          if(muProbe->IsoR03SumPt()/ muProbe->Pt() < 0.10 && muProbe->Pt() > 20.0 && muProbe->Pt() < 40.0){
    	    double sumPt = 0.0; int nTracks = 0;
            Double_t zLepton = 0.0;
            if(muProbe->BestTrk()) zLepton = muProbe->BestTrk()->DzCorrected(*fVertices->At(0));
    	    for(unsigned int i = 0; i < fTracks->GetEntries(); i++) {
    	      const mithep::Track* pTrack = fTracks->At(i);
    	      if(MathUtils::DeltaR(pTrack->Mom(), muProbe->Mom()) < 0.05 ||
    	  	 MathUtils::DeltaR(pTrack->Mom(), muProbe->Mom()) > 0.3) continue;
              Double_t deltaZ = TMath::Abs(pTrack->DzCorrected(*fVertices->At(0)) - zLepton);
	      if(deltaZ > 0.1) continue;
    	      sumPt = sumPt + pTrack->Pt(); nTracks++;
    	    }
            hDZllMuSel[6]->Fill(TMath::Min(sumPt,9.999),NNLOWeight->GetVal());
            hDZllMuSel[7]->Fill(TMath::Min((double)nTracks,9.499),NNLOWeight->GetVal());
    	  }
        }
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
          const Muon *cleanmu = CleanMuons->At(j);
	  if(cleanmu == muProbe){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[24+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[25+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[26+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[30+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[31+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[32+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[3] > 0){ // trigger muon efficiency
        const Muon *muProbe = CleanMuons->At(indexInter[3]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
      	bool isTriggerLepton = kFALSE;
      	for(UInt_t nt=0; nt<objs->GetEntries(); nt++) {
      	  const TriggerObject* to = objs->At(nt);
      	  if(muProbe->AbsEta() < 2.1 && (to->Type() == TriggerObject::TriggerMuon || 
	                                 to->Type() == TriggerObject::TriggerCluster) &&
      	     MathUtils::DeltaR(muProbe->Phi(), muProbe->Eta(), to->Phi(),to->Eta()) < 0.1){
	    isTriggerLepton = kTRUE;
	    break;
	  }
      	}
	if(isTriggerLepton){
          hDZllTP[36+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[37+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[38+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[42+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[43+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[44+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[4] > 0){ // d0 muon efficiency
        const Muon *muProbe = CleanMuonsNoD0->At(indexInter[4]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
          const Muon *cleanmu = CleanMuons->At(j);
	  if(cleanmu == muProbe){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[48+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[49+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[50+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[54+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[55+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[56+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[5] > 0){ // id muon efficiency
        const Muon *muProbe = CleanMuonsNoId->At(indexInter[5]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
          const Muon *cleanmu = CleanMuons->At(j);
	  if(cleanmu == muProbe){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[60+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[61+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[62+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[66+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[67+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[68+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[6] > 0){ // iso muon efficiency
        const Muon *muProbe = CleanMuonsNoIso->At(indexInter[6]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(mu);
        dilepton.AddDaughter(muProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodMuon = kFALSE;
        for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
          const Muon *cleanmu = CleanMuons->At(j);
	  if(cleanmu == muProbe){
	    isGoodMuon = kTRUE;
	    break;
	  }
	}
	if(isGoodMuon){
          hDZllTP[72+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[73+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[74+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[78+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[79+LeptonType*100+theHisto]->Fill(TMath::Min(muProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[80+LeptonType*100+theHisto]->Fill(TMath::Abs(muProbe->Eta()),NNLOWeight->GetVal());
	}
      }
    }

    // Electrons
    LeptonType = 1;
    for (UInt_t i=0; i<CleanElectrons->GetEntries(); i++) {
      const Electron *el = CleanElectrons->At(i);
      if(el->HasGsfTrk() == kFALSE ||
         el->HasTrackerTrk() == kFALSE ||
         el->HasSuperCluster() == kFALSE) continue;
      if(el->Pt() <= 20) continue;
      bool isTriggerLepton = kFALSE;
      for(UInt_t nt=0; nt<objs->GetEntries(); nt++) {
        const TriggerObject* to = objs->At(nt);
        if((to->Type() == TriggerObject::TriggerElectron || 
	    to->Type() == TriggerObject::TriggerCluster  ||
	    to->Type() == TriggerObject::TriggerPhoton) &&
           MathUtils::DeltaR(el->Phi(), el->Eta(), to->Phi(),to->Eta()) < 0.1) {
	  isTriggerLepton = kTRUE;
	  break;
	}
      }
      if(isTriggerLepton == kFALSE) continue;

      // Looking at GEN level
      for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
        MCParticle *gen = GenLeptons->At(j);
        if(gen->AbsEta() > 2.4) continue;
        if(gen->Pt() < 10.0) continue;
	if(gen->Is(MCParticle::kEl) == kFALSE) continue;
        if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), el->Phi(), el->Eta()) < 0.1) continue;
        hDZllTPGen[6]->Fill(TMath::Min(gen->Pt(),199.999));
        hDZllTPGen[7]->Fill(gen->AbsEta());
        for (UInt_t je=0; je<fElectrons->GetEntries(); je++) {
          const Electron *cleanel = fElectrons->At(je);
	  if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), cleanel->Phi(), cleanel->Eta()) < 0.05 &&
	     gen->Charge() * cleanel->Charge() > 0){
            hDZllTPGen[8]->Fill(TMath::Min(gen->Pt(),199.999));
            hDZllTPGen[9]->Fill(gen->AbsEta());
	    break;
	  }
	}
        for (UInt_t je=0; je<CleanElectrons->GetEntries(); je++) {
          const Electron *cleanel = CleanElectrons->At(je);
	  if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), cleanel->Phi(), cleanel->Eta()) < 0.05 &&
	     gen->Charge() * cleanel->Charge() > 0){
            hDZllTPGen[10]->Fill(TMath::Min(gen->Pt(),199.999));
            hDZllTPGen[11]->Fill(gen->AbsEta());
	    break;
	  }
	}
      }

      // Looking at RECO level
      int nInter[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      int indexInter[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
      double ptInter[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
      for (UInt_t j=0; j<GoodTracks->GetEntries(); j++) {
        const Track *track = GoodTracks->At(j);
        if(track->Pt() <= 10) continue;
        if(TMath::Abs(track->Eta()) >= 2.4) continue;
        if(MathUtils::DeltaR(track->Phi(), track->Eta(), el->TrackerTrk()->Phi(), el->TrackerTrk()->Eta()) < 0.1 ||
           track == el->TrackerTrk()) continue;
	double ratio = 1.0;
        GenericParticle *p = new GenericParticle(track->Px()*ratio, track->Py()*ratio, 
	                                         track->Pz()*ratio, track->P()*ratio, track->Charge());
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(p);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[0]++;
	  if(track->Pt() > ptInter[0]){
	    ptInter[0]    = track->Pt();
	    indexInter[0] = j;
	  }
	}
	delete p;
      }
      for (UInt_t j=0; j<GoodSC->GetEntries(); j++) {
        const SuperCluster *sc = GoodSC->At(j);
        if(MathUtils::DeltaR(sc->Phi(), sc->Eta(), el->SCluster()->Phi(), el->SCluster()->Eta()) < 0.1 ||
           sc == el->SCluster()) continue;
        double theP  = sqrt(sc->Point().X() * sc->Point().X() +
	                    sc->Point().Y() * sc->Point().Y() +
			    sc->Point().Z() * sc->Point().Z());
        double thePt = sc->Energy() * sqrt(sc->Point().X() * sc->Point().X() + sc->Point().Y() * sc->Point().Y())/theP;
        if(thePt <= 10) continue;
        if(TMath::Abs(sc->Eta()) >= 2.4) continue;
        GenericParticle *p = new GenericParticle(sc->Point().X()*sc->Energy()/theP,sc->Point().Y()*sc->Energy()/theP,
	                                         sc->Point().Z()*sc->Energy()/theP,sc->Energy(),0);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(p);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[1]++;
	  if(p->Pt() > ptInter[1]){
	    ptInter[1]    = p->Pt();
	    indexInter[1] = j;
	  }
	}
	delete p;
      }
      for (UInt_t j=0; j<fElectrons->GetEntries(); j++) {
        const Electron *elProbe = fElectrons->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(MathUtils::DeltaR(elProbe->Phi(), elProbe->Eta(), el->GsfTrk()->Phi(), el->GsfTrk()->Eta()) < 0.1 ||
           elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[2]++;
	  if(elProbe->Pt() > ptInter[2]){
	    ptInter[2]    = elProbe->Pt();
	    indexInter[2] = j;
	  }
	}
      }
      for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
        const Electron *elProbe = CleanElectrons->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[3]++;
	  if(elProbe->Pt() > ptInter[3]){
	    ptInter[3]    = elProbe->Pt();
	    indexInter[3] = j;
	  }
	}
      }
      for (UInt_t j=0; j<CleanElectronsNoD0->GetEntries(); j++) {
        const Electron *elProbe = CleanElectronsNoD0->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[4]++;
	  if(elProbe->Pt() > ptInter[4]){
	    ptInter[4]    = elProbe->Pt();
	    indexInter[4] = j;
	  }
	}
     }
     for (UInt_t j=0; j<CleanElectronsNoId->GetEntries(); j++) {
        const Electron *elProbe = CleanElectronsNoId->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[5]++;
	  if(elProbe->Pt() > ptInter[5]){
	    ptInter[5]    = elProbe->Pt();
	    indexInter[5] = j;
	  }
	}
      }
     for (UInt_t j=0; j<CleanElectronsNoIso->GetEntries(); j++) {
        const Electron *elProbe = CleanElectronsNoIso->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[6]++;
	  if(elProbe->Pt() > ptInter[6]){
	    ptInter[6]    = elProbe->Pt();
	    indexInter[6] = j;
	  }
	}
      }
     for (UInt_t j=0; j<CleanElectronsNoConvF->GetEntries(); j++) {
        const Electron *elProbe = CleanElectronsNoConvF->At(j);
        if(elProbe->Pt() <= 10) continue;
        if(elProbe->AbsEta() >= 2.4) continue;
        if(elProbe == el) continue;
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Mass() > 50 && dilepton.Mass() < 150) {
	  nInter[7]++;
	  if(elProbe->Pt() > ptInter[7]){
	    ptInter[7]    = elProbe->Pt();
	    indexInter[7] = j;
	  }
	}
      }

      hDZllTPIni[ 0+LeptonType*100]->Fill(nInter[0],NNLOWeight->GetVal());
      hDZllTPIni[ 1+LeptonType*100]->Fill(nInter[1],NNLOWeight->GetVal());
      hDZllTPIni[ 2+LeptonType*100]->Fill(nInter[2],NNLOWeight->GetVal());
      hDZllTPIni[ 3+LeptonType*100]->Fill(nInter[3],NNLOWeight->GetVal());
      hDZllTPIni[ 4+LeptonType*100]->Fill(nInter[4],NNLOWeight->GetVal());
      hDZllTPIni[ 5+LeptonType*100]->Fill(nInter[5],NNLOWeight->GetVal());
      hDZllTPIni[ 6+LeptonType*100]->Fill(nInter[6],NNLOWeight->GetVal());
      if(nInter[0] > 0){ // supercluster electron efficiency
        const Track *track = GoodTracks->At(indexInter[0]);
        GenericParticle *p = new GenericParticle(track->Px(), track->Py(), 
	                                         track->Pz(), track->P(), track->Charge());
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(p);
	int theHisto = 0;
	if(el->Charge()+track->Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
	if(track->SCluster()){
          for (UInt_t j=0; j<GoodSC->GetEntries(); j++) {
            const SuperCluster *sc = GoodSC->At(j);
            if(track->SCluster() == sc){
	      isGoodElectron = kTRUE;
	      break;
	    }
	  }
	}
	if(isGoodElectron){
          hDZllTP[ 0+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[ 1+LeptonType*100+theHisto]->Fill(TMath::Min(track->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[ 2+LeptonType*100+theHisto]->Fill(TMath::Abs(track->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[ 6+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[ 7+LeptonType*100+theHisto]->Fill(TMath::Min(track->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[ 8+LeptonType*100+theHisto]->Fill(TMath::Abs(track->Eta()),NNLOWeight->GetVal());
	}
	delete p;
      }
      if(nInter[1] > 0){ // track electron efficiency
        const SuperCluster *sc = GoodSC->At(indexInter[1]);
        double theP  = sqrt(sc->Point().X() * sc->Point().X() +
	                    sc->Point().Y() * sc->Point().Y() +
			    sc->Point().Z() * sc->Point().Z());
        GenericParticle *p = new GenericParticle(sc->Point().X()*sc->Energy()/theP,sc->Point().Y()*sc->Energy()/theP,
	                                         sc->Point().Z()*sc->Energy()/theP,sc->Energy(),0);

        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(p);
	int theHisto = 0;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<fElectrons->GetEntries(); j++) {
          const Electron *cleanel = fElectrons->At(j);
          if(cleanel->SCluster() == sc){
	    isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[12+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[13+LeptonType*100+theHisto]->Fill(TMath::Min(p->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[14+LeptonType*100+theHisto]->Fill(TMath::Abs(p->Eta()),NNLOWeight->GetVal());
        }
	else {
          hDZllTP[18+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[19+LeptonType*100+theHisto]->Fill(TMath::Min(p->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[20+LeptonType*100+theHisto]->Fill(TMath::Abs(p->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[2] > 0){ // full electron id efficiency
        const Electron *elProbe = fElectrons->At(indexInter[2]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	if(dilepton.Charge() == 0 && TMath::Abs(dilepton.Mass()-91.1876) < 15.0){
          Double_t hOverE	= elProbe->HadronicOverEm();
          Double_t sigmaee	= elProbe->CoviEtaiEta();
          Double_t deltaPhiIn   = TMath::Abs(elProbe->DeltaPhiSuperClusterTrackAtVtx());
          Double_t deltaEtaIn   = TMath::Abs(elProbe->DeltaEtaSuperClusterTrackAtVtx());
          //Double_t eSeedOverPin = elProbe->ESeedClusterOverPIn(); 
          Double_t combinedIso = elProbe->TrackIsolationDr03() + elProbe->EcalRecHitIsoDr03() + elProbe->HcalTowerSumEtDr03();
          if(elProbe->IsEB()) combinedIso = elProbe->TrackIsolationDr03() + TMath::Max(elProbe->EcalRecHitIsoDr03() - 1.0, 0.0) + elProbe->HcalTowerSumEtDr03();
          combinedIso = combinedIso / elProbe->Pt();

          hDZllElSel[ 0+10*elProbe->IsEE()]->Fill(TMath::Min(sigmaee,0.04999),NNLOWeight->GetVal());
          if((elProbe->IsEB() && sigmaee < 0.012) ||
	     (elProbe->IsEE() && sigmaee < 0.032)){
            hDZllElSel[ 2+10*elProbe->IsEE()]->Fill(TMath::Min(hOverE,0.1999),NNLOWeight->GetVal());
            hDZllElSel[ 3+10*elProbe->IsEE()]->Fill(TMath::Min(deltaPhiIn,0.09999),NNLOWeight->GetVal());
            hDZllElSel[ 4+10*elProbe->IsEE()]->Fill(TMath::Min(deltaEtaIn,0.01999),NNLOWeight->GetVal());
            hDZllElSel[ 6+10*elProbe->IsEE()]->Fill(TMath::Min((double)elProbe->BestTrk()->NExpectedHitsInner(),9.49999),NNLOWeight->GetVal());
            hDZllElSel[ 7+10*elProbe->IsEE()]->Fill(TMath::Min(combinedIso,0.99999),NNLOWeight->GetVal());
            hDZllElSel[ 8+10*elProbe->IsEE()]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
            hDZllElSel[ 9+10*elProbe->IsEE()]->Fill(TMath::Min(TMath::Abs(elProbe->Eta()),2.4999),NNLOWeight->GetVal());
    	    if(elProbe->TrackIsolationDr03()/ elProbe->Pt() < 0.10 && el->Pt() > 20.0 && elProbe->Pt() < 40.0){
    	      double sumPt = 0.0; int nTracks = 0;
              Double_t zLepton = 0.0;
              if(elProbe->BestTrk()) zLepton = elProbe->BestTrk()->DzCorrected(*fVertices->At(0));
    	      for(unsigned int i = 0; i < fTracks->GetEntries(); i++) {
    	    	const mithep::Track* pTrack = fTracks->At(i);
    	    	if(MathUtils::DeltaR(pTrack->Mom(), elProbe->Mom()) < 0.05 ||
    	    	   MathUtils::DeltaR(pTrack->Mom(), elProbe->Mom()) > 0.3) continue;
                Double_t deltaZ = TMath::Abs(pTrack->DzCorrected(*fVertices->At(0)) - zLepton);
	        if(deltaZ > 0.1) continue;
    	    	sumPt = sumPt + pTrack->Pt(); nTracks++;
    	      }
              hDZllElSel[ 1+10*elProbe->IsEE()]->Fill(TMath::Min(sumPt,9.999),NNLOWeight->GetVal());
              hDZllElSel[ 5+10*elProbe->IsEE()]->Fill(TMath::Min((double)nTracks,9.499),NNLOWeight->GetVal());
    	    }
	  }
	}
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
          const Electron *cleanel = CleanElectrons->At(j);
	  if(cleanel == elProbe){
            isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[24+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[25+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[26+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[30+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[31+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[32+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[3] > 0){ // trigger electron efficiency
        const Electron *elProbe = CleanElectrons->At(indexInter[3]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
      	bool isTriggerLepton = kFALSE;
      	for(UInt_t nt=0; nt<objs->GetEntries(); nt++) {
      	  const TriggerObject* to = objs->At(nt);
      	  if((to->Type() == TriggerObject::TriggerElectron || 
	      to->Type() == TriggerObject::TriggerCluster  ||
	      to->Type() == TriggerObject::TriggerPhoton) && 
	      MathUtils::DeltaR(elProbe->Phi(), elProbe->Eta(), to->Phi(),to->Eta()) < 0.1){
	    isTriggerLepton = kTRUE;
	    break;
	  }
      	}
	if(isTriggerLepton){
          hDZllTP[36+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[37+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[38+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[42+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[43+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[44+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[4] > 0){ // d0 electron efficiency
        const Electron *elProbe = CleanElectronsNoD0->At(indexInter[4]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
          const Electron *cleanel = CleanElectrons->At(j);
	  if(cleanel == elProbe){
            isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[48+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[49+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[50+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());	 
        }
	else {
          hDZllTP[54+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[55+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[56+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[5] > 0){ // id electron efficiency
        const Electron *elProbe = CleanElectronsNoId->At(indexInter[5]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
          const Electron *cleanel = CleanElectrons->At(j);
	  if(cleanel == elProbe){
            isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[60+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[61+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[62+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[66+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[67+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[68+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[6] > 0){ // iso electron efficiency
        const Electron *elProbe = CleanElectronsNoIso->At(indexInter[6]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
          const Electron *cleanel = CleanElectrons->At(j);
	  if(cleanel == elProbe){
            isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[72+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[73+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[74+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[78+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[79+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[80+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
      if(nInter[7] > 0){ // conversion electron efficiency
        const Electron *elProbe = CleanElectronsNoConvF->At(indexInter[7]);
        CompositeParticle dilepton;
        dilepton.AddDaughter(el);
        dilepton.AddDaughter(elProbe);
	int theHisto = 0;
	if(dilepton.Charge() != 0) theHisto = 3;
	bool isGoodElectron = kFALSE;
        for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {
          const Electron *cleanel = CleanElectrons->At(j);
	  if(cleanel == elProbe){
            isGoodElectron = kTRUE;
	    break;
	  }
	}
	if(isGoodElectron){
          hDZllTP[84+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[85+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[86+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());    
        }
	else {
          hDZllTP[90+LeptonType*100+theHisto]->Fill(TMath::Min(dilepton.Mass(),149.999),NNLOWeight->GetVal());
          hDZllTP[91+LeptonType*100+theHisto]->Fill(TMath::Min(elProbe->Pt(),199.999),NNLOWeight->GetVal());
          hDZllTP[92+LeptonType*100+theHisto]->Fill(TMath::Abs(elProbe->Eta()),NNLOWeight->GetVal());
	}
      }
    }
    delete GoodTracks;
    delete GoodGsfTracks;
  } // End tag and probe study
}
//--------------------------------------------------------------------------------------------------
void ZllEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fEvtHdrName,      fEventHeader);
  ReqBranch(fMuonName,        fMuons);
  ReqBranch(fElectronName,    fElectrons);
  ReqBranch(fPFMetName,   fPFMetStd);
  ReqBranch(fTCMetName,   fTCMetStd);
  ReqBranch(fTrackName,       fTracks);
  ReqBranch(fGsfTrackColName, fGsfTracks);
  ReqBranch(fBarrelSCName,    fBarrelSC);
  ReqBranch(fEndcapSCName,    fEndcapSC);
  ReqBranch(fBeamSpotName,    fBeamSpot);

  char sb[200];
  sprintf(sb,"hDZllZee_%d", 0);  hDZllZee[0]  = new TH1D(sb,sb,50,70,110);
  sprintf(sb,"hDZllZee_%d", 1);  hDZllZee[1]  = new TH1D(sb,sb,50,70,110);
  sprintf(sb,"hDZllZee_%d", 2);  hDZllZee[2]  = new TH1D(sb,sb,50,70,110);
  for(int j=0; j<3; j++){
    AddOutput(hDZllZee[j]);
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZllMET_%d",ind+0);  hDZllMET[ind+0]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllMET_%d",ind+1);  hDZllMET[ind+1]  = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+2);  hDZllMET[ind+2]  = new TH1D(sb,sb,100,-50.0,50.);
    sprintf(sb,"hDZllMET_%d",ind+3);  hDZllMET[ind+3]  = new TH1D(sb,sb,100,-50.0,50.);
    sprintf(sb,"hDZllMET_%d",ind+4);  hDZllMET[ind+4]  = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+5);  hDZllMET[ind+5]  = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+6);  hDZllMET[ind+6]  = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+7);  hDZllMET[ind+7]  = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+8);  hDZllMET[ind+8]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDZllMET_%d",ind+9);  hDZllMET[ind+9]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDZllMET_%d",ind+10); hDZllMET[ind+10] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+11); hDZllMET[ind+11] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+12); hDZllMET[ind+12] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+13); hDZllMET[ind+13] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+14); hDZllMET[ind+14] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+15); hDZllMET[ind+15] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+16); hDZllMET[ind+16] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllMET_%d",ind+17); hDZllMET[ind+17] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDZllMET_%d",ind+18); hDZllMET[ind+18] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+19); hDZllMET[ind+19] = new TH1D(sb,sb,150,0.0,30.);
    sprintf(sb,"hDZllMET_%d",ind+20); hDZllMET[ind+20] = new TH1D(sb,sb,200,0.0,800.);
    sprintf(sb,"hDZllMET_%d",ind+21); hDZllMET[ind+21] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+22); hDZllMET[ind+22] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+23); hDZllMET[ind+23] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+24); hDZllMET[ind+24] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+25); hDZllMET[ind+25] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+26); hDZllMET[ind+26] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+27); hDZllMET[ind+27] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+28); hDZllMET[ind+28] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+29); hDZllMET[ind+29] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+30); hDZllMET[ind+30] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+31); hDZllMET[ind+31] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+32); hDZllMET[ind+32] = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDZllMET_%d",ind+33); hDZllMET[ind+33] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+34); hDZllMET[ind+34] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+35); hDZllMET[ind+35] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+36); hDZllMET[ind+36] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+37); hDZllMET[ind+37] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+38); hDZllMET[ind+38] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+39); hDZllMET[ind+39] = new TH1D(sb,sb,90,0.,180.); 
    sprintf(sb,"hDZllMET_%d",ind+40); hDZllMET[ind+40] = new TH1D(sb,sb,50,0.,50.); 
    sprintf(sb,"hDZllMET_%d",ind+41); hDZllMET[ind+41] = new TH1D(sb,sb,50,0.,50.); 
    sprintf(sb,"hDZllMET_%d",ind+42); hDZllMET[ind+42] = new TH1D(sb,sb,50,0.,50.); 
    sprintf(sb,"hDZllMET_%d",ind+43); hDZllMET[ind+43] = new TH1D(sb,sb,50,0.,50.); 
  }

  for(int i=0; i<44; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZllMET[i+j*100]);
    }
  }

  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZllHLT_%d",ind+0);  hDZllHLT[ind+0]  = new TH1D(sb,sb,9,-0.5,8.5);
    sprintf(sb,"hDZllHLT_%d",ind+1);  hDZllHLT[ind+1]  = new TH1D(sb,sb,201,-100.5,100.5);
    sprintf(sb,"hDZllHLT_%d",ind+2);  hDZllHLT[ind+2]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZllHLT_%d",ind+3);  hDZllHLT[ind+3]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZllHLT_%d",ind+4);  hDZllHLT[ind+4]  = new TH1D(sb,sb,3,-0.5,2.5);
    sprintf(sb,"hDZllHLT_%d",ind+5);  hDZllHLT[ind+5]  = new TH1D(sb,sb,4,-0.5,3.5);
    sprintf(sb,"hDZllHLT_%d",ind+6);  hDZllHLT[ind+6]  = new TH1D(sb,sb,4,-0.5,3.5);
  }

  for(int i=0; i<7; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDZllHLT[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZllSS_%d",ind+ 0); hDZllSS[ind+ 0] = new TH1D(sb,sb,15,-0.5,14.5);
    sprintf(sb,"hDZllSS_%d",ind+ 1); hDZllSS[ind+ 1] = new TH1D(sb,sb,15,-0.5,14.5);
    sprintf(sb,"hDZllSS_%d",ind+ 2); hDZllSS[ind+ 2] = new TH1D(sb,sb,50,0.0,2.5);
    sprintf(sb,"hDZllSS_%d",ind+ 3); hDZllSS[ind+ 3] = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZllSS_%d",ind+ 4); hDZllSS[ind+ 4] = new TH1D(sb,sb,50,0.0,2.5);
    sprintf(sb,"hDZllSS_%d",ind+ 5); hDZllSS[ind+ 5] = new TH1D(sb,sb,200,0.0,200.0);
  }

  for(int i=0; i<6; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZllSS[i+j*100]);
    }
  }

  sprintf(sb,"hDZllSS2D");   hDZllSS2D = new TH2D(sb,sb,25,0.0,2.5,25,0.0,2.5); 

  AddOutput(hDZllSS2D);

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZllTPIni_%d",ind+ 0); hDZllTPIni[ind+ 0] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 1); hDZllTPIni[ind+ 1] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 2); hDZllTPIni[ind+ 2] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 3); hDZllTPIni[ind+ 3] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 4); hDZllTPIni[ind+ 4] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 5); hDZllTPIni[ind+ 5] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllTPIni_%d",ind+ 6); hDZllTPIni[ind+ 6] = new TH1D(sb,sb,10,-0.5,9.5);
  }
  for(int i=0; i<7; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZllTPIni[i+j*100]);
    }
  }
  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDZllTP_%d",ind+ 0); hDZllTP[ind+ 0] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+ 1); hDZllTP[ind+ 1] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+ 2); hDZllTP[ind+ 2] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+ 3); hDZllTP[ind+ 3] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+ 4); hDZllTP[ind+ 4] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+ 5); hDZllTP[ind+ 5] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+ 6); hDZllTP[ind+ 6] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+ 7); hDZllTP[ind+ 7] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+ 8); hDZllTP[ind+ 8] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+ 9); hDZllTP[ind+ 9] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+10); hDZllTP[ind+10] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+11); hDZllTP[ind+11] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+12); hDZllTP[ind+12] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+13); hDZllTP[ind+13] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+14); hDZllTP[ind+14] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+15); hDZllTP[ind+15] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+16); hDZllTP[ind+16] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+17); hDZllTP[ind+17] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+18); hDZllTP[ind+18] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+19); hDZllTP[ind+19] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+20); hDZllTP[ind+20] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+21); hDZllTP[ind+21] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+22); hDZllTP[ind+22] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+23); hDZllTP[ind+23] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+24); hDZllTP[ind+24] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+25); hDZllTP[ind+25] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+26); hDZllTP[ind+26] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+27); hDZllTP[ind+27] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+28); hDZllTP[ind+28] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+29); hDZllTP[ind+29] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+30); hDZllTP[ind+30] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+31); hDZllTP[ind+31] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+32); hDZllTP[ind+32] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+33); hDZllTP[ind+33] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+34); hDZllTP[ind+34] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+35); hDZllTP[ind+35] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+36); hDZllTP[ind+36] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+37); hDZllTP[ind+37] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+38); hDZllTP[ind+38] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+39); hDZllTP[ind+39] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+40); hDZllTP[ind+40] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+41); hDZllTP[ind+41] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+42); hDZllTP[ind+42] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+43); hDZllTP[ind+43] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+44); hDZllTP[ind+44] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+45); hDZllTP[ind+45] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+46); hDZllTP[ind+46] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+47); hDZllTP[ind+47] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+48); hDZllTP[ind+48] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+49); hDZllTP[ind+49] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+50); hDZllTP[ind+50] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+51); hDZllTP[ind+51] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+52); hDZllTP[ind+52] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+53); hDZllTP[ind+53] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+54); hDZllTP[ind+54] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+55); hDZllTP[ind+55] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+56); hDZllTP[ind+56] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+57); hDZllTP[ind+57] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+58); hDZllTP[ind+58] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+59); hDZllTP[ind+59] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+60); hDZllTP[ind+60] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+61); hDZllTP[ind+61] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+62); hDZllTP[ind+62] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+63); hDZllTP[ind+63] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+64); hDZllTP[ind+64] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+65); hDZllTP[ind+65] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+66); hDZllTP[ind+66] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+67); hDZllTP[ind+67] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+68); hDZllTP[ind+68] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+69); hDZllTP[ind+69] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+70); hDZllTP[ind+70] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+71); hDZllTP[ind+71] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+72); hDZllTP[ind+72] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+73); hDZllTP[ind+73] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+74); hDZllTP[ind+74] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+75); hDZllTP[ind+75] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+76); hDZllTP[ind+76] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+77); hDZllTP[ind+77] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+78); hDZllTP[ind+78] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+79); hDZllTP[ind+79] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+80); hDZllTP[ind+80] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+81); hDZllTP[ind+81] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+82); hDZllTP[ind+82] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+83); hDZllTP[ind+83] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+84); hDZllTP[ind+84] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+85); hDZllTP[ind+85] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+86); hDZllTP[ind+86] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+87); hDZllTP[ind+87] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+88); hDZllTP[ind+88] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+89); hDZllTP[ind+89] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+90); hDZllTP[ind+90] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+91); hDZllTP[ind+91] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+92); hDZllTP[ind+92] = new TH1D(sb,sb,48,0,2.4);
    sprintf(sb,"hDZllTP_%d",ind+93); hDZllTP[ind+93] = new TH1D(sb,sb,100,50,150);
    sprintf(sb,"hDZllTP_%d",ind+94); hDZllTP[ind+94] = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDZllTP_%d",ind+95); hDZllTP[ind+95] = new TH1D(sb,sb,48,0,2.4);
  }
  for(int i=0; i<96; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZllTP[i+j*100]);
    }
  }

  sprintf(sb,"hDZllTPGen_%d", 0); hDZllTPGen[ 0] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d", 1); hDZllTPGen[ 1] = new TH1D(sb,sb, 48,0,2.4);
  sprintf(sb,"hDZllTPGen_%d", 2); hDZllTPGen[ 2] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d", 3); hDZllTPGen[ 3] = new TH1D(sb,sb, 48,0,2.4);
  sprintf(sb,"hDZllTPGen_%d", 4); hDZllTPGen[ 4] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d", 5); hDZllTPGen[ 5] = new TH1D(sb,sb, 48,0,2.4);
  sprintf(sb,"hDZllTPGen_%d", 6); hDZllTPGen[ 6] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d", 7); hDZllTPGen[ 7] = new TH1D(sb,sb, 48,0,2.4);
  sprintf(sb,"hDZllTPGen_%d", 8); hDZllTPGen[ 8] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d", 9); hDZllTPGen[ 9] = new TH1D(sb,sb, 48,0,2.4);
  sprintf(sb,"hDZllTPGen_%d",10); hDZllTPGen[10] = new TH1D(sb,sb, 200,0,200);
  sprintf(sb,"hDZllTPGen_%d",11); hDZllTPGen[11] = new TH1D(sb,sb, 48,0,2.4);
  for(int i=0; i<12; i++){
    AddOutput(hDZllTPGen[i]);
  }

  sprintf(sb,"hDZllMuon_%d", 0); hDZllMuon[ 0] = new TH1D(sb,sb,5,-0.5,4.5);
  sprintf(sb,"hDZllMuon_%d", 1); hDZllMuon[ 1] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d", 2); hDZllMuon[ 2] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d", 3); hDZllMuon[ 3] = new TH1D(sb,sb,200, -5.0, 5.0);
  sprintf(sb,"hDZllMuon_%d", 4); hDZllMuon[ 4] = new TH1D(sb,sb,45, 0.0,180.0);
  sprintf(sb,"hDZllMuon_%d", 5); hDZllMuon[ 5] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d", 6); hDZllMuon[ 6] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d", 7); hDZllMuon[ 7] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d", 8); hDZllMuon[ 8] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d", 9); hDZllMuon[ 9] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d",10); hDZllMuon[10] = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDZllMuon_%d",11); hDZllMuon[11] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZllMuon_%d",12); hDZllMuon[12] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDZllMuon_%d",13); hDZllMuon[13] = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDZllMuon_%d",14); hDZllMuon[14] = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDZllMuon_%d",15); hDZllMuon[15] = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDZllMuon_%d",16); hDZllMuon[16] = new TH1D(sb,sb,30,-0.5,29.5);
  sprintf(sb,"hDZllMuon_%d",17); hDZllMuon[17] = new TH1D(sb,sb,100,0.0,20.0);
  sprintf(sb,"hDZllMuon_%d",18); hDZllMuon[18] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",19); hDZllMuon[19] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",20); hDZllMuon[20] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",21); hDZllMuon[21] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",22); hDZllMuon[22] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",23); hDZllMuon[23] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",24); hDZllMuon[24] = new TH1D(sb,sb,200,0.0,100.);
  sprintf(sb,"hDZllMuon_%d",25); hDZllMuon[25] = new TH1D(sb,sb,50,0.0,0.05);
  sprintf(sb,"hDZllMuon_%d",26); hDZllMuon[26] = new TH1D(sb,sb,50,0.0,0.05);

  for(int i=0; i<27; i++){
    AddOutput(hDZllMuon[i]);
  }

  for(int j=0; j<2; j++){
    int ind = 10 * j;
    sprintf(sb,"hDZllElSel_%d",ind+0);  hDZllElSel[ind+0]  = new TH1D(sb,sb,100,0.0,0.05);
    sprintf(sb,"hDZllElSel_%d",ind+1);  hDZllElSel[ind+1]  = new TH1D(sb,sb,100,0.0,10.0);
    sprintf(sb,"hDZllElSel_%d",ind+2);  hDZllElSel[ind+2]  = new TH1D(sb,sb,100,0.0,0.2);
    sprintf(sb,"hDZllElSel_%d",ind+3);  hDZllElSel[ind+3]  = new TH1D(sb,sb,200,0.0,0.1);
    sprintf(sb,"hDZllElSel_%d",ind+4);  hDZllElSel[ind+4]  = new TH1D(sb,sb,100,0.0,0.02);
    sprintf(sb,"hDZllElSel_%d",ind+5);  hDZllElSel[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllElSel_%d",ind+6);  hDZllElSel[ind+6]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllElSel_%d",ind+7);  hDZllElSel[ind+7]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZllElSel_%d",ind+8);  hDZllElSel[ind+8]  = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZllElSel_%d",ind+9);  hDZllElSel[ind+9]  = new TH1D(sb,sb,100,0.0,2.5);
  }

  for(int i=0; i<10; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDZllElSel[i+j*10]);
    }
  }

  for(int j=0; j<1; j++){
    int ind = 10 * j;
    sprintf(sb,"hDZllMuSel_%d",ind+0);  hDZllMuSel[ind+0]  = new TH1D(sb,sb,30,-0.5,29.5);
    sprintf(sb,"hDZllMuSel_%d",ind+1);  hDZllMuSel[ind+1]  = new TH1D(sb,sb,200,0.0,20.0);
    sprintf(sb,"hDZllMuSel_%d",ind+2);  hDZllMuSel[ind+2]  = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZllMuSel_%d",ind+3);  hDZllMuSel[ind+3]  = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDZllMuSel_%d",ind+4);  hDZllMuSel[ind+4]  = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDZllMuSel_%d",ind+5);  hDZllMuSel[ind+5]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDZllMuSel_%d",ind+6);  hDZllMuSel[ind+6]  = new TH1D(sb,sb,100,0.0,10.0);
    sprintf(sb,"hDZllMuSel_%d",ind+7);  hDZllMuSel[ind+7]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDZllMuSel_%d",ind+8);  hDZllMuSel[ind+8]  = new TH1D(sb,sb,200,0.0,20.0);
    sprintf(sb,"hDZllMuSel_%d",ind+9);  hDZllMuSel[ind+9]  = new TH1D(sb,sb,200,0.0,20.0);
  }

  for(int i=0; i<10; i++){
    for(int j=0; j<1; j++){
      AddOutput(hDZllMuSel[i+j*10]);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void ZllEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ZllEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
