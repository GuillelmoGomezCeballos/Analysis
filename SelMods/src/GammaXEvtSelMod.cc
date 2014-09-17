// $Id: GammaXEvtSelMod.cc,v 1.6 2012/01/26 13:17:56 ceballos Exp $

#include "Analysis/SelMods/interface/GammaXEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"

using namespace mithep;
ClassImp(mithep::GammaXEvtSelMod)

//--------------------------------------------------------------------------------------------------
GammaXEvtSelMod::GammaXEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fGoodJetsName(ModNames::gkGoodJetsName),        
  fMCPhotonsName(ModNames::gkMCPhotonsName),
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),  
  fPhotonBranchName(Names::gkPhotonBrn),
  fPhotons(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GammaXEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void GammaXEvtSelMod::Process()
{

  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Get Generator Level information for matching
  MCParticleOArr *GenPhotons = GetObjThisEvt<MCParticleOArr>(fMCPhotonsName);

  //Obtain all the good objects from the event cleaning module
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  JetOArr *GoodJets            = GetObjThisEvt<JetOArr>(fGoodJetsName);
  PhotonOArr *CleanPhotons     = GetObjThisEvt<PhotonOArr>(fCleanPhotonsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  LoadBranch(fPhotonBranchName);

  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);

  hDBeforeHLTSel[0]->Fill(TMath::Min((double)CleanPhotons->GetEntries(),9.499),NNLOWeight->GetVal());

  for (UInt_t j=0; j<GenPhotons->GetEntries(); ++j) {
    MCParticle *gen = GenPhotons->At(j);
    if(gen->AbsEta() > 2.5) continue;
    if(gen->Pt() < 15.0) continue;
    hDGMC[0]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
    hDGMC[1]->Fill(gen->AbsEta(),NNLOWeight->GetVal());

    bool isRecoPhoton = kFALSE;
    for (UInt_t i=0; i<CleanPhotons->GetEntries(); ++i) {
      Photon *ph = CleanPhotons->At(i);
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), ph->Phi(), ph->Eta()) < 0.1){
        isRecoPhoton = kTRUE;
        break;
      }
    }
    if(isRecoPhoton == kTRUE){
      hDGMC[2]->Fill(TMath::Min(gen->Pt(),199.999),NNLOWeight->GetVal());
      hDGMC[3]->Fill(gen->AbsEta(),NNLOWeight->GetVal());
    }
  }

  for(UInt_t i=0; i<fPhotons->GetEntries(); i++){
    const Photon *ph = fPhotons->At(i);
    if(ph->SCluster() == 0) continue;
    if(ph->Pt() <= 15) continue;
    if(ph->AbsEta() >= 2.5) continue;
    //if(ph->IsTightPhoton() == false) continue;
    bool isGoodMCPhoton = kFALSE;
    for (UInt_t j=0; j<GenPhotons->GetEntries(); ++j) {
      MCParticle *gen = GenPhotons->At(j);
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), ph->Phi(), ph->Eta()) < 0.1){
        isGoodMCPhoton = kTRUE;
        break;
      }
    }
    int nHisto = 0;
    if(isGoodMCPhoton == kFALSE) nHisto = 40;

    Bool_t cutPh[7] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

    Double_t totalIso = ph->HollowConeTrkIsoDr04()+
                        ph->EcalRecHitIsoDr04() +
                        ph->HcalTowerSumEtDr04();

    cutPh[0] = ph->HasPixelSeed() == false;

    cutPh[1] = ph->HadOverEm() < 0.05;

    cutPh[2] = ph->R9() > 0.5;

    cutPh[3] = totalIso/ph->Pt() < 0.25;

    cutPh[4] = ph->IsEB() || ph->IsEE();

    cutPh[5] = ph->SCluster()->EtaWidth() < 0.031;
    if(ph->IsEB()) cutPh[5] = ph->SCluster()->EtaWidth() < 0.013;

    cutPh[6] = ph->AbsEta() < 2.5;

    if(ph->IsEB()){
      hDGSel[ 0+nHisto]->Fill(ph->IsTightPhoton(),NNLOWeight->GetVal());	
      hDGSel[ 1+nHisto]->Fill(ph->IsLoosePhoton(),NNLOWeight->GetVal());	
    }
    else if(ph->IsEE()) {
      hDGSel[ 2+nHisto]->Fill(ph->IsTightPhoton(),NNLOWeight->GetVal());	  
      hDGSel[ 3+nHisto]->Fill(ph->IsLoosePhoton(),NNLOWeight->GetVal());
    }
    if(cutPh[1] && cutPh[2] && cutPh[3] && cutPh[4] && cutPh[5] && cutPh[6]){
      hDGSel[ 4+nHisto]->Fill(ph->HasPixelSeed(),NNLOWeight->GetVal());
    }
    if(cutPh[0] && cutPh[2] && cutPh[3] && cutPh[4] && cutPh[5] && cutPh[6]){
      hDGSel[ 5+nHisto]->Fill(TMath::Min(ph->HadOverEm(),0.0999),NNLOWeight->GetVal());
    }
    if(cutPh[0] && cutPh[1] && cutPh[3] && cutPh[4] && cutPh[5] && cutPh[6]){
      hDGSel[ 6+nHisto]->Fill(TMath::Max(TMath::Min(ph->R9(),0.9999),0.50001),NNLOWeight->GetVal());
      if(ph->IsConverted()) hDGSel[ 7+nHisto]->Fill(TMath::Max(TMath::Min(ph->R9(),0.9999),0.50001),NNLOWeight->GetVal());
      else                  hDGSel[ 8+nHisto]->Fill(TMath::Max(TMath::Min(ph->R9(),0.9999),0.50001),NNLOWeight->GetVal());
    }
    if(cutPh[0] && cutPh[1] && cutPh[2] && cutPh[4] && cutPh[5] && cutPh[6]){
      hDGSel[ 9+nHisto]->Fill(TMath::Max(TMath::Min(totalIso/ph->Pt(),0.999),0.0),NNLOWeight->GetVal());
      if     (ph->Pt() > 25.0)
        hDGSel[10+nHisto]->Fill(TMath::Max(TMath::Min(totalIso/ph->Pt(),0.999),0.0),NNLOWeight->GetVal());
      else if(ph->Pt() > 20.0)
        hDGSel[11+nHisto]->Fill(TMath::Max(TMath::Min(totalIso/ph->Pt(),0.999),0.0),NNLOWeight->GetVal());
      else
        hDGSel[12+nHisto]->Fill(TMath::Max(TMath::Min(totalIso/ph->Pt(),0.999),0.0),NNLOWeight->GetVal());

      if(isGoodMCPhoton == kTRUE) hDGSel2D[0]->Fill(TMath::Max(totalIso,0.0),ph->Pt());
      else			  hDGSel2D[1]->Fill(TMath::Max(totalIso,0.0),ph->Pt());
    }
    if(cutPh[0] && cutPh[1] && cutPh[2] && cutPh[3] && cutPh[5] && cutPh[6]){
      if     (ph->IsEB())      hDGSel[13+nHisto]->Fill(0.0,NNLOWeight->GetVal());
      else if(ph->IsEE())      hDGSel[13+nHisto]->Fill(1.0,NNLOWeight->GetVal());
      else if(ph->IsEBGap())   hDGSel[13+nHisto]->Fill(2.0,NNLOWeight->GetVal());
      else if(ph->IsEEGap())   hDGSel[13+nHisto]->Fill(3.0,NNLOWeight->GetVal());
      else if(ph->IsEBEEGap()) hDGSel[13+nHisto]->Fill(4.0,NNLOWeight->GetVal());
      else                     hDGSel[13+nHisto]->Fill(5.0,NNLOWeight->GetVal());
    }
    if(cutPh[0] && cutPh[1] && cutPh[2] && cutPh[3] && cutPh[4] && cutPh[6]){
      if(ph->IsEB()){
	hDGSel[14+nHisto]->Fill(TMath::Min(ph->SCluster()->EtaWidth(),0.049999),NNLOWeight->GetVal());
      }
      else if(ph->IsEE()) {
	hDGSel[15+nHisto]->Fill(TMath::Min(ph->SCluster()->EtaWidth(),0.049999),NNLOWeight->GetVal());
      }
    }
    if(cutPh[0] && cutPh[1] && cutPh[2] && cutPh[3] && cutPh[4] && cutPh[5]){
      hDGSel[16+nHisto]->Fill(TMath::Min(ph->AbsEta(),4.999),NNLOWeight->GetVal());
    }

    if(cutPh[0] && cutPh[1] && cutPh[2] && cutPh[3] && cutPh[4] && cutPh[5] && cutPh[6]){
      hDGSel[17+nHisto]->Fill(TMath::Min(ph->Pt(),199.999),NNLOWeight->GetVal());
      hDGSel[18+nHisto]->Fill(TMath::Min(ph->HollowConeTrkIsoDr04(),39.999),NNLOWeight->GetVal());
      hDGSel[19+nHisto]->Fill(TMath::Max(TMath::Min(ph->EcalRecHitIsoDr04(),39.999),-1.999),NNLOWeight->GetVal());
      hDGSel[20+nHisto]->Fill(TMath::Min(ph->HcalTowerSumEtDr04(),39.999),NNLOWeight->GetVal());
      hDGSel[21+nHisto]->Fill(TMath::Min(ph->SolidConeTrkIsoDr04(),39.999),NNLOWeight->GetVal());
      hDGSel[22+nHisto]->Fill(TMath::Min(ph->HcalDepth1TowerSumEtDr04(),39.999),NNLOWeight->GetVal());
      hDGSel[23+nHisto]->Fill(TMath::Min(ph->HcalDepth2TowerSumEtDr04(),39.999),NNLOWeight->GetVal());
      hDGSel[24+nHisto]->Fill(ph->SolidConeNTrkDr04(),NNLOWeight->GetVal());
      hDGSel[25+nHisto]->Fill(ph->HollowConeNTrkDr04(),NNLOWeight->GetVal());
      hDGSel[26+nHisto]->Fill(ph->IsConverted(),NNLOWeight->GetVal());
      if(ph->IsEB()){
        hDGSel[27+nHisto]->Fill(ph->IsTightPhoton(),NNLOWeight->GetVal());       
        hDGSel[28+nHisto]->Fill(ph->IsLoosePhoton(),NNLOWeight->GetVal());
      }
      else if(ph->IsEE()) {
        hDGSel[29+nHisto]->Fill(ph->IsTightPhoton(),NNLOWeight->GetVal());       
        hDGSel[30+nHisto]->Fill(ph->IsLoosePhoton(),NNLOWeight->GetVal());
      }
    }
  }

  if(CleanPhotons->GetEntries() <= 0) return;

  bool isGoodMCPhoton[2] = {kFALSE, kFALSE};
  Photon *ph0 = CleanPhotons->At(0);
  Photon *ph1 = CleanPhotons->At(0);
  if(CleanPhotons->GetEntries() >= 2) ph1 = CleanPhotons->At(1);
  for (UInt_t j=0; j<GenPhotons->GetEntries(); ++j) {
    MCParticle *gen = GenPhotons->At(j);
    if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), ph0->Phi(), ph0->Eta()) < 0.1){
      isGoodMCPhoton[0] = kTRUE;
    }
    if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), ph1->Phi(), ph1->Eta()) < 0.1){
      isGoodMCPhoton[1] = kTRUE;
    }
  }

  // Sort and count the number of central Jets for vetoing
  int nCentralJets = 0;
  for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
    if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
       CleanJets->At(i)->Pt() > fPtJetCut){
      nCentralJets++;
    }
  }
  
  // di photon analysis
  int nTotal = TMath::Min((int)CleanPhotons->GetEntries(),9);
  nTotal = nTotal + 10 * TMath::Min((int)leptons->GetEntries(),9);
  nTotal = nTotal + 100 * TMath::Min(nCentralJets,5);
  hDGGSel[0]->Fill(nTotal,NNLOWeight->GetVal());
  if(CleanPhotons->GetEntries() >= 2){
    int nHisto = 0;
    if(isGoodMCPhoton[0] == kFALSE || isGoodMCPhoton[1] == kFALSE) nHisto = 20;
    CompositeParticle diphoton;
    diphoton.AddDaughter(CleanPhotons->At(0));
    diphoton.AddDaughter(CleanPhotons->At(1));
    hDGGSel[1+nHisto]->Fill(TMath::Min(diphoton.Mass(),199.999),NNLOWeight->GetVal());
    hDGGSel[2+nHisto]->Fill(TMath::Min(CleanPhotons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
    hDGGSel[3+nHisto]->Fill(TMath::Min(CleanPhotons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
  }
  
  // 1 gamma + 1 lepton
  if(CleanPhotons->GetEntries() >= 1 && leptons->GetEntries() == 1){
    int nHisto = 0;
    if(leptons->At(0)->ObjType() == kElectron) nHisto = 40;
    //if(isGoodMCPhoton[0] == kFALSE) nHisto += 20;
    double deltaRlG =
       MathUtils::DeltaR(CleanPhotons->At(0)->Phi(), CleanPhotons->At(0)->Eta(),
     	                 leptons->At(0)->Phi(), leptons->At(0)->Eta());
    double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), 
                                                   leptons->At(0)->Phi()));

    double mTW = TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
   			    (1.0 - cos(deltaPhiMetLepton)));

    CompositeParticle totalPart;
    totalPart.AddDaughter(leptons->At(0));
    totalPart.AddDaughter(caloMet);
    totalPart.AddDaughter(CleanPhotons->At(0));
    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(CleanPhotons->At(0));
    hDGLSel[ 0+nHisto]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
    if(leptons->At(0)->Pt() > 20){
      hDGLSel[ 1+nHisto]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      if(caloMet->Pt() > 20){
        hDGLSel[ 2+nHisto]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
        if(mTW > 50){
          hDGLSel[ 3+nHisto]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
	  if(leptons->At(0)->ObjType() == kMuon ||
	     TMath::Abs(dilepton.Mass()-91.1886) > 10.0){
            hDGLSel[ 4+nHisto]->Fill(TMath::Min(deltaRlG,4.999),NNLOWeight->GetVal());
	    if(deltaRlG > 0.7){
	      hDGLSel[ 5+nHisto]->Fill(TMath::Min(CleanPhotons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
              hDGLSel[ 6+nHisto]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
              hDGLSel[ 7+nHisto]->Fill(nCentralJets,NNLOWeight->GetVal());
              hDGLSel[ 8+nHisto]->Fill(TMath::Min(totalPart.Pt(),199.999),NNLOWeight->GetVal());
	      hDGLSel[ 9+nHisto]->Fill(TMath::Min(CleanPhotons->At(0)->AbsEta(),2.499),NNLOWeight->GetVal());
	    } // DeltaRLG cut
	  } // mass cut
        } // mtw cut
      } // met cut
    } // pt lepton cut
  }

  //  1 gamma + 2 leptons
  if(CleanPhotons->GetEntries() >= 1 && leptons->GetEntries() == 2 &&
     leptons->At(0)->Pt() > 10 && leptons->At(1)->Pt() > 10){
    CompositeParticle dilepton;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));
    CompositeParticle triparticle;
    triparticle.AddDaughter(leptons->At(0));
    triparticle.AddDaughter(leptons->At(1));
    triparticle.AddDaughter(CleanPhotons->At(0));
    CompositeParticle GLSystem0;
    GLSystem0.AddDaughter(leptons->At(0));
    GLSystem0.AddDaughter(CleanPhotons->At(0));
    CompositeParticle GLSystem1;
    GLSystem1.AddDaughter(leptons->At(1));
    GLSystem1.AddDaughter(CleanPhotons->At(0));
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

    if(dilepton.Charge() == 0 &&
       (pairType == 0 || pairType == 1)){
      int nHisto = 0;
      if(pairType == 1) nHisto = 40;
      //if(isGoodMCPhoton[0] == kFALSE) nHisto += 20;
      double deltaRlG[2] = {
         MathUtils::DeltaR(CleanPhotons->At(0)->Phi(), CleanPhotons->At(0)->Eta(),
    			   leptons->At(0)->Phi(), leptons->At(0)->Eta()),
         MathUtils::DeltaR(CleanPhotons->At(0)->Phi(), CleanPhotons->At(0)->Eta(),
    			   leptons->At(1)->Phi(), leptons->At(1)->Eta())};
      hDGLLSel[ 0+nHisto]->Fill(TMath::Min(TMath::Min(deltaRlG[0],deltaRlG[1]),4.999),NNLOWeight->GetVal());
      if(TMath::Min(deltaRlG[0],deltaRlG[1]) > 0.4){
        hDGLLSel[ 1+nHisto]->Fill(TMath::Min(TMath::Max(deltaRlG[0],deltaRlG[1]),4.999),NNLOWeight->GetVal());
        hDGLLSel[ 2+nHisto]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 3+nHisto]->Fill(TMath::Min(triparticle.Mass(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 4+nHisto]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 5+nHisto]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 6+nHisto]->Fill(TMath::Min(CleanPhotons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 7+nHisto]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
        hDGLLSel[ 8+nHisto]->Fill(nCentralJets,NNLOWeight->GetVal());
        hDGLLSel[ 9+nHisto]->Fill(TMath::Min(CleanPhotons->At(0)->AbsEta(),2.499),NNLOWeight->GetVal());
        hDGLLSel[10+nHisto]->Fill(TMath::Min(GLSystem0.Mass(),499.999),NNLOWeight->GetVal());
        hDGLLSel[10+nHisto]->Fill(TMath::Min(GLSystem1.Mass(),499.999),NNLOWeight->GetVal());
      }
    }
  }

  //  1 gamma + jets
  if(CleanPhotons->GetEntries() >= 1){
    //int nHisto = 0;
    //if(isGoodMCPhoton[0] == kFALSE) nHisto = 20;

    JetOArr *NoGammaJets = new JetOArr;
    NoGammaJets->SetName("NoGammaJets");

    for (UInt_t j=0; j<GoodJets->GetEntries(); ++j) {
      Jet *jet = GoodJets->At(j);
      double deltaRJG = MathUtils::DeltaR(CleanPhotons->At(0)->Phi(), CleanPhotons->At(0)->Eta(),
     	                                  jet->Phi(), jet->Eta());
      hDGJetSel[0]->Fill(TMath::Min(deltaRJG,9.999),NNLOWeight->GetVal());
      if(deltaRJG > 0.3){
        NoGammaJets->Add(jet);
      }
    }
    NoGammaJets->Sort();
    hDGJetSel[1]->Fill(TMath::Min((double)NoGammaJets->GetEntries(),19.499),NNLOWeight->GetVal());
    if(NoGammaJets->GetEntries() >= 1){
      Jet *jet = NoGammaJets->At(0);
      double deltaPhi = fabs(MathUtils::DeltaPhi(jet->Phi(), CleanPhotons->At(0)->Phi()))* 180./TMath::Pi();
      hDGJetSel[2]->Fill(deltaPhi,NNLOWeight->GetVal());
      if(deltaPhi > 160.){
        if(NoGammaJets->GetEntries() >= 2){
          hDGJetSel[3]->Fill(TMath::Min(NoGammaJets->At(1)->Pt()/CleanPhotons->At(0)->Pt(),1.999),NNLOWeight->GetVal());
          hDGJetSel[4]->Fill(TMath::Min(NoGammaJets->At(1)->Pt()/NoGammaJets->At(0)->Pt(),0.999),NNLOWeight->GetVal());
        }
        hDGJetSel[5]->Fill(jet->AbsEta(),NNLOWeight->GetVal());
        hDGJetSel[6]->Fill(TMath::Min(jet->Pt(), 199.999),NNLOWeight->GetVal());
        hDGJetSel[7]->Fill(TMath::Min(CleanPhotons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        hDGJetSel[8]->Fill((jet->Pt()-CleanPhotons->At(0)->Pt())/CleanPhotons->At(0)->Pt(),NNLOWeight->GetVal());
        hDGJetSel[9]->Fill(TMath::Min(caloMet->Pt(), 199.999),NNLOWeight->GetVal());
        hDGJetSel[10]->Fill(TMath::Min(CleanPhotons->At(0)->AbsEta(),2.499),NNLOWeight->GetVal());
      }
    }
    delete NoGammaJets;
  }

}
//--------------------------------------------------------------------------------------------------
void GammaXEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fPhotonBranchName, fPhotons);
  char sb[200];
  // G MC
  sprintf(sb,"hDGMC_%d",0);  hDGMC[0]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDGMC_%d",1);  hDGMC[1]  = new TH1D(sb,sb,100,0,2.5);
  sprintf(sb,"hDGMC_%d",2);  hDGMC[2]  = new TH1D(sb,sb,200,0,200);
  sprintf(sb,"hDGMC_%d",3);  hDGMC[3]  = new TH1D(sb,sb,100,0,2.5);

  for(int i=0; i<4; i++){
    AddOutput(hDGMC[i]);
  }

  // G selection
  for(int j=0; j<2; j++){
    int ind = 40 * j;
    sprintf(sb,"hDGSel_%d",ind+0);  hDGSel[ind+0]  = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDGSel_%d",ind+1);  hDGSel[ind+1]  = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDGSel_%d",ind+2);  hDGSel[ind+2]  = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDGSel_%d",ind+3);  hDGSel[ind+3]  = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDGSel_%d",ind+4);  hDGSel[ind+4]  = new TH1D(sb,sb,2,-0.5,1.5);  
    sprintf(sb,"hDGSel_%d",ind+5);  hDGSel[ind+5]  = new TH1D(sb,sb,100,0.0,0.1); 
    sprintf(sb,"hDGSel_%d",ind+6);  hDGSel[ind+6]  = new TH1D(sb,sb,50,0.5,1.0);
    sprintf(sb,"hDGSel_%d",ind+7);  hDGSel[ind+7]  = new TH1D(sb,sb,50,0.5,1.0);
    sprintf(sb,"hDGSel_%d",ind+8);  hDGSel[ind+8]  = new TH1D(sb,sb,50,0.5,1.0);
    sprintf(sb,"hDGSel_%d",ind+9);  hDGSel[ind+9]  = new TH1D(sb,sb,100,0.0,1.0); 
    sprintf(sb,"hDGSel_%d",ind+10); hDGSel[ind+10] = new TH1D(sb,sb,100,0.0,1.0); 
    sprintf(sb,"hDGSel_%d",ind+11); hDGSel[ind+11] = new TH1D(sb,sb,100,0.0,1.0); 
    sprintf(sb,"hDGSel_%d",ind+12); hDGSel[ind+12] = new TH1D(sb,sb,100,0.0,1.0); 
    sprintf(sb,"hDGSel_%d",ind+13); hDGSel[ind+13] = new TH1D(sb,sb,6,-0.5,5.5); 
    sprintf(sb,"hDGSel_%d",ind+14); hDGSel[ind+14] = new TH1D(sb,sb,100,0.00,0.05); 
    sprintf(sb,"hDGSel_%d",ind+15); hDGSel[ind+15] = new TH1D(sb,sb,100,0.00,0.05); 
    sprintf(sb,"hDGSel_%d",ind+16); hDGSel[ind+16] = new TH1D(sb,sb,200,0.0,5.0); 
    sprintf(sb,"hDGSel_%d",ind+17); hDGSel[ind+17] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDGSel_%d",ind+18); hDGSel[ind+18] = new TH1D(sb,sb,200,0.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+19); hDGSel[ind+19] = new TH1D(sb,sb,210,-2.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+20); hDGSel[ind+20] = new TH1D(sb,sb,200,0.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+21); hDGSel[ind+21] = new TH1D(sb,sb,200,0.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+22); hDGSel[ind+22] = new TH1D(sb,sb,200,0.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+23); hDGSel[ind+23] = new TH1D(sb,sb,200,0.0,40.); 
    sprintf(sb,"hDGSel_%d",ind+24); hDGSel[ind+24] = new TH1D(sb,sb,20,-0.5,19.5); 
    sprintf(sb,"hDGSel_%d",ind+25); hDGSel[ind+25] = new TH1D(sb,sb,20,-0.5,19.5);
    sprintf(sb,"hDGSel_%d",ind+26); hDGSel[ind+26] = new TH1D(sb,sb,2,-0.5,1.5); 
    sprintf(sb,"hDGSel_%d",ind+27); hDGSel[ind+27] = new TH1D(sb,sb,2,-0.5,1.5); 
    sprintf(sb,"hDGSel_%d",ind+28); hDGSel[ind+28] = new TH1D(sb,sb,2,-0.5,1.5); 
    sprintf(sb,"hDGSel_%d",ind+29); hDGSel[ind+29] = new TH1D(sb,sb,2,-0.5,1.5); 
    sprintf(sb,"hDGSel_%d",ind+30); hDGSel[ind+30] = new TH1D(sb,sb,2,-0.5,1.5); 
  }

  for(int i=0; i<31; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDGSel[i+j*40]);
    }
  }

  sprintf(sb,"hDGSel2D_%d",0);   hDGSel2D[0] = new TH2D(sb,sb,50,0.0,20.0,50,0,100);
  sprintf(sb,"hDGSel2D_%d",1);   hDGSel2D[1] = new TH2D(sb,sb,50,0.0,20.0,50,0,100); 
  for(int i=0; i<2; i++){
    AddOutput(hDGSel2D[i]);
  }

  // GG
  for(int j=0; j<2; j++){
    int ind = 20 * j;
    sprintf(sb,"hDGGSel_%d",ind+0);  hDGGSel[ind+0]  = new TH1D(sb,sb,600,0,600);
    sprintf(sb,"hDGGSel_%d",ind+1);  hDGGSel[ind+1]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGGSel_%d",ind+2);  hDGGSel[ind+2]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGGSel_%d",ind+3);  hDGGSel[ind+3]  = new TH1D(sb,sb,200,0,200);
  }

  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDGGSel[i+j*20]);
    }
  }

  // GL
  for(int j=0; j<4; j++){
    int ind = 20 * j;
    sprintf(sb,"hDGLSel_%d",ind+0);  hDGLSel[ind+0]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+1);  hDGLSel[ind+1]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+2);  hDGLSel[ind+2]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+3);  hDGLSel[ind+3]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+4);  hDGLSel[ind+4]  = new TH1D(sb,sb,100,0,5);
    sprintf(sb,"hDGLSel_%d",ind+5);  hDGLSel[ind+5]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+6);  hDGLSel[ind+6]  = new TH1D(sb,sb,180,0,180);
    sprintf(sb,"hDGLSel_%d",ind+7);  hDGLSel[ind+7]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDGLSel_%d",ind+8);  hDGLSel[ind+8]  = new TH1D(sb,sb,200,0,200);
    sprintf(sb,"hDGLSel_%d",ind+9);  hDGLSel[ind+9]  = new TH1D(sb,sb,100,0.0,2.5);
  }

  for(int i=0; i<10; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDGLSel[i+j*20]);
    }
  }

  // GLL
  for(int j=0; j<4; j++){
    int ind = 20 * j;
    sprintf(sb,"hDGLLSel_%d",ind+0); hDGLLSel[ind+0]  = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGLLSel_%d",ind+1); hDGLLSel[ind+1]  = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGLLSel_%d",ind+2); hDGLLSel[ind+2]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+3); hDGLLSel[ind+3]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+4); hDGLLSel[ind+4]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+5); hDGLLSel[ind+5]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+6); hDGLLSel[ind+6]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+7); hDGLLSel[ind+7]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGLLSel_%d",ind+8); hDGLLSel[ind+8]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDGLLSel_%d",ind+9); hDGLLSel[ind+9]  = new TH1D(sb,sb,100,0.0,2.5);
    sprintf(sb,"hDGLLSel_%d",ind+10);hDGLLSel[ind+10] = new TH1D(sb,sb,100,0.0,500); 
  }

  for(int i=0; i<11; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDGLLSel[i+j*20]);
    }
  }

  // G + jets
  for(int j=0; j<2; j++){
    int ind = 20 * j;
    sprintf(sb,"hDGJetSel_%d",ind+0); hDGJetSel[ind+0]  = new TH1D(sb,sb,100,0,10); 
    sprintf(sb,"hDGJetSel_%d",ind+1); hDGJetSel[ind+1]  = new TH1D(sb,sb,20,-0.5,19.5); 
    sprintf(sb,"hDGJetSel_%d",ind+2); hDGJetSel[ind+2]  = new TH1D(sb,sb,90,0,180); 
    sprintf(sb,"hDGJetSel_%d",ind+3); hDGJetSel[ind+3]  = new TH1D(sb,sb,100,0,2); 
    sprintf(sb,"hDGJetSel_%d",ind+4); hDGJetSel[ind+4]  = new TH1D(sb,sb,100,0,1); 
    sprintf(sb,"hDGJetSel_%d",ind+5); hDGJetSel[ind+5]  = new TH1D(sb,sb,100,0.,5); 
    sprintf(sb,"hDGJetSel_%d",ind+6); hDGJetSel[ind+6]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGJetSel_%d",ind+7); hDGJetSel[ind+7]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGJetSel_%d",ind+8); hDGJetSel[ind+8]  = new TH1D(sb,sb,200,-2.0,2.0); 
    sprintf(sb,"hDGJetSel_%d",ind+9); hDGJetSel[ind+9]  = new TH1D(sb,sb,200,0.0,200); 
    sprintf(sb,"hDGJetSel_%d",ind+10);hDGJetSel[ind+10] = new TH1D(sb,sb,100,0.0,2.5); 
  }

  for(int i=0; i<11; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDGJetSel[i+j*20]);
    }
  }

  // Before HLT plot
  sprintf(sb,"hDBeforeHLTSel_%d",0); hDBeforeHLTSel[0]  = new TH1D(sb,sb,10,-0.5,9.5); 

  AddOutput(hDBeforeHLTSel[0]);
}
//--------------------------------------------------------------------------------------------------
void GammaXEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void GammaXEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
