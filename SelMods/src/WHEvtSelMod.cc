// $Id: WHEvtSelMod.cc,v 1.8 2012/01/26 13:17:57 ceballos Exp $

#include "Analysis/SelMods/interface/WHEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"

using namespace mithep;
ClassImp(mithep::WHEvtSelMod)

//-----------------------------------------------------------------------------
WHEvtSelMod::WHEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fCaloJetName0("AKt5Jets"),
  fCaloJet0(0),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fIsData(kFALSE),
  fNEventsProcessed(0)
{
  // Constructor.
}

//-----------------------------------------------------------------------------
void WHEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//-----------------------------------------------------------------------------
void WHEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if(fIsData == kFALSE) LoadBranch(fMCEvInfoName);

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);
  MCParticleOArr *GenLeptons   = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  hDWHSel[0]->Fill((double)leptons->GetEntries(),NNLOWeight->GetVal());

  if(leptons->GetEntries() >= 2){
    // Sort and count the number of central Jets for vetoing
    LoadBranch(fCaloJetName0);
    double sumJetPt = 0.0;
    double btagMax = 0.0;
    double theDijetMass = 999999999.0;
    if(CleanJets->GetEntries() <= 1) theDijetMass = 0.001;
    vector<Jet*> sortedJets;
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(i<CleanJets->GetEntries() - 1){
        for(UInt_t j=i+1; j<sortedJets.size(); j++){
          CompositeParticle dijet;
          dijet.AddDaughter(sortedJets[i]);
          dijet.AddDaughter(sortedJets[j]);
          if(TMath::Abs(dijet.Mass() - 80.40) < TMath::Abs(theDijetMass - 80.40)){
            theDijetMass = dijet.Mass();
          }
	}
      }

      if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	 CleanJets->At(i)->Pt() > fPtJetCut){
        nCentralJets++;
    	Jet* jet_f = new Jet(CleanJets->At(i)->Px(),
    			     CleanJets->At(i)->Py(),
    			     CleanJets->At(i)->Pz(),
    			     CleanJets->At(i)->E() );
	sumJetPt = sumJetPt + CleanJets->At(i)->Pt();
        jet_f->SetMatchedMCFlavor(CleanJets->At(i)->MatchedMCFlavor());
        jet_f->SetCombinedSecondaryVertexBJetTagsDisc(CleanJets->At(i)->CombinedSecondaryVertexBJetTagsDisc());
        jet_f->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJets->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
        jet_f->SetJetProbabilityBJetTagsDisc(CleanJets->At(i)->JetProbabilityBJetTagsDisc());
        jet_f->SetJetBProbabilityBJetTagsDisc(CleanJets->At(i)->JetBProbabilityBJetTagsDisc());
        jet_f->SetTrackCountingHighEffBJetTagsDisc(CleanJets->At(i)->TrackCountingHighEffBJetTagsDisc());
	if(btagMax < jet_f->TrackCountingHighEffBJetTagsDisc()) btagMax = jet_f->TrackCountingHighEffBJetTagsDisc();
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

    if     (leptons->GetEntries() == 2){ // WH -> l+/- l+/- X
      CompositeParticle dilepton;
      dilepton.AddDaughter(leptons->At(0));
      dilepton.AddDaughter(leptons->At(1));

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

      double deltaPhiln = 
          fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), caloMet->Phi()));
      if(deltaPhiln > fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(), caloMet->Phi())))
          deltaPhiln = fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(), caloMet->Phi()));
      for(UInt_t in=0; in<sortedJets.size(); in++){
        if(deltaPhiln > fabs(MathUtils::DeltaPhi(sortedJets[in]->Phi(), caloMet->Phi())))
            deltaPhiln = fabs(MathUtils::DeltaPhi(sortedJets[in]->Phi(), caloMet->Phi()));
      }
      double METdeltaPhilEt = caloMet->Pt();
      if(deltaPhiln < TMath::Pi()/2.) 
        METdeltaPhilEt = METdeltaPhilEt * sin(deltaPhiln);

      double deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(caloMet->Phi(), 
                                   dilepton.Phi()));

      hDWH2lXSel[ 0+100*pairType]->Fill(dilepton.Charge(),NNLOWeight->GetVal());
      if(dilepton.Charge() != 0 &&
         (pairType != 1 || TMath::Abs(dilepton.Mass()-91.1876) >= 10.0)){
        hDWH2lXSel[ 1+100*pairType]->Fill(sortedJets.size(),NNLOWeight->GetVal());
        hDWH2lXSel[ 2+100*pairType]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
        hDWH2lXSel[ 3+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
        hDWH2lXSel[ 4+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
        hDWH2lXSel[ 5+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
        hDWH2lXSel[ 6+100*pairType]->Fill(deltaPhiln * 180./TMath::Pi(),NNLOWeight->GetVal());
        hDWH2lXSel[ 7+100*pairType]->Fill(deltaPhiDileptonMet * 180./TMath::Pi(),NNLOWeight->GetVal());
        hDWH2lXSel[ 8+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
      }
    }
    else if(leptons->GetEntries() == 3){ // WH -> 3l3nu
      CompositeParticle dilepton01;
      dilepton01.AddDaughter(leptons->At(0));
      dilepton01.AddDaughter(leptons->At(1));
      CompositeParticle dilepton02;
      dilepton02.AddDaughter(leptons->At(0));
      dilepton02.AddDaughter(leptons->At(2));
      CompositeParticle dilepton12;
      dilepton12.AddDaughter(leptons->At(1));
      dilepton12.AddDaughter(leptons->At(2));
      CompositeParticle trilepton;
      trilepton.AddDaughter(leptons->At(0));
      trilepton.AddDaughter(leptons->At(1));
      trilepton.AddDaughter(leptons->At(2));
      int pairType = -1;
      if (CleanMuons->GetEntries() == 3)
        pairType = 0;
      else if(CleanElectrons->GetEntries() == 3)
        pairType = 1;
      else if(CleanMuons->GetEntries() == 2)
        pairType = 2;
      else if(CleanElectrons->GetEntries() == 2)
        pairType = 3;
      else {
        cout << "Hey, this is not possible" << endl;
      }

      double deltaPhiln = 
          fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), caloMet->Phi()));
      if(deltaPhiln > fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(), caloMet->Phi())))
          deltaPhiln = fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(), caloMet->Phi()));
      if(deltaPhiln > fabs(MathUtils::DeltaPhi(leptons->At(2)->Phi(), caloMet->Phi())))
          deltaPhiln = fabs(MathUtils::DeltaPhi(leptons->At(2)->Phi(), caloMet->Phi()));
      for(UInt_t i=0; i<sortedJets.size(); i++){
        if(deltaPhiln > fabs(MathUtils::DeltaPhi(sortedJets[i]->Phi(), caloMet->Phi())))
            deltaPhiln = fabs(MathUtils::DeltaPhi(sortedJets[i]->Phi(), caloMet->Phi()));
      }
      double METdeltaPhilEt = caloMet->Pt();
      if(deltaPhiln < TMath::Pi()/2.) 
        METdeltaPhilEt = METdeltaPhilEt * sin(deltaPhiln);

      hDWH3lSel[ 0+100*pairType]->Fill(trilepton.Charge(),NNLOWeight->GetVal());
      
      if(TMath::Abs(trilepton.Charge()) == 1){
	bool isZ = false;
	if(leptons->At(0)->ObjType() == leptons->At(1)->ObjType() &&
           dilepton01.Charge() == 0){
          hDWH3lSel[ 1+100*pairType]->Fill(TMath::Min(dilepton01.Mass(),199.999),NNLOWeight->GetVal());
	  if(TMath::Abs(dilepton01.Mass()-91.1876) < 15.0) isZ = true;
	}
	if(leptons->At(0)->ObjType() == leptons->At(2)->ObjType() &&
           dilepton02.Charge() == 0){
          hDWH3lSel[ 1+100*pairType]->Fill(TMath::Min(dilepton02.Mass(),199.999),NNLOWeight->GetVal());
	  if(TMath::Abs(dilepton02.Mass()-91.1876) < 15.0) isZ = true;
	}
	if(leptons->At(1)->ObjType() == leptons->At(2)->ObjType() &&
           dilepton12.Charge() == 0){
          hDWH3lSel[ 1+100*pairType]->Fill(TMath::Min(dilepton12.Mass(),199.999),NNLOWeight->GetVal());
	  if(TMath::Abs(dilepton12.Mass()-91.1876) < 15.0) isZ = true;
	}

	if(leptons->At(0)->ObjType() == leptons->At(1)->ObjType() &&
           dilepton01.Charge() != 0){
          hDWH3lSel[ 2+100*pairType]->Fill(TMath::Min(dilepton01.Mass(),199.999),NNLOWeight->GetVal());
	  if(leptons->At(0)->ObjType() == kElectron &&
	     TMath::Abs(dilepton01.Mass()-91.1876) < 15.0) isZ = true;
	}
	if(leptons->At(0)->ObjType() == leptons->At(2)->ObjType() &&
           dilepton02.Charge() != 0){
          hDWH3lSel[ 2+100*pairType]->Fill(TMath::Min(dilepton02.Mass(),199.999),NNLOWeight->GetVal());
	  if(leptons->At(0)->ObjType() == kElectron &&
	     TMath::Abs(dilepton02.Mass()-91.1876) < 15.0) isZ = true;
	}
	if(leptons->At(1)->ObjType() == leptons->At(2)->ObjType() &&
           dilepton12.Charge() != 0){
          hDWH3lSel[ 2+100*pairType]->Fill(TMath::Min(dilepton12.Mass(),199.999),NNLOWeight->GetVal());
	  if(leptons->At(2)->ObjType() == kElectron &&
	     TMath::Abs(dilepton12.Mass()-91.1876) < 15.0) isZ = true;
	}
        hDWH3lSel[ 3+100*pairType]->Fill((double)isZ,NNLOWeight->GetVal());

	if(isZ == false){
          hDWH3lSel[ 4+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
	  if(caloMet->Pt() > 40.0){
            hDWH3lSel[ 5+100*pairType]->Fill(sortedJets.size(),NNLOWeight->GetVal());
            hDWH3lSel[ 6+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
            hDWH3lSel[ 7+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
            hDWH3lSel[ 8+100*pairType]->Fill(TMath::Min(leptons->At(2)->Pt(),199.999),NNLOWeight->GetVal());
            hDWH3lSel[ 9+100*pairType]->Fill(deltaPhiln * 180./TMath::Pi(),NNLOWeight->GetVal());
            hDWH3lSel[10+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
	  } // MET > 20
	} // isZ == false
      } // abs(qtot) == 1
    } // Nlep == 3

    // same-sign analysis
    double massCloseZ = 999999.;
    int candSS[2] = {-1, -1};
    for(UInt_t i=0; i<leptons->GetEntries()-1; i++){
      for(UInt_t j=i+1; j<leptons->GetEntries(); j++){
        if(leptons->At(i)->Charge() == leptons->At(j)->Charge() &&
	   candSS[0] == -1){
	  candSS[0] = i;
	  candSS[1] = j;
	}
	else if(leptons->At(i)->Charge()  != leptons->At(j)->Charge() &&
	        leptons->At(i)->ObjType() == leptons->At(j)->ObjType()){
          CompositeParticle dileptonZ;
          dileptonZ.AddDaughter(leptons->At(i));
          dileptonZ.AddDaughter(leptons->At(j));
	  if(TMath::Abs(dileptonZ.Mass()-91.1876) < TMath::Abs(massCloseZ-91.1876))
	    massCloseZ = dileptonZ.Mass();
        }
      }
    }
    if(candSS[0] != -1){
      int pairType = -1;
      if (leptons->At(candSS[0])->ObjType() == kMuon && leptons->At(candSS[1])->ObjType() == kMuon )
	pairType = 0;
      else if(leptons->At(candSS[0])->ObjType() == kElectron && leptons->At(candSS[1])->ObjType() == kElectron)
	pairType = 1;
      else if((leptons->At(candSS[0])->ObjType() == kElectron && leptons->At(candSS[1])->ObjType() == kMuon) || 
              (leptons->At(candSS[0])->ObjType() == kMuon && leptons->At(candSS[1])->ObjType() == kElectron))
	pairType = 2;
      else {
	cout << "Hey, this is not possible, leptonTypes: "
    	     << leptons->At(candSS[0])->ObjType() << " - " 
             << leptons->At(candSS[1])->ObjType() << endl;
      }
      CompositeParticle dilepton;
      dilepton.AddDaughter(leptons->At(candSS[0]));
      dilepton.AddDaughter(leptons->At(candSS[1]));

      if(caloMet->Pt() > 30.0 && leptons->At(candSS[0])->Pt() > 20.0){
        hDWHSSSel[ 0+100*pairType]->Fill(leptons->GetEntries(),NNLOWeight->GetVal());
        hDWHSSSel[ 1+100*pairType]->Fill(TMath::Min(dilepton.Mass(),399.999),NNLOWeight->GetVal());
        hDWHSSSel[ 2+100*pairType]->Fill(TMath::Min(massCloseZ,399.999),NNLOWeight->GetVal());
        int theProcess = 0;
	if(fIsData == kFALSE){
	  if     (fMCEventInfo->ProcessId() ==  24) theProcess = 1;
	  else if(fMCEventInfo->ProcessId() ==  26) theProcess = 2;
	  else if(fMCEventInfo->ProcessId() == 121 ||
	          fMCEventInfo->ProcessId() == 122) theProcess = 3;
	}
	hDWHSSSel[ 3+100*pairType]->Fill(theProcess,NNLOWeight->GetVal());
        if(dilepton.Mass() > 12.0 && 
           TMath::Abs(massCloseZ-91.1876) > 15.0 &&
	   (pairType != 1 || TMath::Abs(dilepton.Mass()-91.1876) > 15.0)){
	  UInt_t leptonGenType[2] = {0, 0};
	  // look for real W/Z -> leptons
	  for(UInt_t i=0; i<2; ++i) {
	    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
              MCParticle *gen = GenLeptons->At(j);
              if(MathUtils::DeltaR(gen->Mom(), leptons->At(candSS[i])->Mom()) < 0.10 &&
	  	 leptonGenType[i] == 0) {
	    	leptonGenType[i] = 1;
	    	if(gen->Charge() != leptons->At(candSS[i])->Charge()) leptonGenType[i] = 2;
	      }
	    }
	  }
	  Int_t eventType = -1;
	  if	 (leptonGenType[0] == 0 &&  leptonGenType[1] == 0)   eventType = 0;
	  else if((leptonGenType[0] == 0 &&  leptonGenType[1] == 1) ||
	  	  (leptonGenType[0] == 1 &&  leptonGenType[1] == 0)) eventType = 1;
	  else if((leptonGenType[0] == 0 &&  leptonGenType[1] == 2) ||
	  	  (leptonGenType[0] == 2 &&  leptonGenType[1] == 0)) eventType = 2;
	  else if(leptonGenType[0] == 1 &&  leptonGenType[1] == 1)   eventType = 3;
 	  else if((leptonGenType[0] == 1 &&  leptonGenType[1] == 2) ||
	  	  (leptonGenType[0] == 2 &&  leptonGenType[1] == 1)) eventType = 4;
	  else if(leptonGenType[0] == 2 &&  leptonGenType[1] == 2)   eventType = 5;
	  else  						     eventType = 6;
          hDWHSSSel[ 4+100*pairType]->Fill(eventType,NNLOWeight->GetVal());
          hDWHSSSel[ 5+100*pairType]->Fill(TMath::Min((double)sortedJets.size(),9.499),NNLOWeight->GetVal());
          hDWHSSSel[ 6+100*pairType]->Fill(theProcess,NNLOWeight->GetVal());
          if(sortedJets.size() >= 1){
            double deltaPhill = 
              fabs(MathUtils::DeltaPhi(leptons->At(candSS[0])->Phi(), leptons->At(candSS[1])->Phi()));
	    hDWHSSSel[ 7+100*pairType]->Fill(TMath::Min(caloMet->Pt(),399.999),NNLOWeight->GetVal());
            hDWHSSSel[ 8+100*pairType]->Fill(TMath::Min(btagMax,0.999),NNLOWeight->GetVal());
            if(btagMax < 0.5) {
              hDWHSSSel[ 9+100*pairType]->Fill(TMath::Min(leptons->At(candSS[0])->Pt(),199.999),NNLOWeight->GetVal());
              hDWHSSSel[10+100*pairType]->Fill(TMath::Min(leptons->At(candSS[1])->Pt(),199.999),NNLOWeight->GetVal());
              hDWHSSSel[11+100*pairType]->Fill(TMath::Min(sortedJets[0]->Pt(),399.999),NNLOWeight->GetVal());
              hDWHSSSel[12+100*pairType]->Fill(TMath::Min(sumJetPt,799.999),NNLOWeight->GetVal());
              hDWHSSSel[13+100*pairType]->Fill(eventType,NNLOWeight->GetVal());
              hDWHSSSel[14+100*pairType]->Fill(deltaPhill * 180./TMath::Pi(),NNLOWeight->GetVal());
              hDWHSSSel[15+100*pairType]->Fill(leptons->GetEntries(),NNLOWeight->GetVal());
              hDWHSSSel[16+100*pairType]->Fill(sortedJets.size(),NNLOWeight->GetVal());
              hDWHSSSel[17+100*pairType]->Fill(theProcess,NNLOWeight->GetVal());	      
	    }
	    if(caloMet->Pt() > 100.0){
              hDWHSSSel[18+100*pairType]->Fill(TMath::Min(leptons->At(candSS[0])->Pt(),199.999),NNLOWeight->GetVal());
              hDWHSSSel[19+100*pairType]->Fill(TMath::Min(leptons->At(candSS[1])->Pt(),199.999),NNLOWeight->GetVal());
              hDWHSSSel[20+100*pairType]->Fill(TMath::Min(sortedJets[0]->Pt(),399.999),NNLOWeight->GetVal());
              hDWHSSSel[21+100*pairType]->Fill(TMath::Min(sumJetPt,799.999),NNLOWeight->GetVal());
              hDWHSSSel[22+100*pairType]->Fill(eventType,NNLOWeight->GetVal());
              hDWHSSSel[23+100*pairType]->Fill(TMath::Min(btagMax,0.999),NNLOWeight->GetVal());
              hDWHSSSel[24+100*pairType]->Fill(leptons->GetEntries(),NNLOWeight->GetVal());
              hDWHSSSel[25+100*pairType]->Fill(sortedJets.size(),NNLOWeight->GetVal());
	    } // MET > 100
          } // Njets >= 1
	} // Preselection
      } // MET > 30
    } // At least one same-sign lepton pair
    for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
  } // Nlep >= 2

}
//--------------------------------------------------------------------------------------------------
void WHEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fCaloJetName0, fCaloJet0);
  if(fIsData == kFALSE) ReqBranch(fMCEvInfoName, fMCEventInfo);

  char sb[200];
  sprintf(sb,"hDWHSel_0"); hDWHSel[0] = new TH1D(sb,sb,10,-0.5,9.5); 
  AddOutput(hDWHSel[0]);

  // WH->2l histograms
  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWH2lXSel_%d",ind+0);  hDWH2lXSel[ind+0]  = new TH1D(sb,sb,5,-2.5,2.5);
    sprintf(sb,"hDWH2lXSel_%d",ind+1);  hDWH2lXSel[ind+1]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWH2lXSel_%d",ind+2);  hDWH2lXSel[ind+2]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH2lXSel_%d",ind+3);  hDWH2lXSel[ind+3]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH2lXSel_%d",ind+4);  hDWH2lXSel[ind+4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH2lXSel_%d",ind+5);  hDWH2lXSel[ind+5]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDWH2lXSel_%d",ind+6);  hDWH2lXSel[ind+6]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDWH2lXSel_%d",ind+7);  hDWH2lXSel[ind+7]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDWH2lXSel_%d",ind+8);  hDWH2lXSel[ind+8]  = new TH1D(sb,sb,200,0.0,200.);
  }

  for(int i=0; i<9; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDWH2lXSel[i+j*100]);
    }
  }

  // WH->2l histograms
  for(int j=0; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWH3lSel_%d",ind+0);  hDWH3lSel[ind+0]  = new TH1D(sb,sb,7,-3.5,3.5);
    sprintf(sb,"hDWH3lSel_%d",ind+1);  hDWH3lSel[ind+1]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH3lSel_%d",ind+2);  hDWH3lSel[ind+2]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH3lSel_%d",ind+3);  hDWH3lSel[ind+3]  = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWH3lSel_%d",ind+4);  hDWH3lSel[ind+4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH3lSel_%d",ind+5);  hDWH3lSel[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWH3lSel_%d",ind+6);  hDWH3lSel[ind+6]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWH3lSel_%d",ind+7);  hDWH3lSel[ind+7]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDWH3lSel_%d",ind+8);  hDWH3lSel[ind+8]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDWH3lSel_%d",ind+9);  hDWH3lSel[ind+9]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDWH3lSel_%d",ind+10); hDWH3lSel[ind+10] = new TH1D(sb,sb,200,0.0,200.); 
  }

  for(int i=0; i<11; i++){
    for(int j=0; j<4; j++){
      AddOutput(hDWH3lSel[i+j*100]);
    }
  }

  // SS histograms
  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWHSSSel_%d",ind+0);  hDWHSSSel[ind+0]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWHSSSel_%d",ind+1);  hDWHSSSel[ind+1]  = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDWHSSSel_%d",ind+2);  hDWHSSSel[ind+2]  = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDWHSSSel_%d",ind+3);  hDWHSSSel[ind+3]  = new TH1D(sb,sb,4,-0.5,3.5);
    sprintf(sb,"hDWHSSSel_%d",ind+4);  hDWHSSSel[ind+4]  = new TH1D(sb,sb,7,-0.5,6.5);
    sprintf(sb,"hDWHSSSel_%d",ind+5);  hDWHSSSel[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWHSSSel_%d",ind+6);  hDWHSSSel[ind+6]  = new TH1D(sb,sb,4,-0.5,3.5);
    sprintf(sb,"hDWHSSSel_%d",ind+7);  hDWHSSSel[ind+7]  = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDWHSSSel_%d",ind+8);  hDWHSSSel[ind+8] = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDWHSSSel_%d",ind+9);  hDWHSSSel[ind+9]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDWHSSSel_%d",ind+10); hDWHSSSel[ind+10] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWHSSSel_%d",ind+11); hDWHSSSel[ind+11] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDWHSSSel_%d",ind+12); hDWHSSSel[ind+12] = new TH1D(sb,sb,200,0.0,800.);
    sprintf(sb,"hDWHSSSel_%d",ind+13); hDWHSSSel[ind+13] = new TH1D(sb,sb,7,-0.5,6.5);
    sprintf(sb,"hDWHSSSel_%d",ind+14); hDWHSSSel[ind+14] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDWHSSSel_%d",ind+15); hDWHSSSel[ind+15] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWHSSSel_%d",ind+16); hDWHSSSel[ind+16] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWHSSSel_%d",ind+17); hDWHSSSel[ind+17] = new TH1D(sb,sb,4,-0.5,3.5);
    sprintf(sb,"hDWHSSSel_%d",ind+18); hDWHSSSel[ind+18] = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDWHSSSel_%d",ind+19); hDWHSSSel[ind+19] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDWHSSSel_%d",ind+20); hDWHSSSel[ind+20] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDWHSSSel_%d",ind+21); hDWHSSSel[ind+21] = new TH1D(sb,sb,200,0.0,800.);
    sprintf(sb,"hDWHSSSel_%d",ind+22); hDWHSSSel[ind+22] = new TH1D(sb,sb,7,-0.5,6.5);
    sprintf(sb,"hDWHSSSel_%d",ind+23); hDWHSSSel[ind+23] = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDWHSSSel_%d",ind+24); hDWHSSSel[ind+24] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWHSSSel_%d",ind+25); hDWHSSSel[ind+25] = new TH1D(sb,sb,10,-0.5,9.5);
  }

  for(int i=0; i<26; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDWHSSSel[i+j*100]);
    }
  }

}

//--------------------------------------------------------------------------------------------------
void WHEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WHEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
