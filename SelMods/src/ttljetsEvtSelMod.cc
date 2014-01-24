// $Id: ttljetsEvtSelMod.cc,v 1.9 2012/01/26 13:17:57 ceballos Exp $

#include "Analysis/SelMods/interface/ttljetsEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitAna/DataTree/interface/ObjTypes.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"

using namespace mithep;
ClassImp(mithep::ttljetsEvtSelMod)

//--------------------------------------------------------------------------------------------------
ttljetsEvtSelMod::ttljetsEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fCaloJetName0("AKt5Jets"),
  fCaloJet0(0),
  fObjsName("randomMet"),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ttljetsEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ttljetsEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);
  if (!objs){
    printf("ZllEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  if(leptons->GetEntries() != 1) return;
  if(leptons->At(0)->Pt() <= 20.0) return;
  if(leptons->At(0)->ObjType() == kMuon && leptons->At(0)->AbsEta() >= 2.1) return;

  int pairType = -1;
  if	 (leptons->At(0)->ObjType() == kMuon	) pairType = 0;
  else if(leptons->At(0)->ObjType() == kElectron) pairType = 1;


  Int_t ents=objs->GetEntries();
  Bool_t passHLT = kFALSE;
  for(Int_t i=0;i<ents;++i) {
    const TriggerObject* to = objs->At(i);
    double DeltaRAux = MathUtils::DeltaR(leptons->At(0)->Mom(), to->Mom());
    if(DeltaRAux < 0.1) {
      passHLT = kTRUE;
      break;
    }
  }

  hDLJHLT[ 0+4*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
  hDLJHLT[ 1+4*pairType]->Fill(TMath::Max(TMath::Min(leptons->At(0)->Eta(),2.499),-2.499),NNLOWeight->GetVal());
  if(passHLT == kTRUE) {
    hDLJHLT[ 2+4*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
    hDLJHLT[ 3+4*pairType]->Fill(TMath::Max(TMath::Min(leptons->At(0)->Eta(),2.499),-2.499),NNLOWeight->GetVal());  
  }

  //if(passHLT == kFALSE) return;

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

  // Sort and count the number of central Jets for vetoing
  LoadBranch(fCaloJetName0);
  vector<Jet*> sortedJets;
  double sumPtJet = 0.0;
  for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
    if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
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
      sumPtJet = sumPtJet + CleanJets->At(i)->Pt();
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

  if(sortedJets.size() >= 1){
    // Angle between MET and closest lepton
    double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi()));

    double METdeltaPhilEt = caloMet->Pt();
    if(deltaPhiMetLepton <  TMath::Pi()/2.0) METdeltaPhilEt = METdeltaPhilEt * sin(deltaPhiMetLepton);

    double mTW = TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
    		            (1.0 - cos(deltaPhiMetLepton)));

    double HT = sumPtJet + leptons->At(0)->Pt() + caloMet->Pt();
    double btagMax[2]  = {-999.0, -999.0};
    UInt_t btagIMax[2] = {99999, 99999};
    for(UInt_t i=0; i<sortedJets.size(); i++){
      if(abs(sortedJets[i]->MatchedMCFlavor()) == 5){
    	hDLJBTag[0]->Fill(TMath::Min(TMath::Max(sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
    	hDLJBTag[1]->Fill(TMath::Min(TMath::Max(sortedJets[i]->CombinedSecondaryVertexMVABJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
    	hDLJBTag[2]->Fill(TMath::Min(TMath::Max(sortedJets[i]->JetProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
    	hDLJBTag[3]->Fill(TMath::Min(TMath::Max(sortedJets[i]->JetBProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
    	hDLJBTag[4]->Fill(TMath::Min(TMath::Max(sortedJets[i]->TrackCountingHighEffBJetTagsDisc(),-19.99999),19.9999),NNLOWeight->GetVal());
      }
      hDLJBTag[5]->Fill(TMath::Min(TMath::Max(sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
      hDLJBTag[6]->Fill(TMath::Min(TMath::Max(sortedJets[i]->CombinedSecondaryVertexMVABJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
      hDLJBTag[7]->Fill(TMath::Min(TMath::Max(sortedJets[i]->JetProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
      hDLJBTag[8]->Fill(TMath::Min(TMath::Max(sortedJets[i]->JetBProbabilityBJetTagsDisc(),0.00001),0.9999),NNLOWeight->GetVal());
      hDLJBTag[9]->Fill(TMath::Min(TMath::Max(sortedJets[i]->TrackCountingHighEffBJetTagsDisc(),-19.99999),19.9999),NNLOWeight->GetVal());

      if     (sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc() > btagMax[0]){
        btagMax[1]  = btagMax[0];
	btagIMax[1] = btagIMax[0];
        btagMax[0]  = sortedJets[i]->TrackCountingHighEffBJetTagsDisc();
	btagIMax[0] = i;	
      }
      else if(sortedJets[i]->CombinedSecondaryVertexBJetTagsDisc() > btagMax[1]){
        btagMax[1]  = sortedJets[i]->TrackCountingHighEffBJetTagsDisc();
	btagIMax[1] = i;	
      }
    }

    hDLJPresel[ 0+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
    hDLJPresel[ 1+100*pairType]->Fill(TMath::Min((double)sortedJets.size(),9.499),NNLOWeight->GetVal());
    hDLJPresel[ 2+100*pairType]->Fill(TMath::Min(TMath::Max(btagMax[0],0.00001),0.9999),NNLOWeight->GetVal());
    hDLJPresel[ 3+100*pairType]->Fill(TMath::Min(TMath::Max(btagMax[1],0.00001),0.9999),NNLOWeight->GetVal());
    hDLJPresel[ 4+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
    hDLJPresel[ 5+100*pairType]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
    hDLJPresel[ 6+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
    hDLJPresel[ 7+100*pairType]->Fill(deltaPhiMetLepton*180/TMath::Pi(),NNLOWeight->GetVal());
    hDLJPresel[ 8+100*pairType]->Fill(TMath::Min(sortedJets[0]->Pt(),199.999),NNLOWeight->GetVal());
    hDLJPresel[ 9+100*pairType]->Fill(TMath::Min(sumPtJet,799.999),NNLOWeight->GetVal());
    if    (sortedJets.size() == 1 && mTW > 50.0)
      hDLJPresel[10+100*pairType]->Fill(TMath::Min(HT,799.999),NNLOWeight->GetVal());
    else if(sortedJets.size() == 2 && mTW > 50.0)
      hDLJPresel[11+100*pairType]->Fill(TMath::Min(HT,799.999),NNLOWeight->GetVal());
    else if(sortedJets.size() == 3 && mTW > 50.0)
      hDLJPresel[12+100*pairType]->Fill(TMath::Min(HT,799.999),NNLOWeight->GetVal());
    else if(sortedJets.size() >= 4 && mTW > 50.0)
      hDLJPresel[13+100*pairType]->Fill(TMath::Min(HT,799.999),NNLOWeight->GetVal());

    Bool_t cuts[3] = {mTW > 50.0,  sortedJets.size() >= 4, btagMax[0] > 2.1};
    
    if(           cuts[1] && cuts[2]){
      hDLJSel[ 0+100*pairType]->Fill(TMath::Min(mTW,199.999),NNLOWeight->GetVal());
    }
    if(cuts[0]            && cuts[2]){
      hDLJSel[ 1+100*pairType]->Fill(TMath::Min((double)sortedJets.size(),9.499),NNLOWeight->GetVal());
    }
    if(cuts[0] && cuts[1]           ){
      hDLJSel[ 2+100*pairType]->Fill(TMath::Min(TMath::Max(btagMax[0],0.00001),0.9999),NNLOWeight->GetVal());
      hDLJSel[ 3+100*pairType]->Fill(TMath::Min(TMath::Max(btagMax[1],0.00001),0.9999),NNLOWeight->GetVal());
    }
    if(cuts[0] && cuts[1] && cuts[2]){
      hDLJSel[ 4+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
      hDLJSel[ 5+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
      hDLJSel[ 6+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999),NNLOWeight->GetVal());
      hDLJSel[ 7+100*pairType]->Fill(deltaPhiMetLepton*180/TMath::Pi(),NNLOWeight->GetVal());
      hDLJSel[ 8+100*pairType]->Fill(TMath::Min(sortedJets[0]->Pt(),199.999),NNLOWeight->GetVal());
      hDLJSel[ 9+100*pairType]->Fill(TMath::Min(sumPtJet,799.999),NNLOWeight->GetVal());
      hDLJSel[10+100*pairType]->Fill(TMath::Min(HT,799.999),NNLOWeight->GetVal());
      if(btagIMax[1] != 99999)
        hDLJSel[11+100*pairType]->Fill(TMath::Min(sortedJets[btagIMax[1]]->Pt(),199.999),NNLOWeight->GetVal());

      double dijetmass[3] = {0.0,  0.0, 999.0};
      UInt_t indexW[2] = {99999, 99999};
      for(UInt_t j0=0; j0<sortedJets.size()-1; j0++){
        if(j0 == btagIMax[0] || j0 == btagIMax[1]) continue;
        for(UInt_t j1=j0+1; j1<sortedJets.size(); j1++){
          if(j1 == btagIMax[0] || j1 == btagIMax[1]) continue;
      	  CompositeParticle dijet;
      	  dijet.AddDaughter(sortedJets[j0]);
      	  dijet.AddDaughter(sortedJets[j1]);
	  if(dijet.Mass() > 0                       && dijetmass[0] <= 0)    dijetmass[0] = dijet.Mass();
	  if(dijet.Mass() > 0 && dijet.Mass() < 200 && dijetmass[1] <= 0)    dijetmass[1] = dijet.Mass();
	  if(TMath::Abs(dijet.Mass()-80.40) < TMath::Abs(dijetmass[2]-80.40)) {
	    dijetmass[2] = dijet.Mass();
	    indexW[0] = j0; indexW[1] = j1;
	  }
	}
      }
      hDLJSel[12+100*pairType]->Fill(TMath::Min(dijetmass[0],199.999),NNLOWeight->GetVal());
      hDLJSel[13+100*pairType]->Fill(TMath::Min(dijetmass[1],199.999),NNLOWeight->GetVal());
      hDLJSel[14+100*pairType]->Fill(TMath::Min(dijetmass[2],199.999),NNLOWeight->GetVal());
      if(btagMax[1] > 2.1){
        hDLJSel[15+100*pairType]->Fill(TMath::Min(dijetmass[1],199.999),NNLOWeight->GetVal());
      }
      if(indexW[1] != 99999 && btagIMax[1] != 99999 && btagMax[1] > 2.1){
      	CompositeParticle *dijetW = new CompositeParticle();
      	dijetW->AddDaughter(sortedJets[indexW[0]]);
      	dijetW->AddDaughter(sortedJets[indexW[1]]);

      	CompositeParticle *topL0 = new CompositeParticle();
      	topL0->AddDaughter(sortedJets[btagIMax[0]]);
      	topL0->AddDaughter(leptons->At(0));
      	topL0->AddDaughter(caloMet);
      	CompositeParticle *topL1 = new CompositeParticle();
      	topL1->AddDaughter(sortedJets[btagIMax[1]]);
      	topL1->AddDaughter(leptons->At(0));
      	topL1->AddDaughter(caloMet);
      	CompositeParticle *topJJ0 = new CompositeParticle();
      	topJJ0->AddDaughter(sortedJets[btagIMax[0]]);
      	topJJ0->AddDaughter(dijetW);
      	CompositeParticle *topJJ1 = new CompositeParticle();
      	topJJ1->AddDaughter(sortedJets[btagIMax[1]]);
      	topJJ1->AddDaughter(dijetW);

	if(topL0->Mt() > topL1->Mt()){
          hDLJSel[16+100*pairType]->Fill(TMath::Min(topL0->Mt(),399.999),NNLOWeight->GetVal());
          hDLJSel[17+100*pairType]->Fill(TMath::Min(topL1->Mt(),399.999),NNLOWeight->GetVal());
          hDLJSel[18+100*pairType]->Fill(TMath::Min(topJJ1->Mass(),399.999),NNLOWeight->GetVal());
          hDLJSel[19+100*pairType]->Fill(TMath::Min(topJJ0->Mass(),399.999),NNLOWeight->GetVal());
	}
	else {
          hDLJSel[16+100*pairType]->Fill(TMath::Min(topL1->Mt(),399.999),NNLOWeight->GetVal());
          hDLJSel[17+100*pairType]->Fill(TMath::Min(topL0->Mt(),399.999),NNLOWeight->GetVal());
          hDLJSel[18+100*pairType]->Fill(TMath::Min(topJJ0->Mass(),399.999),NNLOWeight->GetVal());
          hDLJSel[19+100*pairType]->Fill(TMath::Min(topJJ1->Mass(),399.999),NNLOWeight->GetVal());
	}
        hDLJSel[20+100*pairType]->Fill(TMath::Min(topJJ0->Mass(),399.999),NNLOWeight->GetVal());
        hDLJSel[20+100*pairType]->Fill(TMath::Min(topJJ1->Mass(),399.999),NNLOWeight->GetVal());
	if(TMath::Abs(topJJ0->Mass()-173.0) < TMath::Abs(topJJ1->Mass()-173.0)){
	  hDLJSel[21+100*pairType]->Fill(TMath::Min(topJJ0->Mass(),399.999),NNLOWeight->GetVal());
	}
	else {
	  hDLJSel[21+100*pairType]->Fill(TMath::Min(topJJ1->Mass(),399.999),NNLOWeight->GetVal());
	}
        delete dijetW;
	delete topL0;
	delete topL1;
	delete topJJ0;
	delete topJJ1;
      }
    }
  } // Njets >= 1
  for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
}
//--------------------------------------------------------------------------------------------------
void ttljetsEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fCaloJetName0, fCaloJet0);

  char sb[200];
  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDLJPresel_%d",ind+ 0);  hDLJPresel[ind+ 0]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJPresel_%d",ind+ 1);  hDLJPresel[ind+ 1]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDLJPresel_%d",ind+ 2);  hDLJPresel[ind+ 2]  = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDLJPresel_%d",ind+ 3);  hDLJPresel[ind+ 3]  = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDLJPresel_%d",ind+ 4);  hDLJPresel[ind+ 4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJPresel_%d",ind+ 5);  hDLJPresel[ind+ 5]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJPresel_%d",ind+ 6);  hDLJPresel[ind+ 6]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJPresel_%d",ind+ 7);  hDLJPresel[ind+ 7]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDLJPresel_%d",ind+ 8);  hDLJPresel[ind+ 8]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJPresel_%d",ind+ 9);  hDLJPresel[ind+ 9]  = new TH1D(sb,sb,200,0.0,800.0); 
    sprintf(sb,"hDLJPresel_%d",ind+10);  hDLJPresel[ind+10]  = new TH1D(sb,sb,200,0.0,800.0);
    sprintf(sb,"hDLJPresel_%d",ind+11);  hDLJPresel[ind+11]  = new TH1D(sb,sb,200,0.0,800.0);
    sprintf(sb,"hDLJPresel_%d",ind+12);  hDLJPresel[ind+12]  = new TH1D(sb,sb,200,0.0,800.0);
    sprintf(sb,"hDLJPresel_%d",ind+13);  hDLJPresel[ind+13]  = new TH1D(sb,sb,200,0.0,800.0);
  }

  for(int i=0; i<14; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDLJPresel[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDLJSel_%d",ind+ 0);  hDLJSel[ind+ 0]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+ 1);  hDLJSel[ind+ 1]  = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDLJSel_%d",ind+ 2);  hDLJSel[ind+ 2]  = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDLJSel_%d",ind+ 3);  hDLJSel[ind+ 3]  = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDLJSel_%d",ind+ 4);  hDLJSel[ind+ 4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJSel_%d",ind+ 5);  hDLJSel[ind+ 5]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJSel_%d",ind+ 6);  hDLJSel[ind+ 6]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDLJSel_%d",ind+ 7);  hDLJSel[ind+ 7]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDLJSel_%d",ind+ 8);  hDLJSel[ind+ 8]  = new TH1D(sb,sb,200,0.0,200.0);
    sprintf(sb,"hDLJSel_%d",ind+ 9);  hDLJSel[ind+ 9]  = new TH1D(sb,sb,200,0.0,800.0);
    sprintf(sb,"hDLJSel_%d",ind+10);  hDLJSel[ind+10]  = new TH1D(sb,sb,200,0.0,800.0);
    sprintf(sb,"hDLJSel_%d",ind+11);  hDLJSel[ind+11]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+12);  hDLJSel[ind+12]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+13);  hDLJSel[ind+13]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+14);  hDLJSel[ind+14]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+15);  hDLJSel[ind+15]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJSel_%d",ind+16);  hDLJSel[ind+16]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDLJSel_%d",ind+17);  hDLJSel[ind+17]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDLJSel_%d",ind+18);  hDLJSel[ind+18]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDLJSel_%d",ind+19);  hDLJSel[ind+19]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDLJSel_%d",ind+20);  hDLJSel[ind+20]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDLJSel_%d",ind+21);  hDLJSel[ind+21]  = new TH1D(sb,sb,200,0.0,400.); 
  }

  for(int i=0; i<22; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDLJSel[i+j*100]);
    }
  }

  for(int j=0; j<2; j++){
    int ind = 4 * j;
    sprintf(sb,"hDLJHLT_%d",ind+0);  hDLJHLT[ind+0]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJHLT_%d",ind+1);  hDLJHLT[ind+1]  = new TH1D(sb,sb,100,-2.5,2.5);
    sprintf(sb,"hDLJHLT_%d",ind+2);  hDLJHLT[ind+2]  = new TH1D(sb,sb,200,0.0,200.); 
    sprintf(sb,"hDLJHLT_%d",ind+3);  hDLJHLT[ind+3]  = new TH1D(sb,sb,100,-2.5,2.5);
  }

  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDLJHLT[i+j*4]);
    }
  }

  sprintf(sb,"hDLJBTag_%d",0);  hDLJBTag[0]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",1);  hDLJBTag[1]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",2);  hDLJBTag[2]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",3);  hDLJBTag[3]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",4);  hDLJBTag[4]  = new TH1D(sb,sb,100,-20.0,20.0);
  sprintf(sb,"hDLJBTag_%d",5);  hDLJBTag[5]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",6);  hDLJBTag[6]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",7);  hDLJBTag[7]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",8);  hDLJBTag[8]  = new TH1D(sb,sb,100,0.0,1.);
  sprintf(sb,"hDLJBTag_%d",9);  hDLJBTag[9]  = new TH1D(sb,sb,100,-20.0,20.0);
  for(int i=0; i<10; i++){
    AddOutput(hDLJBTag[i]);
  }

}

//--------------------------------------------------------------------------------------------------
void ttljetsEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ttljetsEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
