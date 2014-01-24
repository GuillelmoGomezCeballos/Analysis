// $Id: QQLLEvtSelMod.cc,v 1.10 2012/01/26 13:17:57 ceballos Exp $

#include "Analysis/SelMods/interface/QQLLEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"

using namespace mithep;
ClassImp(mithep::QQLLEvtSelMod)

//--------------------------------------------------------------------------------------------------
QQLLEvtSelMod::QQLLEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsFastSim(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fCleanJetsName("randomJet"),
  fCaloJetName0("AKt5Jets"),
  fCaloJet0(0),
  fJetScaleSyst(0.0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void QQLLEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void QQLLEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons  = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons         = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  ParticleOArr *leptons         = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet             = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet            = CleanMet->At(0);
  JetOArr *CleanJets            = GetObjThisEvt<JetOArr>(fCleanJetsName);

  LoadBranch(fCaloJetName0);

  if(leptons->GetEntries() != 2) return;
  if(leptons->At(0)->Pt() <= 20) return;
  if(leptons->At(1)->Pt() <= 20) return;

  CompositeParticle *dilepton = new CompositeParticle();
  dilepton->AddDaughter(leptons->At(0));
  dilepton->AddDaughter(leptons->At(1));
  if(dilepton->Mass() > 12.0 && dilepton->Mass()-91.1876 < 20.0 && dilepton->Charge() == 0) {
    double zAverage = 0.0;
    double maxD0Sig = -10.0;
    double maxIsoM = -10.0;
    double maxIsoE = -10.0;

    fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
    // Make lepton vector from muons and electrons
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      zAverage = zAverage + CleanMuons->At(j)->BestTrk()->Z0();
      double isoAux = 1.0 * CleanMuons->At(j)->IsoR03SumPt() + 
                      1.0 * CleanMuons->At(j)->IsoR03EmEt()  +
		      1.0 * CleanMuons->At(j)->IsoR03HadEt();
      isoAux = isoAux/(isoAux+CleanMuons->At(j)->Pt());
      if(isoAux > maxIsoM) maxIsoM = isoAux;
      double d0RealSig = 99999;
      for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        if(fVertices->At(i0)->NTracks() > 0){
	  double pD0 = CleanMuons->At(j)->BestTrk()->D0Corrected(*fVertices->At(i0));
	  d0RealSig = TMath::Abs(pD0);
          break;
	}
      }
      if(d0RealSig > maxD0Sig) maxD0Sig = d0RealSig;
    }
    for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
      zAverage = zAverage + CleanElectrons->At(j)->BestTrk()->Z0();
      double isoAux = 1.0 * CleanElectrons->At(j)->TrackIsolationDr03() + 
                      1.0 * CleanElectrons->At(j)->EcalRecHitIsoDr04();
      isoAux = isoAux/(isoAux+CleanElectrons->At(j)->Pt());
      if(isoAux > maxIsoE) maxIsoE = isoAux;

      double d0RealSig = 99999;
      for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
	if(fVertices->At(i0)->NTracks() > 0){
	  double pD0 = CleanElectrons->At(j)->GsfTrk()->D0Corrected(*fVertices->At(i0));      
	  d0RealSig = TMath::Abs(pD0);
          break;
	}
      }
      if(d0RealSig > maxD0Sig) maxD0Sig = d0RealSig;
    }

    // Computing Z average (our primary vertex)
    if(leptons->GetEntries() > 0) zAverage = zAverage / leptons->GetEntries();

    int pairType = -1;
    if     (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon && dilepton->Mass()-91.1876 < -20.0)
      pairType = 0;
    else if(leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon)
      pairType = 1;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron && dilepton->Mass()-91.1876 < -20.0)
      pairType = 2;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron)
      pairType = 3;
    else if((leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon) || 
 	    (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kElectron))
      pairType = 4;
    else {
      cout << "Hey, this is not possible, leptonTypes: "
 	   << leptons->At(0)->ObjType() << " - " 
 	   << leptons->At(1)->ObjType() << endl;
      assert(0);
    }

    LoadBranch(fCaloJetName0);
    vector<Jet*> sortedJets;
    double sumJetPt[4] = {0.0, 0.0, 0.0, 0.0};
    double bTagMax = -5.0;
    for(UInt_t nj=0; nj<CleanJets->GetEntries(); nj++){
      const Jet *jet = CleanJets->At(nj);        
      if(TMath::Abs(jet->Eta()) >= fEtaJetCut ||
	 jet->Pt()*(1.0+fJetScaleSyst) <= fPtJetCut) continue;
      Jet* jet_f = new Jet(jet->Px()*(1.0+fJetScaleSyst),
        		   jet->Py()*(1.0+fJetScaleSyst),
        		   jet->Pz()*(1.0+fJetScaleSyst),
        		   jet->E()*(1.0+fJetScaleSyst) );
      jet_f->SetMatchedMCFlavor(jet->MatchedMCFlavor());
      jet_f->SetCombinedSecondaryVertexBJetTagsDisc(jet->CombinedSecondaryVertexBJetTagsDisc());
      jet_f->SetCombinedSecondaryVertexMVABJetTagsDisc(jet->CombinedSecondaryVertexMVABJetTagsDisc());
      jet_f->SetJetProbabilityBJetTagsDisc(jet->JetProbabilityBJetTagsDisc());
      jet_f->SetJetBProbabilityBJetTagsDisc(jet->JetBProbabilityBJetTagsDisc());
      jet_f->SetTrackCountingHighEffBJetTagsDisc(jet->TrackCountingHighEffBJetTagsDisc());
      sumJetPt[0] = sumJetPt[0] + jet_f->Px();
      sumJetPt[1] = sumJetPt[1] + jet_f->Py();
      sumJetPt[2] = sumJetPt[2] + jet_f->Pz();
      sumJetPt[3] = sumJetPt[3] + jet_f->P();
      if(bTagMax < jet_f->TrackCountingHighEffBJetTagsDisc()) bTagMax = jet_f->TrackCountingHighEffBJetTagsDisc();
      sortedJets.push_back(jet_f);
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

    hDQQLLSel[ 0+100*pairType]->Fill(TMath::Min((double)sortedJets.size(),9.499),NNLOWeight->GetVal());
    if(sortedJets.size() >= 2){
      double theMass = 999999999.0;
      int theJets[2] = {-1, -1};
      for(UInt_t i=0; i<sortedJets.size()-1; i++){
        for(UInt_t j=i+1; j<sortedJets.size(); j++){
          CompositeParticle *dijet = new CompositeParticle();
          dijet->AddDaughter(sortedJets[i]);
          dijet->AddDaughter(sortedJets[j]);
	  if(i == 0 && j == 1) hDQQLLSel[ 1+100*pairType]->Fill(TMath::Min(dijet->Mass(),399.999),NNLOWeight->GetVal());
	  if(TMath::Abs(dijet->Mass() - 91.1876) < TMath::Abs(theMass - 91.1876)){
            theMass = dijet->Mass();
	    theJets[0] = i;
	    theJets[1] = j;
          }
	  delete dijet;
	}
      }
      if(theJets[0] == -1 || theJets[0] == -1) {printf("Problem in theJets: %d - %d\n",theJets[0],theJets[1]);}

      double deltaPhiMetLepton[2] = {fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi())),
    				     fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(1)->Phi()))};
      double mTW[2] = {TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
        			   (1.0 - cos(deltaPhiMetLepton[0]))),
        	       TMath::Sqrt(2.0*leptons->At(1)->Pt()*caloMet->Pt()*
        			   (1.0 - cos(deltaPhiMetLepton[1])))};

      hDQQLLSel[ 2+100*pairType]->Fill(TMath::Min(theMass,399.999),NNLOWeight->GetVal());
      theJets[0] = 0;
      theJets[1] = 1;
      CompositeParticle *dijet = new CompositeParticle();
      dijet->AddDaughter(sortedJets[theJets[0]]);
      dijet->AddDaughter(sortedJets[theJets[1]]);
      CompositeParticle *higgs = new CompositeParticle();
      higgs->AddDaughter(dilepton);
      higgs->AddDaughter(dijet);

      if(fabs(MathUtils::DeltaPhi(sortedJets[theJets[0]]->Phi(), 
    	  		     sortedJets[theJets[1]]->Phi()))* 180./TMath::Pi() < 160.0 &&
	bTagMax < 2.1){    
          hDQQLLSel[16+100*pairType]->Fill(TMath::Min(dijet->Mass(),399.999),NNLOWeight->GetVal());
      }

      if(TMath::Abs(dijet->Mass()-91.1876) < 20){
	hDQQLLSel[ 6+100*pairType]->Fill(fabs(MathUtils::DeltaPhi(sortedJets[theJets[0]]->Phi(), 
        						     sortedJets[theJets[1]]->Phi()))* 180./TMath::Pi(),NNLOWeight->GetVal());

        hDQQLLSel[ 3+100*pairType]->Fill(TMath::Max(TMath::Min(bTagMax,4.999),-4.999),NNLOWeight->GetVal());
        if(bTagMax < 2.1){
          hDQQLLSel[21+100*pairType]->Fill(TMath::Min(higgs->Mass(),399.999),NNLOWeight->GetVal());
 	  //if(higgs->Mass() > 100.0 && higgs->Mass() < 180.0){
            hDQQLLSel[ 7+100*pairType]->Fill(TMath::Min(MathUtils::DeltaR(sortedJets[theJets[0]]->Mom(),
 						                          sortedJets[theJets[1]]->Mom()),9.999),NNLOWeight->GetVal());
            hDQQLLSel[ 4+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[10+100*pairType]->Fill(TMath::Min(sortedJets[theJets[0]]->Pt(),199.999),NNLOWeight->GetVal());
	    hDQQLLSel[15+100*pairType]->Fill(TMath::Min(dilepton->Mass(),399.999),NNLOWeight->GetVal());
            hDQQLLSel[20+100*pairType]->Fill(TMath::Min(higgs->Mass(),399.999),NNLOWeight->GetVal());
	    hDQQLLSel[ 5+100*pairType]->Fill(TMath::Min(sortedJets[theJets[1]]->Pt()/theMass,4.999),NNLOWeight->GetVal());
            hDQQLLSel[ 8+100*pairType]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
            hDQQLLSel[ 9+100*pairType]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999),NNLOWeight->GetVal());
            hDQQLLSel[11+100*pairType]->Fill(TMath::Min(sortedJets[theJets[1]]->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[12+100*pairType]->Fill(TMath::Sqrt(sumJetPt[0]*sumJetPt[0]+sumJetPt[1]*sumJetPt[1]+sumJetPt[2]*sumJetPt[2])/sumJetPt[3],NNLOWeight->GetVal());
            hDQQLLSel[13+100*pairType]->Fill(TMath::Min(leptons->At(0)->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[14+100*pairType]->Fill(TMath::Min(leptons->At(1)->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[17+100*pairType]->Fill(TMath::Max(sortedJets[theJets[0]]->AbsEta(), 
    	    						sortedJets[theJets[1]]->AbsEta()),NNLOWeight->GetVal());
            hDQQLLSel[18+100*pairType]->Fill(TMath::Min(sortedJets[theJets[0]]->AbsEta(), 
    	    						sortedJets[theJets[1]]->AbsEta()),NNLOWeight->GetVal());
            hDQQLLSel[19+100*pairType]->Fill(TMath::Min(MathUtils::DeltaR(dilepton->Mom(), 
 	    								  dijet->Mom()),9.999),NNLOWeight->GetVal());
	    hDQQLLSel[22+100*pairType]->Fill(TMath::Min(dilepton->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[23+100*pairType]->Fill(TMath::Min(dijet->Pt(),199.999),NNLOWeight->GetVal());
            hDQQLLSel[24+100*pairType]->Fill(TMath::Min(higgs->Pt(),199.999),NNLOWeight->GetVal());
          } // bTagMax cut
	//} // Higgs mass cut
      } // Mjj cut
      delete dijet;
      delete higgs;
    } // Njets >= 2
    for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
  } // Preselection
  delete dilepton ;
}
//--------------------------------------------------------------------------------------------------
void QQLLEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fCaloJetName0, fCaloJet0);

  char sb[200];
  for(int j=0; j<5; j++){
    int ind = 100 * j;
    sprintf(sb,"hDQQLLSel_%d",ind+ 0); hDQQLLSel[ind+ 0] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDQQLLSel_%d",ind+ 1); hDQQLLSel[ind+ 1] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 2); hDQQLLSel[ind+ 2] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 3); hDQQLLSel[ind+ 3] = new TH1D(sb,sb,100,-5.0,5.0);
    sprintf(sb,"hDQQLLSel_%d",ind+ 4); hDQQLLSel[ind+ 4] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 5); hDQQLLSel[ind+ 5] = new TH1D(sb,sb,200,0.0,4.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 6); hDQQLLSel[ind+ 6] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 7); hDQQLLSel[ind+ 7] = new TH1D(sb,sb,100,0.0,10.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 8); hDQQLLSel[ind+ 8] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+ 9); hDQQLLSel[ind+ 9] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+10); hDQQLLSel[ind+10] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+11); hDQQLLSel[ind+11] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+12); hDQQLLSel[ind+12] = new TH1D(sb,sb,100,0.0,1.);
    sprintf(sb,"hDQQLLSel_%d",ind+13); hDQQLLSel[ind+13] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+14); hDQQLLSel[ind+14] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+15); hDQQLLSel[ind+15] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+16); hDQQLLSel[ind+16] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+17); hDQQLLSel[ind+17] = new TH1D(sb,sb,60,0.0,3.0);
    sprintf(sb,"hDQQLLSel_%d",ind+18); hDQQLLSel[ind+18] = new TH1D(sb,sb,60,0.0,3.0);
    sprintf(sb,"hDQQLLSel_%d",ind+19); hDQQLLSel[ind+19] = new TH1D(sb,sb,100,0.0,10.);
    sprintf(sb,"hDQQLLSel_%d",ind+20); hDQQLLSel[ind+20] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+21); hDQQLLSel[ind+21] = new TH1D(sb,sb,400,0.0,400.);
    sprintf(sb,"hDQQLLSel_%d",ind+22); hDQQLLSel[ind+22] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+23); hDQQLLSel[ind+23] = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDQQLLSel_%d",ind+24); hDQQLLSel[ind+24] = new TH1D(sb,sb,200,0.0,200.);
 }

  for(int i=0; i<25; i++){
    for(int j=0; j<5; j++){
      AddOutput(hDQQLLSel[i+j*100]);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void QQLLEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void QQLLEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
