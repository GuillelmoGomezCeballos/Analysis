// $Id: ttEvtSelMod.cc,v 1.10 2012/09/05 09:35:40 ceballos Exp $

#include "Analysis/SelMods/interface/ttEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitAna/DataTree/interface/ObjTypes.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/MetTools.h"

using namespace mithep;
ClassImp(mithep::ttEvtSelMod)

//--------------------------------------------------------------------------------------------------
ttEvtSelMod::ttEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("random"),
  fMuonName(Names::gkMuonBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCAllLeptonsName(ModNames::gkMCAllLeptonsName),
  fMuons(0),
  fPFJetName0("AKt5PFJets"),
  fPFJet0(0),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  //Get Generator Level information for matching
  ObjArray<MCParticle> *GenLeptons = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);

  //Obtain all the good objects from the event cleaning module
  ElectronOArr *CleanElectrons  = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons         = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  ParticleOArr *leptonsFakeable = GetObjThisEvt<ParticleOArr>("MergedLeptonsFakeable");
  ParticleOArr *leptons         = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MetOArr *CleanMet             = GetObjThisEvt<MetOArr>(fMetName);
  const Met *pfMet              = CleanMet->At(0);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  LoadBranch(fPFJetName0);
  LoadBranch(fPFCandidatesName);

  LoadBranch(fMuonName);
  ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(!MuonTools::PassSoftMuonCut(mu, fVertices, 0.2)) continue;
    
    bool isCleanMuon = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      if(fMuons->At(i) == CleanMuons->At(j) &&
  	 CleanMuons->At(j)->Pt() > 10) isCleanMuon = kTRUE;
    }
    if(isCleanMuon == kFALSE) DirtyMuons->Add(mu);
  }

  // Minimun Pt, Nleptons>=2 requirements
  if (leptonsFakeable->GetEntries() >= 2 &&
      leptonsFakeable->At(0)->Pt() > 20 && leptonsFakeable->At(1)->Pt() > 10){

    if (leptons->GetEntries() == 2 &&
        leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10){
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
      if(dilepton.Charge() != 0) hDttPresel[ 0+10*pairType]->Fill(TMath::Min(dilepton.Mass(),299.999),NNLOWeight->GetVal());
      else			 hDttPresel[ 1+10*pairType]->Fill(TMath::Min(dilepton.Mass(),299.999),NNLOWeight->GetVal());
    }

    GenericParticle *thePhoton = new GenericParticle(0,0,0,0);
    MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, fVertices->At(0), 0.1, 8.0, 5.0, 0.0,thePhoton);
    delete thePhoton;
    double pMET[2] = {metTools.GetProjectedMet(leptons,pfMet),
    		      metTools.GetProjectedTrackMet(leptons)};
    bool passCuts = TMath::Min(pMET[0],pMET[1]) > 20.0;

    if(passCuts == true){
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

      for(UInt_t i=0; i<leptons->GetEntries(); i++) {
        Particle *lep = leptons->At(i);
	UInt_t leptonGenType = 0;
        for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
          if(MathUtils::DeltaR(GenLeptons->At(j)->Mom(), lep->Mom()) < 0.10) {
            leptonGenType = GenLeptons->At(j)->PdgId();
            break;
          }
        }
        if(leptonGenType == 0) {
	  int theJet = -1; double deltaRMin = 999;
          for(UInt_t nj=0; nj<fPFJet0->GetEntries(); nj++){
	    const PFJet *jet = fPFJet0->At(nj);
	    Double_t deltaR = MathUtils::DeltaR(jet->Mom(),lep->Mom());
	    if(deltaR < deltaRMin) {deltaRMin = deltaR; theJet = nj;}
	  }
	  if(theJet >= 0) {
	    if     (leptonsFakeable->GetEntries() == 2 && lep->ObjType() == kMuon) {
	      hDttPresel[2]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[3]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[4]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[5]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 2 && lep->ObjType() == kElectron) {
	      hDttPresel[6]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[7]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[8]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[9]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 3 && lep->ObjType() == kMuon) {
	      hDttPresel[12]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[13]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[14]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[15]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 3 && lep->ObjType() == kElectron) {
	      hDttPresel[16]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[17]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[18]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[19]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	  }
	}
      }

      for(UInt_t i=0; i<leptonsOnlyFake->GetEntries(); i++) {
        Particle *lep = leptonsOnlyFake->At(i);
	UInt_t leptonGenType = 0;
        for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
          if(MathUtils::DeltaR(GenLeptons->At(j)->Mom(), lep->Mom()) < 0.10) {
            leptonGenType = GenLeptons->At(j)->PdgId();
            break;
          }
        }
        if(leptonGenType == 0) {
	  int theJet = -1; double deltaRMin = 999;
          for(UInt_t nj=0; nj<fPFJet0->GetEntries(); nj++){
	    const PFJet *jet = fPFJet0->At(nj);
	    Double_t deltaR = MathUtils::DeltaR(jet->Mom(),lep->Mom());
	    if(deltaR < deltaRMin) {deltaRMin = deltaR; theJet = nj;}
	  }
	  if(theJet >= 0) {
	    if     (leptonsFakeable->GetEntries() == 2 && lep->ObjType() == kMuon) {
	      hDttPresel[22]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[23]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[24]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[25]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 2 && lep->ObjType() == kElectron) {
	      hDttPresel[26]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[27]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[28]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[29]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 3 && lep->ObjType() == kMuon) {
	      hDttPresel[32]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[33]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[34]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[35]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	    else if(leptonsFakeable->GetEntries() == 3 && lep->ObjType() == kElectron) {
	      hDttPresel[36]->Fill(TMath::Min(deltaRMin,0.999),NNLOWeight->GetVal());
	      if(deltaRMin < 0.5) {
	        hDttPresel[37]->Fill(TMath::Min(fPFJet0->At(theJet)->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[38]->Fill(TMath::Min(lep->Pt(),199.999),NNLOWeight->GetVal());
	        hDttPresel[39]->Fill(TMath::Min(lep->Pt()/fPFJet0->At(theJet)->Pt(),1.999),NNLOWeight->GetVal());
	      }
	    }
	  }
	}
      }
      delete leptonsOnlyFake;
    } // passCuts
  } // preselection
  delete DirtyMuons;
}
//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMuonName, fMuons);
  ReqBranch(fPFJetName0, fPFJet0);
  ReqBranch(fPFCandidatesName, fPFCandidates);

  char sb[200];

  for(int j=0; j<4; j++){
    int ind = 10 * j;
    if(j <= 2) {
      sprintf(sb,"hDttPresel_%d",ind+0);  hDttPresel[ind+0]  = new TH1D(sb,sb,300,0.0,300.); 
      sprintf(sb,"hDttPresel_%d",ind+1);  hDttPresel[ind+1]  = new TH1D(sb,sb,300,0.0,300.);
    }
    sprintf(sb,"hDttPresel_%d",ind+2);  hDttPresel[ind+2]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDttPresel_%d",ind+3);  hDttPresel[ind+3]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+4);  hDttPresel[ind+4]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+5);  hDttPresel[ind+5]  = new TH1D(sb,sb,200,0.0,2.0); 
    sprintf(sb,"hDttPresel_%d",ind+6);  hDttPresel[ind+6]  = new TH1D(sb,sb,100,0.0,1.0);
    sprintf(sb,"hDttPresel_%d",ind+7);  hDttPresel[ind+7]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+8);  hDttPresel[ind+8]  = new TH1D(sb,sb,200,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+9);  hDttPresel[ind+9]  = new TH1D(sb,sb,200,0.0,2.0); 
  }

  for(int i=0; i<40; i++){
    if(i!=30&&i!=31) AddOutput(hDttPresel[i]);
  }

}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
