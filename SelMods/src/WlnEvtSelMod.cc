// $Id: WlnEvtSelMod.cc,v 1.9 2012/04/18 14:59:41 ceballos Exp $

#include "Ana/SelMods/interface/WlnEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"

using namespace mithep;
ClassImp(mithep::WlnEvtSelMod)

//--------------------------------------------------------------------------------------------------
WlnEvtSelMod::WlnEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName(Names::gkCaloMetBrn),
  fTrackName(Names::gkTrackBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFMetName("PFMet"),
  fPFMet(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "WlnEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  ElectronOArr *CleanElectrons      = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons             = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  JetOArr *CleanJets                = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons             = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  LoadBranch(fPFMetName);
  const PFMet *PFMetStd        = fPFMet->At(0);

  double zAverage = 0.0;

  // Make lepton vector from muons and electrons
  for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
    if(CleanMuons->At(j)->Pt() <= 10) continue;
    zAverage = zAverage + CleanMuons->At(j)->BestTrk()->Z0();
  }
  for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
    const Electron *el = CleanElectrons->At(j);
    if(el->Pt() <= 10) continue;
    zAverage = zAverage + el->BestTrk()->Z0();
  }

  // Computing Z average (our primary vertex)
  if(leptons->GetEntries() > 0) zAverage = zAverage / leptons->GetEntries();

  if(leptons->GetEntries() >= 1){
    hDWlnSel[0]->Fill(TMath::Min((double)leptons->GetEntries(),9.499),NNLOWeight->GetVal());
    hDWmnSel[0]->Fill(TMath::Min((double)CleanMuons->GetEntries(),9.499),NNLOWeight->GetVal());
    hDWenSel[0]->Fill(TMath::Min((double)CleanElectrons->GetEntries(),9.499),NNLOWeight->GetVal());

    if(leptons->GetEntries() == 1){
      // Sort and count the number of central Jets for vetoing
      int nCentralJets = 0;
      for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
        if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	   CleanJets->At(i)->Pt() > fPtJetCut){
          nCentralJets++;
        }
      }
      MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
      const Met *caloMet           = CleanMet->At(0);
      fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

      double deltaPhiMetLepton = fabs(MathUtils::DeltaPhi(caloMet->Phi(), leptons->At(0)->Phi()));

      double mTW = TMath::Sqrt(2.0*leptons->At(0)->Pt()*caloMet->Pt()*
        		      (1.0 - cos(deltaPhiMetLepton)));

      CompositeParticle totalPart;
      totalPart.AddDaughter(leptons->At(0));
      totalPart.AddDaughter(caloMet);
      if(CleanMuons->GetEntries() == 1){
	hDWmnSel[1]->Fill(TMath::Min(leptons->At(0)->Pt(),99.999),NNLOWeight->GetVal());
        if(leptons->At(0)->Pt() > 25.0){
    	  // d0 cut
    	  double d0_real = 99999;
    	  for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
	    if(fVertices->At(i0)->NTracks() > 0){
    	      double pD0 = CleanMuons->At(0)->BestTrk()->D0Corrected(*fVertices->At(i0));
    	      d0_real = TMath::Abs(pD0);
	      break;
	    }
    	  }
          hDWmnSel[2]->Fill(TMath::Min(TMath::Abs(d0_real),0.999),NNLOWeight->GetVal());
          if(TMath::Abs(d0_real) < 0.05){
            hDWmnSel[4]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
	    //if(mTW > 50){
	    if(caloMet->Pt() > 30){
              hDWmnSel[3]->Fill(TMath::Min(mTW, 199.999),NNLOWeight->GetVal());
              hDWmnSel[5]->Fill(TMath::Min(TMath::Abs(d0_real)/CleanMuons->At(0)->BestTrk()->D0Err(),19.999),NNLOWeight->GetVal());
              hDWmnSel[6]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
              hDWmnSel[7]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	      hDWmnSel[8]->Fill(nCentralJets,NNLOWeight->GetVal());
              hDWmnSel[10]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
              hDWmnSel[11]->Fill(TMath::Min(totalPart.Pt(),199.999),NNLOWeight->GetVal());
            }
	  }
	}
      }
      if(CleanElectrons->GetEntries() == 1){
        hDWenSel[1]->Fill(TMath::Min(leptons->At(0)->Pt(),99.999),NNLOWeight->GetVal());
        if(leptons->At(0)->Pt() > 30.0){
    	  // d0 cut
    	  double d0_real = 99999;
    	  for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
    	      if(fVertices->At(i0)->NTracks() > 0){
	      double pD0 = CleanElectrons->At(0)->GsfTrk()->D0Corrected(*fVertices->At(i0));
    	      d0_real = TMath::Abs(pD0);
    	      break;
	    }
	  }
          hDWenSel[2]->Fill(TMath::Min(TMath::Abs(d0_real),0.999),NNLOWeight->GetVal());
          if(TMath::Abs(d0_real) < 0.05){
            hDWenSel[4]->Fill(TMath::Min(caloMet->Pt(), 199.999),NNLOWeight->GetVal());
	    //if(mTW > 50){
	    if(caloMet->Pt() > 30){
              hDWenSel[3]->Fill(TMath::Min(mTW, 199.999),NNLOWeight->GetVal());
	      hDWenSel[5]->Fill(TMath::Min(TMath::Abs(d0_real)/CleanElectrons->At(0)->BestTrk()->D0Err(),19.999),NNLOWeight->GetVal());
              hDWenSel[6]->Fill(TMath::Min(caloMet->Pt(),199.999),NNLOWeight->GetVal());
              hDWenSel[7]->Fill(deltaPhiMetLepton * 180./TMath::Pi(),NNLOWeight->GetVal());
	      hDWenSel[8]->Fill(nCentralJets,NNLOWeight->GetVal());
              hDWenSel[10]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
              hDWenSel[11]->Fill(TMath::Min(totalPart.Pt(),199.999),NNLOWeight->GetVal());
	    }
          }
	}
      }
    } // Nlep == 1
  } // Nlep >= 1
}
//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fPFMetName,     fPFMet);
  ReqBranch(fTrackName,     fTracks);

  char sb[200];
  // Event
  sprintf(sb,"hDWlnSel_%d",0);  hDWlnSel[0]  = new TH1D(sb,sb,10,-0.5,9.5);
  for(int i=0; i<1; i++){
    AddOutput(hDWlnSel[i]);
  }

  sprintf(sb,"hDWmnSel_%d", 0);  hDWmnSel[ 0]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWmnSel_%d", 1);  hDWmnSel[ 1]  = new TH1D(sb,sb,100,0.0,100.);
  sprintf(sb,"hDWmnSel_%d", 2);  hDWmnSel[ 2]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDWmnSel_%d", 3);  hDWmnSel[ 3]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWmnSel_%d", 4);  hDWmnSel[ 4]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWmnSel_%d", 5);  hDWmnSel[ 5]  = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDWmnSel_%d", 6);  hDWmnSel[ 6]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWmnSel_%d", 7);  hDWmnSel[ 7]  = new TH1D(sb,sb,90,0.0,180.);
  sprintf(sb,"hDWmnSel_%d", 8);  hDWmnSel[ 8]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWmnSel_%d", 9);  hDWmnSel[ 9]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWmnSel_%d",10);  hDWmnSel[10]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWmnSel_%d",11);  hDWmnSel[11]  = new TH1D(sb,sb,200,0.0,200.);

  for(int i=0; i<12; i++){
    AddOutput(hDWmnSel[i]);
  }

  sprintf(sb,"hDWenSel_%d", 0);  hDWenSel[ 0]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWenSel_%d", 1);  hDWenSel[ 1]  = new TH1D(sb,sb,100,0.0,100.);
  sprintf(sb,"hDWenSel_%d", 2);  hDWenSel[ 2]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDWenSel_%d", 3);  hDWenSel[ 3]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWenSel_%d", 4);  hDWenSel[ 4]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWenSel_%d", 5);  hDWenSel[ 5]  = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDWenSel_%d", 6);  hDWenSel[ 6]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWenSel_%d", 7);  hDWenSel[ 7]  = new TH1D(sb,sb,90,0.0,180.);
  sprintf(sb,"hDWenSel_%d", 8);  hDWenSel[ 8]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWenSel_%d", 9);  hDWenSel[ 9]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWenSel_%d",10);  hDWenSel[10]  = new TH1D(sb,sb,200,0.0,200.);
  sprintf(sb,"hDWenSel_%d",11);  hDWenSel[11]  = new TH1D(sb,sb,200,0.0,200.);

  for(int i=0; i<12; i++){
    AddOutput(hDWenSel[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
