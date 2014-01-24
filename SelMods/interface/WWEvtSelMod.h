//------------------------------------------------------------------------------
// $Id: WWEvtSelMod.h,v 1.21 2012/05/05 08:49:29 ceballos Exp $
//
// WWEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WWEVTSELMOD_H
#define HWWMODS_WWEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WWEvtSelMod : public BaseMod
  {
    public:
    WWEvtSelMod(const char *name="WWEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WWEvtSelMod() {}
      void      SetPrintDebug(bool b)	           { fPrintDebug = b;          }
      void      SetIsFastSim(bool b)	           { fIsFastSim = b;           }
      void      SetIsData(bool b)	           { fIsData = b;              }
      void      SetCleanJetsNoPtCutName (TString s){ fCleanJetsNoPtCutName = s;}
      void      SetMetName(TString s)              { fMetName = s;             }   
      void      SetJetScaleSyst(double x)          { fJetScaleSyst = x;        }   
      void      SetPtJetCut(double x)              { fPtJetCut = x;            }
      void      SetEtaJetCut(double x)             { fEtaJetCut = x;           }
      void      SetUsePDFs(bool b)                 { fUsePDFs = b;             }
      void      SetAllVertexName(TString s)        { fAllVertexName = s;       }   
      void      SetIntRadius(Double_t dr)          { fIntRadius = dr;          }
      void      SetIsOldSelection(Bool_t b)        { fIsOldSelection = b;      }

    protected:
      bool      fPrintDebug;
      bool      fIsFastSim;
      bool      fIsData;
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fPlotType;	 // Type of histograms to make
      TString   fMetName;	 // name of met collection
      TString   fMuonName;	 // name of muon collection
      TString   fElectronName;	 // name of electron collection
      TString   fTrackName;	 // name of track collection
      TString   fAllVertexName;	 // name of all vertex collection
      TString   fVertexName;	 // name of vertex collection
      TString   fCleanJetsNoPtCutName; // name of clean central jets collection with no pt cut
      const MuonCol   *fMuons;	 // Muon branch
      const ElectronCol   *fElectrons;	 // Electron branch
      const TrackCol  *fTracks;	 // Track branch     
      const VertexCol *fAllVertices;	 // Vertices branches
      const VertexCol *fVertices;	 // Vertices branches
      TString             fPFMetName;
      const PFMetCol     *fPFMetStd;
      TString             fTCMetName;
      const MetCol       *fTCMetStd;
      TString              fPileupInfoName;
      const PileupInfoCol *fPileupInfo;

      TString   fCaloJetName0;
      TString   fCaloJetName1;
      TString   fCaloJetName2;
      TString   fCaloJetName3;
      TString   fTrackJetName0;
      TString   fPFJetName0;
      TString   fPFJetName1;
      TString   fPFJetName2;
      TString   fPFJetName3;
      const CaloJetCol *fCaloJet0;
      const CaloJetCol *fCaloJet1;
      const CaloJetCol *fCaloJet2;
      const CaloJetCol *fCaloJet3;
      const TrackJetCol *fTrackJet0;
      const PFJetCol   *fPFJet0;
      const PFJetCol   *fPFJet1;
      const PFJetCol   *fPFJet2;
      const PFJetCol   *fPFJet3;
      Double_t  fJetScaleSyst;
      bool fUsePDFs;
      TString           fEvtHdrName;                   // name of event header branch
      const EventHeader *fEventHeader;                 // event header for current event
      TString           fPFCandidatesName;
      const PFCandidateCol    *fPFCandidates;
      TString                   fBeamSpotName;           //name of beamspot collection
      const BeamSpotCol        *fBeamSpot;               //!beamspot branch
      TString                   fConversionBranchName;   //name of electron collection (input)
      const DecayParticleCol   *fConversions;            //!conversion collection
      TString              fMCPartName;         //name of MCParticle branch
      const MCParticleCol *fParticles;	        //!MCParticle branch
      TString              fPileupEnergyDensityName;
      const PileupEnergyDensityCol *fPileupEnergyDensity;
      Double_t  fIntRadius;
      Bool_t  fIsOldSelection;
 
      TH1D      *hDwwBTagCheck[2];
      TH1D      *hDwwXS[5];
      TH1D      *hDwwSkim[5];
      TH1D      *hDwwSel[400];
      TH1D      *hDwwMassSel[50];
      TH1D      *hDwwMET[90];
      TH1D      *hDwwPresel[30];
      TH2D      *hDwwSelAlphaEP0;
      TH2D      *hDwwSelAlphaEP1;
      TH2D      *hDwwDeltaPhiMetLeptonMet;
      TH2D      *hDwwSelD0Phi;
      TH1D      *hDwwJet[400];
      TH1D      *hDwwFake[400];
      TH1D      *hDwwJetSel[400];
      TH1D      *hDWWPDF[50];
      TH1D      *hDwwJetVar[15];
      TH1D      *hDwwBTag[90];
      TH1D      *hDwwOverlap[10];
      TH1D      *hDGenZwwSel[10];

      int       fNEventsProcessed;

      MuonTools myMuonTools;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(WWEvtSelMod,1) // TAM example analysis module
  };
}
#endif
