//------------------------------------------------------------------------------
// $Id: LeptonEvtSelMod.h,v 1.10 2012/04/18 14:59:44 ceballos Exp $
//
// LeptonEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_LeptonEvtSelMod_H
#define HWWMODS_LeptonEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class LeptonEvtSelMod : public BaseMod
  {
    public:
    LeptonEvtSelMod(const char *name="LeptonEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~LeptonEvtSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetIsData(bool b)             { fIsData = b;       }
      void         SetMCLeptonsName(TString s)   { fMCLeptonsName = s;}     
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetAllVertexName(TString s)   { fAllVertexName = s;}   

    protected:
      bool         fPrintDebug;     // Debug info
      bool         fIsData;         // isData boolean
      TString      fMetName;	    // name of met collection
      TString      fTrackName;      // name of track collection
      TString      fTrack0Name;     // name of track collection
      TString      fTrack1Name;     // name of track collection
      TString      fTrack2Name;     // name of track collection
      TString      fMuonName;	    // name of muon collection
      TString      fElectronName;   // name of electron collection
      TString      fMCLeptonsName;  // new lepton coll
      TString      fAllVertexName;	 // name of all vertex collection
      TString      fVertexName;	    // name of vertex collection
      TString      fConversionName;	    // name of vertex collection
      const TrackCol     *fTracks;        // Track branch     
      const TrackCol     *fTracks0;        // Track branch     
      const TrackCol     *fTracks1;        // Track branch     
      const TrackCol     *fTracks2;        // Track branch     
      const MuonCol      *fMuons;     	    
      const ElectronCol  *fElectrons; 		
      const JetCol       *fJets;      	    
      const VertexCol    *fAllVertices;	 // Vertices branches
      const VertexCol    *fVertices;     // Vertices branches
      const DecayParticleCol *fConversions; //!conversion collection
      TString           fPFCandidatesName;
      const PFCandidateCol    *fPFCandidates;
      TString           fPileupEnergyDensityName;
      const PileupEnergyDensityCol *fPileupEnergyDensity;
      TString                   fBeamSpotName;           //name of beamspot collection
      const BeamSpotCol        *fBeamSpot;               //!beamspot branch

      TH1D         *hDLepSel[300];
      TH2D         *hDLepSel2D[10];
      Double_t      fCuts[6][8];
      TH1D         *hDCutEleSel[400];
      TH1D         *hDEleConvSel[80];
      TH1D         *hDD0LepSel[40];

      TH1D         *hDIsoMLepSel0[80];
      TH1D         *hDIsoMLepSel1[80];
      TH1D         *hDIsoELepSel0[80];
      TH1D         *hDIsoELepSel1[80];

      int          fNEventsProcessed;

      MuonTools    myMuonTools;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TRandom RND;

      ClassDef(LeptonEvtSelMod,1) // TAM example analysis module
  };
}
#endif
