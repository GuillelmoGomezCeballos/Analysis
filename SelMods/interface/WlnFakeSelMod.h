//------------------------------------------------------------------------------
// $Id: WlnFakeSelMod.h,v 1.4 2012/04/18 14:59:44 ceballos Exp $
//
// WlnFakeSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WlnFakeSelMod_H
#define HWWMODS_WlnFakeSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WlnFakeSelMod : public BaseMod
  {
    public:
    WlnFakeSelMod(const char *name="WlnFakeSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WlnFakeSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetMCLeptonsName(TString s)   { fMCLeptonsName = s;}     
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetPtJetCut(double x)         { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)        { fEtaJetCut = x;    }
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}

    protected:
      bool         fPrintDebug;     // Debug info
      Double_t     fPtJetCut;
      Double_t     fEtaJetCut;
      TString      fPlotType;	    // Type of histograms to make
      TString      fMetName;	    // name of met collection
      TString      fTrackName;      // name of track collection
      TString      fCleanJetsName;  // name of clean central jets collection
      TString      fMCLeptonsName;  // new lepton coll
      TString      fMuonName;       // new muon coll
      TString      fElectronName;   // new electron coll
      TString      fVertexName;	    // name of vertex collection
      TString      fConversionBranchName;   //name of electron collection (input)
      const TrackCol     *fTracks;        // Track branch     
      const MuonCol      *fMuons;     	    
      const ElectronCol  *fElectrons; 		
      const JetCol       *fJets;      	    
      const VertexCol    *fVertices;	    // Vertices branches
      Bool_t fApplyConversionFilter;  //!whether remove conversions
      const DecayParticleCol *fConversions;            //!conversion collection

      TH1D         *hDWlnFakeSel[400];
      TH1D         *hDWlnFakeGenType[50];
    
      int          fNEventsProcessed;

      MuonTools    myMuonTools;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TRandom RND;

      ClassDef(WlnFakeSelMod,1) // TAM example analysis module
  };
}
#endif
