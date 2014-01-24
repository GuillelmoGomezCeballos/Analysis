//------------------------------------------------------------------------------
// $Id: WlnEvtSelMod.h,v 1.4 2011/11/27 06:17:05 ceballos Exp $
//
// WlnEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WlnEvtSelMod_H
#define HWWMODS_WlnEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WlnEvtSelMod : public BaseMod
  {
    public:
    WlnEvtSelMod(const char *name="WlnEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WlnEvtSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}	
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetPtJetCut(double x)         { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)        { fEtaJetCut = x;    }

    protected:
      bool         fPrintDebug;     // Debug info
      Double_t     fPtJetCut;
      Double_t     fEtaJetCut;
      TString      fPlotType;	    // Type of histograms to make
      TString      fMetName;	    // name of met collection
      TString      fTrackName;      // name of track collection
      TString      fCaloTowerName;  // name of calotower collection
      TString      fCleanJetsName;  // name of clean central jets collection
      TString      fVertexName;	    // name of vertex collection
      const TrackCol     *fTracks;        // Track branch     
      const CaloTowerCol *fCaloTowers;    // CaloTowers branch     
      const JetCol       *fJets;      	    
      const VertexCol    *fVertices;	    // Vertices branches
      TString             fPFMetName;	 // name of pfmet collection
      const PFMetCol     *fPFMet; // pfMet branch

      TH1D         *hDWlnSel[10];
      TH1D         *hDWmnSel[50];
      TH1D         *hDWenSel[50];

      int          fNEventsProcessed;

      MuonTools    myMuonTools;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(WlnEvtSelMod,1) // TAM example analysis module
  };
}
#endif
