//------------------------------------------------------------------------------
// $Id: ttEvtSelMod.h,v 1.4 2012/09/04 15:18:38 ceballos Exp $
//
// ttEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_TTEVTSELMOD_H
#define HWWMODS_TTEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ttEvtSelMod : public BaseMod
  {
    public:
    ttEvtSelMod(const char *name="ttEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ttEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }   
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }   
      void      SetMCLeptonsName(TString s)    { fMCLeptonsName = s;   }     
      void      SetMCAllLeptonsName(TString s) { fMCAllLeptonsName = s;}     
      void      SetMetName(TString s)          { fMetName = s;         }   
      void      SetPtJetCut(double x)          { fPtJetCut = x;        }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;       }

    protected:
      bool      fPrintDebug;	   // Debug info
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;	   // name of met collection
      TString   fMuonName;	   // name of muon collection
      TString   fVertexName;	   // name of vertex collection
      TString   fCleanJetsName;    // name of clean central jets collection
      TString   fMCLeptonsName;    // new lepton coll (from W)
      TString   fMCAllLeptonsName; // new lepton coll (all)
      const MuonCol	*fMuons;	   // Muon branch
      const VertexCol *fVertices;	   // Vertices branches
      TString   fPFJetName0;
      const PFJetCol   *fPFJet0;
      TString           fPFCandidatesName;
      const PFCandidateCol *fPFCandidates;

      TH1D      *hDttPresel[50];

      int       fNEventsProcessed;

      MuonTools myMuonTools;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(ttEvtSelMod,1) // TAM example analysis module
  };
}
#endif
