//------------------------------------------------------------------------------
// $Id: LowEvtSelMod.h,v 1.4 2012/01/19 08:19:37 ceballos Exp $
//
// LowEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_LowEvtSelMod_H
#define HWWMODS_LowEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class LowEvtSelMod : public BaseMod
  {
    public:
    LowEvtSelMod(const char *name="LowEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~LowEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;	   }
      void      SetIsFastSim(bool b)	       { fIsFastSim = b;	   }
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;	   }
      void      SetMetName(TString s)          { fMetName = s;  	   }   
      void      SetJetScaleSyst(double x)      { fJetScaleSyst = x;	   }   
      void      SetPtJetCut(double x)          { fPtJetCut = x; 	   }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;           }
      void      SetTrigObjsName(const char *n) { fObjsName = n; 	   }

    protected:
       // Trigger info
      TString      fObjsName;   //name of trigger objects
 
      bool      fPrintDebug;	 // Debug info
      bool      fIsFastSim;      // isFastSim
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fPlotType;	 // Type of histograms to make
      TString   fMetName;	 // name of met collection
      TString   fMuonName;	 // name of muon collection
      TString   fVertexName;	 // name of vertex collection
      TString   fCleanJetsName;  // name of clean central jets collection
      const MuonCol   *fMuons;	 // Muon branch
      const VertexCol *fVertices;	 // Vertices branches
      TString             fCaloMetName;	 // name of calomet collection
      const CaloMetCol   *fCaloMet; // caloMet branch
      Double_t  fJetScaleSyst;

      TH1D      *hDLowPresel[100];
      TH1D      *hDLowSel[300];

      int       fNEventsProcessed;

      MuonTools myMuonTools;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(LowEvtSelMod,1) // TAM example analysis module
  };
}
#endif
