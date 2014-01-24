//------------------------------------------------------------------------------
// $Id: ZttEvtSelMod.h,v 1.3 2010/10/19 22:09:43 ceballos Exp $
//
// ZttEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_ZttEvtSelMod_H
#define HWWMODS_ZttEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ZttEvtSelMod : public BaseMod
  {
    public:
    ZttEvtSelMod(const char *name="ZttEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ZttEvtSelMod() {}
      void      SetPrintDebug(bool b)	      { fPrintDebug = b;   } 
      void      SetCleanJetsName (TString s)  { fCleanJetsName = s;} 
      void      SetTrigObjsName(const char *n){ fObjsName = n;     }
      void      SetMetName(TString s)	      { fMetName = s;	   }
      void      SetPtJetCut(double x)	      { fPtJetCut = x;     }
      void      SetEtaJetCut(double x)        { fEtaJetCut = x;    }
      void      SetJetScaleSyst(double x)     { fJetScaleSyst = x; }   

    protected:
      bool         fPrintDebug;     // Debug info
      Double_t     fPtJetCut;
      Double_t     fEtaJetCut;
      Double_t     fJetScaleSyst;
      TString      fMetName;	    // name of met collection
      TString      fCleanJetsName;  // name of clean central jets collection
      TString      fCaloJetName0;
      const CaloJetCol *fCaloJet0;

      // Trigger info
      TString      fObjsName;   //name of trigger objects

      TH1D         *hDZttSel[200];
      TH1D         *hDZttMCSel[50];
      TH1D         *hDZttHLT[200];
      TH1D         *hDZemSel[30];
    
      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TRandom RND;

      ClassDef(ZttEvtSelMod,1) // TAM example analysis module
  };
}
#endif
