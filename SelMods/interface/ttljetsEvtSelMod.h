//------------------------------------------------------------------------------
// $Id: ttljetsEvtSelMod.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// ttljetsEvtSelMod
//
// A Module to select tt->ljets events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_TTLJETSEVTSELMOD_H
#define HWWMODS_TTLJETSEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ttljetsEvtSelMod : public BaseMod
  {
    public:
    ttljetsEvtSelMod(const char *name="ttljetsEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ttljetsEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }   
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }   
      void      SetMetName(TString s)          { fMetName = s;         }   
      void      SetPtJetCut(double x)          { fPtJetCut = x;        }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;       }
      void      SetTrigObjsName(const char *n) { fObjsName = n;        }

    protected:
      bool      fPrintDebug;	   // Debug info
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;	   // name of met collection
      TString   fCleanJetsName;    // name of clean central jets collection
      TString   fCaloJetName0;
      const CaloJetCol *fCaloJet0; // ItrCone5Jets
      TString      fObjsName;   //name of trigger objects

      TH1D      *hDLJPresel[200];
      TH1D      *hDLJSel[200];
      TH1D      *hDLJBTag[10];
      TH1D      *hDLJHLT[8];

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(ttljetsEvtSelMod,1) // TAM example analysis module
  };
}
#endif
