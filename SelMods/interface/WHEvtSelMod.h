//-----------------------------------------------------------------------------
// $Id: WHEvtSelMod.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// WHEvtSelMod
//
// A Module for Selecting WH -> >=2 leptons
// and produces some distributions
//
//
// Authors: ceballos
//-----------------------------------------------------------------------------

#ifndef HWWMODS_WHEvtSelMod_H
#define HWWMODS_WHEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WHEvtSelMod : public BaseMod
  {
    public:
    WHEvtSelMod(const char *name="WHEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WHEvtSelMod() {}
      void         SetPrintDebug(bool b)            { fPrintDebug = b;   }   
      void         SetCleanJetsName (TString s)     { fCleanJetsName = s;}   
      void         SetMetName(TString s)            { fMetName = s;      }   
      void         SetPtJetCut(double x)            { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)           { fEtaJetCut = x;    }
      void         SetMCEventInfoName(const char *s){ fMCEvInfoName  = s;}
      void         SetIsData(bool b)                { fIsData  = b;      }

    protected:
      bool         fPrintDebug;      // Debug info
      Double_t     fPtJetCut;
      Double_t     fEtaJetCut;
      TString      fMetName;	     // name of met collection
      TString      fCleanJetsName;   // name of clean central jets collection
      TString      fMCLeptonsName ;  // new lepton coll
      TString      fCaloJetName0;
      const CaloJetCol *fCaloJet0; // ItrCone5Jets
      TString      fMCEvInfoName;     //event info branch name
      const MCEventInfo *fMCEventInfo;      //!event info branch pointer
      bool         fIsData;

      TH1D         *hDWHSel[1];
      TH1D         *hDWH2lXSel[300];
      TH1D         *hDWH3lSel[400];
      TH1D         *hDWHSSSel[300];

      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(WHEvtSelMod,1) // TAM example analysis module
  };
}
#endif
