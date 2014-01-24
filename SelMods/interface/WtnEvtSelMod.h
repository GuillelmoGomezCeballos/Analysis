//------------------------------------------------------------------------------
// $Id: WtnEvtSelMod.h,v 1.4 2011/03/07 12:48:06 ceballos Exp $
//
// WtnEvtSelMod
//
// A Module for Selecting tau events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WtnEvtSelMod_H
#define HWWMODS_WtnEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WtnEvtSelMod : public BaseMod
  {
    public:
    WtnEvtSelMod(const char *name="WtnEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WtnEvtSelMod() {}
      void         SetPrintDebug(bool b)          { fPrintDebug        = b;}	
      void         SetPFTaus0Name(TString s)      { fPFTaus0Name       = s;}    
      void         SetPFTaus1Name(TString s)      { fPFTaus1Name       = s;}    
      void         SetCleanJetsName (TString s)   { fCleanJetsName     = s;}	
      void         SetMetName(TString s)          { fMetName           = s;}   
      void         SetCleanPFTausName(TString s)  { fCleanPFTausName   = s;}	

    protected:
      bool              fPrintDebug;
      TString           fMCTausName;
      TString           fCleanJetsName;
      TString           fMetName;
      TString           fPFTaus0Name;
      const PFTauCol   *fPFTaus0; 	       
      TString           fPFTaus1Name;
      const PFTauCol   *fPFTaus1; 	       
      TString           fCleanPFTausName;

      TH1D             *hDWtnSel[400];
 
      int               fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(WtnEvtSelMod,1) // TAM example analysis module
  };
}
#endif
