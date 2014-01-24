//------------------------------------------------------------------------------
// $Id: SkimEvtSelMod.h,v 1.4 2012/04/24 11:26:27 ceballos Exp $
//
// SkimEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_SKIMEVTSELMOD_H
#define HWWMODS_SKIMEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class SkimEvtSelMod : public BaseMod
  {
    public:
    SkimEvtSelMod(const char *name="SkimEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~SkimEvtSelMod() {}
      void      SetTrigObjsName(const char *n)     { fObjsName = n;     }

    protected:
      TH1D      *hDskimPresel[5];
      TH1D      *hDskimSel[40];

      TString   fObjsName;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(SkimEvtSelMod,1) // TAM example analysis module
  };
}
#endif
