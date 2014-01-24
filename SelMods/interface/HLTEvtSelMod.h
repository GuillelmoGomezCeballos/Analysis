//------------------------------------------------------------------------------
// $Id: HLTEvtSelMod.h,v 1.1 2010/10/23 04:44:54 ceballos Exp $
//
// HLTEvtSelMod
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_HLTEvtSelMod_H
#define HWWMODS_HLTEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class HLTEvtSelMod : public BaseMod
  {
    public:
    HLTEvtSelMod(const char *name="HLTEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~HLTEvtSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}	
      void         SetTrigObjsName(const char *n){ fObjsName = n;     }
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetPtJetCut(double x)         { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)        { fEtaJetCut = x;    }

    protected:
      bool      fPrintDebug;
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;
      TString   fCleanJetsName;
      TString             fEvtHdrName;
      const EventHeader  *fEventHeader;
      TString             fVertexName;
      const VertexCol    *fVertices;
      TString      fObjsName;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      float             fTreeVariables[42];
      TTree            *fTree;

      ClassDef(HLTEvtSelMod,1) // TAM example analysis module
  };
}
#endif
