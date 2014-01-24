//------------------------------------------------------------------------------
// $Id: ScaleFactorMod.h,v 1.2 2011/04/18 22:28:42 ceballos Exp $
//
// Module to apply scale factor weights
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_ScaleFactorMod_H
#define HWWMODS_ScaleFactorMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ScaleFactorMod : public BaseMod
  {
    public:
    ScaleFactorMod(const char *name="ScaleFactorMod", 
		 const char *title="Example analysis module with all branches");
      ~ScaleFactorMod() {}
      void         SetTrigObjsName(const char *n){ fObjsName    = n; }
      void         SetIsData(Bool_t b)           { fIsData      = b; }
      void         SetApplyVertex(Bool_t b)      { fApplyVertex = b; }
      void         SetVertexName(TString s)      { fVertexName  = s; }   

    protected:
      TString          fObjsName;
      Bool_t           fIsData;
      Bool_t           fApplyVertex;
      TString          fVertexName;
      const VertexCol *fVertices;

      double       FactorHLT(Bool_t isMu, Double_t pt, Double_t eta);
      double       FactorEff(Bool_t isMu, Double_t pt, Double_t eta);
      double       FactorVrt(Int_t vrt);

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(ScaleFactorMod,1) // TAM example analysis module
  };
}
#endif
