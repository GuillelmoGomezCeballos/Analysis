//------------------------------------------------------------------------------
// $Id: AAWWEvtSelMod.h,v 1.2 2013/01/09 10:20:13 ceballos Exp $
//
// AAWWEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_AAWWEvtSelMod_H
#define HWWMODS_AAWWEvtSelMod_H

#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class AAWWEvtSelMod : public BaseMod
  {
    public:
    AAWWEvtSelMod(const char *name="AAWWEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~AAWWEvtSelMod() {}
      void      SetPrintDebug(bool b)	           { fPrintDebug = b; }
      void      SetIsData(bool b)	           { fIsData = b;     }

    protected:
      bool                  fPrintDebug;
      bool                  fIsData;
      TString               fVertexName;
      const VertexCol      *fVertices;
      TString               fPFMetName;
      const PFMetCol       *fPFMetStd;
      TString               fPileupInfoName;
      const PileupInfoCol  *fPileupInfo;
      TString               fEvtHdrName;
      const EventHeader     *fEventHeader;
      TString               fPFCandidatesName;
      const PFCandidateCol *fPFCandidates;

      TH1D      *hDaawwPresel[20];

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(AAWWEvtSelMod,1) // TAM example analysis module
  };
}
#endif
