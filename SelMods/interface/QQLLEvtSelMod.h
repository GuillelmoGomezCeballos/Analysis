//------------------------------------------------------------------------------
// $Id: QQLLEvtSelMod.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// QQLLEvtSelMod
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_QQLLEVTSELMOD_H
#define HWWMODS_QQLLEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class QQLLEvtSelMod : public BaseMod
  {
    public:
    QQLLEvtSelMod(const char *name="QQLLEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~QQLLEvtSelMod() {}
      void      SetPrintDebug(bool b)	           { fPrintDebug = b;          }
      void      SetIsFastSim(bool b)	           { fIsFastSim = b;           }
      void      SetMetName(TString s)              { fMetName = s;             }   
      void      SetJetScaleSyst(double x)          { fJetScaleSyst = x;        }   
      void      SetCleanJetsName (TString s)       { fCleanJetsName = s;       }
      void      SetPtJetCut(double x)              { fPtJetCut = x;            }
      void      SetEtaJetCut(double x)             { fEtaJetCut = x;           }

    protected:
      bool      fPrintDebug;	 // Debug info
      bool      fIsFastSim;      // isFastSim
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;	 // name of met collection
      TString   fVertexName;	 // name of vertex collection
      TString   fCleanJetsName;  // name of clean central jets collection
      const VertexCol *fVertices; // Vertices branches

      TString   fCaloJetName0;
      const CaloJetCol *fCaloJet0; // ICCone5Jets
      Double_t  fJetScaleSyst;

      TH1D      *hDQQLLSel[500];

      int       fNEventsProcessed;

      MuonTools myMuonTools;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(QQLLEvtSelMod,1) // TAM example analysis module
  };
}
#endif
