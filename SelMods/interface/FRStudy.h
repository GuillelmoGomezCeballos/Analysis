//------------------------------------------------------------------------------
// $Id: FRStudy.h,v 1.2 2011/05/08 11:50:38 ceballos Exp $
//
// FRStudy
//
// A Module for Selecting lepton(s)+lepton fakes(s) events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_FRStudy_H
#define HWWMODS_FRStudy_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class FRStudy : public BaseMod
  {
    public:
    FRStudy(const char *name="FRStudy", 
		 const char *title="Example analysis module with all branches");
      ~FRStudy() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }
      void      SetIsFastSim(bool b)	       { fIsFastSim = b;       }
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }
      void      SetMetName(TString s)          { fMetName = s;         }   
      void      SetPtJetCut(double x)          { fPtJetCut = x;        }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;       }
      void      SetIsData(bool b)              { fIsData = b;          }
      void      SetSelectGenLeptons(bool b)    { fSelectGenLeptons = b;}

    protected:
      bool      fPrintDebug;
      bool      fIsFastSim;
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;
      TString   fCleanJetsName;
      TString   fVertexName;
      const VertexCol   *fVertices;
      bool      fIsData;
      bool      fSelectGenLeptons;

      TH1D      *hDFRStudyPresel[20];
      TH1D      *hDFRStudySel[80];

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(FRStudy,1) // TAM example analysis module
  };
}
#endif
