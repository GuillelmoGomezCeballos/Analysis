//------------------------------------------------------------------------------
// $Id: IsoStudy.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// IsoStudy
//
// A Module for Selecting lepton(s)+lepton fakes(s) events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_IsoStudy_H
#define HWWMODS_IsoStudy_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class IsoStudy : public BaseMod
  {
    public:
    IsoStudy(const char *name="IsoStudy", 
		 const char *title="Example analysis module with all branches");
      ~IsoStudy() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }
      void      SetIsFastSim(bool b)	       { fIsFastSim = b;       }
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }
      void      SetMetName(TString s)          { fMetName = s;         }   
      void      SetPtJetCut(double x)          { fPtJetCut = x;        }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;       }

    protected:
      bool      fPrintDebug;
      bool      fIsFastSim;
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;
      TString   fCleanJetsName;
      TString   fVertexName;
      const VertexCol   *fVertices;
      TString   fMuonName;
      const MuonCol     *fMuons;
      TString   fElectronName;
      const ElectronCol *fElectrons;
      TString           fEvtHdrName;   // name of event header branch
      const EventHeader *fEventHeader; // event header for current event

      TH1D      *hDisoPresel[20];
      TH1D      *hDiso2lSel[70];
      TH1D      *hDiso3lSel[60];

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(IsoStudy,1) // TAM example analysis module
  };
}
#endif
