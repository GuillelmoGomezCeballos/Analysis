//------------------------------------------------------------------------------
// $Id: GammaXEvtSelMod.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// GammaXEvtSelMod
//
// A Module for Selecting GammaX events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_GammaXEvtSelMod_H
#define HWWMODS_GammaXEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class GammaXEvtSelMod : public BaseMod
  {
    public:
    GammaXEvtSelMod(const char *name="GammaXEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~GammaXEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }   
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }   
      void      SetGoodJetsName (TString s)    { fGoodJetsName = s;    }   
      void      SetMetName(TString s)          { fMetName = s;         }   
      void      SetPtJetCut(double x)          { fPtJetCut = x;        }
      void      SetEtaJetCut(double x)         { fEtaJetCut = x;       }

    protected:
      bool      fPrintDebug;	  
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fMetName;	  
      TString   fCleanJetsName;   
      TString   fGoodJetsName;   
      TString   fMCPhotonsName;   
      TString   fCleanPhotonsName;
      TString   fPhotonBranchName;
      const PhotonCol *fPhotons;  	  
      TString              fMCPartName;
      const MCParticleCol *fParticles;

      TH1D      *hDGMC[4];
      TH1D      *hDGSel[80];
      TH2D      *hDGSel2D[2];
      TH1D      *hDGGSel[40];
      TH1D      *hDGLSel[80];
      TH1D      *hDGLLSel[80];
      TH1D      *hDGJetSel[40];
      TH1D      *hDBeforeHLTSel[1];

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(GammaXEvtSelMod,1) // TAM example analysis module
  };
}
#endif
