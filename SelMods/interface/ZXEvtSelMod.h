//-----------------------------------------------------------------------------
// $Id: ZXEvtSelMod.h,v 1.6 2011/11/01 10:15:50 ceballos Exp $
//
// ZXEvtSelMod
//
// A Module for Selecting WZ/ZZ-> >=3 leptons
// and produces some distributions
//
//
// Authors: ceballos
//-----------------------------------------------------------------------------

#ifndef HWWMODS_ZXEvtSelMod_H
#define HWWMODS_ZXEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ZXEvtSelMod : public BaseMod
  {
    public:
    ZXEvtSelMod(const char *name="ZXEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ZXEvtSelMod() {}
      void         SetPrintDebug(bool b)          { fPrintDebug = b;   }   
      void         SetCleanJetsName (TString s)   { fCleanJetsName = s;}   
      void         SetMCLeptonsName(TString s)    { fMCLeptonsName = s;}     
      void         SetMetName(TString s)          { fMetName = s;      }   
      void         SetPtJetCut(double x)          { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)         { fEtaJetCut = x;    }
      void         SetUsePDFs(bool b)             { fUsePDFs = b;      }

    protected:
      bool         fPrintDebug;      // Debug info
      Double_t     fPtJetCut;
      Double_t     fEtaJetCut;
      TString      fMetName;	     // name of met collection
      TString      fCleanJetsName;   // name of clean central jets collection
      TString      fMCLeptonsName ;  // new lepton coll
      TString           fCaloJetName0;
      const CaloJetCol *fCaloJet0;
      bool         fUsePDFs;
      TString               fPFCandidatesName;
      const PFCandidateCol *fPFCandidates;
      TString                   fPileupEnergyDensityName;
      const PileupEnergyDensityCol *fPileupEnergyDensity;
      Double_t fRho;
      TString   fTrackName;	  // name of track collection
      TString   fGsfTrackName; //name of input lepton collection
      const TrackCol     *fTracks;        // Track branch     
      const TrackCol     *fGsfTracks;     //!pointer to collection 
  
      TH1D         *hDZZGenSel[20];
      TH1D         *hDWGGenSel[10];
      TH1D         *hDZXSel[10];
      TH1D         *hDWZSel[400];
      TH1D         *hDZZSel[500];
      TH1D         *hDZZtt0Sel[200];
      TH1D         *hDZZtt1Sel[700];
      TH1D         *hDWZPDF[50];
      TH1D         *hDZZPDF[50];
      TH1D         *hDZHSel[30];
      TH1D         *hDHZZForwardSel[30];
      TH1D         *hDHZZForwardPSel[30];
      TH1D         *hDHZZForwardMSel[30];
      void          HZZForward(const ParticleOArr *leptons, const PFCandidateCol *fPFCandidates, double weight);
      void          HZZForward(const ParticleOArr *leptons, const PhotonCol *fPhotons, double weight);
      void          HZZForward(const ParticleOArr *leptons, const ParticleOArr *fakes, double weight);

      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(ZXEvtSelMod,1) // TAM example analysis module
  };
}
#endif
