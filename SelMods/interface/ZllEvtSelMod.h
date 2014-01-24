//------------------------------------------------------------------------------
// $Id: ZllEvtSelMod.h,v 1.6 2012/04/18 14:59:44 ceballos Exp $
//
// ZllEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_ZLLEVTSELMOD_H
#define HWWMODS_ZLLEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ZllEvtSelMod : public BaseMod
  {
    public:
    ZllEvtSelMod(const char *name="ZllEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ZllEvtSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}	
      void         SetTrigObjsName(const char *n){ fObjsName = n;     }
      void         SetMCLeptonsName(TString s)   { fMCLeptonsName = s;}     
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetPtJetCut(double x)         { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)        { fEtaJetCut = x;    }

    protected:
      bool      fPrintDebug;	  // Debug info
      Double_t  fPtJetCut;
      Double_t  fEtaJetCut;
      TString   fPlotType;	  // Type of histograms to make
      TString   fMetName;	  // name of met collection
      TString   fTrackName;	  // name of track collection
      TString   fGsfTrackColName; //name of input lepton collection
      TString   fCleanJetsName;   // name of clean central jets collection
      TString   fMCLeptonsName;   // new lepton coll
      TString   fMuonName;	  // name of luon collection
      TString   fElectronName;    // name of electron collection
      const TrackCol     *fTracks;        // Track branch     
      const TrackCol     *fGsfTracks;     //!pointer to collection 
      const MuonCol      *fMuons;     	    
      const ElectronCol  *fElectrons; 		
      const JetCol       *fJets;
      TString             fPFMetName;
      const PFMetCol     *fPFMetStd;
      TString             fTCMetName;
      const MetCol       *fTCMetStd;
      TString             fBarrelSCName;
      TString             fEndcapSCName;
      const SuperClusterCol *fBarrelSC;
      const SuperClusterCol *fEndcapSC;
      TString             fEvtHdrName;   // name of event header branch
      const EventHeader  *fEventHeader;  // event header for current event
      TString             fVertexName;	 // name of vertex collection
      const VertexCol    *fVertices;	 // Vertices branches
      TString             fBeamSpotName;
      const BeamSpotCol  *fBeamSpot;     // BeamSpot branches

      // Trigger info
      TString      fObjsName;   //name of trigger objects

      TH1D         *hDZllZee[3];
      TH1D         *hDZllMET[200];
      TH2D         *hDZllSS2D;
      TH1D         *hDZllSS[200];
      TH1D         *hDZllHLT[300];
      TH1D         *hDZllTPIni[200];
      TH1D         *hDZllTP[200];
      TH1D         *hDZllTPGen[100];
      TH1D         *hDZllMuon[30];
      TH1D         *hDZllElSel[20];
      TH1D         *hDZllMuSel[10];

      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TRandom RND;

      void         getTrackIsolation(double theEta, double thePhi,
                                     ParticleOArr* leptons,
                                     const mithep::TrackCol *iTracks,
				     double &sumPt, int &nTracks);
      void         getCaloIsolation(double theEta, double thePhi,
                                    ParticleOArr* leptons,
                                    const mithep::CaloTowerCol *iTowers,
				    double lIso[4]);

      ClassDef(ZllEvtSelMod,1) // TAM example analysis module
  };
}
#endif
