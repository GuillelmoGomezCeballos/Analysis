//-----------------------------------------------------------------------------
// $Id: ZZEvtSelMod.h,v 1.6 2011/11/01 10:15:50 ceballos Exp $
//
// ZZEvtSelMod
//
// A Module for Selecting WZ/ZZ-> >=3 leptons
// and produces some distributions
//
//
// Authors: ceballos
//-----------------------------------------------------------------------------

#ifndef HWWMODS_ZZEvtSelMod_H
#define HWWMODS_ZZEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/EventHeader.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ZZEvtSelMod : public BaseMod
  {
    public:
    ZZEvtSelMod(const char *name="ZZEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ZZEvtSelMod() {}
      void         SetPrintDebug       (bool b)       { fPrintDebug       = b;}
      void         SetMuonName         (TString s)    { fMuonName         = s;}
      void         SetElectronName     (TString s)    { fElectronName     = s;}
      void         SetLeptonName       (TString s)    { fLeptonName       = s;}
      void         SetJetName          (TString s)    { fJetName          = s;}
      void         SetMuonFakeName     (TString s)    { fMuonFakeName     = s;}
      void         SetElectronFakeName (TString s)    { fElectronFakeName = s;}
      void         SetLeptonFakeName   (TString s)    { fLeptonFakeName   = s;}
      void         SetIsData           (bool b)       { fIsData           = b;}

    protected:
      bool               fPrintDebug;
      TString            fMuonName;
      TString            fElectronName;
      TString            fLeptonName;
      TString            fJetName;
      TString            fMuonFakeName;
      TString            fElectronFakeName;
      TString            fLeptonFakeName;
      TString            fPFMetName;
      const PFMetCol    *fPFMet;
      bool               fIsData;
      TString            fEvtHdrName;
      const EventHeader *fEventHeader;
  
      TH1D         *hDZZXSel[10];
      TH1D         *hDZZXGenSel[10];

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(ZZEvtSelMod,1) // TAM example analysis module
  };
}
#endif
