// $Id: runHWWSelection_ZH.C,v 1.1 2013/11/06 14:34:33 ceballos Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <exception>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "Ana/SelMods/interface/LeptonEvtSelMod.h"
#include "Ana/SelMods/interface/ZXEvtSelMod.h"
#include "Ana/SelMods/interface/ttEvtSelMod.h"
#include "Ana/SelMods/interface/ttljetsEvtSelMod.h"
#include "Ana/SelMods/interface/ZllEvtSelMod.h"
#include "Ana/SelMods/interface/ZttEvtSelMod.h"
#include "Ana/SelMods/interface/WlnEvtSelMod.h"
#include "Ana/SelMods/interface/WWEvtSelMod.h"
#include "Ana/SelMods/interface/AAWWEvtSelMod.h"
#include "MitHiggs/HwwMods/interface/HwwMakeNtupleMod.h"
#include "Ana/SelMods/interface/GammaXEvtSelMod.h"
#include "Ana/SelMods/interface/WlnFakeSelMod.h"
#include "Ana/SelMods/interface/WtnEvtSelMod.h"
#include "Ana/SelMods/interface/IsoStudy.h"
#include "Ana/SelMods/interface/FRStudy.h"
#include "Ana/SelMods/interface/SkimEvtSelMod.h"
#include "Ana/SelMods/interface/QQLLEvtSelMod.h"
#include "Ana/SelMods/interface/LowEvtSelMod.h"
#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runHWWSelection(const char *catalogDir   = "~/scratch0/catalog",
		     const char *book	      = "cern/filefi/025",
                     const char *dataset      = "f11-wwj-v14b-bp",
                     const char *fileset      = "0001",
                     int nsel = 620, int NEvents = 100)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;
  if(nsel == 299) gDebugLevel = 3;

  double theMH     = 0.0;
  double theWidth  = 0.0;
  int    theBWflag = 0;

  TString theOutputDir    = "test/";
  bool usePDFProducer     = false;
  bool isData             = false;
  bool isDataMuonElectron = false;
  bool isDataDMuon        = false;
  bool isDataSMuon        = false;
  bool isDataDElectron    = false;
  bool isDataSElectron    = false;
  bool isDataPhoton       = false;
  bool isPhotonMCSel      = false;
  bool applyMllGenCut     = false;
  bool applyZVetoGenCut   = false;
  bool is2011Corr         = false;
  TString theWWMuId  = "WWMuIdV4";
  TString theWWMuIso = "PFIsoBetaPUCorrected";
  RhoUtilities::RhoType theRhoType = RhoUtilities::CMS_RHO_RHOKT6PFJETS;
  double muIsoCut = -999.0;
  bool applyMVACut = true;

  if(nsel > 400){ // 42x
    theOutputDir = "test1/";
    theWWMuId    = "WWMuIdV3";
    theWWMuIso   = "TrackCaloSliding";
    theRhoType   = RhoUtilities::DEFAULT;
    is2011Corr   = true;
    muIsoCut     = -999.0;
    applyMVACut  = false;
    if(nsel != 299) nsel = nsel - 400;
  }

  double fIntRadius = 0.0;
  double ptJetCut = 30.0;
  double etaJetCut = 4.7;
  const int Nfiles = 25;
  TString files[Nfiles] = {"", "", "", "", "",
                           "", "", "", "", "",
		           "", "", "", "", "",
		           "", "", "", "", "",
		           "", "", "", "", ""};
  int processid = -999999999;
  TString fInputFilenameKF = string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/HWW_KFactors_PowhegToNNLL_160_7TeV.dat"));
  // /home/mitprod/cms/jobs/filler/014/Productions
  TString MCType = "kMCTypeUndef";
  Bool_t applyPartonFlavorFilter = kFALSE;
  Bool_t applyISRFilter          = kFALSE;
  Bool_t applyWWFilter           = kFALSE;
  Bool_t applyVGFilter           = kFALSE;
  Bool_t applyFilterBTEvents     = kFALSE;
  int fDecay = -1;
  Bool_t useSelectGenLeptons   = kTRUE;
  TString  myRootFile;
  if(TString(fileset) == TString("all")) myRootFile = TString("histo_test");
  else                                   myRootFile = TString("histo_") + dataset;
  files[0] = "*.root";

  if     (nsel ==  1){
    fDecay = 110; theMH = 110; theWidth =  2.82e-03; theBWflag = -1;
  }
  else if(nsel ==  2){
    fDecay = 1110; theMH = 110; theWidth =  2.82e-03; theBWflag = -1;
  }
  else if(nsel ==  3){
    fDecay = 115; theMH = 115; theWidth =  3.09e-03; theBWflag = -1;
  }
  else if(nsel ==  4){
    fDecay = 1115; theMH = 115; theWidth =  3.09e-03; theBWflag = -1;
  }
  else if(nsel ==  5){
    fDecay = 120; theMH = 120; theWidth =  3.47e-03; theBWflag = -1;
  }
  else if(nsel ==  6){
    fDecay = 1120; theMH = 120; theWidth =  3.47e-03; theBWflag = -1;
  }
  else if(nsel ==  7){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; theBWflag = -1;
  }
  else if(nsel ==  8){
    fDecay = 1125; theMH = 125; theWidth =  4.03e-03; theBWflag = -1;
  }
  else if(nsel ==  9){
    fDecay = 130; theMH = 130; theWidth =  4.87e-03; theBWflag = -1;
  }
  else if(nsel == 10){
    fDecay = 1130; theMH = 130; theWidth =  4.87e-03; theBWflag = -1;
  }
  else if(nsel == 11){
    fDecay = 135; theMH = 135; theWidth =  6.14e-03; theBWflag = -1;
  }
  else if(nsel == 12){
    fDecay = 1135; theMH = 135; theWidth =  6.14e-03; theBWflag = -1;
  }
  else if(nsel == 13){
    fDecay = 140; theMH = 140; theWidth =  8.12e-03; theBWflag = -1;
  }
  else if(nsel == 14){
    fDecay = 1140; theMH = 140; theWidth =  8.12e-03; theBWflag = -1;
  }
  else if(nsel == 15){
    fDecay = 145; theMH = 145; theWidth =  1.14e-02; theBWflag = -1;
  }
  else if(nsel == 16){
    fDecay = 1145; theMH = 145; theWidth =  1.14e-02; theBWflag = -1;
  }
  else if(nsel == 17){
    fDecay = 150; theMH = 150; theWidth =  1.73e-02; theBWflag = -1;
  }
  else if(nsel == 18){
    fDecay = 1150; theMH = 150; theWidth =  1.73e-02; theBWflag = -1;
  }
  else if(nsel == 19){
    fDecay = 155; theMH = 155; theWidth =  3.02e-02; theBWflag = -1;
  }
  else if(nsel == 20){
    fDecay = 1155; theMH = 155; theWidth =  3.02e-02; theBWflag = -1;
  }
  else if(nsel == 21){
    fDecay = 160; theMH = 160; theWidth =  8.29e-02; theBWflag = -1;
  }
  else if(nsel == 22){
    fDecay = 1160; theMH = 160; theWidth =  8.29e-02; theBWflag = -1;
  }
  else if(nsel == 23){
    fDecay = 170; theMH = 170; theWidth =  3.80e-01; theBWflag = -1;
  }
  else if(nsel == 24){
    fDecay = 1170; theMH = 170; theWidth =  3.80e-01; theBWflag = -1;
  }
  else if(nsel == 25){
    fDecay = 180; theMH = 180; theWidth =  6.31e-01; theBWflag = -1;
  }
  else if(nsel == 26){
    fDecay = 1180; theMH = 180; theWidth =  6.31e-01; theBWflag = -1;
  }
  else if(nsel == 27){
    fDecay = 190; theMH = 190; theWidth =  1.04e+00; theBWflag = -1;
  }
  else if(nsel == 28){
    fDecay = 1190; theMH = 190; theWidth =  1.04e+00; theBWflag = -1;
  }
  else if(nsel == 29){
    fDecay = 200; theMH = 200; theWidth =  1.43e+00; theBWflag = -1;
  }
  else if(nsel == 30){
    fDecay = 1200; theMH = 200; theWidth =  1.43e+00; theBWflag = -1;
  }
  else if(nsel == 31){
    fDecay = 250; theMH = 250; theWidth =  4.04e+00; theBWflag = 1;
  }
  else if(nsel == 32){
    fDecay = 1250; theMH = 250; theWidth =  4.04e+00; theBWflag = 0;
  }
  else if(nsel == 33){
    fDecay = 300; theMH = 300; theWidth =  8.43e+00; theBWflag = 1;
  }
  else if(nsel == 34){
    fDecay = 1300; theMH = 300; theWidth =  8.43e+00; theBWflag = 0;
  }
  else if(nsel == 35){
    fDecay = 350; theMH = 350; theWidth =  1.52e+01; theBWflag = 1;
  }
  else if(nsel == 36){
    fDecay = 1350; theMH = 350; theWidth =  1.52e+01; theBWflag = 0;
  }
  else if(nsel == 37){
    fDecay = 400; theMH = 400; theWidth =  2.92e+01; theBWflag = 1;
  }
  else if(nsel == 38){
    fDecay = 1400; theMH = 400; theWidth =  2.92e+01; theBWflag = 0;
  }
  else if(nsel == 39){
    fDecay = 450; theMH = 450; theWidth =  4.69e+01; theBWflag = 1;
  }
  else if(nsel == 40){
    fDecay = 1450; theMH = 450; theWidth =  4.69e+01; theBWflag = 0;
  }
  else if(nsel == 41){
    fDecay = 500; theMH = 500; theWidth =  6.80e+01; theBWflag = 1;
  }
  else if(nsel == 42){
    fDecay = 1500; theMH = 500; theWidth =  6.80e+01; theBWflag = 0;
  }
  else if(nsel == 43){
    fDecay = 550; theMH = 550; theWidth =  9.31e+01; theBWflag = 1;
  }
  else if(nsel == 44){
    fDecay = 1550; theMH = 550; theWidth =  9.31e+01; theBWflag = 0;
  }
  else if(nsel == 45){
    fDecay = 600; theMH = 600; theWidth =  1.23e+02; theBWflag = 1;
  }
  else if(nsel == 46){
    fDecay = 1600; theMH = 600; theWidth =  1.23e+02; theBWflag = 0;
  }
  else if(nsel == 47){
    fDecay = 700; theMH = 700; theWidth =  1.99E+02; theBWflag = 1;
  }
  else if(nsel == 48){
    fDecay = 1700; theMH = 700; theWidth =  1.99E+02; theBWflag = 0;
  }
  else if(nsel == 49){
    fDecay = 800; theMH = 800; theWidth =  3.04E+02; theBWflag = 1;
  }
  else if(nsel == 50){
    fDecay = 1800; theMH = 800; theWidth =  3.04E+02; theBWflag = 0;
  }
  else if(nsel == 51){
    fDecay = 900; theMH = 900; theWidth =  4.49E+02; theBWflag = 1;
  }
  else if(nsel == 52){
    fDecay = 1900; theMH = 900; theWidth =  4.49E+02; theBWflag = 0;
  }
  else if(nsel == 53){
    fDecay = 1000; theMH = 1000; theWidth =  6.47E+02; theBWflag = 1;
  }
  else if(nsel == 54){
    fDecay = 2000; theMH = 1000; theWidth =  6.47E+02; theBWflag = 0;
  }
  else if(nsel == 55){
    fDecay = 2110;
  }
  else if(nsel == 56){
    fDecay = 2115;
  }
  else if(nsel == 57){
    fDecay = 2120;
  }
  else if(nsel == 58){
    fDecay = 2125;
  }
  else if(nsel == 59){
    fDecay = 2130;
  }
  else if(nsel == 60){
    fDecay = 2135;
  }
  else if(nsel == 61){
    fDecay = 2140;
  }
  else if(nsel == 62){
    fDecay = 2145;
  }
  else if(nsel == 63){
    fDecay = 2150;
  }
  else if(nsel == 64){
    fDecay = 2155;
  }
  else if(nsel == 65){
    fDecay = 2160;
  }
  else if(nsel == 66){
    fDecay = 2170;
  }
  else if(nsel == 67){
    fDecay = 2180;
  }
  else if(nsel == 68){
    fDecay = 2190;
  }
  else if(nsel == 69){
    fDecay = 2200;
  }
  else if(nsel == 70){
    fDecay = 2250;
  }
  else if(nsel == 71){
    fDecay = 2300;
  }
  else if(nsel == 72){
    fDecay = 118;
  }
  else if(nsel == 73){
    fDecay = 1118;
  }
  else if(nsel == 74){
    fDecay = 122;
  }
  else if(nsel == 75){
    fDecay = 1122;
  }
  else if(nsel == 76){
    fDecay = 124;
  }
  else if(nsel == 77){
    fDecay = 1124;
  }
  else if(nsel == 78){
    fDecay = 126;
  }
  else if(nsel == 79){
    fDecay = 1126;
  }
  else if(nsel == 80){
    fDecay = 128;
  }
  else if(nsel == 81){
    fDecay = 1128;
  }
  else if(nsel ==  91){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-0m-7tev-dkr"; theBWflag = -1;
  }
  else if(nsel ==  92){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-0p-7tev-dkr"; theBWflag = -1;
  }
  else if(nsel ==  93){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-2p-7tev-dkr"; theBWflag = -1;
  }
  else if(nsel ==  94){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-0m-8tev-dkr"; theBWflag = -1;
  }
  else if(nsel ==  95){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-0p-8tev-dkr"; theBWflag = -1;
  }
  else if(nsel ==  96){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; myRootFile = "histo_s12-x125ww-2p-8tev-dkr"; theBWflag = -1;
  }
  else if(nsel == 200){
    fDecay = 0;
  }
  else if(nsel == 201){
    fDecay = 1;
  }
  else if(nsel == 202){
    fDecay = 2;
  }
  else if(nsel == 203){
    fDecay = 3;
    isPhotonMCSel = true;
  }
  else if(nsel == 204){
    fDecay = 5;
  }
  else if(nsel == 205){
    fDecay = 5;
  }
  else if(nsel == 206){
    fDecay = 19;
    isPhotonMCSel = true;
  }
  else if(nsel == 207){
    fDecay = 19;
  }
  else if(nsel == 208){
    fDecay = 19;
  }
  else if(nsel == 209){
    fDecay = 19;
  }
  else if(nsel == 210){
    fDecay = 19;
  }
  else if(nsel == 211){
    fDecay = 19;
  }
  else if(nsel == 215){
    fDecay = 15;
  }
  else if(nsel == 212){
    fDecay = 25;
  }
  else if(nsel == 218){ // ww powheg
    fDecay = 32;
  }
  else if(nsel == 219){
    fDecay = 28;
  }
  else if(nsel == 220){
    fDecay = 29;
  }
  else if(nsel == 221){
    fDecay = 30;
  }
  else if(nsel == 222){
    fDecay = 27;
  }
  else if(nsel == 223){ // ww 2j
    fDecay = 33;
    applyFilterBTEvents = kTRUE;
  }
  else if(nsel == 225){ // ww mcatnlo
    fDecay = 29;
    processid = 999;
  }
  else if(nsel == 226){ // wwup mcatnlo
    fDecay = 30; // this is a hack
    processid = 999;
  }
  else if(nsel == 227){ // wwdown mcatnlo
    fDecay = 31; // this is a hack
    processid = 999;
  }
  else if(nsel == 228){
    fDecay = 28;
  }
  else if(nsel == 229){
    fDecay = 28;
  }
  else if(nsel == 230){
    fDecay = 28;
  }
  else if(nsel == 231){
    fDecay = 6;
  }
  else if(nsel == 232){
    fDecay = 6;
  }
  else if(nsel == 233){
    fDecay = 7;
  }
  else if(nsel == 234){
    fDecay = 7;
  }
  else if(nsel == 235){
    fDecay = 8;
  }
  else if(nsel == 236){
    fDecay = 8;
  }
  else if(nsel == 237){
    fDecay = 9;
  }
  else if(nsel == 238){
    fDecay = 19;
  }
  else if(nsel == 239){
    fDecay = 125; theMH = 125; theWidth =  4.03e-03; theBWflag = -1;
    processid = 998;
  }
  else if(nsel == 240){
    fDecay = 11;
  }
  else if(nsel == 241){
    fDecay = 11;
  }
  else if(nsel == 242){
    fDecay = 12;
  }
  else if(nsel == 243){
    fDecay = 12;
  }
  else if(nsel == 244){
    fDecay = 13;
  }
  else if(nsel == 245){
    fDecay = 13;
  }
  else if(nsel == 246){
    fDecay = 13;
  }
  else if(nsel == 247){
    fDecay = 13;
  }
  else if(nsel == 248){
    fDecay = 14;
    //applyWWFilter = kTRUE;
  }
  else if(nsel == 250){
    fDecay = 10;
    fIntRadius = 0.05;
    processid = 997;
  }
  else if(nsel == 251){
    fDecay = 10;
    fIntRadius = 0.05;
    processid = 997;
  }
  else if(nsel == 252){
    fDecay = 10;
    fIntRadius = 0.05;
    processid = 997;
  }
  else if(nsel == 253){
    applyVGFilter = kTRUE;
    fDecay = 20;
  }
  else if(nsel == 255){
    fDecay = 19;
  }
  else if(nsel == 299){
    myRootFile = "histo_ww_new";
    //files[0]   = "/data/blue/cmsprod/test2_sync_Summer12/WWIncl_py.root";
    //files[0]   = "/castor/cern.ch/user/p/paus/filefi/025/f11-wwj-v14b-bp/622DF074-A0F1-E011-BF1C-E0CB4EA0A91E.root";
    files[0]   = "/data/blue/cmsprod/test3_sync_Summer12/A6DE4085-8191-E111-BF4E-001E67396D5B.root";
    //myRootFile = "histo_1ldata";
    //isData = true;
    //isDataMuonElectron = true;
  }
  else if(nsel == 300){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataMuonElectron = true;
  }
  else if(nsel == 301){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataDMuon = true;
  }
  else if(nsel == 302){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataSMuon = true;
  }
  else if(nsel == 303){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataDElectron = true;
  }
  else if(nsel == 304){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataSElectron = true;
  }
  else if(nsel == 305){
    files[0]   = "*.root";
    is2011Corr = true;
    isData = true;
    isDataPhoton = true;
  }
  else if(nsel == 330){
    files[0]   = "*.root";
    isData = true; isDataMuonElectron = true; isDataDMuon = true;isDataSMuon = true;isDataDElectron = true;isDataSElectron = true;
    fDecay = 10;
    fIntRadius = 0.05;
    processid = 997;
  }
  else if(nsel == 350){
    files[0]   = "*.root";
    isData = true;
    isDataMuonElectron = true;
  }
  else if(nsel == 351){
    files[0]   = "*.root";
    isData = true;
    isDataDMuon = true;
  }
  else if(nsel == 352){
    files[0]   = "*.root";
    isData = true;
    isDataSMuon = true;
  }
  else if(nsel == 353){
    files[0]   = "*.root";
    isData = true;
    isDataDElectron = true;
  }
  else if(nsel == 354){
    files[0]   = "*.root";
    isData = true;
    isDataSElectron = true;
  }
  else if(nsel == 355){
    files[0]   = "*.root";
    isData = true;
    isDataPhoton = true;
  }
  else {
    printf("ERROR, WRONG OPTION: %d\n",nsel);
    return;
  }

  // Generator info
  GeneratorMod *GeneratorMod1 = new GeneratorMod;
  GeneratorMod1->SetPrintDebug(kFALSE);
  GeneratorMod1->SetPtLeptonMin(0.0);
  GeneratorMod1->SetEtaLeptonMax(2.7);
  GeneratorMod1->SetPtPhotonMin(15.0);
  GeneratorMod1->SetEtaPhotonMax(2.7);
  GeneratorMod1->SetPtRadPhotonMin(10.0);
  GeneratorMod1->SetEtaRadPhotonMax(2.7);
  GeneratorMod1->SetIsData(isData);
  GeneratorMod1->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMaxCut(50.);
  }
  else if(applyZVetoGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMinCut(20000.);
    GeneratorMod1->SetMassMaxCut(20000.);
  }
  GeneratorMod1->SetApplyISRFilter(applyISRFilter);
  GeneratorMod1->SetApplyVVFilter(applyWWFilter);
  GeneratorMod1->SetApplyVGFilter(applyVGFilter);
  GeneratorMod1->SetFilterBTEvents(applyFilterBTEvents);

  PartonFlavorHistoryMod *PartonFlavorHistoryMod1 = new PartonFlavorHistoryMod;
  PartonFlavorHistoryMod1->SetMCSampleType(MCType);
  PartonFlavorHistoryMod1->SetApplyPartonFlavorFilter(applyPartonFlavorFilter);

  // HLT info
  HLTMod *hltmod = new HLTMod;
  if(isData == true && isDataMuonElectron == true) {
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",163262,164237);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",163262,164237);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",165085,165888);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",165085,165888);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v5",165900,166967);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v5",165900,166967);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v6",166968,170053);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v6",166968,170053);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v8",170054,173198);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3",170054,173198);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4",173199,178380);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4",173199,178380);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7",178381,179889);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7",178381,179889);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8",179890,180000);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8",179890,180000);

    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",190456,190738);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",190456,190738);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",190739,191419);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",190739,191419);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",191420,193686);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",191420,193686);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",193687,197669);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",193687,197669);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",197770,199631);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",197770,199631);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",199632,999999);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",199632,999999);
  }
  if(isData == true && isDataDMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&HLT_DoubleMu7_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&HLT_DoubleMu7_v1",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&HLT_DoubleMu7_v2",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&HLT_Mu13_Mu8_v2" ,165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&HLT_Mu13_Mu8_v2" ,165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&HLT_Mu13_Mu8_v4" ,167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&HLT_Mu13_Mu8_v6" ,170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&HLT_Mu13_Mu8_v7" ,173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&HLT_Mu17_Mu8_v10" ,178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&HLT_Mu17_TkMu8_v3" ,178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&HLT_Mu17_Mu8_v11"  ,179890,180000);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&HLT_Mu17_TkMu8_v4" ,179890,180000);

    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&HLT_Mu17_Mu8_v16"  ,190456,190738);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&HLT_Mu17_TkMu8_v9" ,190456,190738);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&HLT_Mu17_Mu8_v16"  ,190739,191419);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&HLT_Mu17_TkMu8_v9" ,190739,191419);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&HLT_Mu17_Mu8_v16"  ,191420,193686);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&HLT_Mu17_TkMu8_v9" ,191420,193686);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Mu17_Mu8_v17"  ,193687,196045);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Mu17_TkMu8_v10",193687,196045);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Mu17_Mu8_v18"  ,196046,197669);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Mu17_TkMu8_v11",196046,197669);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&HLT_Mu17_Mu8_v19"  ,197770,199631);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&HLT_Mu17_TkMu8_v12",197770,199631);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Mu17_Mu8_v21"  ,199632,205256);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Mu17_TkMu8_v13",199632,205256);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Mu17_Mu8_v22"  ,205257,999999);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Mu17_TkMu8_v14",205257,999999);
  }
  if(isData == true && isDataSMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&HLT_Mu15_v2",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&HLT_Mu15_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_Mu24_v2",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_IsoMu17_v6",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&HLT_Mu30_v3",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&HLT_IsoMu17_v8",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&HLT_Mu30_v3",165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&HLT_IsoMu17_v9",165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v4&HLT_Mu30_v5",167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v4&HLT_IsoMu17_eta2p1_v1",167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&HLT_Mu40_v5",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&HLT_IsoMu24_v8",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_Mu40_eta2p1_v1",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_IsoMu30_eta2p1_v3",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_IsoMu24_eta2p1_v3",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_Mu40_eta2p1_v4",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_IsoMu30_eta2p1_v6",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_IsoMu24_eta2p1_v6",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_Mu40_eta2p1_v5",179890,180000);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_IsoMu30_eta2p1_v7",179890,180000);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_IsoMu24_eta2p1_v7",179890,180000);

    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&HLT_IsoMu24_eta2p1_v11" ,190456,190738);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&HLT_IsoMu24_eta2p1_v12" ,190739,191419);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&HLT_IsoMu24_eta2p1_v12" ,191420,193686);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v17&!HLT_Mu17_TkMu8_v10&HLT_IsoMu24_eta2p1_v13",193687,196045);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v18&!HLT_Mu17_TkMu8_v11&HLT_IsoMu24_eta2p1_v13",196046,197669);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu17_Mu8_v19&!HLT_Mu17_TkMu8_v12&HLT_IsoMu24_eta2p1_v14",197770,199631);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v21&!HLT_Mu17_TkMu8_v13&HLT_IsoMu24_eta2p1_v15",199632,205256);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v22&!HLT_Mu17_TkMu8_v14&HLT_IsoMu24_eta2p1_v15",205257,999999);
  }
  if(isData == true && isDataDElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&!HLT_Mu30_v3&!HLT_IsoMu17_v8&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&!HLT_Mu40_eta2p1_v1&!HLT_IsoMu30_eta2p1_v3&!HLT_IsoMu24_eta2p1_v3&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&!HLT_Mu40_eta2p1_v4&!HLT_IsoMu30_eta2p1_v6&!HLT_IsoMu24_eta2p1_v6&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&!HLT_Mu40_eta2p1_v5&!HLT_IsoMu30_eta2p1_v7&!HLT_IsoMu24_eta2p1_v7&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,180000);

    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v11&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15" ,190456,190738);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v12&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16" ,190739,191419);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v12&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17" ,191420,193686);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v17&!HLT_Mu17_TkMu8_v10&!HLT_IsoMu24_eta2p1_v13&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",193687,196045);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v18&!HLT_Mu17_TkMu8_v11&!HLT_IsoMu24_eta2p1_v13&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",196046,197669);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu17_Mu8_v19&!HLT_Mu17_TkMu8_v12&!HLT_IsoMu24_eta2p1_v14&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",197770,199631);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v21&!HLT_Mu17_TkMu8_v13&!HLT_IsoMu24_eta2p1_v15&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",199632,205256);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v22&!HLT_Mu17_TkMu8_v14&!HLT_IsoMu24_eta2p1_v15&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19",205257,999999);
  }
  if(isData == true && isDataSElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&!HLT_Mu30_v3&!HLT_IsoMu17_v8&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4&HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5&HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4",165900,166967);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6&HLT_Ele52_CaloIdVT_TrkIdT_v3",166968,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Ele65_CaloIdVT_TrkIdT_v3",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&!HLT_Mu40_eta2p1_v1&!HLT_IsoMu30_eta2p1_v3&!HLT_IsoMu24_eta2p1_v3&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&HLT_Ele65_CaloIdVT_TrkIdT_v4",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&!HLT_Mu40_eta2p1_v4&!HLT_IsoMu30_eta2p1_v6&!HLT_IsoMu24_eta2p1_v6&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Ele80_CaloIdVT_TrkIdT_v2",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&!HLT_Mu40_eta2p1_v5&!HLT_IsoMu30_eta2p1_v7&!HLT_IsoMu24_eta2p1_v7&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10&HLT_Ele80_CaloIdVT_TrkIdT_v3",179890,180000);

    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v11&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15&HLT_Ele27_WP80_v8",  190456,190738);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v12&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16&HLT_Ele27_WP80_v9",  190739,191419);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Mu17_Mu8_v16&!HLT_Mu17_TkMu8_v9&!HLT_IsoMu24_eta2p1_v12&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17&HLT_Ele27_WP80_v10", 191420,193686);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v17&!HLT_Mu17_TkMu8_v10&!HLT_IsoMu24_eta2p1_v13&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17&HLT_Ele27_WP80_v10",193687,196045);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&!HLT_Mu17_Mu8_v18&!HLT_Mu17_TkMu8_v11&!HLT_IsoMu24_eta2p1_v13&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17&HLT_Ele27_WP80_v10",196046,197669);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&!HLT_Mu17_Mu8_v19&!HLT_Mu17_TkMu8_v12&!HLT_IsoMu24_eta2p1_v14&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18&HLT_Ele27_WP80_v11",197770,199631);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v21&!HLT_Mu17_TkMu8_v13&!HLT_IsoMu24_eta2p1_v15&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18&HLT_Ele27_WP80_v11",199632,205256);
    hltmod->AddTrigger("!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&!HLT_Mu17_Mu8_v22&!HLT_Mu17_TkMu8_v14&!HLT_IsoMu24_eta2p1_v15&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19&HLT_Ele27_WP80_v11",205257,999999);
  }
  if(isData == true && isDataPhoton == true) {
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_v*", 0,999999);
  }
  if(isData == false){
    if(nsel == 299) {
      hltmod->AddTrigger("HLT_Mu15_v9");
      hltmod->AddTrigger("!HLT_Mu15_v9");
      hltmod->AddTrigger("HLT_Mu15_v2");
      hltmod->AddTrigger("!HLT_Mu15_v2");
      hltmod->AddTrigger("HLT_Mu9");
      hltmod->AddTrigger("!HLT_Mu9");
      hltmod->AddTrigger("HLT_Mu12_v13");
      hltmod->AddTrigger("!HLT_Mu12_v13");
      hltmod->AddTrigger("HLT_Mu12_v14");
      hltmod->AddTrigger("!HLT_Mu12_v14");
      hltmod->AddTrigger("HLT_Mu12_v16");
      hltmod->AddTrigger("!HLT_Mu12_v16");
    }
    else if(nsel != 999) {
      hltmod->AddTrigger("HLT_Mu15_v9");
      hltmod->AddTrigger("!HLT_Mu15_v9");
      hltmod->AddTrigger("HLT_Mu15_v1");
      hltmod->AddTrigger("!HLT_Mu15_v1");
      hltmod->AddTrigger("HLT_Mu15_v2");
      hltmod->AddTrigger("!HLT_Mu15_v2");
      hltmod->AddTrigger("HLT_Mu12_v13");
      hltmod->AddTrigger("!HLT_Mu12_v13");
      hltmod->AddTrigger("HLT_Mu12_v14");
      hltmod->AddTrigger("!HLT_Mu12_v14");
      hltmod->AddTrigger("HLT_Mu12_v16");
      hltmod->AddTrigger("!HLT_Mu12_v16");
      /*
      hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*");
      hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*");
      hltmod->AddTrigger("HLT_Mu17_Mu8_v*");
      hltmod->AddTrigger("HLT_Mu17_TkMu8_v*");
      hltmod->AddTrigger("HLT_IsoMu24_eta2p1_v*");
      hltmod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      hltmod->AddTrigger("HLT_Ele27_WP80_v*");
      */
    }
    else {
      hltmod->AddTrigger("HLT_Mu9");
      hltmod->AddTrigger("HLT_Mu11");
      hltmod->AddTrigger("HLT_Mu15_v1");
      hltmod->AddTrigger("HLT_Ele10_LW_L1R");
      hltmod->AddTrigger("HLT_Ele15_SW_L1R");
      hltmod->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R");
      hltmod->AddTrigger("HLT_Ele17_SW_CaloEleId_L1R");
      hltmod->AddTrigger("HLT_Ele17_SW_L1R");
      hltmod->AddTrigger("HLT_Ele17_SW_TightEleId_L1R");
      hltmod->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
      hltmod->AddTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
  }
  hltmod->SetTrigObjsName("myhltobjs");

  HLTMod *hltmodFR = new HLTMod;
  if(isData == false) {
    hltmodFR->AddTrigger("HLT_Mu17_v*");
    hltmodFR->AddTrigger("HLT_Mu8_v*");
    hltmodFR->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*");
    hltmodFR->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    hltmodFR->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*");
    hltmodFR->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  }
  else if(isData == true && isDataMuonElectron == true) {
    hltmodFR->AddTrigger("dummy");
  }
  else if(isData == true && isDataDMuon == true) {
    hltmodFR->AddTrigger("HLT_Mu17_v*");
    hltmodFR->AddTrigger("HLT_Mu8_v*");
  }
  else if(isData == true && isDataSMuon == true) {
    hltmodFR->AddTrigger("dummy");
  }
  else if(isData == true && isDataDElectron == true) {
    hltmodFR->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*");
    hltmodFR->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    hltmodFR->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*");
    hltmodFR->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  }
  else if(isData == true && isDataSElectron == true) {
    hltmodFR->AddTrigger("dummy");
  }
  hltmodFR->SetTrigObjsName("myhltobjsFR");

  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelection = new RunLumiSelectionMod;      
  runLumiSelection->SetAcceptMC(!isData);
  runLumiSelection->SetAcceptAll(kTRUE);
  if(nsel >= 350 || nsel <= 360){
    runLumiSelection->AddJSONFile(string(getenv("CMSSW_BASE")+string("/src/json/Cert_Current_JSON.txt")));
  }
  else {
    runLumiSelection->AddJSONFile(string(getenv("CMSSW_BASE")+string("/src/json/hww.Full2011.json")));
  }
  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  Bool_t isFastSim = kFALSE;

  TString rootFileHKF = TString(theOutputDir);
  rootFileHKF += TString(myRootFile);
  rootFileHKF += TString("_pdf_") + TString(fileset);
  rootFileHKF += TString("_") + TString("noskim");
  rootFileHKF += TString(".root");
  HKFactorProducer *HKFactorProducer1 = new HKFactorProducer;
  HKFactorProducer1->SetProcessID(processid);
  HKFactorProducer1->SetInputFilename(fInputFilenameKF);
  HKFactorProducer1->SetIsData(isData);
  HKFactorProducer1->SetMakePDFNtuple(!isData);
  HKFactorProducer1->SetFillHist(!isData);
  HKFactorProducer1->SetOutputName(rootFileHKF);
  HKFactorProducer1->SetMh(theMH);
  HKFactorProducer1->SetWidth(theWidth);
  HKFactorProducer1->SetBWflag(theBWflag);

  // PFNoPU
  SeparatePileUpMod *separatePileUpMod1 = new SeparatePileUpMod;

  // Special lepton id
  MuonIDMod *muonIDBS = new MuonIDMod;
  muonIDBS->SetClassType("GlobalTracker");
  muonIDBS->SetIDType(theWWMuId);
  muonIDBS->SetIsoType(theWWMuIso);
  muonIDBS->SetApplyD0Cut(kTRUE);
  muonIDBS->SetApplyDZCut(kFALSE);
  muonIDBS->SetWhichVertex(-2);
  muonIDBS->SetCleanMuonsName("CleanMuonsBS");
  muonIDBS->SetIntRadius(fIntRadius);
  muonIDBS->SetRhoType(theRhoType);

  ElectronIDMod *electronIDBS = new ElectronIDMod;
  electronIDBS->SetIDType("VBTFWorkingPointLowPtId");
  electronIDBS->SetIsoType("PFIso");
  electronIDBS->SetApplyConversionFilterType1(kTRUE);
  electronIDBS->SetApplyConversionFilterType2(kFALSE);
  electronIDBS->SetChargeFilter(kFALSE);
  electronIDBS->SetApplyD0Cut(kTRUE);
  electronIDBS->SetApplyDZCut(kFALSE);
  electronIDBS->SetWhichVertex(-2);
  electronIDBS->SetNExpectedHitsInnerCut(0);
  electronIDBS->SetGoodElectronsName("GoodElectronsBS");
  electronIDBS->SetIntRadius(fIntRadius);
  electronIDBS->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaningBS = new ElectronCleaningMod;
  electronCleaningBS->SetCleanMuonsName("CleanMuonsBS");
  electronCleaningBS->SetGoodElectronsName("GoodElectronsBS");
  electronCleaningBS->SetCleanElectronsName("CleanElectronsBS");

  MergeLeptonsMod *mergerBS = new MergeLeptonsMod;
  mergerBS->SetMuonsName(muonIDBS->GetOutputName());
  mergerBS->SetElectronsName(electronCleaningBS->GetOutputName());
  mergerBS->SetMergedName("MergedLeptonsBS");

  // Object ID and Cleaning Sequence
  MuonIDMod *muonID1 = new MuonIDMod;
  muonID1->SetIntRadius(fIntRadius);
  muonID1->SetClassType("GlobalTracker");
  //muonID1->SetIDType("MVA_BDTG_IDIso");
  //muonID1->SetIsoType("MVA_BDTG_IDIso");
  muonID1->SetIDType(theWWMuId);
  muonID1->SetIsoType(theWWMuIso);
  muonID1->SetApplyD0Cut(kTRUE);
  muonID1->SetApplyDZCut(kTRUE);
  if(nsel > 400){ // 42x
    muonID1->SetD0Cut(0.20);
    muonID1->SetDZCut(0.50);
  }
  muonID1->SetWhichVertex(0);
  muonID1->SetRhoType(theRhoType);

  ElectronIDMod *electronID1 = new ElectronIDMod;
  electronID1->SetIntRadius(fIntRadius);
  //electronID1->SetIDType("VBTFWorkingPointLowPtId");
  //electronID1->SetIsoType("PFIso");
  //electronID1->SetIDType("MVA_BDTG_IDIsoCombined");
  //electronID1->SetIsoType("MVA_BDTG_IDIsoCombined");
  //electronID1->SetElectronMVAWeightsSubdet0Pt10To20(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml")));
  //electronID1->SetElectronMVAWeightsSubdet1Pt10To20(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml")));
  //electronID1->SetElectronMVAWeightsSubdet2Pt10To20(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml")));
  //electronID1->SetElectronMVAWeightsSubdet0Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml")));
  //electronID1->SetElectronMVAWeightsSubdet1Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml")));
  //electronID1->SetElectronMVAWeightsSubdet2Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml")));
  if(is2011Corr == true){
    electronID1->SetIDType("VBTFWorkingPoint80Id");
    electronID1->SetIsoType("CombinedRelativeConeAreaCorrected");
  } else {
    electronID1->SetIDType("CustomLoose");
    electronID1->SetIsoType("PFIso_HWW2012TrigV0");
  }
  electronID1->SetChargeFilter(kFALSE);
  electronID1->SetApplyD0Cut(kTRUE);
  electronID1->SetApplyDZCut(kTRUE);
  electronID1->SetWhichVertex(0);
  electronID1->SetApplyConversionFilterType2(kFALSE);
  electronID1->SetApplyConversionFilterType1(kTRUE);
  if(nsel > 400){ // 42x
    electronID1->SetNExpectedHitsInnerCut(0);
  } else {
    electronID1->SetNExpectedHitsInnerCut(1);
  }
  electronID1->SetRhoType(theRhoType);
  electronID1->SetCombRelativeIsoCut(0.10);

  PhotonIDMod *photonIDMod1 = new PhotonIDMod;
  photonIDMod1->SetIDType("BaseLineCiC");
  photonIDMod1->SetIsoType("MITPUCorrected");
  photonIDMod1->SetApplyPixelSeed(kFALSE);
  photonIDMod1->SetApplyElectronVetoConvRecovery(kTRUE);
  photonIDMod1->SetApplyConversionId(kTRUE);
  photonIDMod1->SetApplyFiduciality(kFALSE);
  photonIDMod1->SetPtMin(10.0);

  PFTauIDMod *pftauIDMod1 = new PFTauIDMod;
  pftauIDMod1->SetPFTausName("HPSTaus");
  pftauIDMod1->SetIsHPSSel(kFALSE);

  const char *jetInput1 = "AKt5PFJets";
  PublisherMod<PFJet,Jet> *pubJet1 = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet1->SetInputName(jetInput1);
  pubJet1->SetOutputName(Form("Pub%s",jetInput1));

  JetCorrectionMod *jetCorr1_ntuple = new JetCorrectionMod;
  if     (is2011Corr == true) {
    jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt"))); 
    jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L2Relative.txt"))); 
    jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt")));
    if(isData == true){
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt")));
    }
  }
  else if(is2011Corr == false){
    if     (isData == true){
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L1FastJet_AK5PF.txt"))); 
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L2Relative_AK5PF.txt"))); 
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L3Absolute_AK5PF.txt")));
      if(isData == true){
        jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L2L3Residual_AK5PF.txt")));
      }
    }
    else if(isData == false){
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L1FastJet_AK5PF.txt"))); 
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L2Relative_AK5PF.txt"))); 
      jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L3Absolute_AK5PF.txt")));
      if(isData == true){
        jetCorr1_ntuple->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L2L3Residual_AK5PF.txt")));
      }
    }
  }
  jetCorr1_ntuple->SetInputName(pubJet1->GetOutputName());
  jetCorr1_ntuple->SetCorrectedName("CorrectedJets_ntuple");
  jetCorr1_ntuple->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaning1 = new ElectronCleaningMod;
  PhotonCleaningMod *photonCleaningMod1 = new PhotonCleaningMod;
  PFTauCleaningMod *pftauCleaningMod1 = new PFTauCleaningMod;

  MergeLeptonsMod *merger1 = new MergeLeptonsMod;
  merger1->SetMuonsName(muonID1->GetOutputName());
  merger1->SetElectronsName(electronCleaning1->GetOutputName());

  JetIDMod *theJetID2_ntuple = new JetIDMod;
  theJetID2_ntuple->SetInputName(jetCorr1_ntuple->GetOutputName());
  theJetID2_ntuple->SetPtCut(0.0);
  theJetID2_ntuple->SetEtaMaxCut(etaJetCut);
  theJetID2_ntuple->SetJetEEMFractionMinCut(0.00);
  theJetID2_ntuple->SetOutputName("GoodJetsNoPtCut_ntuple");
  theJetID2_ntuple->SetApplyBetaCut(kFALSE);
  theJetID2_ntuple->SetApplyMVACut(applyMVACut);

  JetCleaningMod *theJetCleaning2_ntuple = new JetCleaningMod;
  theJetCleaning2_ntuple->SetGoodJetsName("GoodJetsNoPtCut_ntuple");
  theJetCleaning2_ntuple->SetCleanJetsName("CleanJetsNoPtCut_ntuple");

  const char *metPFInput = "PFMet";
  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName(metPFInput);
  pubPFMet->SetOutputName(Form("Pub%s",metPFInput));

  // Analyses modules
  ZXEvtSelMod *zxEvtSelMod1 = new ZXEvtSelMod;
  zxEvtSelMod1->SetPrintDebug(kFALSE);
  zxEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  zxEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  zxEvtSelMod1->SetPtJetCut(ptJetCut);
  zxEvtSelMod1->SetEtaJetCut(etaJetCut);
  zxEvtSelMod1->SetUsePDFs(usePDFProducer);

  ttEvtSelMod *ttEvtSelMod1 = new ttEvtSelMod;
  ttEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  ttEvtSelMod1->SetPrintDebug(kFALSE);
  ttEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  ttEvtSelMod1->SetPtJetCut(ptJetCut);
  ttEvtSelMod1->SetEtaJetCut(etaJetCut);

  ttljetsEvtSelMod *ttljetsEvtSelMod1 = new ttljetsEvtSelMod;
  ttljetsEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  ttljetsEvtSelMod1->SetPrintDebug(kFALSE);
  ttljetsEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  ttljetsEvtSelMod1->SetPtJetCut(30.0);
  ttljetsEvtSelMod1->SetEtaJetCut(etaJetCut);
  ttljetsEvtSelMod1->SetTrigObjsName("myhltobjs");

  WlnEvtSelMod *WlnEvtSelMod1 = new WlnEvtSelMod;
  WlnEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  WlnEvtSelMod1->SetPrintDebug(kFALSE);
  WlnEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  WlnEvtSelMod1->SetPtJetCut(ptJetCut);
  WlnEvtSelMod1->SetEtaJetCut(etaJetCut);

  LeptonEvtSelMod *LeptonEvtSelMod1 = new LeptonEvtSelMod;
  LeptonEvtSelMod1->SetPrintDebug(kTRUE);
  LeptonEvtSelMod1->SetIsData(isData);
  LeptonEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  LeptonEvtSelMod1->SetAllVertexName("PrimaryVertexes");

  WtnEvtSelMod *WtnEvtSelMod1 = new WtnEvtSelMod;
  WtnEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  WtnEvtSelMod1->SetPrintDebug(kFALSE);
  WtnEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  WtnEvtSelMod1->SetPFTaus0Name("ShrinkingConePFTaus");
  WtnEvtSelMod1->SetPFTaus1Name("HPSTaus");

  ZllEvtSelMod *ZllEvtSelMod1 = new ZllEvtSelMod;
  ZllEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  ZllEvtSelMod1->SetTrigObjsName("myhltobjs");
  ZllEvtSelMod1->SetPrintDebug(kFALSE);
  ZllEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  ZllEvtSelMod1->SetPtJetCut(ptJetCut);
  ZllEvtSelMod1->SetEtaJetCut(etaJetCut);

  ZttEvtSelMod *ZttEvtSelMod1 = new ZttEvtSelMod;
  ZttEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  ZttEvtSelMod1->SetTrigObjsName("myhltobjs");
  ZttEvtSelMod1->SetPrintDebug(kFALSE);
  ZttEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  ZttEvtSelMod1->SetPtJetCut(ptJetCut);
  ZttEvtSelMod1->SetEtaJetCut(etaJetCut);

  HwwExampleAnalysisMod *HwwExampleAnalysisMod1 = new HwwExampleAnalysisMod;
  HwwExampleAnalysisMod1->SetMetName(pubPFMet->GetOutputName());
  HwwExampleAnalysisMod1->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");

  WWEvtSelMod *WWEvtSelMod1 = new WWEvtSelMod;
  WWEvtSelMod1->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  if(nsel == 25 || nsel == 299)
    WWEvtSelMod1->SetPrintDebug(kTRUE);
  else
    WWEvtSelMod1->SetPrintDebug(kFALSE);
  WWEvtSelMod1->SetIsFastSim(isFastSim);
  WWEvtSelMod1->SetIsData(isData);
  WWEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  WWEvtSelMod1->SetJetScaleSyst(0.0);
  WWEvtSelMod1->SetPtJetCut(ptJetCut);
  WWEvtSelMod1->SetEtaJetCut(etaJetCut);
  WWEvtSelMod1->SetUsePDFs(usePDFProducer);
  WWEvtSelMod1->SetAllVertexName("PrimaryVertexes");
  WWEvtSelMod1->SetIntRadius(fIntRadius);
  WWEvtSelMod1->SetIsOldSelection(is2011Corr);

  AAWWEvtSelMod *AAWWEvtSelMod1 = new AAWWEvtSelMod;
  AAWWEvtSelMod1->SetPrintDebug(kFALSE);
  AAWWEvtSelMod1->SetIsData(isData);

  SkimEvtSelMod *SkimEvtSelMod1 = new SkimEvtSelMod;
  SkimEvtSelMod1->SetTrigObjsName("myhltobjs");

  GammaXEvtSelMod *gammaXEvtSelMod1 = new GammaXEvtSelMod;
  gammaXEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  gammaXEvtSelMod1->SetGoodJetsName("GoodJetsNoPtCut_ntuple");
  gammaXEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  gammaXEvtSelMod1->SetPtJetCut(ptJetCut);
  gammaXEvtSelMod1->SetEtaJetCut(etaJetCut);

  WlnFakeSelMod *WlnFakeSelMod1 = new WlnFakeSelMod;
  WlnFakeSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  WlnFakeSelMod1->SetMetName(pubPFMet->GetOutputName());
  WlnFakeSelMod1->SetPtJetCut(ptJetCut);
  WlnFakeSelMod1->SetEtaJetCut(etaJetCut);

  TString rootFileHwwMake0 = TString(theOutputDir);
  rootFileHwwMake0 += TString(myRootFile);
  rootFileHwwMake0 += TString("_smurf0_") + TString(fileset);
  rootFileHwwMake0 += TString("_") + TString("noskim");
  rootFileHwwMake0 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod0 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod0->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod0->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod0->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod0->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod0->SetProcessID(fDecay);
  HwwMakeNtupleMod0->SetIsData(isData);
  HwwMakeNtupleMod0->SetFakeRatePredictionType(0);
  HwwMakeNtupleMod0->SetFillNtupleType(0);
  HwwMakeNtupleMod0->SetOutputName(rootFileHwwMake0);
  HwwMakeNtupleMod0->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod0->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod0->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod0->SetIntRadius(fIntRadius);
  HwwMakeNtupleMod0->SetIs42x(is2011Corr);
  HwwMakeNtupleMod0->SetMVAElVersion(1);
  HwwMakeNtupleMod0->SetMVAMuVersion(1);
  HwwMakeNtupleMod0->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  HwwMakeNtupleMod0->SetPFTauName(pftauCleaningMod1->GetOutputName());

  TString rootFileHwwMake1 = TString(theOutputDir);
  rootFileHwwMake1 += TString(myRootFile);
  rootFileHwwMake1 += TString("_smurf1_") + TString(fileset);
  rootFileHwwMake1 += TString("_") + TString("noskim");
  rootFileHwwMake1 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod1 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod1->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod1->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod1->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod1->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod1->SetProcessID(fDecay);
  HwwMakeNtupleMod1->SetIsData(isData);
  HwwMakeNtupleMod1->SetFakeRatePredictionType(1);
  HwwMakeNtupleMod1->SetFillNtupleType(1);
  HwwMakeNtupleMod1->SetOutputName(rootFileHwwMake1);
  HwwMakeNtupleMod1->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod1->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod1->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod1->SetIntRadius(fIntRadius);
  HwwMakeNtupleMod1->SetIs42x(is2011Corr);
  HwwMakeNtupleMod1->SetMVAElVersion(1);
  HwwMakeNtupleMod1->SetMVAMuVersion(1);
  HwwMakeNtupleMod1->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  HwwMakeNtupleMod1->SetPFTauName(pftauCleaningMod1->GetOutputName());

  TString rootFileHwwMake2 = TString(theOutputDir);
  rootFileHwwMake2 += TString(myRootFile);
  rootFileHwwMake2 += TString("_smurf2_") + TString(fileset);
  rootFileHwwMake2 += TString("_") + TString("noskim");
  rootFileHwwMake2 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod2 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod2->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod2->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod2->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod2->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod2->SetProcessID(fDecay);
  HwwMakeNtupleMod2->SetIsData(isData);
  HwwMakeNtupleMod2->SetFakeRatePredictionType(2);
  HwwMakeNtupleMod2->SetFillNtupleType(2);
  HwwMakeNtupleMod2->SetOutputName(rootFileHwwMake2);
  HwwMakeNtupleMod2->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod2->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod2->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod2->SetIntRadius(fIntRadius);
  HwwMakeNtupleMod2->SetIs42x(is2011Corr);
  HwwMakeNtupleMod2->SetMVAElVersion(1);
  HwwMakeNtupleMod2->SetMVAMuVersion(1);
  HwwMakeNtupleMod2->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  HwwMakeNtupleMod2->SetPFTauName(pftauCleaningMod1->GetOutputName());

  TString rootFileHwwMake3 = TString(theOutputDir);
  rootFileHwwMake3 += TString(myRootFile);
  rootFileHwwMake3 += TString("_smurf3_") + TString(fileset);
  rootFileHwwMake3 += TString("_") + TString("noskim");
  rootFileHwwMake3 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod3 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod3->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod3->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod3->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod3->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod3->SetProcessID(fDecay);
  HwwMakeNtupleMod3->SetIsData(isData);
  HwwMakeNtupleMod3->SetFakeRatePredictionType(3);
  HwwMakeNtupleMod3->SetFillNtupleType(3);
  HwwMakeNtupleMod3->SetOutputName(rootFileHwwMake3);
  HwwMakeNtupleMod3->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod3->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod3->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod3->SetIntRadius(fIntRadius);
  HwwMakeNtupleMod3->SetIs42x(is2011Corr);
  HwwMakeNtupleMod3->SetMVAElVersion(1);
  HwwMakeNtupleMod3->SetMVAMuVersion(1);
  HwwMakeNtupleMod3->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  HwwMakeNtupleMod3->SetPFTauName(pftauCleaningMod1->GetOutputName());

  TString rootFileHwwMake4 = TString(theOutputDir);
  rootFileHwwMake4 += TString(myRootFile);
  rootFileHwwMake4 += TString("_smurf4_") + TString(fileset);
  rootFileHwwMake4 += TString("_") + TString("noskim");
  rootFileHwwMake4 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod4 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod4->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod4->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod4->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod4->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod4->SetProcessID(fDecay);
  HwwMakeNtupleMod4->SetIsData(isData);
  HwwMakeNtupleMod4->SetFakeRatePredictionType(0);
  HwwMakeNtupleMod4->SetFillNtupleType(4);
  HwwMakeNtupleMod4->SetOutputName(rootFileHwwMake4);
  HwwMakeNtupleMod4->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod4->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod4->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod4->SetIntRadius(fIntRadius);
  HwwMakeNtupleMod4->SetIs42x(is2011Corr);
  HwwMakeNtupleMod4->SetMVAElVersion(1);
  HwwMakeNtupleMod4->SetMVAMuVersion(1);
  HwwMakeNtupleMod4->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  HwwMakeNtupleMod4->SetFillPhotonTemplate(kTRUE);
  HwwMakeNtupleMod4->SetPFTauName(pftauCleaningMod1->GetOutputName());

  // Lepton ID with no D0
  MuonIDMod *muonID2 = new MuonIDMod;
  muonID2->SetClassType("GlobalTracker");
  muonID2->SetIDType(theWWMuId);
  muonID2->SetIsoType("PFIso");
  muonID2->SetApplyD0Cut(kFALSE);
  muonID2->SetCleanMuonsName("CleanMuonsNoD0");
  muonID2->SetWhichVertex(0);
  muonID2->SetIntRadius(fIntRadius);
  muonID2->SetRhoType(theRhoType);

  ElectronIDMod *electronID2 = new ElectronIDMod;
  electronID2->SetIDType("VBTFWorkingPointLowPtId");
  electronID2->SetIsoType("PFIso");
  electronID2->SetApplyConversionFilterType1(kTRUE);
  electronID2->SetApplyConversionFilterType2(kFALSE);
  electronID2->SetChargeFilter(kFALSE);
  electronID2->SetApplyD0Cut(kFALSE);
  electronID2->SetNExpectedHitsInnerCut(0);
  electronID2->SetGoodElectronsName("GoodElectronsNoD0");
  electronID2->SetWhichVertex(0);
  electronID2->SetIntRadius(fIntRadius);
  electronID2->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaning2 = new ElectronCleaningMod;
  electronCleaning2->SetCleanMuonsName("CleanMuonsNoD0");
  electronCleaning2->SetGoodElectronsName("GoodElectronsNoD0");
  electronCleaning2->SetCleanElectronsName("CleanElectronsNoD0");

  MergeLeptonsMod *merger2 = new MergeLeptonsMod;
  merger2->SetMuonsName(muonID2->GetOutputName());
  merger2->SetElectronsName(electronCleaning2->GetOutputName());
  merger2->SetMergedName("MergedLeptonsNoD0");

  // Lepton ID with no Id
  MuonIDMod *muonID3 = new MuonIDMod;
  muonID3->SetClassType("GlobalTracker");
  muonID3->SetIDType("NoId");
  muonID3->SetIsoType("PFIso");
  muonID3->SetApplyD0Cut(kTRUE);
  muonID3->SetCleanMuonsName("CleanMuonsNoId");
  muonID3->SetWhichVertex(0);
  muonID3->SetIntRadius(fIntRadius);
  muonID3->SetRhoType(theRhoType);

  ElectronIDMod *electronID3 = new ElectronIDMod;
  electronID3->SetIDType("NoId");
  electronID3->SetIsoType("PFIso");
  electronID3->SetApplyConversionFilterType1(kTRUE);
  electronID3->SetApplyConversionFilterType2(kFALSE);
  electronID3->SetChargeFilter(kFALSE);
  electronID3->SetApplyD0Cut(kTRUE);
  electronID3->SetNExpectedHitsInnerCut(0);
  electronID3->SetGoodElectronsName("GoodElectronsNoId");
  electronID3->SetWhichVertex(0);
  electronID3->SetIntRadius(fIntRadius);
  electronID3->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaning3 = new ElectronCleaningMod;
  electronCleaning3->SetCleanMuonsName("CleanMuonsNoId");
  electronCleaning3->SetGoodElectronsName("GoodElectronsNoId");
  electronCleaning3->SetCleanElectronsName("CleanElectronsNoId");

  MergeLeptonsMod *merger3 = new MergeLeptonsMod;
  merger3->SetMuonsName(muonID3->GetOutputName());
  merger3->SetElectronsName(electronCleaning3->GetOutputName());
  merger3->SetMergedName("MergedLeptonsNoId");

  // Lepton ID with no Iso
  MuonIDMod *muonID4 = new MuonIDMod;
  muonID4->SetClassType("GlobalTracker");
  muonID4->SetIDType(theWWMuId);
  muonID4->SetIsoType("NoIso");
  muonID4->SetApplyD0Cut(kTRUE);
  muonID4->SetCleanMuonsName("CleanMuonsNoIso");
  muonID4->SetWhichVertex(0);
  muonID4->SetIntRadius(fIntRadius);
  muonID4->SetRhoType(theRhoType);

  ElectronIDMod *electronID4 = new ElectronIDMod;
  electronID4->SetIDType("VBTFWorkingPointLowPtId");
  electronID4->SetIsoType("NoIso");
  electronID4->SetApplyConversionFilterType1(kTRUE);
  electronID4->SetApplyConversionFilterType2(kFALSE);
  electronID4->SetChargeFilter(kFALSE);
  electronID4->SetApplyD0Cut(kTRUE);
  electronID4->SetNExpectedHitsInnerCut(0);
  electronID4->SetGoodElectronsName("GoodElectronsNoIso");
  electronID4->SetWhichVertex(0);
  electronID4->SetIntRadius(fIntRadius);
  electronID4->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaning4 = new ElectronCleaningMod;
  electronCleaning4->SetCleanMuonsName("CleanMuonsNoIso");
  electronCleaning4->SetGoodElectronsName("GoodElectronsNoIso");
  electronCleaning4->SetCleanElectronsName("CleanElectronsNoIso");

  MergeLeptonsMod *merger4 = new MergeLeptonsMod;
  merger4->SetMuonsName(muonID4->GetOutputName());
  merger4->SetElectronsName(electronCleaning4->GetOutputName());
  merger4->SetMergedName("MergedLeptonsNoIso");

  // Lepton ID with no Conversion Filter
  ElectronIDMod *electronID5 = new ElectronIDMod;
  electronID5->SetIDType("VBTFWorkingPointLowPtId");
  electronID5->SetIsoType("PFIso");
  electronID5->SetApplyConversionFilterType1(kFALSE);
  electronID5->SetApplyConversionFilterType2(kFALSE);
  electronID5->SetChargeFilter(kFALSE);
  electronID5->SetApplyD0Cut(kTRUE);
  electronID5->SetNExpectedHitsInnerCut(999);
  electronID5->SetGoodElectronsName("GoodElectronsNoConvF");
  electronID5->SetWhichVertex(0);
  electronID5->SetIntRadius(fIntRadius);
  electronID5->SetRhoType(theRhoType);

  ElectronCleaningMod *electronCleaning5 = new ElectronCleaningMod;
  electronCleaning5->SetCleanMuonsName("CleanMuons");
  electronCleaning5->SetGoodElectronsName("GoodElectronsNoConvF");
  electronCleaning5->SetCleanElectronsName("CleanElectronsNoConvF");

  MergeLeptonsMod *merger5 = new MergeLeptonsMod;
  merger5->SetMuonsName(muonID1->GetOutputName());
  merger5->SetElectronsName(electronCleaning5->GetOutputName());
  merger5->SetMergedName("MergedLeptonsNoConvF");

  // Lepton ID with loose requirements
  MuonIDMod *muonIDFakeable = new MuonIDMod;
  muonIDFakeable->SetClassType("GlobalTracker");
  muonIDFakeable->SetIDType(theWWMuId);
  muonIDFakeable->SetIsoType(theWWMuIso);
  muonIDFakeable->SetApplyD0Cut(kTRUE);
  muonIDFakeable->SetApplyDZCut(kTRUE);
  muonIDFakeable->SetD0Cut(0.20);
  muonIDFakeable->SetPFIsoCut(muIsoCut);
  muonIDFakeable->SetCleanMuonsName("CleanMuonsFakeable");
  muonIDFakeable->SetWhichVertex(0);
  muonIDFakeable->SetIntRadius(fIntRadius);
  muonIDFakeable->SetRhoType(theRhoType);

  ElectronIDMod *electronIDFakeable = new ElectronIDMod;
  electronIDFakeable->SetIntRadius(fIntRadius);
  electronIDFakeable->SetIDType("VBTFWorkingPointFakeableId");
  electronIDFakeable->SetIsoType("NoIso");
  electronIDFakeable->SetChargeFilter(kFALSE);
  electronIDFakeable->SetApplyD0Cut(kTRUE);
  electronIDFakeable->SetApplyDZCut(kTRUE);
  electronIDFakeable->SetWhichVertex(0);
  electronIDFakeable->SetApplyConversionFilterType2(kFALSE);
  electronIDFakeable->SetApplyConversionFilterType1(kTRUE);
  if(nsel > 400){ // 42x
    electronIDFakeable->SetNExpectedHitsInnerCut(0);
  } else {
    electronIDFakeable->SetNExpectedHitsInnerCut(1);
  }
  electronIDFakeable->SetElectronMVAWeightsSubdet0Pt10To20 (string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat1.weights.xml")));
  electronIDFakeable->SetElectronMVAWeightsSubdet1Pt10To20 (string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat2.weights.xml")));
  electronIDFakeable->SetElectronMVAWeightsSubdet2Pt10To20 (string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat3.weights.xml")));
  electronIDFakeable->SetElectronMVAWeightsSubdet0Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat4.weights.xml")));
  electronIDFakeable->SetElectronMVAWeightsSubdet1Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat5.weights.xml")));
  electronIDFakeable->SetElectronMVAWeightsSubdet2Pt20ToInf(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat6.weights.xml")));
  electronIDFakeable->SetRhoType(theRhoType);
  electronIDFakeable->SetGoodElectronsName("GoodElectronsFakeable");

  ElectronCleaningMod *electronCleaningFakeable = new ElectronCleaningMod;
  electronCleaningFakeable->SetCleanMuonsName(muonIDFakeable->GetOutputName());
  electronCleaningFakeable->SetGoodElectronsName(electronIDFakeable->GetOutputName());
  electronCleaningFakeable->SetCleanElectronsName("CleanElectronsFakeable");

  MergeLeptonsMod *mergerFakeable = new MergeLeptonsMod;
  mergerFakeable->SetMuonsName(muonIDFakeable->GetOutputName());
  mergerFakeable->SetElectronsName(electronCleaningFakeable->GetOutputName());
  mergerFakeable->SetMergedName("MergedLeptonsFakeable");

  FRStudy *FRStudy1 = new FRStudy;
  FRStudy1->SetCleanJetsName(theJetID2_ntuple->GetOutputName());
  FRStudy1->SetPrintDebug(kFALSE);
  FRStudy1->SetIsFastSim(isFastSim);
  FRStudy1->SetMetName(pubPFMet->GetOutputName());
  FRStudy1->SetPtJetCut(30.0);
  FRStudy1->SetEtaJetCut(etaJetCut);
  FRStudy1->SetIsData(isData);
  FRStudy1->SetSelectGenLeptons(useSelectGenLeptons);

  IsoStudy *isoStudy1 = new IsoStudy;
  isoStudy1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  isoStudy1->SetPrintDebug(kFALSE);
  isoStudy1->SetIsFastSim(isFastSim);
  isoStudy1->SetMetName(pubPFMet->GetOutputName());
  isoStudy1->SetPtJetCut(ptJetCut);
  isoStudy1->SetEtaJetCut(etaJetCut);

  QQLLEvtSelMod *QQLLEvtSelMod1 = new QQLLEvtSelMod;
  QQLLEvtSelMod1->SetPrintDebug(kFALSE);
  QQLLEvtSelMod1->SetIsFastSim(isFastSim);
  QQLLEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  QQLLEvtSelMod1->SetJetScaleSyst(0.0);
  QQLLEvtSelMod1->SetPtJetCut(30.0);
  QQLLEvtSelMod1->SetEtaJetCut(etaJetCut);
  QQLLEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");

  LowEvtSelMod *LowEvtSelMod1 = new LowEvtSelMod;
  LowEvtSelMod1->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  LowEvtSelMod1->SetPrintDebug(kFALSE);
  LowEvtSelMod1->SetIsFastSim(isFastSim);
  LowEvtSelMod1->SetMetName(pubPFMet->GetOutputName());
  LowEvtSelMod1->SetJetScaleSyst(0.0);
  LowEvtSelMod1->SetPtJetCut(ptJetCut);
  LowEvtSelMod1->SetEtaJetCut(etaJetCut);
  LowEvtSelMod1->SetTrigObjsName("myhltobjs");

  // Chain modules together
  if(isData == false){
    GeneratorMod1->Add(HKFactorProducer1);
    HKFactorProducer1->Add(PartonFlavorHistoryMod1);
    PartonFlavorHistoryMod1->Add(goodPVFilterMod);
  }
  else if(nsel == 299){
    GeneratorMod1->Add(HKFactorProducer1);
    HKFactorProducer1->Add(goodPVFilterMod);
  }
  else {
    GeneratorMod1->Add(HKFactorProducer1);
    HKFactorProducer1->Add(runLumiSelection);
    runLumiSelection->Add(goodPVFilterMod);
  }
  goodPVFilterMod->Add(separatePileUpMod1);
  separatePileUpMod1->Add(muonIDBS);
  muonIDBS->Add(electronIDBS);
  electronIDBS->Add(electronCleaningBS);
  electronCleaningBS->Add(mergerBS);
  mergerBS->Add(muonID1);
  muonID1->Add(electronID1);
  electronID1->Add(photonIDMod1);
  photonIDMod1->Add(pftauIDMod1);
  pftauIDMod1->Add(pubJet1);
  pubJet1->Add(jetCorr1_ntuple);
  jetCorr1_ntuple->Add(electronCleaning1);
  electronCleaning1->Add(photonCleaningMod1);
  photonCleaningMod1->Add(pftauCleaningMod1);
  pftauCleaningMod1->Add(theJetID2_ntuple);
  theJetID2_ntuple->Add(theJetCleaning2_ntuple);
  theJetCleaning2_ntuple->Add(merger1);
  merger1->Add(pubPFMet);
  pubPFMet->Add(WtnEvtSelMod1);
  WtnEvtSelMod1->Add(muonIDFakeable);
  muonIDFakeable->Add(electronIDFakeable);
  electronIDFakeable->Add(electronCleaningFakeable);
  electronCleaningFakeable->Add(mergerFakeable);
  mergerFakeable->Add(hltmod);
  hltmod->Add(gammaXEvtSelMod1);
  gammaXEvtSelMod1->Add(WWEvtSelMod1);
  WWEvtSelMod1->Add(AAWWEvtSelMod1);
  AAWWEvtSelMod1->Add(HwwExampleAnalysisMod1);
  HwwExampleAnalysisMod1->Add(ttEvtSelMod1);
  ttEvtSelMod1->Add(ttljetsEvtSelMod1);
  ttljetsEvtSelMod1->Add(muonID2);
  muonID2->Add(electronID2);
  electronID2->Add(electronCleaning2);
  electronCleaning2->Add(merger2);
  merger2->Add(muonID3);
  muonID3->Add(electronID3);
  electronID3->Add(electronCleaning3);
  electronCleaning3->Add(merger3);
  merger3->Add(muonID4);
  muonID4->Add(electronID4);
  electronID4->Add(electronCleaning4);
  electronCleaning4->Add(merger4);
  merger4->Add(electronID5);
  electronID5->Add(electronCleaning5);
  electronCleaning5->Add(merger5);
  merger5->Add(zxEvtSelMod1);
  zxEvtSelMod1->Add(ZllEvtSelMod1);
  ZllEvtSelMod1->Add(ZttEvtSelMod1);
  ZttEvtSelMod1->Add(SkimEvtSelMod1);
  SkimEvtSelMod1->Add(QQLLEvtSelMod1);
  QQLLEvtSelMod1->Add(LowEvtSelMod1);
  LowEvtSelMod1->Add(isoStudy1);
  isoStudy1->Add(WlnFakeSelMod1);
  WlnFakeSelMod1->Add(LeptonEvtSelMod1);
  LeptonEvtSelMod1->Add(WlnEvtSelMod1);
  if     (isDataPhoton == true){
    WlnEvtSelMod1    ->Add(HwwMakeNtupleMod4);
    HwwMakeNtupleMod4->Add(HwwMakeNtupleMod3);
  }
  else {
    WlnEvtSelMod1    ->Add(HwwMakeNtupleMod0);
    HwwMakeNtupleMod0->Add(HwwMakeNtupleMod1);
    if(nsel == 203 || nsel == 206 || nsel >= 300){
      HwwMakeNtupleMod1->Add(HwwMakeNtupleMod2);
      HwwMakeNtupleMod2->Add(HwwMakeNtupleMod3);
      if(isPhotonMCSel == true){
         cout << "isPhotonMCSel is true" << endl;
        HwwMakeNtupleMod3->Add(HwwMakeNtupleMod4);
      }
    }
  }

  mergerFakeable->Add(hltmodFR);
  hltmodFR->Add(FRStudy1);

  // Set up analysis
  Analysis *ana1 = new Analysis;
  ana1->SetUseHLT(kTRUE);
  ana1->SetKeepHierarchy(kFALSE);
  ana1->SetProcessNEvents(NEvents);
  ana1->SetSuperModule(GeneratorMod1);
  ana1->SetPrintScale(100);
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  if(TString(fileset) != TString("all")){
    printf("Rely on Catalog: %s\n",catalogDir);
    printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
    Catalog *c = new Catalog(catalogDir);
    Dataset *d = c->FindDataset(book,dataset,fileset);
    ana1->AddDataset(d);
  }
  else {
    printf("Reading all files at once\n");
    ana1->AddFile(files[0]);
    for(int i=1; i<Nfiles; i++){
      if(TString(files[i]) != TString("")) ana1->AddFile(files[i]);
    }
  }
  //ana1->SetCacheSize(0);
  ana1->SetCacheSize(64*1024*1024);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(theOutputDir);
  //TString rootFile = TString("");
  rootFile += TString(myRootFile);
  rootFile += TString("_") + TString(fileset);
  rootFile += TString("_") + TString("noskim");
  rootFile += TString(".root");

  printf("Root output: %s\n",rootFile.Data());  
  ana1->SetOutputName(rootFile.Data());

  // run the analysis after successful initialisation
  ana1->Run(!gROOT->IsBatch());  

  return;
}
