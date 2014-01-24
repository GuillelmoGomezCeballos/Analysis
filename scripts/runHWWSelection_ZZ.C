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
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "Ana/SelMods/interface/ZZEvtSelMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
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
    theWWMuIso   = "PFIsoBetaPUCorrected";
    theRhoType   = RhoUtilities::DEFAULT;
    is2011Corr   = true;
    muIsoCut     = -999.0;
    applyMVACut  = false;
    if(nsel != 299) nsel = nsel - 400;
  }

  double fIntRadius = 0.0;
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

  // full lepton id
  MuonIDMod *muonID1 = new MuonIDMod;
  muonID1->SetIntRadius(fIntRadius);
  muonID1->SetClassType("GlobalTracker");
  muonID1->SetIDType(theWWMuId);
  muonID1->SetIsoType(theWWMuIso);
  muonID1->SetApplyD0Cut(kTRUE);
  muonID1->SetApplyDZCut(kTRUE);
  muonID1->SetD0Cut(0.05);
  muonID1->SetDZCut(0.20);
  muonID1->SetPFIsoCut(0.25);
  muonID1->SetWhichVertex(0);
  muonID1->SetRhoType(theRhoType);
  muonID1->SetCleanMuonsName("CleanMuons");

  ElectronIDMod *electronID1 = new ElectronIDMod;
  electronID1->SetIntRadius(fIntRadius);
  electronID1->SetIDType("VBTFWorkingPoint95Id");
  electronID1->SetIsoType("CombinedRelativeConeAreaCorrected");
  electronID1->SetChargeFilter(kFALSE);
  electronID1->SetApplyD0Cut(kTRUE);
  electronID1->SetApplyDZCut(kTRUE);
  electronID1->SetD0Cut(0.05);
  electronID1->SetDZCut(0.20);
  electronID1->SetWhichVertex(0);
  electronID1->SetApplyConversionFilterType2(kFALSE);
  electronID1->SetApplyConversionFilterType1(kTRUE);
  electronID1->SetNExpectedHitsInnerCut(1);
  electronID1->SetRhoType(theRhoType);
  electronID1->SetCombRelativeIsoCut(0.25);

  ElectronCleaningMod *electronCleaning1 = new ElectronCleaningMod;
  electronCleaning1->SetCleanMuonsName(muonID1->GetOutputName());
  electronCleaning1->SetGoodElectronsName(electronID1->GetOutputName());
  electronCleaning1->SetCleanElectronsName("CleanElectrons");

  MergeLeptonsMod *merger1 = new MergeLeptonsMod;
  merger1->SetMuonsName(muonID1->GetOutputName());
  merger1->SetElectronsName(electronCleaning1->GetOutputName());
  merger1->SetMergedName("MergedLeptons");

  PhotonIDMod *photonIDMod1 = new PhotonIDMod;
  photonIDMod1->SetIDType("BaseLineCiC");
  photonIDMod1->SetIsoType("MITPUCorrected");
  photonIDMod1->SetApplyPixelSeed(kFALSE);
  photonIDMod1->SetApplyElectronVetoConvRecovery(kTRUE);
  photonIDMod1->SetApplyConversionId(kTRUE);
  photonIDMod1->SetApplyFiduciality(kFALSE);
  photonIDMod1->SetPtMin(10.0);

  PhotonCleaningMod *photonCleaningMod1 = new PhotonCleaningMod;

  // fakeable lepton id
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
  electronIDFakeable->SetNExpectedHitsInnerCut(1);
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

  const char *jetInput1 = "AKt5PFJets";
  PublisherMod<PFJet,Jet> *pubJet1 = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet1->SetInputName(jetInput1);
  pubJet1->SetOutputName(Form("Pub%s",jetInput1));

  JetCorrectionMod *jetCorr1 = new JetCorrectionMod;
  if     (is2011Corr == true) {
    jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt"))); 
    jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L2Relative.txt"))); 
    jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt")));
    if(isData == true){
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt")));
    }
  }
  else if(is2011Corr == false){
    if     (isData == true){
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L1FastJet_AK5PF.txt"))); 
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L2Relative_AK5PF.txt"))); 
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L3Absolute_AK5PF.txt")));
      if(isData == true){
        jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/GR_P_V42_AN3_L2L3Residual_AK5PF.txt")));
      }
    }
    else if(isData == false){
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L1FastJet_AK5PF.txt"))); 
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L2Relative_AK5PF.txt"))); 
      jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L3Absolute_AK5PF.txt")));
      if(isData == true){
        jetCorr1->AddCorrectionFromFile(string(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/START53_V15_L2L3Residual_AK5PF.txt")));
      }
    }
  }
  jetCorr1->SetInputName(pubJet1->GetOutputName());
  jetCorr1->SetCorrectedName("CorrectedJets");
  jetCorr1->SetRhoType(theRhoType);

  JetIDMod *theJetID1 = new JetIDMod;
  theJetID1->SetInputName(jetCorr1->GetOutputName());
  theJetID1->SetPtCut(0.0);
  theJetID1->SetEtaMaxCut(4.7);
  theJetID1->SetJetEEMFractionMinCut(0.00);
  theJetID1->SetOutputName("GoodJets");
  theJetID1->SetApplyBetaCut(kFALSE);
  theJetID1->SetApplyMVACut(applyMVACut);

  JetCleaningMod *theJetCleaning1 = new JetCleaningMod;
  theJetCleaning1->SetGoodJetsName(theJetID1->GetOutputName());
  theJetCleaning1->SetCleanJetsName("CleanJets");

  // analysis module
  ZZEvtSelMod *ZZEvtSelMod1 = new ZZEvtSelMod;
  ZZEvtSelMod1->SetPrintDebug(kFALSE);
  ZZEvtSelMod1->SetMuonName        (muonID1->GetOutputName());
  ZZEvtSelMod1->SetElectronName    (electronCleaning1->GetOutputName());
  ZZEvtSelMod1->SetLeptonName      (merger1->GetOutputName());
  ZZEvtSelMod1->SetJetName         (theJetCleaning1->GetOutputName());
  ZZEvtSelMod1->SetMuonFakeName    (muonIDFakeable->GetOutputName());
  ZZEvtSelMod1->SetElectronFakeName(electronCleaningFakeable->GetOutputName());
  ZZEvtSelMod1->SetLeptonFakeName  (mergerFakeable->GetOutputName());
  ZZEvtSelMod1->SetIsData          (isData);

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
  goodPVFilterMod->Add(hltmod);
  hltmod->Add(separatePileUpMod1);
  separatePileUpMod1->Add(muonID1);
  muonID1->Add(electronID1);
  electronID1->Add(electronCleaning1);
  electronCleaning1->Add(merger1);
  merger1->Add(photonIDMod1);
  photonIDMod1->Add(photonCleaningMod1);
  photonCleaningMod1->Add(pubJet1);
  pubJet1->Add(jetCorr1);
  jetCorr1->Add(theJetID1);
  theJetID1->Add(theJetCleaning1);
  theJetCleaning1->Add(muonIDFakeable);
  muonIDFakeable->Add(electronIDFakeable);
  electronIDFakeable->Add(electronCleaningFakeable);
  electronCleaningFakeable->Add(mergerFakeable);
  mergerFakeable->Add(ZZEvtSelMod1);

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
