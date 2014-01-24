//****************************************************************************************************************
//Test a file 
//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"rfio:/castor/cern.ch/user/p/paus/filler/013/c10-minb-goodcoll-v7/c10-minb-goodcoll-v7_000_25.root\",\"TestOutput.root\",10000,1\) 
//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"rfio:/castor/cern.ch/user/s/sixie/BAMBU/013/s10-minb7-su3v26a/outfilename_000_70.root\",\"TestOutput.root\",10000,1\) 
//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"/uscms/home/sxie/resilient/BAMBU/013/c10-minb-rra01/outfilename_000_1_1.root\",\"TestOutput\",10000,1\)
//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"/home/sixie/download/BAMBU/SingleLeptonSkims/c10-minb-prv8/SingleLeptonDataSkim_SingleLeptonSkim_c10-minb-prv8_148.root_000.root\",\"TestOutput\",10000,-1\)


//****************************************************************************************************************

//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"0000\",\"noskim\",\"c10-minb-goodcoll-v7\",\"cern/filler/013/\",\"/home/mitprod/catalog\",\"Commissioning\",10000,-1\)
//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runCommissioning/runCommissioning.C+\(\"0000\",\"noskim\",\"p10-minb7000-st3\",\"cern/filler/013/\",\"/home/mitprod/catalog\",\"Commissioning\",10000,1\)





 



#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <exception>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/L1Mod.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h" 
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/FakeMods/interface/GenFakeableObjsMod.h"
#include "MitPhysics/FakeMods/interface/GenFakesMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/PhysicsMod/interface/RunSelectionMod.h"
#include "MitHiggs/Commissioning/interface/ElectronCommissioning.h"
#include "MitHiggs/Commissioning/interface/MuonCommissioning.h"
#include "MitHiggs/Commissioning/interface/JetCommissioning.h"
#include "MitHiggs/Commissioning/interface/DileptonCommissioning.h"
#include "MitHiggs/Commissioning/interface/WSelectionCommissioning.h"
#include "MitHiggs/FakeRateMods/interface/ComputeElectronFakeRateMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void executeCommissioning(const char *fileset,
                        const char *skim     ,
                        const char *dataset  ,
                        const char *book     ,
                        const char *catalogDir,
                        const char *outputName ,
                        int   sampleID        ,
                        int   nEvents        )
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
//   gDebugMask  = Debug::kAnalysis;
  gDebugMask  = Debug::kAll;
  gDebugLevel = 50;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  Bool_t isFastSim = kFALSE;
  if(sampleID >= 100) isFastSim = kTRUE;

  Bool_t isData = kFALSE;
  if (sampleID == -1) {
    isData = kTRUE;
    cout << "RUNNING ON DATA\n";
  }

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;

  //------------------------------------------------------------------------------------------------
  // L1 information
  //------------------------------------------------------------------------------------------------
  L1Mod *l1mod = new L1Mod;
  l1mod->SetBitsName("L1TechBitsBeforeMask");
  l1mod->SetPrintTable(kTRUE);
  l1mod->AddTrigger("L1Tech_BSC_minBias_threshold1.v0&!L1Tech_BSC_halo_beam1_inner.v0&!L1Tech_BSC_halo_beam1_outer.v0&!L1Tech_BSC_halo_beam2_inner.v0&!L1Tech_BSC_halo_beam2_outer.v0");  
  l1mod->AddTrigger("L1Tech_BSC_minBias_threshold2.v0&!L1Tech_BSC_halo_beam1_inner.v0&!L1Tech_BSC_halo_beam1_outer.v0&!L1Tech_BSC_halo_beam2_inner.v0&!L1Tech_BSC_halo_beam2_outer.v0");

  L1Mod *l1algomod = new L1Mod;
  l1algomod->SetBitsName("L1AlgoBitsBeforeMask");
  l1algomod->AddTrigger("L1_BptxMinus&L1_BptxPlus");


  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;
  hltmod->SetPrintTable(kTRUE);
  hltmod->AddTrigger("HLT_PhysicsDeclared");
//   hltmod->AddTrigger("HLT_Ele15_SW_L1R");
//   hltmod->AddTrigger("HLT_Mu15");
  hltmod->SetTrigObjsName("myhltobjs");

  //------------------------------------------------------------------------------------------------
  // Run SelectionMod
  //------------------------------------------------------------------------------------------------
  RunSelectionMod * runSelection = new RunSelectionMod;
  runSelection->SetDefaultAccept(kTRUE);
//   runSelection->ExcludeRun(132442); //ECAL Scan
//   runSelection->ExcludeRun(132476); //Pixel/HF delay scan
//   runSelection->ExcludeRun(132477); //Pixel/HF delay scan
//   runSelection->ExcludeRun(132478); //CSC/ECAL delay scan


  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<CaloJet,Jet> *pubJet = new PublisherMod<CaloJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5Jets");
  pubJet->SetOutputName("AKt5Jets");

  PublisherMod<CaloMet,Met> *pubMet = new PublisherMod<CaloMet,Met>("MetPub");
  pubMet->SetInputName("CorMuonMet");
  pubMet->SetOutputName(ModNames::gkCleanCaloMetName);

  PublisherMod<Met,Met> *pubTCMet = new PublisherMod<Met,Met>("TCMetPub");
  pubTCMet->SetInputName("TCMet");
  pubTCMet->SetOutputName("pubTCMet");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("PFMetPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("pubPFMet");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonID           = new MuonIDMod;  
  muonID->SetIDType("Loose");
  muonID->SetIsoType("TrackCaloSliding"); 
  ElectronIDMod       *electronID       = new ElectronIDMod;
  electronID->SetIDType(TString("CustomTight"));
  electronID->SetIsoType(TString("TrackJuraSliding"));
  PhotonIDMod         *photonID       = new PhotonIDMod;
  photonID->SetIDType(TString("Custom"));
  photonID->SetPtMin(20.0);
  TauIDMod *tauID = new TauIDMod;
  JetIDMod            *jetID            = new JetIDMod;
  jetID->SetInputName(pubJet->GetOutputName());
  jetID->SetUseCorrection(kTRUE);
  jetID->SetPtCut(35.0);
  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;

  //------------------------------------------------------------------------------------------------
  // Commissioning modules
  //------------------------------------------------------------------------------------------------
  ElectronCommissioning *electronCommissioningAnalysis = new ElectronCommissioning;
  electronCommissioningAnalysis->SetName("ElectronCommissioningMod");
  electronCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  electronCommissioningAnalysis->SetSampleType(sampleID);
  
  MuonCommissioning *muonCommissioningAnalysis = new MuonCommissioning;
  muonCommissioningAnalysis->SetName("MuonCommissioningMod");
  muonCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  muonCommissioningAnalysis->SetSampleType(sampleID);

  JetCommissioning *jetCommissioningAnalysis = new JetCommissioning;
  jetCommissioningAnalysis->SetName("JetCommissioningMod");
  jetCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  jetCommissioningAnalysis->SetSampleType(sampleID);

  DileptonCommissioning *dileptonCommissioningAnalysis = new DileptonCommissioning;
  dileptonCommissioningAnalysis->SetName("DileptonCommissioningMod");
  dileptonCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  dileptonCommissioningAnalysis->SetSampleType(sampleID);

  WSelectionCommissioning *wSelectionCommissioningAnalysis = new WSelectionCommissioning;
  wSelectionCommissioningAnalysis->SetName("WSelectionCommissioning");
  wSelectionCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  wSelectionCommissioningAnalysis->SetSampleType(sampleID);

  //------------------------------------------------------------------------------------------------
  // Compute Fake Rate
  //------------------------------------------------------------------------------------------------
  ComputeElectronFakeRateMod *computeElectronFakeRate = new ComputeElectronFakeRateMod;
  computeElectronFakeRate->SetName("ComputeElectronFakeRateMod");
  computeElectronFakeRate->SetMetName(pubMet->GetOutputName());
  computeElectronFakeRate->SetTriggerType(0);
  computeElectronFakeRate->SetTriggerName("");
  computeElectronFakeRate->SetTriggerObjectsName("myhltobjs");
  computeElectronFakeRate->SetIsData(isData);


  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  l1mod->Add(l1algomod);
  l1algomod->Add(hltmod);
  hltmod->Add(runSelection); 
  runSelection->Add(pubJet);
  pubJet->Add(pubMet);
  pubMet->Add(pubTCMet);
  pubTCMet->Add(pubPFMet);
  pubPFMet->Add(muonID);
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(jetID);
  jetID->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning->Add(tauCleaning);
  tauCleaning->Add(jetCleaning);
  jetCleaning->Add(mergeLeptonsMod);
  mergeLeptonsMod->Add(electronCommissioningAnalysis);
  mergeLeptonsMod->Add(muonCommissioningAnalysis);
  mergeLeptonsMod->Add(jetCommissioningAnalysis);
  mergeLeptonsMod->Add(dileptonCommissioningAnalysis);
  mergeLeptonsMod->Add(wSelectionCommissioningAnalysis);
  mergeLeptonsMod->Add(computeElectronFakeRate);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  if (isData) {
    ana->SetSuperModule(l1mod);
  } else {
    ana->SetSuperModule(pubJet);
  }
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);

  cout << book << " " << skimdataset.Data() << " " << fileset << endl;
  

  //ana->AddFile("castor:/castor/cern.ch/cms/store/caf/user/frankma/MinimumBias-Express/132440-veryloosecuts/bambu_000_10.root");


  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  ana->SetOutputName(rootFile.Data());

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());  

  return;
}



//--------------------------------------------------------------------------------------------------
void executeCommissioning(const char *inputfile, 
                          const char *outputName,
                          int   sampleID       ,
                          int   nEvents       )
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  Bool_t isFastSim = kFALSE;
  if(sampleID >= 100) isFastSim = kTRUE;
  Bool_t isData = kFALSE;
  if (sampleID == -1) {
    isData =kTRUE;
    cout << "RUNNING ON DATA\n";
  }

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;

  //------------------------------------------------------------------------------------------------
  // L1 information
  //------------------------------------------------------------------------------------------------
  L1Mod *l1mod = new L1Mod;
  l1mod->SetBitsName("L1TechBitsBeforeMask");
  l1mod->SetPrintTable(kTRUE);
  l1mod->AddTrigger("L1Tech_BSC_minBias_threshold1.v0&!L1Tech_BSC_halo_beam1_inner.v0&!L1Tech_BSC_halo_beam1_outer.v0&!L1Tech_BSC_halo_beam2_inner.v0&!L1Tech_BSC_halo_beam2_outer.v0");  
  l1mod->AddTrigger("L1Tech_BSC_minBias_threshold2.v0&!L1Tech_BSC_halo_beam1_inner.v0&!L1Tech_BSC_halo_beam1_outer.v0&!L1Tech_BSC_halo_beam2_inner.v0&!L1Tech_BSC_halo_beam2_outer.v0");

  L1Mod *l1algomod = new L1Mod;
  l1algomod->SetBitsName("L1AlgoBitsBeforeMask");
  l1algomod->AddTrigger("L1_BptxMinus&L1_BptxPlus");


  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;
  hltmod->SetPrintTable(kTRUE);
  hltmod->AddTrigger("HLT_Ele15_SW_L1R");
  hltmod->AddTrigger("HLT_Mu15");
  hltmod->SetTrigObjsName("myhltobjs");

  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<CaloJet,Jet> *pubJet = new PublisherMod<CaloJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5Jets");
  pubJet->SetOutputName("AKt5Jets");

  PublisherMod<CaloMet,Met> *pubMet = new PublisherMod<CaloMet,Met>("MetPub");
  pubMet->SetInputName("CorMuonMet");
  pubMet->SetOutputName(ModNames::gkCleanCaloMetName);

  PublisherMod<Met,Met> *pubTCMet = new PublisherMod<Met,Met>("TCMetPub");
  pubTCMet->SetInputName("TCMet");
  pubTCMet->SetOutputName("pubTCMet");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("PFMetPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("pubPFMet");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonID           = new MuonIDMod;  
  muonID->SetIDType("Loose");
  muonID->SetIsoType("TrackCaloSliding"); 
  ElectronIDMod       *electronID       = new ElectronIDMod;
  electronID->SetIDType(TString("CustomTight"));
  electronID->SetIsoType(TString("TrackJuraSliding"));
  PhotonIDMod         *photonID       = new PhotonIDMod;
  photonID->SetIDType(TString("Custom"));
  photonID->SetPtMin(20.0);
  TauIDMod *tauID = new TauIDMod;
  JetIDMod            *jetID            = new JetIDMod;
  jetID->SetInputName(pubJet->GetOutputName());
  jetID->SetUseCorrection(kTRUE);
  jetID->SetPtCut(35.0);
  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;

  //------------------------------------------------------------------------------------------------
  // analyses modules
  //------------------------------------------------------------------------------------------------
  ElectronCommissioning *electronCommissioningAnalysis = new ElectronCommissioning;
  electronCommissioningAnalysis->SetName("ElectronCommissioningMod");
  electronCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  electronCommissioningAnalysis->SetSampleType(sampleID);
  
  MuonCommissioning *muonCommissioningAnalysis = new MuonCommissioning;
  muonCommissioningAnalysis->SetName("MuonCommissioningMod");
  muonCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  muonCommissioningAnalysis->SetSampleType(sampleID);

  JetCommissioning *jetCommissioningAnalysis = new JetCommissioning;
  jetCommissioningAnalysis->SetName("JetCommissioningMod");
  jetCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  jetCommissioningAnalysis->SetSampleType(sampleID);

  DileptonCommissioning *dileptonCommissioningAnalysis = new DileptonCommissioning;
  dileptonCommissioningAnalysis->SetName("DileptonCommissioningMod");
  dileptonCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  dileptonCommissioningAnalysis->SetSampleType(sampleID);

  WSelectionCommissioning *wSelectionCommissioningAnalysis = new WSelectionCommissioning;
  wSelectionCommissioningAnalysis->SetName("WSelectionCommissioning");
  wSelectionCommissioningAnalysis->SetMetName(pubMet->GetOutputName());
  wSelectionCommissioningAnalysis->SetSampleType(sampleID);

  //------------------------------------------------------------------------------------------------
  // Compute Fake Rate
  //------------------------------------------------------------------------------------------------
  ComputeElectronFakeRateMod *computeElectronFakeRate = new ComputeElectronFakeRateMod;
  computeElectronFakeRate->SetName("ComputeElectronFakeRateMod");
  computeElectronFakeRate->SetMetName(pubMet->GetOutputName());
  computeElectronFakeRate->SetTriggerType(0);
  computeElectronFakeRate->SetTriggerName("");
  computeElectronFakeRate->SetTriggerObjectsName("myhltobjs");
  computeElectronFakeRate->SetIsData(isData);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  l1mod->Add(l1algomod);
  l1algomod->Add(pubJet);
  pubJet->Add(pubMet);
  pubMet->Add(pubTCMet);
  pubTCMet->Add(pubPFMet);
  pubPFMet->Add(muonID);
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(jetID);
  jetID->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning->Add(tauCleaning);
  tauCleaning->Add(jetCleaning);
  jetCleaning->Add(mergeLeptonsMod);
  mergeLeptonsMod->Add(electronCommissioningAnalysis);
  mergeLeptonsMod->Add(muonCommissioningAnalysis);
  mergeLeptonsMod->Add(jetCommissioningAnalysis);
  mergeLeptonsMod->Add(dileptonCommissioningAnalysis);
  mergeLeptonsMod->Add(wSelectionCommissioningAnalysis);
  mergeLeptonsMod->Add(computeElectronFakeRate);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  if (isData) {
    ana->SetSuperModule(l1mod);
  } else {
    ana->SetSuperModule(pubJet);
  }
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------

  ana->AddFile(inputfile);
//   ana->AddFile("/home/sixie/download/BAMBU/SingleLeptonSkims/c10-minb-prv8/SingleLeptonDataSkim_SingleLeptonSkim_c10-minb-prv8_99.root_000.root");




  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());  

  return;
}



//--------------------------------------------------------------------------------------------------
void runCommissioning(const char *fileset     ,
                       const char *skim        ,
                       const char *dataset    ,
                       const char *book        ,
                       const char *catalogDir  ,
                       const char *outputName ,
                       int         nEvents    , 
                       int         runTypeIndex)
{
  TString outfileName = TString(outputName);
  executeCommissioning(fileset,skim,dataset,book,catalogDir,outfileName,runTypeIndex,nEvents); 
}

//--------------------------------------------------------------------------------------------------
void runCommissioning(const char *inputfiles, 
                            const char *outputName,
                            int         nEvents,
                            int         runTypeIndex )
{
  TString outfileName = TString(outputName);
  executeCommissioning(inputfiles,outputName, runTypeIndex, nEvents);
}
