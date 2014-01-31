//
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/LP2011/mitf_noweights/data_2fake.root\",\"/data/smurf/ceballos/test/data_2fake.root\",-1,\"/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_EPS2011.txt\"\) 
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/LP2011/mitf_noweights/data_2l.root\",\"/data/smurf/ceballos/test/data_2l.root\",-1,\"/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_EPS2011.txt\"\) 
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/LP2011/mitf_noweights/data_lfake.root\",\"/data/smurf/ceballos/test/data_lfake.root\",-1,\"/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_EPS2011.txt\"\) 
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/LP2011/mitf_noweights/background42x.root\",\"/data/smurf/ceballos/test/background42x.root\",-1,\"/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_EPS2011.txt\"\) 
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/LP2011/mitf_noweights/hww150.root\",\"/data/smurf/ceballos/test/test.root\",-1,\"/home/ceballos/releases/CMSSW_4_2_2/src/json/Cert_EPS2011.txt\"\) 
//root -l -b -q Analysis/nt_scripts/SkimGoodLumiList.C+\(\"/data/smurf/data/Run2011_Spring11_SmurfV7_42X/tas/data.root\",\"data.root\",-1,\"/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json\"\) 

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Smurf/Core/SmurfTree.h"
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname)
{
  bool verbose = true;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  //TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwMakeNtupleMod");
  //assert(dir);

  TTree* t;
  t = (TTree*) inf->Get("HwwTree0");
  if(!t) t = (TTree*) inf->Get("HwwTree1");
  if(!t) t = (TTree*) inf->Get("HwwTree2");
  if(!t) t = (TTree*) inf->Get("tree");
  assert(t);
  t->SetName("tree");

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void SkimGoodLumiList(const string InputFilename, const string OutputFilename, Int_t type,
                      const string GoodLumiFile) {

  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(GoodLumiFile.c_str()); 

  TTree* HwwTree = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree);
  SmurfTree sigEvent;
  sigEvent.LoadTree(InputFilename.c_str(),type);
  sigEvent.InitTree(0);
  sigEvent.tree_->SetName("tree");

  //*************************************************************************************************
  //Create new normalized tree
  //*************************************************************************************************
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = sigEvent.tree_->CloneTree(0);
  
  cout << "Events in the ntuple: " << sigEvent.tree_->GetEntries() << endl;
  
  for (int n=0;n<sigEvent.tree_->GetEntries();n++) { 
    sigEvent.tree_->GetEntry(n);
    if (n%100000==0) cout << "Processed Event " << n << endl;    

    Bool_t passSkim = kFALSE;
    mithep::RunLumiRangeMap::RunLumiPairType rl(sigEvent.run_, sigEvent.lumi_);      
    if( sigEvent.dstype_ !=0 || rlrm.HasRunLumi(rl) ) passSkim = kTRUE;

    if (passSkim) {
      normalizedTree->Fill(); 
    }
  }
  cout << "Events in output ntuple: " << normalizedTree->GetEntries() << endl;
  normalizedTree->Write();
  outputFile->Close();
}
