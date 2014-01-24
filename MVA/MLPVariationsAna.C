/**********************************************************************************
* Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
* Package   : TMVA                                                               *
* Root Macro: HiggsAna2, developed from TMVAnalysis                              *
*                                                                                *
* This macro provides examples for the training and testing of all the           *
* TMVA classifiers.                                                              *
*                                                                                *
* As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
* and linearly correlated input variables.                                       *
*                                                                                *
* The methods to be used can be switched on and off by means of booleans, or     *
* via the prompt command, for example:                                           *
*                                                                                *
*    root -l TMVAnalysis.C\(\"Fisher,Likelihood\"\)                              *
*                                                                                *
* (note that the backslashes are mandatory)                                      *
*                                                                                *
* The output file "TMVA.root" can be analysed with the use of dedicated          *
* macros (simply say: root -l <macro.C>), which can be conveniently              *
* invoked through a GUI that will appear at the end of the run of this macro.    *
**********************************************************************************/

#include <iostream> 

#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVAGui.C"

// -----------------------------------------------------------------------------
// Choose MVA methods to be trained + tested
// -----------------------------------------------------------------------------
Bool_t Use_BDT             = 1;
Bool_t Use_MLP_0           = 1;
Bool_t Use_MLP_1           = 1;
Bool_t Use_MLP_2           = 1;
Bool_t Use_MLP_3           = 1;
Bool_t Use_MLP_4           = 1;
Bool_t Use_MLP_5           = 1;
Bool_t Use_MLP_6           = 1;
// -----------------------------------------------------------------------------
Bool_t Use_dim01           = 0;
Bool_t Use_pt2             = 0;
Bool_t Use_pt1             = 0;
Bool_t Use_met             = 0;
Bool_t Use_delphil         = 0;
Bool_t Use_mtw2            = 0;
Bool_t Use_mtw1            = 0;
Bool_t Use_eta             = 0;
Bool_t Use_ltype           = 0;
Bool_t Use_njets           = 0;

// Read input data file with ascii format (otherwise ROOT)
Bool_t ReadDataFromAsciiIFormat = kFALSE;

//------------------------------------------------------------------------------
// MLPVariations
//------------------------------------------------------------------------------
void MLPVariationsAna
(
 int     njets           = 0,			   
 TString signalInputFile = "data/inputNtuple-train-data-standard-histo_H150_3W2l_all.root",
 TString bgdInputFile    = "data/inputNtuple-train-data-standard-histo_HBCK2.root",
 TString outTag          = "mlp-test",
 TString myMethodList    = "BDT,MLP_0,MLP_1,MLP_2,MLP_3,MLP_4,MLP_5,MLP_6",
 TString myVarList       = "dim01,pt2,pt1,met,delphil,ltype"
 ) 
{
  TString mycut("njets==");
  mycut += njets;
	
  printf("\n==> Start TMVAnalysis\n");
  printf("\n");
  printf("           njets: %d\n",njets                 );
  printf("           mycut: %s\n",mycut.Data()          );
  printf(" signalInputFile: %s\n",signalInputFile.Data());
  printf("    bgdInputFile: %s\n",bgdInputFile.Data()   );
  printf("          outTag: %s\n",outTag.Data()         );
  printf("    myMethodList: %s\n",myMethodList.Data()   );
  printf("       myVarList: %s\n",myVarList.Data()      );
  printf("\n");

  parseOptions   (myMethodList);
  parseVarOptions(myVarList   );

  if (njets == 0) {
    outTag += ".0jets";
  }
  else if (njets == 1) {
    outTag += ".1jets";
  }
  else {
    cout << "==> Unknown cut " << mycut.Data() << endl;
    return;
  }

  TString outfileName = outTag + ".root";
  TFile*  outputFile  = TFile::Open(outfileName.Data(),"recreate");
	
  // The 1st argument is the user-defined job name that will reappear in the
  // name of the weight files containing the training results.
  // The 2nd argument is the pointer to a writable TFile target file created
  // by the user, where control and performance histograms are stored.
  TMVA::Factory *factory =
    new TMVA::Factory(outTag.Data(),
		      outputFile,
		      Form("!V:%sColor",gROOT->IsBatch()?"!":""));
		
  // Get inputs (can use more than one tree per sample)
  TTree* signal = getTreeFromFile(signalInputFile.Data());
  assert(signal);
  factory->AddSignalTree(signal,1.0);
	
  // All backgrounds stored in a single tree
  TTree* bgdTree = getTreeFromFile(bgdInputFile.Data());
  assert(bgdTree);
  factory->AddBackgroundTree(bgdTree);

  // The weights for each event are stored in the tree under variable "weight"
  factory->SetWeightExpression("abs(weight)");
	
  // Note that you may also use expressions such as "3*var1/var2*abs(var3)"

  if (Use_dim01   ) factory->AddVariable("dim01",    'F');
  if (Use_pt2  ) factory->AddVariable("pt2",   'F');
  if (Use_pt1  ) factory->AddVariable("pt1",   'F');
  if (Use_met    ) factory->AddVariable("met",     'F');
  if (Use_delphil) factory->AddVariable("delphil", 'F');
  if (Use_mtw2 ) factory->AddVariable("mtw2",  'F');
  if (Use_mtw1 ) factory->AddVariable("mtw1",  'F');
  if (Use_eta    ) factory->AddVariable("eta",     'F');
  if (Use_ltype  ) factory->AddVariable("ltype",   'I');
  if (Use_njets  ) factory->AddVariable("njets",   'I');
	
  factory->PrepareTrainingAndTestTree(TCut(mycut.Data()),
				      "SplitMode=Random:NormMode=NumEvents:!V");

  // Book MVA methods
  //----------------------------------------------------------------------------
  if (Use_BDT)
    factory->BookMethod(TMVA::Types::kBDT,"BDT",
			"!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=4.5");
  
//  if (Use_MLP_0)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_0",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N,N:TestRate=5");
//  if (Use_MLP_1)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_1",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+1,N:TestRate=5");
//  if (Use_MLP_2)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_2",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+2,N:TestRate=5");
//  if (Use_MLP_3)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_3",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+2,N+1,N:TestRate=5");
//  if (Use_MLP_4)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_4",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N:TestRate=5");
//  if (Use_MLP_5)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_5",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+1:TestRate=5");
//  if (Use_MLP_6)
//    factory->BookMethod(TMVA::Types::kMLP,"MLP_6",
//			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+2:TestRate=5");


  if (Use_MLP_0)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_0",
			"Normalise:!H:!V:NCycles=50:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_1)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_1",
			"Normalise:!H:!V:NCycles=100:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_2)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_2",
			"Normalise:!H:!V:NCycles=200:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_3)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_3",
			"Normalise:!H:!V:NCycles=500:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_4)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_4",
			"Normalise:!H:!V:NCycles=700:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_5)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_5",
			"Normalise:!H:!V:NCycles=1000:HiddenLayers=N+1,N:TestRate=5");
  if (Use_MLP_6)
    factory->BookMethod(TMVA::Types::kMLP,"MLP_6",
			"Normalise:!H:!V:NCycles=2000:HiddenLayers=N+1,N:TestRate=5");
	
  // ---- Train MVAs using the set of training events
  factory->TrainAllMethods();
	
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
	
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();    
	
  //----------------------------------------------------------------------------
	
  // Save the output
  outputFile->Close();
	
  std::cout << "==> wrote root file " << outfileName.Data() << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;      
	
  // Clean up
  delete factory;
	
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVAGui(outfileName);
}

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname)
{
  cout << "Open file:\t" << infname << "\tfor reading" << endl;

  TFile* inf = new TFile(infname,"read");
  assert(inf);
  TTree* t = (TTree*)inf->Get("all");
  assert(t);
  cout << "Recovered tree:\t" << t->GetName()
       << "\twith:\t"         << t->GetEntries()
       << "\tentries"         <<endl;

  return t;
}

//------------------------------------------------------------------------------
// parseOptions
//------------------------------------------------------------------------------
void parseOptions(TString myMethodList)
{
  TList* mlist = TMVA::Tools::ParseFormatLine(myMethodList,":,");
	
  if (mlist->GetSize() > 0) {
    Use_BDT
      = Use_MLP_0
      = Use_MLP_1
      = Use_MLP_2
      = Use_MLP_3
      = Use_MLP_4
      = Use_MLP_5
      = Use_MLP_6
      = 0;
		
    if (mlist->FindObject("BDT"  ) != 0) Use_BDT   = 1; 
    if (mlist->FindObject("MLP_0") != 0) Use_MLP_0 = 1; 
    if (mlist->FindObject("MLP_1") != 0) Use_MLP_1 = 1; 
    if (mlist->FindObject("MLP_2") != 0) Use_MLP_2 = 1; 
    if (mlist->FindObject("MLP_3") != 0) Use_MLP_3 = 1; 
    if (mlist->FindObject("MLP_4") != 0) Use_MLP_4 = 1; 
    if (mlist->FindObject("MLP_5") != 0) Use_MLP_5 = 1; 
    if (mlist->FindObject("MLP_6") != 0) Use_MLP_6 = 1; 
    
    delete mlist;
  }
}

//------------------------------------------------------------------------------
// parseVarOptions
//------------------------------------------------------------------------------
void parseVarOptions(TString myMethodList)
{
  TList* mlist = TMVA::Tools::ParseFormatLine(myMethodList,":,");
	
  if (mlist->GetSize() > 0) {
    Use_dim01   
      = Use_pt2  
      = Use_pt1  
      = Use_met    
      = Use_delphil
      = Use_mtw2 
      = Use_mtw1 
      = Use_eta  
      = Use_ltype  
      = Use_njets
      = 0;
		
    if (mlist->FindObject("dim01"   ) != 0) Use_dim01    = 1; 
    if (mlist->FindObject("pt2"  ) != 0) Use_pt2   = 1; 
    if (mlist->FindObject("pt1"  ) != 0) Use_pt1   = 1; 
    if (mlist->FindObject("met"    ) != 0) Use_met     = 1; 
    if (mlist->FindObject("delphil") != 0) Use_delphil = 1; 
    if (mlist->FindObject("mtw2" ) != 0) Use_mtw2  = 1; 
    if (mlist->FindObject("mtw1" ) != 0) Use_mtw1  = 1; 
    if (mlist->FindObject("eta"    ) != 0) Use_eta     = 1; 
    if (mlist->FindObject("ltype"  ) != 0) Use_ltype   = 1; 
    if (mlist->FindObject("njets"  ) != 0) Use_njets   = 1; 
    
    delete mlist;
  }
}
