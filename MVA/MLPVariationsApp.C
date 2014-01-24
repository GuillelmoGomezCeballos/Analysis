#include "MitNtupleEvent.C" // Created by TTree::MakeClass()

//------------------------------------------------------------------------------
// Choose MVA methods to be applied
//------------------------------------------------------------------------------
Bool_t Use_BDT           = 1;
Bool_t Use_MLP_0         = 1;
Bool_t Use_MLP_1         = 1;
Bool_t Use_MLP_2         = 1;
Bool_t Use_MLP_3         = 1;
Bool_t Use_MLP_4         = 1;
Bool_t Use_MLP_5         = 1;
Bool_t Use_MLP_6         = 1;
// -----------------------------------------------------------------------------
Bool_t Use_dim01          = 0;
Bool_t Use_pt2         = 0;
Bool_t Use_pt1         = 0;
Bool_t Use_met           = 0;
Bool_t Use_delphil       = 0;
Bool_t Use_mtw2        = 0;
Bool_t Use_mtw1        = 0;
Bool_t Use_eta           = 0;
Bool_t Use_ltype         = 0;
Bool_t Use_njets         = 0;

//------------------------------------------------------------------------------
// Get Tree from file
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname)
{
  cout << "*** Open file " << infname << " for reading" << endl;
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get("all");
  assert(t);
  
  cout << "*** Recovered Tree "<< t->GetName()
       << " with " << t->GetEntries() << " entries" << endl;

  return t;
}

//------------------------------------------------------------------------------
// MLPVariationsApp
//------------------------------------------------------------------------------
void MLPVariationsApp
(
 const char* inputFile,
 const char* directory,
 TString     tag          = "150.0jets",
 TString     myMethodList = "MLP_0,MLP_1,MLP_2,MLP_3,MLP_4,MLP_5,MLP_6",
 TString     myVarList    = "dim01,pt2,pt1,met,delphil,ltype",
 TString     test         = ""
 )
{
  printf("\n==> Start TMVApplication\n");
  printf("\n");
  printf("    inputFile: %s\n",inputFile          );
  printf("    directory: %s\n",directory          );
  printf("          tag: %s\n",tag.Data()         );
  printf(" myMethodList: %s\n",myMethodList.Data());
  printf("    myVarList: %s\n",myVarList.Data()   );
  printf("\n");

  parseOptions   (myMethodList);
  parseVarOptions(myVarList   );

  // Create the Reader object
  TMVA::Reader *reader = new TMVA::Reader("!Color");  

  // Create a set of variables and declare them to the reader.
  // Must be the same as used for training (HiggsAna2.C)
  //  float dim01, pt1, pt2, met,, delphil, ltype;
  float dim01, pt2, pt1, met, delphil, mtw2, mtw1, eta;
  int   ltype, njets;

  if (Use_dim01  ) reader->AddVariable("dim01",   &dim01  );
  if (Use_pt2    ) reader->AddVariable("pt2",     &pt2    );
  if (Use_pt1    ) reader->AddVariable("pt1",     &pt1    );
  if (Use_met    ) reader->AddVariable("met",     &met    );
  if (Use_delphil) reader->AddVariable("delphil", &delphil);
  if (Use_mtw2 ) reader->AddVariable("mtw2",  &mtw2 );
  if (Use_mtw1 ) reader->AddVariable("mtw1",  &mtw1 );
  if (Use_eta    ) reader->AddVariable("eta",     &eta    );
  if (Use_ltype  ) reader->AddVariable("ltype",   &ltype  );
  if (Use_njets  ) reader->AddVariable("njets",   &njets  );

  TString dirarg(directory);

  dir = dirarg + "_BDT.weights.txt";
  method = "BDT method";
  if (Use_BDT) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_0.weights.txt";
  method = "MLP_0 method";
  if (Use_MLP_0) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_1.weights.txt";
  method = "MLP_1 method";
  if (Use_MLP_1) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_2.weights.txt";
  method = "MLP_2 method";
  if (Use_MLP_2) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_3.weights.txt";
  method = "MLP_3 method";
  if (Use_MLP_3) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_4.weights.txt";
  method = "MLP_4 method";
  if (Use_MLP_4) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_5.weights.txt";
  method = "MLP_5 method";
  if (Use_MLP_5) reader->BookMVA(method,dir);

  dir = dirarg + "_MLP_6.weights.txt";
  method = "MLP_6 method";
  if (Use_MLP_6) reader->BookMVA(method,dir);
 
  // Prepare the input tree
  TTree* theTree = getTreeFromFile(inputFile);
  assert(theTree);
  MitNtupleEvent myEvent(theTree);
	
  // Open output file
  TString ofn(inputFile);
  ofn.ReplaceAll("root",test + "-" + tag + ".tmva.root");
  ofn.ReplaceAll("data/","proc/");
	
  cout << "*** Open file " << ofn << " for writing" << endl;

  TFile*   target = new TFile(ofn.Data(),"recreate");
  TNtuple* nt     = new TNtuple("nt","Classified Ntuple",
				"run:event:bdt:mlp_0:mlp_1:mlp_2:mlp_3:mlp_4:mlp_5:mlp_6");
		
  cout << "*** Processing " << theTree->GetEntries() << " events" << endl;
	
  TStopwatch sw;
  sw.Start();
	
  int nTotal = theTree->GetEntries();
  for (Long64_t ievt=0; ievt<theTree->GetEntries(); ievt++) {
		
    if (ievt%1000 == 0)
      cout << "*** --- event " << ievt << " of " << nTotal << endl;
		
    myEvent.GetEntry(ievt);
		
    // Now set variables
    if (Use_dim01   ) dim01    = myEvent.H_dim01;
    if (Use_pt2  ) pt2   = myEvent.H_pt2;
    if (Use_pt1  ) pt1   = myEvent.H_pt1;
    if (Use_met    ) met     = myEvent.H_met;
    if (Use_delphil) delphil = myEvent.H_delphil;
    if (Use_mtw2 ) mtw2  = myEvent.H_mtw2;
    if (Use_mtw1 ) mtw1  = myEvent.H_mtw1;
    if (Use_eta    ) eta     = myEvent.H_eta;
    if (Use_ltype  ) ltype   = myEvent.H_ltype;
    if (Use_njets  ) njets   = myEvent.H_njets;

    double bdt, mlp_0, mlp_1, mlp_2, mlp_3, mlp_4, mlp_5, mlp_6;
    bdt = mlp_0 = mlp_1 = mlp_2 = mlp_3 = mlp_4 = mlp_5 = mlp_6 = 0.0;

    if (Use_BDT)
      bdt = reader->EvaluateMVA("BDT method");
    if (Use_MLP_0)
      mlp_0 = reader->EvaluateMVA("MLP_0 method");
    if (Use_MLP_1)
      mlp_1 = reader->EvaluateMVA("MLP_1 method");
    if (Use_MLP_2)
      mlp_2 = reader->EvaluateMVA("MLP_2 method");
    if (Use_MLP_3)
      mlp_3 = reader->EvaluateMVA("MLP_3 method");
    if (Use_MLP_4)
      mlp_4 = reader->EvaluateMVA("MLP_4 method");
    if (Use_MLP_5)
      mlp_5 = reader->EvaluateMVA("MLP_5 method");
    if (Use_MLP_6)
      mlp_6 = reader->EvaluateMVA("MLP_6 method");
				
    nt->Fill(myEvent.H_run,myEvent.H_event,bdt,mlp_0,mlp_1,mlp_2,mlp_3,mlp_4,mlp_5,mlp_6);
  }
	
  //  nt->Scan("*");	
  //  nt->Draw("fisher");
	
  sw.Stop();
  cout << "*** End of event loop: ";
  sw.Print();

  // Write output file, close
  target->Write();
  target->Close();
	
  delete reader;
    
  cout << "==> End TMVApplication" << endl << endl;
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
