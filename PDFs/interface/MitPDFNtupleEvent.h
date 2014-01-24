//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 13 19:39:28 2010 by ROOT version 5.26/00c
// from TTree PDFTree/PDFTree
// found on file: histo_f10-ww2l-z2-v12_all_noskim.root
//////////////////////////////////////////////////////////

#ifndef MitPDFNtupleEvent_h
#define MitPDFNtupleEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MitPDFNtupleEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         H_weight;
   Float_t         H_lq;
   Float_t         H_lid1;
   Float_t         H_lx1;
   Float_t         H_lpdf1;
   Float_t         H_lid2;
   Float_t         H_lx2;
   Float_t         H_lpdf2;

   // List of branches
   TBranch        *b_H;   //!

   MitPDFNtupleEvent(TTree *tree=0);
   virtual ~MitPDFNtupleEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MitPDFNtupleEvent_cxx
MitPDFNtupleEvent::MitPDFNtupleEvent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histo_f10-ww2l-z2-v12_all_noskim.root");
      if (!f) {
         f = new TFile("histo_f10-ww2l-z2-v12_all_noskim.root");
      }
      tree = (TTree*)gDirectory->Get("PDFTree");

   }
   Init(tree);
}

MitPDFNtupleEvent::~MitPDFNtupleEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MitPDFNtupleEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MitPDFNtupleEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MitPDFNtupleEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("H", &H_weight, &b_H);
   Notify();
}

Bool_t MitPDFNtupleEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MitPDFNtupleEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MitPDFNtupleEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MitPDFNtupleEvent_cxx
