//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 13 16:58:13 2010 by ROOT version 5.22/00d
// from TTree HwwTree/HwwTree
// found on file: histo_f10-ww2l-z2-v12_all_noskim.root
//////////////////////////////////////////////////////////

#ifndef MitNtupleEvent_h
#define MitNtupleEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MitNtupleEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         H_weight;
   Float_t         H_ltype;
   Float_t         H_pt1;
   Float_t         H_pt2;
   Float_t         H_pt3;
   Float_t         H_eta1;
   Float_t         H_eta2;
   Float_t         H_eta3;
   Float_t         H_met;
   Float_t         H_dilpt;
   Float_t         H_dily;
   Float_t         H_delphil;
   Float_t         H_deltarll;
   Float_t         H_deltaphihmet;
   Float_t         H_deltaphilmet;
   Float_t         H_deltaphil3met;
   Float_t         H_metdeltaphilmet;
   Float_t         H_mtw3;
   Float_t         H_mtw1;
   Float_t         H_mtw2;
   Float_t         H_mthiggs;
   Float_t         H_deltaphijetmet;
   Float_t         H_nlep;
   Float_t         H_sumptj;
   Float_t         H_ptj1;
   Float_t         H_ptj2;
   Float_t         H_ptj3;
   Float_t         H_etaj1;
   Float_t         H_etaj2;
   Float_t         H_etaj3;
   Float_t         H_phij1;
   Float_t         H_phij2;
   Float_t         H_phij3;
   Float_t         H_njets;
   Float_t         H_njettiness;
   Float_t         H_nmuons;
   Float_t         H_muonpt;
   Float_t         H_isotot1;
   Float_t         H_isotot2;
   Float_t         H_isotot3;
   Float_t         H_d0sig1;
   Float_t         H_d0sig2;
   Float_t         H_d0sig3;
   Float_t         H_mindr;
   Float_t         H_massz;
   Float_t         H_dim01;
   Float_t         H_dim02;
   Float_t         H_dim12;
   Float_t         H_type1;
   Float_t         H_type2;
   Float_t         H_type3;
   Float_t         H_q1;
   Float_t         H_q2;
   Float_t         H_q3;
   Float_t         H_qqptj1;
   Float_t         H_qqptj2;
   Float_t         H_qqdeta;
   Float_t         H_qqdphi;
   Float_t         H_qqmass;
   Float_t         H_qqeta12;
   Float_t         H_qqnjet;
   Float_t         H_qqgen;
   Float_t         H_qqbtag1;
   Float_t         H_qqbtag2;
   Float_t         H_recotaum;
   Float_t         H_btag1;
   Float_t         H_btag2;
   Float_t         H_maxmbtag;
   Float_t         H_lepgentype1;
   Float_t         H_lepgentype2;
   Float_t         H_lepgentype3;
   Float_t         H_dijetmass;
   Float_t         H_mr;
   Float_t         H_razor;
   Float_t         H_cosomega;
   Float_t         H_ptoth;
   Float_t         H_ptotqqh;
   Float_t         H_deltaphilljet;
   Float_t         H_lq;
   Float_t         H_lid1;
   Float_t         H_lx1;
   Float_t         H_lpdf1;
   Float_t         H_lid2;
   Float_t         H_lx2;
   Float_t         H_lpdf2;
   Float_t         H_zmax;
   Float_t         H_nvert;
   Float_t         H_run;
   Float_t         H_event;
   Float_t         H_decay;

   // List of branches
   TBranch        *b_H;   //!

   MitNtupleEvent(TTree *tree=0);
   virtual ~MitNtupleEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MitNtupleEvent_cxx
MitNtupleEvent::MitNtupleEvent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("histo_f10-ww2l-z2-v12_all_noskim.root");
      if (!f) {
         f = new TFile("histo_f10-ww2l-z2-v12_all_noskim.root");
      }
      tree = (TTree*)gDirectory->Get("HwwTree");

   }
   Init(tree);
}

MitNtupleEvent::~MitNtupleEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MitNtupleEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MitNtupleEvent::LoadTree(Long64_t entry)
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

void MitNtupleEvent::Init(TTree *tree)
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

Bool_t MitNtupleEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MitNtupleEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MitNtupleEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MitNtupleEvent_cxx
