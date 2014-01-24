#include "TFile.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TRandom.h"

TFile *_file0;
TFile *_file1;
TTree* oldtree;
TTree *newtree;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

// "jetLowBtag>2.1||nSoftMuons>0||jet1Btag>2.1||jet2Btag>2.1"
// "auxVar0>0.65||nSoftMuons>0||jet1ProbBtag>0.65||jet2ProbBtag>0.65"
void topFakeAna(TString inputFile = "/home/ceballos/condor/old_42x_pu/histo_f11-wwj-v14b-pu_all_noskim.root",
                int nVtxMin = 1, int nVtxMax = 20, bool isZPeak = false, int bFilterType = 0,
		char bTaggingString[50] = "jetLowBtag>2.1||nSoftMuons>0"){

const int NjetsMax = 2;

_file0 = TFile::Open(inputFile);
if(!_file0) return;
_file0->cd();
TTree* tree = (TTree*)_file0->Get("HwwTree0");
oldtree = tree;

_file1 = TFile::Open(Form("trash_%d.root",gRandom->Integer(1000)),"recreate");
_file1->cd();

char bFilterTypeString[200];
sprintf(bFilterTypeString," ");
if     (bFilterType == 1) sprintf(bFilterTypeString ,"&&(jet1McId==5||jet2McId==5||jet3McId==5||jet4McId==5)");
else if(bFilterType == 2) sprintf(bFilterTypeString,"&&!(jet1McId==5||jet2McId==5||jet3McId==5||jet4McId==5)");

if(isZPeak == true){

  char cutZSelString[200];
  sprintf(cutZSelString,"TMath::Abs(dilep->mass()-91.1876)<15&&lid3==0&&(type==0||type==3)&&met<30&&lep1->pt()>20&&lep2->pt()>20&&nvtx>=%d&&nvtx<=%d%s",nVtxMin,nVtxMax,bFilterTypeString);
  cout << "cutZSelString: " << cutZSelString << endl;

  double topTagZ[NjetsMax+1][6],topTagZIni[NjetsMax+1];
  for(int i=0; i<=NjetsMax; i++){
    newtree = oldtree->CopyTree(Form("njets==%d&&%s",i,cutZSelString));
    topTagZIni[i] = newtree->GetEntries();
    if(topTagZIni[i] == 0) topTagZIni[i] = 1;

    newtree = oldtree->CopyTree(Form("njets==%d&&%s&&(%s)",i,cutZSelString,bTaggingString));
    topTagZ[i][0] = newtree->GetEntries()/topTagZIni[i];
    if(topTagZ[i][0] == 0) topTagZ[i][0] = 1/topTagZIni[i];

    printf(" ZSel(%d) & %7d & %6.4f $\\pm$ %6.4f\\\\\n",
           i,(int)topTagZIni[i],topTagZ[i][0],sqrt(topTagZ[i][0]*(1.0-topTagZ[i][0])/topTagZIni[i]));
  }

} else {

  char cutWWSelString[200];
  sprintf(cutWWSelString,"min(pmet,pTrackMet)>20&&(TMath::Abs(dilep->mass()-91.1876)>15||type==1||type==2)&&lid3==0&&lep1->pt()>20&&lep2->pt()>10&&nvtx>=%d&&nvtx<=%d%s",nVtxMin,nVtxMax,bFilterTypeString);
  cout << "cutWWSelString: " << cutWWSelString << endl;

  double topTagWW[NjetsMax+1][6],topTagWWIni[NjetsMax+1];
  for(int i=0; i<=NjetsMax; i++){
    newtree = oldtree->CopyTree(Form("njets==%d&&%s",i,cutWWSelString));
    topTagWWIni[i] = newtree->GetEntries();
    if(topTagWWIni[i] == 0) topTagWWIni[i] = 1;

    newtree = oldtree->CopyTree(Form("njets==%d&&%s&&(%s)",i,cutWWSelString,bTaggingString));
    topTagWW[i][0] = newtree->GetEntries()/topTagWWIni[i];
    if(topTagWW[i][0] == 0) topTagWW[i][0] = 1/topTagWWIni[i];

    printf("WWSel(%d) & %7d & %6.4f $\\pm$ %6.4f\\\\\n",
           i,(int)topTagWWIni[i],topTagWW[i][0],sqrt(topTagWW[i][0]*(1.0-topTagWW[i][0])/topTagWWIni[i]));
  }

}
/*
if(isZPeak == true){

  char cutZSelStringMuon[200];
  sprintf(cutZSelStringMuon,"TMath::Abs(dilep->mass()-91.1876)<15&&lid3==0&&(type==0||type==3)&&met<30&&lep1->pt()>20&&lep2->pt()>20&&nvtx>=%d&&nvtx<=%d%s",nVtxMin,nVtxMax,bFilterTypeString);
  cout << "cutZSelStringMuon: " << cutZSelStringMuon << endl;

  double topTagZMuon[NjetsMax+1][6],topTagZMuonIni[NjetsMax+1];
  for(int i=0; i<=NjetsMax; i++){
    newtree = oldtree->CopyTree(Form("njets==%d&&%s",i,cutZSelStringMuon));
    topTagZMuonIni[i] = newtree->GetEntries();
    if(topTagZMuonIni[i] == 0) topTagZMuonIni[i] = 1;

    newtree = oldtree->CopyTree(Form("nSoftMuons!=0&&njets==%d&&%s",i,cutZSelStringMuon));
    topTagZMuon[i][0] = newtree->GetEntries()/topTagZMuonIni[i];
    if(topTagZMuon[i][0] == 0) topTagZMuon[i][0] = 1/topTagZMuonIni[i];

    printf(" ZMuonSel(%d) & %7d & %6.4f $\\pm$ %6.4f \\\\\n",
           i,(int)topTagZMuonIni[i],topTagZMuon[i][0],sqrt(topTagZMuon[i][0]*(1.0-topTagZMuon[i][0])/topTagZMuonIni[i]));
  }

}
else {

  char cutWWSelStringMuon[200];
  sprintf(cutWWSelStringMuon,"min(pmet,pTrackMet)>20&&(TMath::Abs(dilep->mass()-91.1876)>15||type==1||type==2)&&lid3==0&&lep1->pt()>20&&lep2->pt()>10&&nvtx>=%d&&nvtx<=%d%s",nVtxMin,nVtxMax,bFilterTypeString);
  cout << "cutWWSelStringMuon: " << cutWWSelStringMuon << endl;

  double topTagWWMuon[NjetsMax+1][6],topTagWWMuonIni[NjetsMax+1];
  for(int i=0; i<=NjetsMax; i++){
    newtree = oldtree->CopyTree(Form("njets==%d&&%s",i,cutWWSelStringMuon));
    topTagWWMuonIni[i] = newtree->GetEntries();
    if(topTagWWMuonIni[i] == 0) topTagWWMuonIni[i] = 1;

    newtree = oldtree->CopyTree(Form("nSoftMuons!=0&&njets==%d&&%s",i,cutWWSelStringMuon));
    topTagWWMuon[i][0] = newtree->GetEntries()/topTagWWMuonIni[i];
    if(topTagWWMuon[i][0] == 0) topTagWWMuon[i][0] = 1/topTagWWMuonIni[i];

    printf(" WWMuonSel(%d) & %7d & %6.4f $\\pm$ %6.4f\\\\\n",
           i,(int)topTagWWMuonIni[i],topTagWWMuon[i][0],sqrt(topTagWWMuon[i][0]*(1.0-topTagWWMuon[i][0])/topTagWWMuonIni[i]));
  }

}
*/

//_file0->Close();
//_file1->Close();
//delete oldtree;
//delete newtree;

}
