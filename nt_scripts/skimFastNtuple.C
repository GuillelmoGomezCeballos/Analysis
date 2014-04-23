#include "TFile.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void skimFastNtuple(TString inputFile = "input.root",TString outputFile = "output.root", int skim = 0){

TFile *_file0 = TFile::Open(inputFile);
_file0->cd();
TTree* tree = (TTree*)_file0->Get("tree");
TTree* oldtree = tree;
cout << oldtree->GetEntries() << endl;

TFile *_file1 = TFile::Open(outputFile,"recreate");
_file1->cd();

char cutString[200];
if     (skim == 0) sprintf(cutString,"lep1->Pt()>10.");
else if(skim == 1) sprintf(cutString,"(type==1||type==2)&&lid3==0&&((cuts & 524288) == 524288)");
else if(skim == 2) sprintf(cutString,"min(pmet,pTrackMet)>20&&(TMath::Abs(dilep->mass()-91.1876)>15||type==1||type==2)&&lid3==0&&dilep->pt()>45&&njets<=3");
else if(skim == 3) sprintf(cutString,"((cuts & 4) == 4) && ((cuts & 512) == 512)");
else if(skim == 4) sprintf(cutString,"((cuts & 4) != 4) || ((cuts & 512) != 512)");
else if(skim == 5) sprintf(cutString,"lid3!=0");
else if(skim == 6) sprintf(cutString,"min(pmet,pTrackMet)>20&&((TMath::Abs(dilep->mass()-91.1876)>15&&dilep->pt()>45)||type==1||type==2)&&lid3==0&&dilep->pt()>30&&njets<=3");
else if(skim == 7) sprintf(cutString,"((cuts & 4) != 4) && ((cuts & 512) != 512)");
else if(skim == 8) sprintf(cutString,"lq1*lq2>0&&lid3==0");
else if(skim == 9) sprintf(cutString,"min(pmet,pTrackMet)>20&&((TMath::Abs(dilep->mass()-91.1876)>15&&((dymva>0.88&&njets==0)||(dymva>0.84&&njets==1)||(met>45&&njets>=2)))||type==1||type==2)&&lid3==0&&dilep->M()>12&&dilep->pt()>0");
else if(skim ==10) sprintf(cutString,"dilep->pt()>60&&met>60&&dilep->mass()>55&&dilep->mass()<200&&njets<=1");
else if(skim ==11) sprintf(cutString,"dilep->pt()>60&&met>60&&njets<=1");
else if(skim ==12) sprintf(cutString,"min(pmet,pTrackMet)>20&&(TMath::Abs(dilep->mass()-91.1876)>15||type==1||type==2)&&lid3==0&&lq1*lq2<0&&dilep->mass()>12&&njets>=2&&((jet1ProbBtag>=0.244&&jet2ProbBtag>=0.244)||jet1ProbBtag>=0.679||jet2ProbBtag>=0.679)");
else if(skim ==13) sprintf(cutString,"dilep->pt()>60");
else if(skim ==14) sprintf(cutString,"dilep->mass()>15&&TMath::Abs(jet1->Eta()-jet2->Eta())>2.5&&jet1->Pt()>25&&jet2->Pt()>25&&met>40");
else {cout << "No good option" << endl; return;}
cout << "cut: " << cutString << endl;

TTree *newtree1;
newtree1 = oldtree->CopyTree(cutString);
cout << newtree1->GetEntries() << endl;

_file1->Write();
_file1->Close();
}
