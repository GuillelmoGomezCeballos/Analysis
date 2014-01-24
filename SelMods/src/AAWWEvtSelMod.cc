// $Id: AAWWEvtSelMod.cc,v 1.4 2013/01/09 10:20:13 ceballos Exp $

#include "Analysis/SelMods/interface/AAWWEvtSelMod.h"

using namespace mithep;
ClassImp(mithep::AAWWEvtSelMod)

//--------------------------------------------------------------------------------------------------
AAWWEvtSelMod::AAWWEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fIsData(kFALSE),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFMetName("PFMet"),
  fPFMetStd(0),
  fPileupInfoName(Names::gkPileupInfoBrn),
  fPileupInfo(0),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void AAWWEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void AAWWEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  ElectronOArr *CleanElectrons  = GetObjThisEvt<ElectronOArr>("CleanElectronsBS");
  MuonOArr  *CleanMuons         = GetObjThisEvt<MuonOArr>("CleanMuonsBS");
  ParticleOArr *leptons         = GetObjThisEvt<ParticleOArr>("MergedLeptonsBS");
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  if(leptons->GetEntries() <  2) return;  
  if(leptons->At(0)->Pt() <= 20) return;
  if(leptons->At(1)->Pt() <= 10) return;

  CompositeParticle dilepton;;
  dilepton.AddDaughter(leptons->At(0));
  dilepton.AddDaughter(leptons->At(1));

  LoadBranch(fPFCandidatesName);

  LoadBranch(fPFMetName);
  const PFMet *PFMetStd = fPFMetStd->At(0);
  LoadBranch(fEvtHdrName);
  //if(fIsData == kFALSE){
  //  LoadBranch(fPileupInfoName);
  //}

  double zAverage = 0.0;
  double zDiffMax = 0.0;
  std::vector<double> leptonsDz;

  // Make lepton vector from muons and electrons
  for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
    double pDz = CleanMuons->At(j)->BestTrk()->DzCorrected(*fVertices->At(0));
    leptonsDz.push_back(pDz);
  }
  for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
    double pDz = CleanElectrons->At(j)->GsfTrk()->DzCorrected(*fVertices->At(0));
    leptonsDz.push_back(pDz);
  }
  for(UInt_t t=0; t<leptonsDz.size(); t++) {
    zAverage = zAverage + leptonsDz[t];
    for(UInt_t i=t+1; i<leptonsDz.size(); i++) {
      if(TMath::Abs(leptonsDz[t]-leptonsDz[i]) > zDiffMax) zDiffMax = TMath::Abs(leptonsDz[t]-leptonsDz[i]);
    }
  }
  leptonsDz.clear();
  zAverage = zAverage / leptons->GetEntries();
  hDaawwPresel[0]->Fill(TMath::Min(zDiffMax,0.999),NNLOWeight->GetVal());
  if(zDiffMax >= 0.1) return;
  hDaawwPresel[1]->Fill(TMath::Min(TMath::Max(TMath::Abs(leptonsDz[0]),TMath::Abs(leptonsDz[1])),0.999),NNLOWeight->GetVal());

  double trackNumeratorX  =0.0, trackNumeratorY  = 0.0;
  int nTracks[3] = {0, 0, 0};
  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pf = fPFCandidates->At(i);
    // charged
    if (pf->HasTrackerTrk() || pf->HasGsfTrk()){
      bool isMuonTrack = false;
      for (UInt_t m = 0; m < CleanMuons->GetEntries(); ++m) {
	if (CleanMuons->At(m)->TrackerTrk() == pf->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
      }      
      if (isMuonTrack) continue;
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < CleanElectrons->GetEntries(); ++m) {
	if ( (CleanElectrons->At(m)->TrackerTrk() == pf->TrackerTrk()) or
	     (CleanElectrons->At(m)->HasGsfTrk() and CleanElectrons->At(m)->GsfTrk() == pf->GsfTrk()) ) {
	  isElectronTrack = true;
	  break;
	}
      }      
      if (isElectronTrack) continue;

      if((pf->HasTrackerTrk() && fVertices->At(0)->HasTrack(pf->TrackerTrk()) && fVertices->At(0)->TrackWeight(pf->TrackerTrk()) > 0) ||
         (pf->HasGsfTrk()     && fVertices->At(0)->HasTrack(pf->GsfTrk())     && fVertices->At(0)->TrackWeight(pf->GsfTrk())     > 0 )) nTracks[0]++;

      if ((pf->HasTrackerTrk() && fabs(pf->TrackerTrk()->DzCorrected(*fVertices->At(0))) < 0.1) ||
          (pf->HasGsfTrk()     && fabs(pf->GsfTrk()->DzCorrected(*fVertices->At(0))    ) < 0.1)) {
	trackNumeratorX =+ pf->Px();
	trackNumeratorY =+ pf->Py();
	nTracks[1]++;
	if(pf->Pt() > 0.5) nTracks[2]++;
      }
    }
  }

  int pairType = -1;
  if (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon )
    pairType = 0;
  else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron)
    pairType = 1;
  else if((leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon) || 
  	  (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kElectron))
    pairType = 2;
  else {
    cout << "Hey, this is not possible, leptonTypes: "
  	 << leptons->At(0)->ObjType() << " - " 
  	 << leptons->At(1)->ObjType() << endl;
  }

  double sumPtTracks = sqrt(trackNumeratorX*trackNumeratorX+trackNumeratorY*trackNumeratorY);
  hDaawwPresel[2]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
  hDaawwPresel[3]->Fill(TMath::Min((double)nTracks[0],99.499),NNLOWeight->GetVal());
  hDaawwPresel[4]->Fill(TMath::Min((double)nTracks[1],99.499),NNLOWeight->GetVal());
  hDaawwPresel[5]->Fill(TMath::Min((double)nTracks[2],99.499),NNLOWeight->GetVal());
  if(nTracks[1] == 0 || sumPtTracks <= 0){
    hDaawwPresel[6]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    hDaawwPresel[7]->Fill(TMath::Min(dilepton.Pt(),199.999),NNLOWeight->GetVal());
    hDaawwPresel[8]->Fill(TMath::Min(PFMetStd->Pt(),199.999),NNLOWeight->GetVal());
    if(dilepton.Charge() == 0)
      hDaawwPresel[9]->Fill(pairType+0,NNLOWeight->GetVal());
    else
      hDaawwPresel[9]->Fill(pairType+3,NNLOWeight->GetVal());
    
    if(fIsData == kTRUE) {
      printf("DATAAAWW: %d %d %d %f\n",fEventHeader->EvtNum(),fEventHeader->RunNum(),fEventHeader->LumiSec(),dilepton.Pt());
    }

    if(dilepton.Pt() > 20.0 && dilepton.Charge() == 0){
      hDaawwPresel[10]->Fill(TMath::Min(dilepton.Mass(),199.999),NNLOWeight->GetVal());
    }
  }
}
//--------------------------------------------------------------------------------------------------
void AAWWEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fEvtHdrName,       fEventHeader);
  ReqBranch(fPFMetName,        fPFMetStd);
  ReqBranch(fPFCandidatesName, fPFCandidates);

  //if(fIsData == kFALSE){
  //  ReqBranch(fPileupInfoName,fPileupInfo);
  //}

  char sb[200];
  sprintf(sb,"hDaawwPresel_%d",0);   hDaawwPresel[ 0]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDaawwPresel_%d",1);   hDaawwPresel[ 1]  = new TH1D(sb,sb,100,0.0,1.0);
  sprintf(sb,"hDaawwPresel_%d",2);   hDaawwPresel[ 2]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDaawwPresel_%d",3);   hDaawwPresel[ 3]  = new TH1D(sb,sb,100,-0.5,99.5);
  sprintf(sb,"hDaawwPresel_%d",4);   hDaawwPresel[ 4]  = new TH1D(sb,sb,100,-0.5,99.5); 
  sprintf(sb,"hDaawwPresel_%d",5);   hDaawwPresel[ 5]  = new TH1D(sb,sb,100,-0.5,99.5); 
  sprintf(sb,"hDaawwPresel_%d",6);   hDaawwPresel[ 6]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDaawwPresel_%d",7);   hDaawwPresel[ 7]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDaawwPresel_%d",8);   hDaawwPresel[ 8]  = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDaawwPresel_%d",9);   hDaawwPresel[ 9]  = new TH1D(sb,sb,6,-0.5,5.5); 
  sprintf(sb,"hDaawwPresel_%d",10);  hDaawwPresel[10]  = new TH1D(sb,sb,200,0.0,200.0);

  for(int i=0; i<=10; i++){
    AddOutput(hDaawwPresel[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void AAWWEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void AAWWEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
