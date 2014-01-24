// $Id: WtnEvtSelMod.cc,v 1.5 2011/03/07 12:48:04 ceballos Exp $

#include "Ana/SelMods/interface/WtnEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

using namespace mithep;
ClassImp(mithep::WtnEvtSelMod)

//--------------------------------------------------------------------------------------------------
WtnEvtSelMod::WtnEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMetName(Names::gkCaloMetBrn),
  fPFTaus0Name("random"),
  fPFTaus0(0),
  fPFTaus1Name("random"),
  fPFTaus1(0),
  fCleanPFTausName(ModNames::gkCleanPFTausName),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WtnEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WtnEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "WtnEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  //Get Generator Level information for matching
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");
  MCParticleOArr *GenTaus      = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCTausName);
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *CleanMuons        = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  //JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  //MetOArr *CleanMet	       = GetObjThisEvt<MetOArr>(fMetName);
  //const Met *caloMet	       = CleanMet->At(0);
  const PFTauCol   *CleanPFTaus   = 0;
  if (!fCleanPFTausName.IsNull()) CleanPFTaus = GetObjThisEvt<PFTauCol>(fCleanPFTausName);

  // fPFTaus0 loop == ShrinkingConePFTaus
  LoadBranch(fPFTaus0Name);
  for (UInt_t i=0; i<fPFTaus0->GetEntries(); ++i) {
    const PFTau *tau = fPFTaus0->At(i);

    bool isLepton = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); ++j) {
      const Muon *mu = CleanMuons->At(j);
      if(MathUtils::DeltaR(mu->Phi(), mu->Eta(), tau->Phi(), tau->Eta()) < 0.3){
        isLepton = kTRUE;
	break;
      }
    }
    if(isLepton == kTRUE) continue;

    for (UInt_t j=0; j<CleanElectrons->GetEntries(); ++j) {
      const Electron *el = CleanElectrons->At(j);
      if(MathUtils::DeltaR(el->Phi(), el->Eta(), tau->Phi(), tau->Eta()) < 0.3){
    	isLepton = kTRUE;
        break;
      }
    }
    if(isLepton == kTRUE) continue;

    Int_t indGenTau = -1;
    for (UInt_t j=0; j<GenTaus->GetEntries(); ++j) {
      MCParticle *gen = GenTaus->At(j);
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), tau->Phi(), tau->Eta()) < 0.3){
        indGenTau = j;
	break;
      }
    }

    Int_t theHisto = 200;
    if(indGenTau == -1) theHisto = 300;

    if(indGenTau != -1 && TMath::Abs(GenTaus->At(indGenTau)->Eta()) < 2.5 &&
       GenTaus->At(indGenTau)->Pt() > 10.0){
      hDWtnSel[151]->Fill(TMath::Min(GenTaus->At(indGenTau)->Pt(),199.999),NNLOWeight->GetVal());
      hDWtnSel[152]->Fill(TMath::Abs(GenTaus->At(indGenTau)->Eta()),NNLOWeight->GetVal());
    }

    hDWtnSel[ 0+theHisto]->Fill(TMath::Min((double)tau->NSignalPFCands(),9.499),NNLOWeight->GetVal());
    if(tau->NSignalPFCands() == 0) continue;

    Int_t nTrk = 0;
    CompositeParticle tauSystem;
    CompositeParticle tauChargedSystem;
    for (UInt_t j=0; j<tau->NSignalPFCands(); ++j) {
      tauSystem.AddDaughter(tau->SignalPFCand(j));
      if(tau->SignalPFCand(j) != 0 &&
         tau->SignalPFCand(j)->Charge() != 0){
        nTrk++;
        tauChargedSystem.AddDaughter(tau->SignalPFCand(j));
      }
    }

    Bool_t cuts[6] = {tauSystem.Pt() > 15, (nTrk == 1 || nTrk == 3) && TMath::Abs(tau->Charge()) == 1, tau->IsoChargedHadronPtSum() < 2.0,
                      tau->IsoGammaEtSum() < 3.0, tauChargedSystem.Mass() > 0.13 && tauChargedSystem.Mass() < 2.00, 
		      tau->LeadChargedHadronPFCand() && tau->LeadChargedHadronPFCand()->Pt() > 5.0};
    
    if(           cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5]){
      hDWtnSel[ 1+theHisto]->Fill(TMath::Min(tauSystem.Pt(),99.999),NNLOWeight->GetVal());
    }

    if(cuts[0] &&            cuts[2] && cuts[3] && cuts[4] && cuts[5]){
      hDWtnSel[ 2+theHisto]->Fill(TMath::Min((double)nTrk,7.499),NNLOWeight->GetVal());
      hDWtnSel[ 3+theHisto]->Fill(TMath::Min(TMath::Abs(tau->Charge()),7.499),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] &&            cuts[3] && cuts[4] && cuts[5]){
      hDWtnSel[ 4+theHisto]->Fill(TMath::Min(tau->IsoChargedHadronPtSum(),19.999),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] &&            cuts[4] && cuts[5]){
      hDWtnSel[ 5+theHisto]->Fill(TMath::Min(tau->IsoGammaEtSum(),99.999),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] &&            cuts[5]){
      hDWtnSel[ 6+theHisto]->Fill(TMath::Min(tauChargedSystem.Mass(),3.999),NNLOWeight->GetVal());
      hDWtnSel[29+theHisto]->Fill(           tauChargedSystem.Mass(),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4]           ){
      if(tau->LeadChargedHadronPFCand()) hDWtnSel[26+theHisto]->Fill(TMath::Min(tau->LeadChargedHadronPFCand()->Pt(),99.999),NNLOWeight->GetVal());
      else                               hDWtnSel[26+theHisto]->Fill(0.0,NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5]){
      hDWtnSel[ 7+theHisto]->Fill(TMath::Min(tau->Pt(),199.999),NNLOWeight->GetVal());
      hDWtnSel[ 8+theHisto]->Fill(TMath::Min(tau->Mass(),3.999),NNLOWeight->GetVal());
      hDWtnSel[ 9+theHisto]->Fill(TMath::Min(tauChargedSystem.Pt(),99.999),NNLOWeight->GetVal());
      hDWtnSel[10+theHisto]->Fill(TMath::Min(tauSystem.Mass(),3.999),NNLOWeight->GetVal());

      if(!TMath::IsNaN(tau->ElectronPreIDDecision())) hDWtnSel[11+theHisto]->Fill(tau->ElectronPreIDDecision(),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->BremRecoveryEOverP()))    hDWtnSel[12+theHisto]->Fill(TMath::Min(tau->BremRecoveryEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->SegmentCompatibility()))  hDWtnSel[13+theHisto]->Fill(TMath::Min(tau->SegmentCompatibility(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->CaloCompatibility()))     hDWtnSel[14+theHisto]->Fill(TMath::Min(tau->CaloCompatibility(),99.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->MuonDecision()))	      hDWtnSel[15+theHisto]->Fill(tau->MuonDecision(),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->ElectronPreIDOutput()))   hDWtnSel[16+theHisto]->Fill(TMath::Min(tau->ElectronPreIDOutput(),0.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->ECalStripSumEOverP()))    hDWtnSel[17+theHisto]->Fill(TMath::Min(tau->ECalStripSumEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->EMFraction()))	      hDWtnSel[18+theHisto]->Fill(TMath::Min(tau->EMFraction(),0.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCal3x3EOverP()))	      hDWtnSel[19+theHisto]->Fill(TMath::Min(tau->HCal3x3EOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCalMaxEOverP()))	      hDWtnSel[20+theHisto]->Fill(TMath::Min(tau->HCalMaxEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCalTotalEOverP()))       hDWtnSel[21+theHisto]->Fill(TMath::Min(tau->HCalTotalEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->LeadPFCandSignD0Sig()))   hDWtnSel[22+theHisto]->Fill(TMath::Min(tau->LeadPFCandSignD0Sig(),9.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->MaxHCalPFClusterEt()))    hDWtnSel[23+theHisto]->Fill(TMath::Min(tau->MaxHCalPFClusterEt(),99.999),NNLOWeight->GetVal());
      if(tau->SourcePFJet())			      hDWtnSel[24+theHisto]->Fill(TMath::Min(tau->SourcePFJet()->Pt(),99.999),NNLOWeight->GetVal());
      if(tau->LeadPFCand())			      hDWtnSel[25+theHisto]->Fill(TMath::Min(tau->LeadPFCand()->Pt(),99.999),NNLOWeight->GetVal());
      hDWtnSel[27+theHisto]->Fill(TMath::Min((double)tau->NIsoPFCandS(),9.499),NNLOWeight->GetVal());
      hDWtnSel[28+theHisto]->Fill(TMath::Abs(tau->Eta()),NNLOWeight->GetVal());

      if(indGenTau != -1 && TMath::Abs(GenTaus->At(indGenTau)->Eta()) < 2.5 &&
         GenTaus->At(indGenTau)->Pt() > 10.0){
        hDWtnSel[153]->Fill(TMath::Min(GenTaus->At(indGenTau)->Pt(),199.999),NNLOWeight->GetVal());
        hDWtnSel[154]->Fill(TMath::Abs(GenTaus->At(indGenTau)->Eta()),NNLOWeight->GetVal());
      }
    }
  } // fPFTaus0 loop

  // fPFTaus1 loop == HPSTaus
  LoadBranch(fPFTaus1Name);
  for (UInt_t i=0; i<fPFTaus1->GetEntries(); ++i) {
    const PFTau *tau = fPFTaus1->At(i);

    bool isLepton = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); ++j) {
      const Muon *mu = CleanMuons->At(j);
      if(MathUtils::DeltaR(mu->Phi(), mu->Eta(), tau->Phi(), tau->Eta()) < 0.3){
        isLepton = kTRUE;
	break;
      }
    }
    if(isLepton == kTRUE) continue;

    for (UInt_t j=0; j<CleanElectrons->GetEntries(); ++j) {
      const Electron *el = CleanElectrons->At(j);
      if(MathUtils::DeltaR(el->Phi(), el->Eta(), tau->Phi(), tau->Eta()) < 0.3){
    	isLepton = kTRUE;
        break;
      }
    }
    if(isLepton == kTRUE) continue;

    Int_t indGenTau = -1;
    for (UInt_t j=0; j<GenTaus->GetEntries(); ++j) {
      MCParticle *gen = GenTaus->At(j);
      if(MathUtils::DeltaR(gen->Phi(), gen->Eta(), tau->Phi(), tau->Eta()) < 0.3){
        indGenTau = j;
	break;
      }
    }

    Int_t theHisto = 0;
    if(indGenTau == -1) theHisto = 100;

    if(indGenTau != -1 && TMath::Abs(GenTaus->At(indGenTau)->Eta()) < 2.5 &&
       GenTaus->At(indGenTau)->Pt() > 10.0){
      hDWtnSel[51]->Fill(TMath::Min(GenTaus->At(indGenTau)->Pt(),199.999),NNLOWeight->GetVal());
      hDWtnSel[52]->Fill(TMath::Abs(GenTaus->At(indGenTau)->Eta()),NNLOWeight->GetVal());
    }

    hDWtnSel[ 0+theHisto]->Fill(TMath::Min((double)tau->NSignalPFCands(),9.499),NNLOWeight->GetVal());
    if(tau->NSignalPFCands() == 0) continue;

    Int_t nTrk = 0;
    CompositeParticle tauSystem;
    CompositeParticle tauChargedSystem;
    for (UInt_t j=0; j<tau->NSignalPFCands(); ++j) {
      tauSystem.AddDaughter(tau->SignalPFCand(j));
      if(tau->SignalPFCand(j) != 0 &&
         tau->SignalPFCand(j)->Charge() != 0){
        nTrk++;
        tauChargedSystem.AddDaughter(tau->SignalPFCand(j));
      }
    }

    Bool_t cuts[7] = {tauSystem.Pt() > 15.0,
                      tau->DiscriminationAgainstElectron() == 1,
		      tau->DiscriminationAgainstMuon() == 1,
                      tau->DiscriminationByDecayModeFinding() == 1,
		      tau->DiscriminationByMediumIsolation() == 1,
		      nTrk == 1 || nTrk == 3,
		      TMath::Abs(tau->Charge()) == 1};
         
    if( 	  cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5] && cuts[6]){
      hDWtnSel[ 1+theHisto]->Fill(TMath::Min(tauSystem.Pt(),99.999),NNLOWeight->GetVal());
    }

    if(cuts[0]            && cuts[2] && cuts[3] && cuts[4] && cuts[5] && cuts[6]){
      hDWtnSel[30+theHisto]->Fill(TMath::Min(tau->DiscriminationAgainstElectron(),0.999),NNLOWeight->GetVal());
    }
    if(cuts[0] && cuts[1]            && cuts[3] && cuts[4] && cuts[5] && cuts[6]){
      hDWtnSel[31+theHisto]->Fill(TMath::Min(tau->DiscriminationAgainstMuon(),0.999),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2]            && cuts[4] && cuts[5] && cuts[6]){
      hDWtnSel[32+theHisto]->Fill(TMath::Min(tau->DiscriminationByDecayModeFinding(),0.999),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3]            && cuts[5] && cuts[6]){
      hDWtnSel[33+theHisto]->Fill(TMath::Min(tau->DiscriminationByLooseIsolation(),0.999),NNLOWeight->GetVal());
      hDWtnSel[34+theHisto]->Fill(TMath::Min(tau->DiscriminationByMediumIsolation(),0.999),NNLOWeight->GetVal());
      hDWtnSel[35+theHisto]->Fill(TMath::Min(tau->DiscriminationByTightIsolation(),0.999),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4]            && cuts[6]){
      hDWtnSel[ 2+theHisto]->Fill(TMath::Min((double)nTrk,7.499),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5]           ){
      hDWtnSel[ 3+theHisto]->Fill(TMath::Min(TMath::Abs(tau->Charge()),7.499),NNLOWeight->GetVal());
    }

    if(cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5] && cuts[6]){
      hDWtnSel[ 4+theHisto]->Fill(TMath::Min(tau->IsoChargedHadronPtSum(),19.999),NNLOWeight->GetVal());
      hDWtnSel[ 5+theHisto]->Fill(TMath::Min(tau->IsoGammaEtSum(),99.999),NNLOWeight->GetVal());
      hDWtnSel[ 6+theHisto]->Fill(TMath::Min(tauChargedSystem.Mass(),3.999),NNLOWeight->GetVal());
      hDWtnSel[29+theHisto]->Fill(           tauChargedSystem.Mass(),NNLOWeight->GetVal());
      if(tau->LeadChargedHadronPFCand()) hDWtnSel[26+theHisto]->Fill(TMath::Min(tau->LeadChargedHadronPFCand()->Pt(),99.999),NNLOWeight->GetVal());
      else                               hDWtnSel[26+theHisto]->Fill(0.0,NNLOWeight->GetVal());
      hDWtnSel[ 7+theHisto]->Fill(TMath::Min(tau->Pt(),199.999),NNLOWeight->GetVal());
      hDWtnSel[ 8+theHisto]->Fill(TMath::Min(tau->Mass(),3.999),NNLOWeight->GetVal());
      hDWtnSel[ 9+theHisto]->Fill(TMath::Min(tauChargedSystem.Pt(),99.999),NNLOWeight->GetVal());
      hDWtnSel[10+theHisto]->Fill(TMath::Min(tauSystem.Mass(),3.999),NNLOWeight->GetVal());

      if(!TMath::IsNaN(tau->ElectronPreIDDecision())) hDWtnSel[11+theHisto]->Fill(tau->ElectronPreIDDecision(),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->BremRecoveryEOverP()))    hDWtnSel[12+theHisto]->Fill(TMath::Min(tau->BremRecoveryEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->SegmentCompatibility()))  hDWtnSel[13+theHisto]->Fill(TMath::Min(tau->SegmentCompatibility(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->CaloCompatibility()))     hDWtnSel[14+theHisto]->Fill(TMath::Min(tau->CaloCompatibility(),99.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->MuonDecision()))	      hDWtnSel[15+theHisto]->Fill(tau->MuonDecision(),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->ElectronPreIDOutput()))   hDWtnSel[16+theHisto]->Fill(TMath::Min(tau->ElectronPreIDOutput(),0.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->ECalStripSumEOverP()))    hDWtnSel[17+theHisto]->Fill(TMath::Min(tau->ECalStripSumEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->EMFraction()))	      hDWtnSel[18+theHisto]->Fill(TMath::Min(tau->EMFraction(),0.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCal3x3EOverP()))	      hDWtnSel[19+theHisto]->Fill(TMath::Min(tau->HCal3x3EOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCalMaxEOverP()))	      hDWtnSel[20+theHisto]->Fill(TMath::Min(tau->HCalMaxEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->HCalTotalEOverP()))       hDWtnSel[21+theHisto]->Fill(TMath::Min(tau->HCalTotalEOverP(),4.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->LeadPFCandSignD0Sig()))   hDWtnSel[22+theHisto]->Fill(TMath::Min(tau->LeadPFCandSignD0Sig(),9.999),NNLOWeight->GetVal());
      if(!TMath::IsNaN(tau->MaxHCalPFClusterEt()))    hDWtnSel[23+theHisto]->Fill(TMath::Min(tau->MaxHCalPFClusterEt(),99.999),NNLOWeight->GetVal());
      if(tau->SourcePFJet())			      hDWtnSel[24+theHisto]->Fill(TMath::Min(tau->SourcePFJet()->Pt(),99.999),NNLOWeight->GetVal());
      if(tau->LeadPFCand())			      hDWtnSel[25+theHisto]->Fill(TMath::Min(tau->LeadPFCand()->Pt(),99.999),NNLOWeight->GetVal());
      hDWtnSel[27+theHisto]->Fill(TMath::Min((double)tau->NIsoPFCandS(),9.499),NNLOWeight->GetVal());
      hDWtnSel[28+theHisto]->Fill(TMath::Abs(tau->Eta()),NNLOWeight->GetVal());

      if(indGenTau != -1 && TMath::Abs(GenTaus->At(indGenTau)->Eta()) < 2.5 &&
         GenTaus->At(indGenTau)->Pt() > 10.0){
        hDWtnSel[53]->Fill(TMath::Min(GenTaus->At(indGenTau)->Pt(),199.999),NNLOWeight->GetVal());
        hDWtnSel[54]->Fill(TMath::Abs(GenTaus->At(indGenTau)->Eta()),NNLOWeight->GetVal());
      }
    }
  } // fPFTaus1 loop

}
//--------------------------------------------------------------------------------------------------
void WtnEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fPFTaus0Name, fPFTaus0);
  ReqBranch(fPFTaus1Name, fPFTaus1);

  char sb[200];
  // PFTaus0
  for(int j=2; j<4; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWtnSel_%d",ind+ 0);  hDWtnSel[ind+ 0] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 1);  hDWtnSel[ind+ 1] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+ 2);  hDWtnSel[ind+ 2] = new TH1D(sb,sb,8,-0.5,7.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 3);  hDWtnSel[ind+ 3] = new TH1D(sb,sb,8,-0.5,7.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 4);  hDWtnSel[ind+ 4] = new TH1D(sb,sb,200,0.,20.);
    sprintf(sb,"hDWtnSel_%d",ind+ 5);  hDWtnSel[ind+ 5] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+ 6);  hDWtnSel[ind+ 6] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+ 7);  hDWtnSel[ind+ 7] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDWtnSel_%d",ind+ 8);  hDWtnSel[ind+ 8] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+ 9);  hDWtnSel[ind+ 9] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+10);  hDWtnSel[ind+10] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+11);  hDWtnSel[ind+11] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+12);  hDWtnSel[ind+12] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+13);  hDWtnSel[ind+13] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+14);  hDWtnSel[ind+14] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+15);  hDWtnSel[ind+15] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+16);  hDWtnSel[ind+16] = new TH1D(sb,sb,200,0.,2.);
    sprintf(sb,"hDWtnSel_%d",ind+17);  hDWtnSel[ind+17] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+18);  hDWtnSel[ind+18] = new TH1D(sb,sb,200,0.,1.);
    sprintf(sb,"hDWtnSel_%d",ind+19);  hDWtnSel[ind+19] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+20);  hDWtnSel[ind+20] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+21);  hDWtnSel[ind+21] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+22);  hDWtnSel[ind+22] = new TH1D(sb,sb,200,0.,10.);
    sprintf(sb,"hDWtnSel_%d",ind+23);  hDWtnSel[ind+23] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+24);  hDWtnSel[ind+24] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+25);  hDWtnSel[ind+25] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+26);  hDWtnSel[ind+26] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+27);  hDWtnSel[ind+27] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWtnSel_%d",ind+28);  hDWtnSel[ind+28] = new TH1D(sb,sb,50,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+29);  hDWtnSel[ind+29] = new TH1D(sb,sb,200,0.0,0.2);
  }
  for(int i=0; i<30; i++){
    for(int j=2; j<4; j++){
      AddOutput(hDWtnSel[i+j*100]);
    }
  }

  sprintf(sb,"hDWtnSel_%d", 151);  hDWtnSel[151] = new TH1D(sb,sb,200,0.,200.);
  sprintf(sb,"hDWtnSel_%d", 152);  hDWtnSel[152] = new TH1D(sb,sb,25,0.,2.5);
  sprintf(sb,"hDWtnSel_%d", 153);  hDWtnSel[153] = new TH1D(sb,sb,200,0.,200.);
  sprintf(sb,"hDWtnSel_%d", 154);  hDWtnSel[154] = new TH1D(sb,sb,25,0.,2.5);
  for(int i=151; i<155; i++){
    AddOutput(hDWtnSel[i]);
  }

  // PFTaus1
  for(int j=0; j<2; j++){
    int ind = 100 * j;
    sprintf(sb,"hDWtnSel_%d",ind+ 0);  hDWtnSel[ind+ 0] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 1);  hDWtnSel[ind+ 1] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+ 2);  hDWtnSel[ind+ 2] = new TH1D(sb,sb,8,-0.5,7.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 3);  hDWtnSel[ind+ 3] = new TH1D(sb,sb,8,-0.5,7.5);
    sprintf(sb,"hDWtnSel_%d",ind+ 4);  hDWtnSel[ind+ 4] = new TH1D(sb,sb,200,0.,20.);
    sprintf(sb,"hDWtnSel_%d",ind+ 5);  hDWtnSel[ind+ 5] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+ 6);  hDWtnSel[ind+ 6] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+ 7);  hDWtnSel[ind+ 7] = new TH1D(sb,sb,200,0.,200.);
    sprintf(sb,"hDWtnSel_%d",ind+ 8);  hDWtnSel[ind+ 8] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+ 9);  hDWtnSel[ind+ 9] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+10);  hDWtnSel[ind+10] = new TH1D(sb,sb,400,0.,4.);
    sprintf(sb,"hDWtnSel_%d",ind+11);  hDWtnSel[ind+11] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+12);  hDWtnSel[ind+12] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+13);  hDWtnSel[ind+13] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+14);  hDWtnSel[ind+14] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+15);  hDWtnSel[ind+15] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+16);  hDWtnSel[ind+16] = new TH1D(sb,sb,200,0.,2.);
    sprintf(sb,"hDWtnSel_%d",ind+17);  hDWtnSel[ind+17] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+18);  hDWtnSel[ind+18] = new TH1D(sb,sb,200,0.,1.);
    sprintf(sb,"hDWtnSel_%d",ind+19);  hDWtnSel[ind+19] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+20);  hDWtnSel[ind+20] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+21);  hDWtnSel[ind+21] = new TH1D(sb,sb,200,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+22);  hDWtnSel[ind+22] = new TH1D(sb,sb,200,0.,10.);
    sprintf(sb,"hDWtnSel_%d",ind+23);  hDWtnSel[ind+23] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+24);  hDWtnSel[ind+24] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+25);  hDWtnSel[ind+25] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+26);  hDWtnSel[ind+26] = new TH1D(sb,sb,200,0.,100.);
    sprintf(sb,"hDWtnSel_%d",ind+27);  hDWtnSel[ind+27] = new TH1D(sb,sb,10,-0.5,9.5);
    sprintf(sb,"hDWtnSel_%d",ind+28);  hDWtnSel[ind+28] = new TH1D(sb,sb,50,0.,5.);
    sprintf(sb,"hDWtnSel_%d",ind+29);  hDWtnSel[ind+29] = new TH1D(sb,sb,200,0.0,0.2);
    sprintf(sb,"hDWtnSel_%d",ind+30);  hDWtnSel[ind+30] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+31);  hDWtnSel[ind+31] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+32);  hDWtnSel[ind+32] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+33);  hDWtnSel[ind+33] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+34);  hDWtnSel[ind+34] = new TH1D(sb,sb,2,-0.5,1.5);
    sprintf(sb,"hDWtnSel_%d",ind+35);  hDWtnSel[ind+35] = new TH1D(sb,sb,2,-0.5,1.5);
  }
  for(int i=0; i<36; i++){
    for(int j=0; j<2; j++){
      AddOutput(hDWtnSel[i+j*100]);
    }
  }

  sprintf(sb,"hDWtnSel_%d", 51);  hDWtnSel[51] = new TH1D(sb,sb,200,0.,200.);
  sprintf(sb,"hDWtnSel_%d", 52);  hDWtnSel[52] = new TH1D(sb,sb,25,0.,2.5);
  sprintf(sb,"hDWtnSel_%d", 53);  hDWtnSel[53] = new TH1D(sb,sb,200,0.,200.);
  sprintf(sb,"hDWtnSel_%d", 54);  hDWtnSel[54] = new TH1D(sb,sb,25,0.,2.5);
  for(int i=51; i<55; i++){
    AddOutput(hDWtnSel[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void WtnEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WtnEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
