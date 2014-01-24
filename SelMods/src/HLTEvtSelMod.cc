// $Id: HLTEvtSelMod.cc,v 1.3 2010/11/12 10:31:40 ceballos Exp $

#include "Analysis/SelMods/interface/HLTEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TTree.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"

using namespace mithep;
ClassImp(mithep::HLTEvtSelMod)

//--------------------------------------------------------------------------------------------------
HLTEvtSelMod::HLTEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fMetName("randomMet"),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fVertexName(ModNames::gkGoodVertexesName),
  fVertices(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  //Obtain all the good objects from the event cleaning modul
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ParticleOArr *leptons        = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  LoadBranch(fEvtHdrName);

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    printf("HLTEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  MetOArr *CleanMet            = GetObjThisEvt<MetOArr>(fMetName);
  const Met *caloMet           = CleanMet->At(0);

  // Minimun Pt, Nleptons>=2 requirements
  if (leptons->GetEntries() >= 2 &&
      leptons->At(0)->Pt() > 20 && leptons->At(1)->Pt() > 10){

    CompositeParticle dilepton;;
    dilepton.AddDaughter(leptons->At(0));
    dilepton.AddDaughter(leptons->At(1));

    int pairType = -1;
    if (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon )
      pairType = 0;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kElectron)
      pairType = 1;
    else if(leptons->At(0)->ObjType() == kElectron && leptons->At(1)->ObjType() == kMuon)
      pairType = 2;
    else if(leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kElectron)
      pairType = 3;
    else {
      cout << "Hey, this is not possible, leptonTypes: "
    	   << leptons->At(0)->ObjType() << " - " 
           << leptons->At(1)->ObjType() << endl;
    }

    if((TMath::Abs(dilepton.Mass()-91.1876) < 30.0 || 
       (pairType >= 2 && dilepton.Mass() > 12.0)) && 
       dilepton.Charge() == 0){
      // Count the number of central Jets for vetoing
      int nCentralJets = 0;
      for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
        if(TMath::Abs(CleanJets->At(i)->Eta()) < fEtaJetCut &&
	   CleanJets->At(i)->Pt() > fPtJetCut){
          nCentralJets++;
        }
      }
      double HLTMu[2][5]  = {{0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.}};
      double HLTEl[2][10] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
      Int_t ents=objs->GetEntries();
      for(Int_t i=0;i<ents;++i) {
        const TriggerObject* to = objs->At(i);
        for(int nl=0; nl<2; nl++){
          if     (leptons->At(nl)->ObjType() == kMuon && (to->Type() == TriggerObject::TriggerMuon || to->Type() == TriggerObject::TriggerCluster) &&
	          MathUtils::DeltaR(leptons->At(nl)->Mom(), to->Mom()) < 0.3){
            if(strcmp(to->TrigName(),"HLT_Mu15_v1")==0)       HLTMu[nl][0] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu11_Ele8_v1")==0)  HLTMu[nl][1] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu5_Ele13_v2")==0)  HLTMu[nl][2] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu5_Ele17_v1")==0)  HLTMu[nl][3] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu8_Ele8_v1")==0)   HLTMu[nl][4] = 1;
	  }
          else if(leptons->At(nl)->ObjType() == kElectron && (to->Type() == TriggerObject::TriggerElectron || 
	                                                      to->Type() == TriggerObject::TriggerCluster  ||
							      to->Type() == TriggerObject::TriggerPhoton) &&
		  MathUtils::DeltaR(leptons->At(nl)->Mom(), to->Mom()) < 0.3){
            if(strcmp(to->TrigName(),"HLT_DoubleEle17_SW_L1R_v1")==0)                 HLTEl[nl][0] = 1;
            if(strcmp(to->TrigName(),"HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1")==0) HLTEl[nl][1] = 1;
            if(strcmp(to->TrigName(),"HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")==0)      HLTEl[nl][2] = 1;
            if(strcmp(to->TrigName(),"HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")==0)      HLTEl[nl][3] = 1;
            if(strcmp(to->TrigName(),"HLT_Ele22_SW_TighterEleId_L1R_v2")==0)          HLTEl[nl][4] = 1;
            if(strcmp(to->TrigName(),"HLT_Ele17_SW_TightEleId_L1R")==0)               HLTEl[nl][5] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu11_Ele8_v1")==0)                          HLTEl[nl][6] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu5_Ele13_v2")==0)                          HLTEl[nl][7] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu5_Ele17_v1")==0)                          HLTEl[nl][8] = 1;
            if(strcmp(to->TrigName(),"HLT_Mu8_Ele8_v1")==0)                           HLTEl[nl][9] = 1;
          }
        }
      } // Loop over HLT objects
      fTreeVariables[ 0] = HLTMu[0][0];
      fTreeVariables[ 1] = HLTMu[0][1];
      fTreeVariables[ 2] = HLTMu[0][2];
      fTreeVariables[ 3] = HLTMu[0][3];
      fTreeVariables[ 4] = HLTMu[0][4];
      fTreeVariables[ 5] = HLTMu[1][0];
      fTreeVariables[ 6] = HLTMu[1][1];
      fTreeVariables[ 7] = HLTMu[1][2];
      fTreeVariables[ 8] = HLTMu[1][3];
      fTreeVariables[ 9] = HLTMu[1][4];
      fTreeVariables[10] = HLTEl[0][0];
      fTreeVariables[11] = HLTEl[0][1];
      fTreeVariables[12] = HLTEl[0][2];
      fTreeVariables[13] = HLTEl[0][3];
      fTreeVariables[14] = HLTEl[0][4];
      fTreeVariables[15] = HLTEl[0][5];
      fTreeVariables[16] = HLTEl[0][6];
      fTreeVariables[17] = HLTEl[0][7];
      fTreeVariables[18] = HLTEl[0][8];
      fTreeVariables[19] = HLTEl[0][9];
      fTreeVariables[20] = HLTEl[1][0];
      fTreeVariables[21] = HLTEl[1][1];
      fTreeVariables[22] = HLTEl[1][2];
      fTreeVariables[23] = HLTEl[1][3];
      fTreeVariables[24] = HLTEl[1][4];
      fTreeVariables[25] = HLTEl[1][5];
      fTreeVariables[26] = HLTEl[1][6];
      fTreeVariables[27] = HLTEl[1][7];
      fTreeVariables[28] = HLTEl[1][8];
      fTreeVariables[29] = HLTEl[1][9];
      fTreeVariables[30] = pairType;
      fTreeVariables[31] = dilepton.Mass();
      fTreeVariables[32] = caloMet->Pt();
      fTreeVariables[33] = nCentralJets;
      fTreeVariables[34] = leptons->At(0)->Pt();
      fTreeVariables[35] = leptons->At(0)->Eta();
      fTreeVariables[36] = leptons->At(0)->Phi();
      fTreeVariables[37] = leptons->At(1)->Pt();
      fTreeVariables[38] = leptons->At(1)->Eta();
      fTreeVariables[39] = leptons->At(1)->Phi();
      fTreeVariables[40] = fVertices->GetEntries();
      fTreeVariables[41] = fEventHeader->RunNum();
      fTree->Fill();
    } // Mass or emu requirements
  } // Minimun Pt and Nleptons >= 2 requirements
}
//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fEvtHdrName,      fEventHeader);

  //***********************************************************************************************
  //Create Ntuple Tree  
  //***********************************************************************************************
  printf("... init HLT ntuple ...\n");
  fTree = new TTree("HLTTree", "HLTTree");
  const char* TreeFormat;
  TreeFormat = "hltm00/F:hltm01:hltm02:hltm03:hltm04:hltm10:hltm11:hltm12:hltm13:hltm14:hlte00:hlte01:hlte02:hlte03:hlte04:hlte05:hlte06:hlte07:hlte08:hlte09:hlte10:hlte11:hlte12:hlte13:hlte14:hlte15:hlte16:hlte17:hlte18:hlte19:ltype:mass:met:njets:pt0:eta0:phi0:pt1:eta1:phi1:vert:run";
  fTree->Branch("H", &fTreeVariables,TreeFormat);
  AddOutput(fTree);
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}
