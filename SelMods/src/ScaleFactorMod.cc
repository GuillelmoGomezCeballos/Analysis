// $Id: ScaleFactorMod.cc,v 1.3 2011/04/18 22:28:42 ceballos Exp $

#include "Analysis/SelMods/interface/ScaleFactorMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

using namespace mithep;
ClassImp(mithep::ScaleFactorMod)

//--------------------------------------------------------------------------------------------------
ScaleFactorMod::ScaleFactorMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fObjsName("randomSet"),
  fIsData(kTRUE),
  fApplyVertex(kFALSE),
  fVertexName(ModNames::gkGoodVertexesName),
  fVertices(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ScaleFactorMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ScaleFactorMod::Process()
{
  
  // Make sure nothing is done for data
  if(fIsData == kTRUE) return;

  // Obtain all the good leptons from the event cleaning module
  ParticleOArr *leptons            = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  MCParticleOArr *GenLeptons       = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);

  // Variable to be reweighted
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  // Good vertexes
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  // Trigger objects
  const TriggerObjectCol *objs     = GetHLTObjects(fObjsName);
  if (!objs){
    printf("ScaleFactorMod::TriggerObjectCol not found\n");
    return;
  }

  for (UInt_t i=0; i<leptons->GetEntries(); ++i) {
    Bool_t isGenLepton = kFALSE; 
    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
      MCParticle *gen = GenLeptons->At(j);
      if(MathUtils::DeltaR(gen->Mom(), leptons->At(i)->Mom()) < 0.10) {
    	isGenLepton = kTRUE;
    	break;
      }
    }
    if(isGenLepton == kFALSE) continue;
    double sf1 = FactorHLT((leptons->At(i)->ObjType() == kMuon),
                            leptons->At(i)->Pt(),
		  	    leptons->At(i)->AbsEta());
    double sf2 = FactorEff((leptons->At(i)->ObjType() == kMuon),
    		            leptons->At(i)->Pt(),
	 	            leptons->At(i)->AbsEta());
    NNLOWeight->SetVal(NNLOWeight->GetVal()*sf1*sf2);
  }
  
  if(fApplyVertex == kTRUE){
    double sf = FactorVrt(fVertices->GetEntries());
    NNLOWeight->SetVal(NNLOWeight->GetVal()*sf);
  }
}

//--------------------------------------------------------------------------------------------------
void ScaleFactorMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.
}

//--------------------------------------------------------------------------------------------------
void ScaleFactorMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ScaleFactorMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
double ScaleFactorMod::FactorHLT(Bool_t isMu, Double_t pt, Double_t eta)
{
  if     (isMu == kTRUE){
    if(eta < 1.479) return 1.0;
    else            return 1.0;
  }
  else {
    if(eta < 1.479) return 1.0;
    else            return 1.0;
  }
  
  return 1.0;
}
//--------------------------------------------------------------------------------------------------
double ScaleFactorMod::FactorEff(Bool_t isMu, Double_t pt, Double_t eta)
{
  if     (isMu == kTRUE){
    if(eta < 1.479) return 1.0;
    else            return 1.0;
  }
  else {
    if(eta < 1.479) return 1.0;
    else            return 1.0;
  }
  
  return 1.0;
}
//--------------------------------------------------------------------------------------------------
double ScaleFactorMod::FactorVrt(Int_t vrt)
{
  if	 (vrt == 0) return 0.00000;
  else if(vrt == 1) return 0.21591;
  else if(vrt == 2) return 0.74238;
  else if(vrt == 3) return 1.41497;
  else if(vrt == 4) return 1.78058;
  else if(vrt == 5) return 1.79655;
  else if(vrt == 6) return 1.48655;
  else if(vrt == 7) return 1.09271;
  else if(vrt == 8) return 0.75290;
  else if(vrt == 9) return 0.55617;
  else if(vrt ==10) return 0.36933;
  else if(vrt ==11) return 0.25578;
  else if(vrt ==12) return 0.16628;
  else if(vrt ==13) return 0.11484;
  else if(vrt ==14) return 0.10664;
  else if(vrt ==15) return 0.05628;
  else if(vrt ==16) return 0.01000;
  else if(vrt >=17) return 0.01000;

  return 1.0;
}
