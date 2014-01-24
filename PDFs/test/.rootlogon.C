// $Id: .rootlogon.C,v 1.1 2010/10/17 12:32:12 ceballos Exp $
{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libAnaSelMods.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libAnaPDFs.so");
}
