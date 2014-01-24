#ifndef MODS_LINKDEF_H
#define MODS_LINKDEF_H
#include "Analysis/PDFs/interface/pdf.h"
#include "Analysis/PDFs/interface/pdf_reweighting.h"
#include "Analysis/PDFs/interface/MitNtupleEvent.h"
#include "Analysis/PDFs/interface/MitPDFNtupleEvent.h"
#endif


#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::pdf+;
#pragma link C++ class mithep::pdf_reweighting+;
#pragma link C++ class MitNtupleEvent+;
#pragma link C++ class MitPDFNtupleEvent+;
#endif

