#ifndef MODS_LINKDEF_H
#define MODS_LINKDEF_H
#include "Ana/PDFs/interface/pdf.h"
#include "Ana/PDFs/interface/pdf_reweighting.h"
#include "Ana/PDFs/interface/MitNtupleEvent.h"
#include "Ana/PDFs/interface/MitPDFNtupleEvent.h"
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

