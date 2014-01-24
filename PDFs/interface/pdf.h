#ifndef MITWLNU_PDF_H
#define MITWLNU_PDF_H
#include <string>
#include "TObject.h"

namespace mithep  {
  class pdf {
  public:
    pdf(std::string iName,int iPDF,std::string iPDFName="MSTW2008nlo90cl.LHgrid", 
        int nsel = 0, unsigned int nJetsType = 0);
    ~pdf() {}
    ClassDef(pdf, 1)
  };
}
#endif
