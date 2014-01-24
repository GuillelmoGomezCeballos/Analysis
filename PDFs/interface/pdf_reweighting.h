#ifndef MITWLNU_PDF_REWEIGHTING_H
#define MITWLNU_PDF_REWEIGHTING_H
#include <string>
#include "TObject.h"

namespace mithep  {
  class pdf_reweighting {
  public:
    pdf_reweighting(std::string iName, int iPDF, std::string iPDFName="cteq66.LHgrid");
    ~pdf_reweighting() {}
    ClassDef(pdf_reweighting, 1)
  };
}
#endif
