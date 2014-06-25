#ifndef MITWLNU_PDFFILLER_H
#define MITWLNU_PDFFILLER_H
#include <string>
#include "TObject.h"

namespace mithep  {
  class pdfFiller {
  public:
    pdfFiller(std::string iName);
    ~pdfFiller() {}
    ClassDef(pdfFiller, 1)
  };
}
#endif
