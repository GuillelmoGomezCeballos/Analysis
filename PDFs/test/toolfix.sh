#!/bin/bash

#To Check LHAPDF do this  --> 5.8.4
scramv1 tool info LHAPDF
#default is /afs/cern.ch/cms/slc5_ia32_gcc434/external/lhapdf/5.6.0-cms2  etc....

#Then do a 
scramv1 setup -i LHAPDF
#
#LHAPDF_BASE=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt/
#LIBDIR=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt//lib
#INCLUDE=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt//include
#LHAPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets
#LHAPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets
#
#Check again 
scramv1 tool info LHAPDF
#Should be
#SCRAM_PROJECT=no
#LHAPDF_BASE=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt/
#LHAPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets
#LIB=LHAPDF
#LIBDIR=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt//lib
#INCLUDE=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/i686-slc5-gcc43-opt//include
#USE=f77compiler
#LHAPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/share/lhapdf/PDFsets

#Finally
scramv1 b ToolUpdated_lhapdf

