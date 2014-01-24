#!/bin/sh

/home/ceballos/releases/CMSSW_5_2_3_patch3/src/LandS/test/lands.exe -M Hybrid -d $1 --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 1000 --nToysForCLb 500 --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 -rMin 0.1 -rMax 50
