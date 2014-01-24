#!/bin/sh

/home/ceballos/releases/CMSSW_5_2_3_patch3/src/LandS/test/lands.exe -tB 1000 -tPB 30 -M Bayesian -d $1 --doExpectation 1 -t 10000
