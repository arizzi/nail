from nail import *
import ROOT
import sys
flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
histosPerSelection={"" : ["MET_pt"]}
flow.printRDFCpp([],debug=False,outname="tmp.C",selections=histosPerSelection)
import os
print("code generated")
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print("code compiled")
os.system("./tmp 4")
print("code run")



