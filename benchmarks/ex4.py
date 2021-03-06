from nail import *
import ROOT
import sys
flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
flow.Selection("twoHighPtJets","Sum(Jet_pt > 40 && abs(Jet_eta) < 1.0) >= 2")
histosPerSelection={
"twoHighPtJets",["MET_pt"]
}

flow.binningRules = [(".*_pt","100,0,500")]


flow.printRDFCpp([],debug=False,outname="tmp.C",selections=histosPerSelection)
import os
print "code generated"
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print "code compiled"
os.system("./tmp 4")
print "code run"



