from nail import *
import ROOT
import sys


flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")

flow.Define("Muon_p4","@p4v(Muon)")
flow.Selection("twoMuons","nMuon>=2") 
flow.Distinct("MuMu","Muon")
flow.Define("OppositeSignMuMuInRange","Nonzero(MuMu0_charge != MuMu1_charge && abs(MemberMap(MuMu0_p4+MuMu1_p4,M())-90) < 30 )",requires=["twoMuons"])
flow.Selection("pairInRange","OppositeSignMuMuInRange.size() > 0")

histosPerSelection={
"pairInRange" : ["MET_pt"],
}

flow.binningRules = [
(".*_pt","100,0,500"),
]
flow.printRDFCpp([],debug=False,outname="tmp.C",selections=histosPerSelection)
import os
print("code generated")
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print("code compiled")
os.system("./tmp 4")
print("code run")



