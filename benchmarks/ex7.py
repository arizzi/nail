from nail import *
import ROOT
import sys


flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")

#- sum of the pt of all jets of pt > 30 GeV that are not within DR 0.4 from a lepton of pt > 10 GeV [looping on two separate collections]
flow.Define("Electron_p4","@p4v(Electron)")
flow.Define("Muon_p4","@p4v(Muon)")
flow.Define("Electron_pid","Electron_pt*0+11")
flow.Define("Muon_pid","Muon_pt*0+13")
flow.MergeCollections("Lepton",["Muon","Electron"])
flow.SubCollection("SelectedLepton","Lepton","Lepton_pt > 10")
flow.MatchDeltaR("Jet","SelectedLepton")
flow.SubCollection("CleanJet","Jet",'''
Jet_pt > 30 &&
(Jet_SelectedLeptonIdx==-1 || Jet_SelectedLeptonDr > 0.4)
''')
flow.Define("SumJet_pt","Sum(CleanJet_pt)")

histosPerSelection={
"" : ["SumJet_pt"],
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



