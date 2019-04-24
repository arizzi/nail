from nail import *
import ROOT
import sys


flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")

#- in events with >=3 leptons and a same-flavour opposite-sign lepton pair, find the best same-flavour opposite-sign lepton pair (mass closest to 91.2 GeV), and plot the transverse mass of the missing energy and the leading other lepton [ something whose formulation in an imperative language is easy, but whose translations to a functional language may be less clear and/or possibly inefficient]
flow.Define("Muon_p4","@p4v(Muon)")
flow.Define("Electron_p4","@p4v(Electron)")
flow.Define("Electron_pid","Electron_pt*0+11")
flow.Define("Muon_pid","Muon_pt*0+13")
flow.MergeCollections("Lepton",["Muon","Electron"])
flow.Define("Lepton_index","Range(nLepton)")
flow.Distinct("LPair","Lepton")
flow.Selection("twoLeptons","nLepton>=2")
flow.Define("isOSSF","LPair0_charge != LPair1_charge && LPair0_pid == LPair1_pid",requires=["twoLeptons"])
flow.Selection("hasOSSF","Sum(isOSSF) > 0")
flow.TakePair("Z","Lepton","OSPair","Argmax(-abs(MemberMap((LPair0_p4+LPair1_p4),M() )*isOSSF-91.2))",requires=["hasOSSF"])
flow.Selection("threeLeptons","nLepton>=3",requires=["twoLeptons","hasOSSF"]) 
flow.SubCollection("ResidualLeptons","Lepton",sel="Lepton_index != LPair0[Z_indices] && Lepton_index != LPair1[Z_indices]",requires=["threeLeptons"])
flow.ObjectAt("ResidualLepton","ResidualLeptons","Argmax(ResidualLeptons_pt)")
flow.Define("MET_eta","0.f")
flow.Define("MET_mass","0.f")
flow.Define("MET_p4","@p4(MET)")
flow.Define("METplusLepton_p4","ResidualLepton_p4+MET_p4")
flow.Define("METplusLepton_Mt","METplusLepton_p4.Mt()")

histosPerSelection={
"threeLeptons" : ["METplusLepton_Mt"]
}

flow.binningRules = [
(".*Mt.*","100,0,1000")
]
flow.printRDFCpp([],debug=False,outname="tmp.C",selections=histosPerSelection)
import os
print "code generated"
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print "code compiled"
os.system("./tmp 4")
print "code run"



