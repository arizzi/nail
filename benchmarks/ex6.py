from nail import *
import ROOT
import sys


flow=SampleProcessing("Just MET","root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")

# pt of the tri-jet system with mass closest to 172.5 GeV, and leading b-tag discriminator among the 3 jets in the triplet [requires looping on combination of objects in the same collection, 4-vector algebra, and extracting properties of a combination other than the key used to sort them]
flow.Define("Jet_p4","@p4v(Jet)")
flow.Define("Jet_btag","Jet_btagCSV")
flow.Selection("threeJets","nJet>=3")
flow.Distinct("TriJets","Jet",n=3,requires=["threeJets"])
flow.TakeTriplet("SelTriJet","Jet","TriJets","Argmax(-abs(MemberMap((TriJets0_p4+TriJets1_p4+TriJets2_p4),M() )-172.5))",requires=["threeJets"])
flow.Define("SelTriJet_btagMax","std::max(std::max(SelTriJet0_btag,SelTriJet1_btag),SelTriJet2_btag)")

histosPerSelection={
"threeJets" : ["SelTriJet_btagMax"],
}

flow.binningRules = [
(".*btag.*","100,0,1"),
]
flow.printRDFCpp([],debug=False,outname="tmp.C",selections=histosPerSelection)
import os
print "code generated"
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print "code compiled"
os.system("./tmp 4")
print "code run"



