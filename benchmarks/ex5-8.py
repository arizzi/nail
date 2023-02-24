import os
from nail import *
import ROOT
import sys


flow = SampleProcessing(
    "Just MET", "root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")

# etmiss for events that have an opposite-sign muon pair of mass 60-120 GeV [looping on combinations of objects in one collection]
flow.Define("Muon_p4", "@p4v(Muon)")
flow.Selection("twoMuons", "nMuon>=2")
flow.Distinct("MuMu", "Muon")
flow.Define("OppositeSignMuMuInRange",
            "Nonzero(MuMu0_charge != MuMu1_charge && abs(MemberMap(MuMu0_p4+MuMu1_p4,M())-90) < 30 )", requires=["twoMuons"])
flow.Selection("pairInRange", "OppositeSignMuMuInRange.size() > 0")

# pt of the tri-jet system with mass closest to 172.5 GeV, and leading b-tag discriminator among the 3 jets in the triplet [requires looping on combination of objects in the same collection, 4-vector algebra, and extracting properties of a combination other than the key used to sort them]
flow.Define("Jet_p4", "@p4v(Jet)")
flow.Define("Jet_btag", "Jet_btagCSV")
flow.Selection("threeJets", "nJet>=3")
flow.Distinct("TriJets", "Jet", n=3, requires=["threeJets"])
flow.TakeTriplet("SelTriJet", "Jet", "TriJets",
                 "Argmax(-abs(MemberMap((TriJets0_p4+TriJets1_p4+TriJets2_p4),M() )-172.5))", requires=["threeJets"])
flow.Define("SelTriJet_btagMax",
            "std::max(std::max(SelTriJet0_btag,SelTriJet1_btag),SelTriJet2_btag)")

# - sum of the pt of all jets of pt > 30 GeV that are not within DR 0.4 from a lepton of pt > 10 GeV [looping on two separate collections]
flow.Define("Electron_p4", "@p4v(Electron)")
flow.Define("Electron_pid", "Electron_pt*0+11")
flow.Define("Muon_pid", "Muon_pt*0+13")
flow.MergeCollections("Lepton", ["Muon", "Electron"])
flow.SubCollection("SelectedLepton", "Lepton", "Lepton_pt > 10")
flow.MatchDeltaR("Jet", "SelectedLepton")
flow.SubCollection("CleanJet", "Jet", '''
Jet_pt > 30 &&
(Jet_SelectedLeptonIdx==-1 || Jet_SelectedLeptonDr > 0.4)
''')
flow.Define("SumJet_pt", "Sum(CleanJet_pt)")
# - in events with >=3 leptons and a same-flavour opposite-sign lepton pair, find the best same-flavour opposite-sign lepton pair (mass closest to 91.2 GeV), and plot the transverse mass of the missing energy and the leading other lepton [ something whose formulation in an imperative language is easy, but whose translations to a functional language may be less clear and/or possibly inefficient]
flow.Define("Lepton_index", "Range(nLepton)")
flow.Distinct("LPair", "Lepton")
flow.Selection("twoLeptons", "nLepton>=2")
flow.Define("isOSSF", "LPair0_charge != LPair1_charge && LPair0_pid == LPair1_pid",
            requires=["twoLeptons"])
flow.Selection("hasOSSF", "Sum(isOSSF) > 0")
flow.TakePair("Z", "Lepton", "LPair",
              "Argmax(-abs(MemberMap((LPair0_p4+LPair1_p4),M() )*isOSSF-91.2))", requires=["hasOSSF"])
flow.Selection("threeLeptons", "nLepton>=3",
               requires=["twoLeptons", "hasOSSF"])
flow.SubCollection("ResidualLeptons", "Lepton",
                   sel="Lepton_index != LPair0[Z_index] && Lepton_index != LPair1[Z_index]", requires=["threeLeptons"])
flow.ObjectAt("ResidualLepton", "ResidualLeptons",
              "Argmax(ResidualLeptons_pt)")
flow.Define("MET_eta", "0.f")
flow.Define("MET_mass", "0.f")
flow.Define("MET_p4", "@p4(MET)")
flow.Define("METplusLepton_p4", "ResidualLepton_p4+MET_p4")
flow.Define("METplusLepton_Mt", "METplusLepton_p4.Mt()")

histosPerSelection = {
    "": ["SumJet_pt"],
    "threeJets": ["SelTriJet_btagMax"],
    "pairInRange": ["MET_pt"],
    "threeLeptons": ["METplusLepton_Mt"]
}

flow.binningRules = [
    (".*_pt", "100,0,500"),
    (".*btag.*", "100,0,1"),
    (".*Mt.*", "100,0,1000")
]
flow.printRDFCpp([], debug=False, outname="tmp.C",
                 selections=histosPerSelection)
print("code generated")
os.system("g++  -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o tmp")
print("code compiled")
os.system("./tmp 4")
print("code run")
