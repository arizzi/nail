import os
from nail import *
import ROOT
import sys
# Create a group of plots (jet pT, eta, phi, N_jets). Now make it for all events, for events with missing et > 20 GeV, and for events with missing et > 20 GeV and 2 jets with 40 GeV and abs(eta) < 1.0.


flow = SampleProcessing(
    "Just MET", "root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
flow.Selection("Met20", "MET_pt>20")
flow.Selection("twoHighPtJets",
               "Sum(Jet_pt > 40 && abs(Jet_eta) < 1.) >= 2", requires=["Met20"])

histos = ["Jet_pt", "Jet_eta", "Jet_phi", "nJet"]
histosPerSelection = {
    "": histos,
    "Met20": histos,
    "twoHighPtJets": histos
}

flow.binningRules = [
    (".*_pt", "100,0,500"),
    (".*_eta", "100,-5,5"),
    (".*_phi", "100,-3.1416,3.1416"),
    ("nJet", "10,0,10"),
]


processor = flow.CreateProcessor("eventProcessor", [], histosPerSelection, [], "", 4)
rdf = ROOT.RDataFrame("Events", "root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
result = processor(rdf)
c1=ROOT.TCanvas()
for h in result.histos:
    h.Draw()
    c1.SaveAs("Histo"+h.GetName()+".png")
