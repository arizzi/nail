import os
from nail import *
import ROOT
import sys
flow = SampleProcessing(
    "Just MET", "root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
histosPerSelection = {"": ["MET_pt", "Jet_pt"]}
flow.binningRules = [(".*_pt", "100,0,500")]
processor = flow.CreateProcessor("eventProcessor", [], histosPerSelection, [], "", 4)
rdf = ROOT.RDataFrame("Events", "root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root")
result = processor(rdf)
for h in result.histos:
        h.Draw()
                                                                                                                                
