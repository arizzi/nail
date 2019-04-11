from nail import *
import ROOT
import sys


#define some event weights
def addBtagWeight(flow):
    flow.Define("SelectedJet_btagWeight","vector_map(btagWeight,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
    flow.Define("btagEventWeight","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
    flow.CentralWeight("genWeight")
    flow.CentralWeight("btagEventWeight")

def addMuEffWeight(flow):
    flow.Define("muEffWeight","effMu2016(LeadMuon_pt,LeadMuon_eta)*effMu2016(SubMuon_pt,SubMuon_eta)")
    flow.CentralWeight("muEffWeight",["twoMuons"])


