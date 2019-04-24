from nail import *
import ROOT
import sys


def addLheScale(flow):
    flow.VariationWeightArray("LHEScaleWeight",8,filt=lambda sname,hname,wname : "__syst__" not in hname and "__syst__" not in sname ) #systematic variations are 1D, let's avoid systematics of systematic
    #this is not obvious as N replicas can change... think about it
    #flow.AddVariationWeightArray("LHEPdfWeight",30,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic


#create btag systematics
#this should be simplified
def addBtag(flow):
    flow.Define("SelectedJet_btagWeight_up","vector_map(btagWeightUp,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
    #flow.Define("btagEventWeightUp","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
    flow.Systematic("BTagUp","SelectedJet_btagWeight","SelectedJet_btagWeight_up")
    flow.createVariationBranch("BTagUp",["btagEventWeight"])
    flow.VariationWeight("btagEventWeight__syst__BTagUP","btagEventWeight")


def addMuScale(flow):
#Define Systematic variations
    flow.Define("Muon_pt_scaleUp","Muon_pt*1.01f") #this should be protected against systematic variations
    flow.Define("Muon_pt_scaleDown","Muon_pt*0.97f")
    flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
    flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

def addCompleteJecs(flow):
    for i in range(2):
        flow.Define("Jet_pt_JEC%s"%i,"Jet_pt+%s/100.f"%i)
        flow.Systematic("JEC%s"%i,"Jet_pt","Jet_pt_JEC%s"%i) #name, target, replacement


