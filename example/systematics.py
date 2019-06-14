from nail import *
import ROOT
import sys


def addLheScale(flow):  

    #we really need only 0,1,3
    flow.VariationWeightArray("LHEScaleWeightSafe",4,filt=lambda sname,hname,wname : "__syst__" not in hname and "__syst__" not in sname ) #systematic variations are 1D, let's avoid systematics of systematic
    #this is not obvious as N replicas can change... think about it
    #flow.AddVariationWeightArray("LHEPdfWeight",30,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic


#create btag systematics
#this should be simplified
def addBtag(flow):
    pass
#    flow.Define("SelectedJet_btagWeight_up","vector_map(btagWeightUp,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
    #flow.Define("btagEventWeightUp","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
#    flow.Systematic("BTagUp","SelectedJet_btagWeight","SelectedJet_btagWeight_up")
#    flow.createVariationBranch("BTagUp",["btagEventWeight"])
#    flow.VariationWeight("btagEventWeight__syst__BTagUP","btagEventWeight")


def addMuScale(flow):
#Define Systematic variations
    flow.Systematic("MuScaleDown","Muon_corrected_pt","Muon_correctedDown_pt") #name, target, replacement
    flow.Systematic("MuScaleUp","Muon_corrected_pt","Muon_correctedUp_pt") #name, target, replacement

def addCompleteJecs(flow):
        flow.Systematic("JERDown","Jet_pt_touse","Jet_pt_jerDown") #name, target, replacement
        flow.Systematic("JERUp","Jet_pt_touse","Jet_pt_jerUp") #name, target, replacement
        flow.Systematic("JESDown","Jet_pt_touse","Jet_pt_jesTotalDown") #name, target, replacement
        flow.Systematic("JESUp","Jet_pt_touse","Jet_pt_jesTotalUp") #name, target, replacement
        flow.Systematic("WithJER","Jet_pt_touse","Jet_pt_nom") #name, target, replacement


