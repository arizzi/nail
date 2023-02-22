from .nail import *
import ROOT
import sys

#flow=SampleProcessing("Example Analysis","root://xrootd-cms.infn.it//store/mc/RunIIAutumn18NanoAODv4/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/270000/320474A7-2A79-E042-BD91-BD48021177A2.root")
flow = SampleProcessing("Example Analysis", "/gpfs/ddn/srm/cms//store/mc/RunIIAutumn18NanoAODv4/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/270000/320474A7-2A79-E042-BD91-BD48021177A2.root")

# cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible
flow.DefaultConfig(muIsoCut=0.25, muIdCut=0,
                   muPtCut=10, dzCut=0.2, dxyCut=0.05)
flow.Define("Muon_id", "Muon_tightId*4+Muon_mediumId*2+Muon_softId")
flow.Define("Muon_iso", "Muon_pfRelIso04_all")
flow.SubCollection("SelectedMuon", "Muon",
                   sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut && abs(Muon_dz) < dzCut && abs(Muon_dxy) < dxyCut")
flow.Define("SelectedMuon_p4", "@p4v(SelectedMuon)")
flow.Selection("oneMuon", "nSelectedMuon > 0")
# here "requires" actually means "implies"... to be adjusted
flow.Selection("twoMuons", "nSelectedMuon>=2", requires=["oneMuon"])
flow.Distinct("MuMu", "SelectedMuon")
flow.Define("OppositeSignMuMu",
            "Nonzero(MuMu0_charge != MuMu1_charge)", requires=["twoMuons"])
flow.Selection("twoOppositeSignMuons", "OppositeSignMuMu.size() > 0")
flow.TakePair("Mu", "SelectedMuon", "MuMu",
              "OppositeSignMuMu[0]", requires=["twoOppositeSignMuons"])
flow.Define("DiMuon", "Mu0_p4+Mu1_p4")

# loose leptons for isolation
flow.DefaultConfig(muIsoCutL=0.25, muIdCutL=0, muPtCutL=10)
flow.DefaultConfig(eIsoCutL=0.25, eIdCutL=3, ePtCutL=10)
flow.SubCollection("LooseMuon", "Muon",
                   sel="Muon_iso < muIsoCutL && Muon_id > muIdCutL && Muon_pt > muPtCutL && abs(Muon_dz) < dzCut && abs(Muon_dxy) < dxyCut")
flow.Define("Electron_id", "Electron_cutBased")
flow.Define("Electron_iso", "Electron_pfRelIso03_all")
flow.SubCollection("LooseEle", "Electron",
                   sel="Electron_iso < eIsoCutL && Electron_id > eIdCutL && Electron_pt > ePtCutL && abs(Electron_dz) < dzCut && abs(Electron_dxy) < dxyCut")
flow.Define("LooseMuon_p4", "@p4v(LooseMuon)")
flow.Define("LooseMuon_pid", "(LooseMuon_pt*0)+13")
flow.Define("LooseEle_p4", "@p4v(LooseEle)")
flow.Define("LooseEle_pid", "(LooseEle_pt*0)+11")
flow.MergeCollections("Lepton", ["LooseMuon", "LooseEle"])

# Match jets to selected loose Leptons
flow.Define("Jet_p4", "@p4v(Jet)")
flow.MatchDeltaR("Jet", "Lepton")

# Jet Selection
flow.DefaultConfig(jetPtCut=25)
# pt cut plus DeltaR cleaning
flow.SubCollection("SelectedJet", "Jet",
                   "Jet_pt > jetPtCut && (Jet_LeptonIdx==-1 || Jet_LeptonDr > 0.3)")
flow.MatchDeltaR("GenJet", "SelectedJet")
flow.Define("GenJet_p4", "@p4v(GenJet)")

# Feature request=> Sum with non T(0) start value
flow.Define("JetSum_p4", "std::accumulate(SelectedJet_p4.begin(), SelectedJet_p4.end(), ROOT::Math::PtEtaPhiMVector())")
flow.Define("HT", "Sum(SelectedJet_pt)")
flow.Define("MHT", "-JetSum_p4.pt()")


flow.Define("Z_pt", "DiMuon.Pt()")
flow.Define("Z_m", "DiMuon.M()")

# define some event weights
flow.CentralWeight("genWeight")
flow.CentralWeight("btagWeight_CSVV2")
flow.ObjectAt("LeadMuon", "SelectedMuon", "0", requires=["oneMuon"])
flow.ObjectAt("SubMuon", "SelectedMuon", "1", requires=["twoMuons"])

# Define Systematic variations via weights
flow.VariationWeightArray("LHEScaleWeight", 8, filt=lambda selname,
                          hname, wname: "__syst__" not in hname and "__syst__" not in selname)

# Define Systematic variations via replacement of inputs
flow.Define("Muon_pt_scaleUp", "Muon_pt*1.01f")
flow.Define("Muon_pt_scaleDown", "Muon_pt*0.97f")
# name, target, replacement
flow.Systematic("MuScaleDown", "Muon_pt", "Muon_pt_scaleDown")
# name, target, replacement
flow.Systematic("MuScaleUp", "Muon_pt", "Muon_pt_scaleUp")

# define hist mapping
genericHistos = ["LeadMuon_pt", "LeadMuon_eta", "HT", "MHT"]
signalRegionHistos = ["Z_pt", "Z_m", "SubMuon_pt"]
histosPerSelection = {
    "oneMuon": genericHistos,
    "twoOppositeSignMuons": genericHistos+signalRegionHistos,
}

systematics = ["MuScaleDown", "MuScaleUp"]
histosWithSystematics = flow.createSystematicBranches(
    systematics, histosPerSelection)

print("Histograms I will produce in various regions")
for sel in histosWithSystematics:
    print((sel, ":", histosWithSystematics[sel]))

snap = genericHistos+signalRegionHistos + \
    ["nMuon"]+[x for x in flow.validCols if x[:5] == "Muon_"]
print(("columns selected for snapshot ntuple", snap))
flow.printRDFCpp(snap, debug=False,
                 outname=sys.argv[1], selections=histosWithSystematics, snap=snap, snapsel="twoOppositeSignMuons")

##
# generate and compile with
# python  helloworld.py test.C && g++  -fPIC -Wall -O3 test.C $(root-config --libs --cflags)  -o test
