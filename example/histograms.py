#define hist mapping
genericHistos=["LeadMuon_pt","LeadMuon_eta","SubMuon_pt","SubMuon_eta","QJet0_pt","QJet1_pt","QJet0_eta","QJet1_eta"]
signalHistos=["Higgs_m","SBClassifier"]
histosPerSelection={
"PreSel" : genericHistos+["Mqq"],
"SignalRegion": genericHistos+signalHistos,
"SideBand" : genericHistos+signalHistos,
"TwoJetsTwoMu" : genericHistos
}

