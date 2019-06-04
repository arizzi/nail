#define hist mapping
genericHistos=["LeadMuon_pt","LeadMuon_eta","SubMuon_pt","SubMuon_eta","QJet0_pt","QJet1_pt","QJet0_eta","QJet1_eta","QJet0_pt_nom","QJet1_pt_nom","PV_npvs","LeadingSAJet_pt","NSoft5New2","NSoft2New2","NSoft10New2","SAHT","SAHT5","nFootprintSAJet","FootHT"]
signalHistos=["Mqq_log","mmjj_pt","qqDeltaEta","NSoft5","ll_zstar","Higgs_pt","theta2","mmjj_pz","MaxJetAbsEta","Higgs_m","SBClassifier","BDTAtan","Higgs_m_uncalib","NSoft5New","ll_zstar_log"]
signalHistosZ=["Higgs_m","SBClassifierZ","BDTAtanZ","Higgs_m_uncalib","Mqq_log","mmjj_pt","qqDeltaEta","NSoft5","ll_zstar","Higgs_pt","theta2","mmjj_pz","MaxJetAbsEta","NSoft5New","ll_zstar_log"]
histosPerSelection={
"PreSel" : genericHistos+["Mqq"],
"SignalRegion": genericHistos+signalHistos,
"ZRegion": genericHistos+signalHistosZ,
"ZRegionNadya": genericHistos+signalHistosZ,
"SideBand" : genericHistos+signalHistos,
#"TwoJetsTwoMu" : genericHistos
}

