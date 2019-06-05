#define hist mapping
genericHistos=["LeadMuon_pt","LeadMuon_eta","SubMuon_pt","SubMuon_eta","QJet0_eta","QJet1_eta","QJet0_pt_nom","QJet1_pt_nom","PV_npvs","LeadingSAJet_pt","NSoft5New2","NSoft2New2","NSoft10New2","SAHT","SAHT5","nFootprintSAJet","FootHT"]
bdtInputHistos=["Mqq_log","mmjj_pt","qqDeltaEta","NSoft5","ll_zstar","Higgs_pt","theta2","mmjj_pz","MaxJetAbsEta","Higgs_m","Higgs_m_uncalib","NSoft5New","ll_zstar_log"]
signalHistosVariations=["SBClassifierNoMass","BDTAtanNoMass"]
signalHistos=["SBClassifier","BDTAtan"]
histosPerSelection={
#"PreSel" : genericHistos+["Mqq"],
"SignalRegion": genericHistos+signalHistos+bdtInputHistos,
#"ZRegion": genericHistos+signalHistosZ,
#"ZRegionNadya": genericHistos+signalHistosZ,
"SideBand" : signalHistos+signalHistosVariations
#"TwoJetsTwoMu" : genericHistos
}

