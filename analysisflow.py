from nail import *
flow.SetAlias("n(.*)", "\\1.size()", defaultPersitency=True)
flow.SetAlias(
    "(.*)_p4", "{TLorentzVector ret; ret.SetPtEtaPhiM(\\1_pt,\\1_eta,\\1_phi,\\1_mass); return ret;}", defaultPersistency=False)
# SubCikkectuib actuib"
flow.SetAlias("SelectedMuon_(.*)([\.*\])", "Muon_\1[SelectedMuon[\2]]")

flow = SampleProcessing("")

# cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible
flow.DefaultConfig(muIsoCut=0.13, muIdCut=3, muPtCut=25)
# Higgs to mumu reconstruction
# Maps to plain RDF VecOps
flow.DefineCollAttr("Muon_id", "Muon_tightId*3+Muon_looseId")
# this should generate some kind of wrapper/ref that can be used as the parent collection
flow.SubCollection("SelectedMuon", "Muon",
                   sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut")

flow.Filter("twoOppositeSignMuons",
            "nSelectedMuon==2 && SelectedMuon_charge[0]*SelectedMuon_charge[1] < 0")
# p4 should be handled somehow ... any syntax is ok such as p4(SelectedMuon[0]) or _p4 or .p4 etc..
flow.Define("Higgs", "p4at(SelectedMuon,0)+p4at(SelectedMuon,1)",
            requires=["twoOppositeSignMuons"])
# the following could work
# define p4at(x,y)  ROOT::Math::PtEtaPhiMVector(x##_pt[y] , x##_eta[y], x##_phi[y], x##_mass[y])
# define p4(x)  ROOT::Math::PtEtaPhiMVector(x##_pt , x##_eta, x##_phi, x##_mass)


# VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet", "Jet",
                   "Jet_pt > jetPtCut && (Jet_muonIdx1 == -1 || Muon_iso[Jet_muonIdx1] > muIsoCut || Muon_id[Jet_muonIdx1] > 0")
flow.Filter("twoJets", "nSelectedJet>=2")
flow.Define("Qjet1", "SelectedJet[0].p4()", requires=["twoJets"])
flow.Define("Qjet2", "SelectedJet[1].p4()", requires=["twoJets"])
flow.Define("qq", "Qjet1+Qjet2")
flow.Define("Mqq", "qq.M()")
flow.Define("qq_pt", "qq.Pt()")
flow.Define("qqDeltaEta", "TMath::Abs(Qjet1.Eta()-Qjet2.Eta())")
flow.Define("qqDeltaPhi", "TMath::Abs(Qjet1.DeltaPhi(Qjet2))")

# QQ vs ll kinematic
flow.Define(
    "ll_ystar", "Higgs.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity())")
flow.Define(
    "ll_zstar", " TMath::Abs( ll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ")
flow.Define("DeltaEtaQQSum",
            "TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta())")
flow.Define("PhiZQ1", "TMath::Abs(Higgs.DeltaPhi(Qjet1))")
flow.Define("PhiZQ2", "TMath::Abs(Higgs.DeltaPhi(Qjet2))")
flow.Define("EtaHQ1", "TMath::Abs(Higgs.Eta() - Qjet1.Eta())")
flow.Define("EtaHQ2", "TMath::Abs(Higgs.Eta() - Qjet2.Eta())")
flow.Define("DeltaRelQQ", "(Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt())")
flow.Define(
    "Rpt", "(Qjet1+Qjet2+ Higgs).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Higgs.Pt())")

flow.DefaultConfig(higgsMassWindowWidth=15, mQQcut=400, nominalHMass=125.03)
flow.Filter("MassWindow", "abs(Higgs_m-nominalHMass)<higgsMassWindowWidth")
flow.Filter("SideBand", "! MassWindow")
flow.Filter("VBFRegion", "Mqq > mQQcut")
flow.Filter("SignalRegion", "VBFRegion && MassWindow")

# flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion")


print(flow.NeededInputs())

# flow.AddSystematic("MuScaleUp","Muon_pt","Muon_pt*1.01") #name, target, replacement
# flow.AddSystematic("HMassUncertainityUp","nominalHMass","125.1") #name, target, replacement
# flow.OptimizationScan("MuPtCutScan","muPtCut","30") #name, target, replacement

#from samples import background,signal,data
