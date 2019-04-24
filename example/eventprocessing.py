from nail import *
import ROOT
import sys

flow=SampleProcessing("VBF Hmumu Analysis","/scratch/arizzi/Hmm/nail/samples/6B8A2AC8-35E6-1146-B8A8-B1BA90E3F3AA.root")
#/scratch/mandorli/Hmumu/samplePerAndrea/GluGlu_HToMuMu_skim_nano2016.root")
#flow.Define("LHEScaleWeight","ROOT::VecOps::RVec<float>(9,1.)") #this result in NOOP if already defined, otherwise it is a failsafe

#Higgs to mumu reconstruction
flow.Selection("hasHiggs","Sum(GenPart_pdgId == 25) > 0")
flow.Define("GenHiggs_idx","Nonzero(GenPart_pdgId == 25)", requires=["hasHiggs"])
flow.SubCollection("QParton","GenPart",sel="GenPart_genPartIdxMother==At(Take(GenPart_genPartIdxMother,GenHiggs_idx),0,-100) && GenPart_pdgId!= 25")
flow.Define("QParton_p4","@p4v(QParton)")
flow.Distinct("QQ","QParton")
flow.Selection("twoQ","nQParton>=2")
flow.Define("QQ_p4","QQ0_p4+QQ1_p4",requires=["twoQ"])
flow.Define("QQ_mass","MemberMap(QQ_p4,M())")
flow.Define("HighestGenQQMass","At(QQ_mass,Argmax(QQ_mass),-99)")

flow.DefaultConfig(muIsoCut=0.25,muIdCut=0,muPtCut=10, dzCut=0.2,dxyCut=0.05) #cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible
flow.Define("Muon_id","Muon_tightId*4+Muon_mediumId*2+Muon_softId") 
flow.Define("Muon_iso","Muon_pfRelIso04_all")
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut && abs(Muon_dz) < dzCut && abs(Muon_dxy) < dxyCut") 
flow.Define("SelectedMuon_p4","@p4v(SelectedMuon)")
flow.Selection("twoUnpreselMuons","nMuon>=2")
flow.Selection("twoMuons","nSelectedMuon>=2") 
flow.Distinct("MuMu","SelectedMuon")
flow.Define("OppositeSignMuMu","Nonzero(MuMu0_charge != MuMu1_charge)",requires=["twoMuons"])
flow.Selection("twoOppositeSignMuons","OppositeSignMuMu.size() > 0")
flow.TakePair("Mu","SelectedMuon","MuMu","At(OppositeSignMuMu,0,-200)",requires=["twoOppositeSignMuons"])
flow.Define("Higgs","Mu0_p4+Mu1_p4")


#loose leptons for isolation
flow.DefaultConfig(muIsoCutL=0.25,muIdCutL=0,muPtCutL=10)
flow.DefaultConfig(eIsoCutL=0.25,eIdCutL=3,ePtCutL=10)
flow.SubCollection("LooseMuon","Muon",sel="Muon_iso < muIsoCutL && Muon_id > muIdCutL && Muon_pt > muPtCutL && abs(Muon_dz) < dzCut && abs(Muon_dxy) < dxyCut")
flow.Define("Electron_id","Electron_cutBased")
flow.Define("Electron_iso","Electron_pfRelIso03_all")
flow.SubCollection("LooseEle","Electron",sel="Electron_iso < eIsoCutL && Electron_id > eIdCutL && Electron_pt > ePtCutL && abs(Electron_dz) < dzCut && abs(Electron_dxy) < dxyCut")
flow.Define("LooseMuon_p4","@p4v(LooseMuon)")
flow.Define("LooseMuon_pid","(LooseMuon_pt*0)+13")
flow.Define("LooseEle_p4","@p4v(LooseEle)")
flow.Define("LooseEle_pid","(LooseEle_pt*0)+11")
flow.MergeCollections("Lepton",["LooseMuon","LooseEle"])

#Match jets to selected loose Leptons
flow.Define("Jet_p4","@p4v(Jet)")
flow.MatchDeltaR("Jet","Lepton",defIdx=-500)

#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet","Jet_pt > jetPtCut && (Jet_LeptonIdx==-1 || Jet_LeptonDr > 0.3)")
flow.Selection("twoJets","nSelectedJet>=2")
flow.Define("GenJet_p4","@p4v(GenJet)")
flow.Define("SelectedJet_p4","@p4v(SelectedJet)")
flow.Distinct("JetPair","SelectedJet")
flow.TakePair("QJet","SelectedJet","JetPair","Argmax(MemberMap((JetPair0_p4+JetPair1_p4),M() ))",requires=["twoJets"])
#flow.ObjectAt("QJet0","SelectedJet","0",requires=["twoJets"])
#flow.ObjectAt("QJet1","SelectedJet","1",requires=["twoJets"])

#compute number of softjets removing signal footprint
flow.Define("SoftActivityJet_mass","SoftActivityJet_pt*0")
flow.Define("SoftActivityJet_p4","@p4v(SoftActivityJet)")
flow.Match("SelectedJet","SoftActivityJet") #associate signal jets
flow.Match("SelectedMuon","SoftActivityJet") #associate signal muons
flow.Define("NSoft2",'''SoftActivityJetNjets2-Sum( 
(	   (SoftActivityJet_SelectedJetDr<0.2 && ( SoftActivityJet_SelectedJetIdx == QJet0 ||  SoftActivityJet_SelectedJetIdx == QJet1)) ||
	   (SoftActivityJet_SelectedMuonDr<0.2 && ( SoftActivityJet_SelectedJetIdx == Mu0 || SoftActivityJet_SelectedJetIdx == Mu1) ) || 
	   (SoftActivityJet_eta > std::max(QJet0_eta, QJet1_eta) || SoftActivityJet_eta < std::min(QJet0_eta, QJet1_eta))
)&&
	   SoftActivityJet_pt > 2. 
)''')


flow.Define("qq","QJet0_p4+QJet1_p4")
flow.Define("Mqq","qq.M()")
flow.Define("MqqGenJet","(QJet0_genJetIdx>=0&&QJet1_genJetIdx>=0&&QJet0_genJetIdx<nGenJet&&QJet1_genJetIdx<nGenJet)?(At(GenJet_p4,QJet0_genJetIdx)+At(GenJet_p4,QJet1_genJetIdx)).M():-99")
flow.Define("qq_pt","qq.Pt()")
flow.Define("qqDeltaEta","abs(QJet0_eta-QJet1_eta)")
flow.Define("qqDeltaPhi","abs(ROOT::Math::VectorUtil::DeltaPhi(QJet0_p4,QJet1_p4))")

#QQ vs ll kinematic
flow.Define("ll_ystar","Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())")
flow.Define("ll_zstar"," abs( ll_ystar/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )) ")
flow.Define("DeltaEtaQQSum","abs(QJet0_eta) +  abs(QJet1_eta)")
flow.Define("PhiZQ1","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet0_p4))")
flow.Define("PhiZQ2","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet1_p4))")
flow.Define("EtaHQ1","abs(Higgs.Eta() - QJet0_eta)")
flow.Define("EtaHQ2","abs(Higgs.Eta() - QJet1_eta)")
flow.Define("DeltaRelQQ","(QJet0_p4+QJet1_p4).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt())")
flow.Define("Rpt","(QJet0_p4+QJet1_p4+ Higgs).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt() + Higgs.Pt())")

flow.DefaultConfig(higgsMassWindowWidth=15,mQQcut=400,nominalHMass=125.03)
flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Selection("VBFRegion","Mqq > mQQcut && QJet0_pt > 35")
flow.Selection("PreSel","VBFRegion && twoOppositeSignMuons",requires=["VBFRegion","twoOppositeSignMuons"])
flow.Selection("SideBand","! MassWindow && VBFRegion",requires=["VBFRegion"])
flow.Selection("SignalRegion","VBFRegion && MassWindow", requires=["VBFRegion","MassWindow"])
flow.Selection("TwoJetsTwoMu","twoJets && twoOppositeSignMuons", requires=["twoJets","twoOppositeSignMuons"])

#flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion") 
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("SBClassifier","Higgs_pt+Higgs_m+Mqq+Rpt+DeltaRelQQ+NSoft2") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])
flow.ObjectAt("LeadMuon","SelectedMuon","0",requires=["twoMuons"])
flow.ObjectAt("SubMuon","SelectedMuon","1",requires=["twoMuons"])


