from nail.nail import *
import ROOT
import sys

#flow=SampleProcessing("VBF Hmumu Analysis","/scratch/arizzi/Hmm/nail/samples/6B8A2AC8-35E6-1146-B8A8-B1BA90E3F3AA.root")
flow=SampleProcessing("VBF Hmumu Analysis","/scratch/mandorli/Hmumu/fileSkimFromNanoAOD/fileSkim2016_Z/VBF_HToMuMu_nano2016.root")
#/scratch/mandorli/Hmumu/samplePerAndrea/GluGlu_HToMuMu_skim_nano2016.root")
#flow.Define("LHEScaleWeight","ROOT::VecOps::RVec<float>(9,1.)") #this result in NOOP if already defined, otherwise it is a failsafe

#variables that we will add file by file before passing the RNode to the event processor
flow.AddExpectedInput("year","int")
flow.AddExpectedInput("isMC","bool")
#flow.AddExpectedInput("Muon_sf","ROOT::VecOps::RVec<float>")

flow.Define("LHEScaleWeightSafe","nLHEScaleWeight>=8?LHEScaleWeight:std::vector<float>(9,1)")
flow.Define("Jet_pt_touse","Jet_pt")

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
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id >= muIdCut && Muon_corrected_pt > muPtCut && abs(Muon_dz) < dzCut && abs(Muon_dxy) < dxyCut") 
flow.Define("SelectedMuon_p4","vector_map_t<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >        >(SelectedMuon_corrected_pt , SelectedMuon_eta, SelectedMuon_phi, SelectedMuon_mass)")
flow.Define("SelectedMuon_p4uncalib","@p4v(SelectedMuon)")
flow.Selection("twoUnpreselMuons","nMuon>=2")
flow.Selection("twoMuons","nSelectedMuon>=2") 
flow.Distinct("MuMu","SelectedMuon")
flow.Define("OppositeSignMuMu","Nonzero(MuMu0_charge != MuMu1_charge)",requires=["twoMuons"])
flow.Selection("twoOppositeSignMuons","OppositeSignMuMu.size() > 0")
flow.TakePair("Mu","SelectedMuon","MuMu","At(OppositeSignMuMu,0,-200)",requires=["twoOppositeSignMuons"])
flow.Define("Higgs","Mu0_p4+Mu1_p4")
flow.Define("HiggsUncalib","Mu0_p4uncalib+Mu1_p4uncalib")

flow.SubCollection("GenLepton","GenPart",sel="(abs(GenPart_pdgId)==13 || abs(GenPart_pdgId)==11 || abs(GenPart_pdgId)==15)")
flow.MatchDeltaR("GenLepton","GenJet") 
flow.SubCollection("GenJetVBFFilter","GenJet",sel="GenJet_GenLeptonDr>0.3 || GenJet_GenLeptonIdx==-1 ")
flow.Define("GenJetVBFFilter_p4","@p4v(GenJetVBFFilter)")
flow.Selection("twoVBFFilterGenJet","nGenJetVBFFilter > 1")
flow.Define("VBFFilterjj_p4","At(GenJetVBFFilter_p4,0)+At(GenJetVBFFilter_p4,1)",requires=["twoVBFFilterGenJet"])
flow.Define("VBFFilterjj_mass","VBFFilterjj_p4.M()")
flow.Selection("VBFFilterFlag", "VBFFilterjj_mass>350")
flow.Selection("VBFFilterAntiFlag", "!VBFFilterFlag")

flow.Define("Jet_p4","vector_map_t<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >        >(Jet_pt_touse , Jet_eta, Jet_phi, Jet_mass)")
#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet",'''
(year != 2017 ||  Jet_pt_touse > 50 || abs(Jet_eta) < 2.7 || abs(Jet_eta) > 3.0 ||  Jet_neEmEF<0.55 ) && 
Jet_pt_touse > jetPtCut && Jet_puId >0 &&  Jet_jetId > 0 && abs(Jet_eta) < 4.7 && (abs(Jet_eta)<2.5 || Jet_puId > 6) && 
(Jet_muonIdx1==-1 || TakeDef(Muon_pfRelIso04_all,Jet_muonIdx1,100) > 0.25 || abs(TakeDef(Muon_dz,Jet_muonIdx1,100)) > 0.2 || abs(TakeDef(Muon_dxy,Jet_muonIdx1,100) > 0.05)) &&
(Jet_muonIdx2==-1 || TakeDef(Muon_pfRelIso04_all,Jet_muonIdx2,100) > 0.25 || abs(TakeDef(Muon_dz,Jet_muonIdx2,100)) > 0.2 || abs(TakeDef(Muon_dxy,Jet_muonIdx2,100) > 0.05)) 
''')
flow.Selection("twoJets","nSelectedJet>=2")
flow.Define("GenJet_p4","@p4v(GenJet)")
#flow.Define("SelectedJet_p4","@p4v(SelectedJet)")
flow.Distinct("JetPair","SelectedJet")
#flow.TakePair("QJet","SelectedJet","JetPair","Argmax(MemberMap((JetPair0_p4+JetPair1_p4),M() ))",requires=["twoJets"])
flow.Define("SortedSelectedJetIndices","Argsort(-SelectedJet_pt_touse)")
flow.ObjectAt("QJet0","SelectedJet","SortedSelectedJetIndices[0]",requires=["twoJets"])
flow.ObjectAt("QJet1","SelectedJet","SortedSelectedJetIndices[1]",requires=["twoJets"])

#compute number of softjets removing signal footprint
flow.Define("SoftActivityJet_mass","SoftActivityJet_pt*0")
flow.Define("SoftActivityJet_p4","@p4v(SoftActivityJet)")
flow.Match("SelectedJet","SoftActivityJet") #associate signal jets
flow.Match("SelectedMuon","SoftActivityJet") #associate signal muons
flow.Define("NSoft5Saved","Nonzero(SoftActivityJet_pt > 5.).size()")
flow.SubCollection("SAJet","SoftActivityJet",'''
! (   (SoftActivityJet_SelectedJetDr<0.4 && ( SoftActivityJet_SelectedJetIdx == QJet0 ||  SoftActivityJet_SelectedJetIdx == QJet1)) ||
           (SoftActivityJet_SelectedMuonDr<0.4 && ( SoftActivityJet_SelectedMuonIdx == Mu0 || SoftActivityJet_SelectedMuonIdx == Mu1) ) || 
           (SoftActivityJet_eta > std::max(QJet0_eta, QJet1_eta) || SoftActivityJet_eta < std::min(QJet0_eta, QJet1_eta)))
''')
flow.SubCollection("FootprintSAJet","SoftActivityJet",'''
 (   (SoftActivityJet_SelectedJetDr<0.4 && ( SoftActivityJet_SelectedJetIdx == QJet0 ||  SoftActivityJet_SelectedJetIdx == QJet1)) ||
           (SoftActivityJet_SelectedMuonDr<0.4 && ( SoftActivityJet_SelectedMuonIdx == Mu0 || SoftActivityJet_SelectedMuonIdx == Mu1) ) || 
           (SoftActivityJet_eta > std::max(QJet0_eta, QJet1_eta) || SoftActivityJet_eta < std::min(QJet0_eta, QJet1_eta)))
''')

flow.Define("NSoft5",'''Sum( 
 (  (SoftActivityJet_SelectedJetDr>0.2 || ( SoftActivityJet_SelectedJetIdx != QJet0 &&  SoftActivityJet_SelectedJetIdx != QJet1)) &&
	   (SoftActivityJet_SelectedMuonDr>0.2 || ( SoftActivityJet_SelectedMuonIdx != Mu0 && SoftActivityJet_SelectedMuonIdx != Mu1) ) && 
	   (SoftActivityJet_eta < std::max(QJet0_eta, QJet1_eta) && SoftActivityJet_eta > std::min(QJet0_eta, QJet1_eta))
)&&
	   SoftActivityJet_pt > 5. 
)''')

flow.Define("LeadingSAJet_pt","At(SAJet_pt,0,-1)")

flow.Define("NSoft5New",'''SoftActivityJetNjets5-Sum( 
(	   (SoftActivityJet_SelectedJetDr<0.4 && ( SoftActivityJet_SelectedJetIdx == QJet0 ||  SoftActivityJet_SelectedJetIdx == QJet1)) ||
	   (SoftActivityJet_SelectedMuonDr<0.4 && ( SoftActivityJet_SelectedMuonIdx == Mu0 || SoftActivityJet_SelectedMuonIdx == Mu1) ) || 
	   (SoftActivityJet_eta > std::max(QJet0_eta, QJet1_eta) || SoftActivityJet_eta < std::min(QJet0_eta, QJet1_eta))
)&&
	   SoftActivityJet_pt > 5. 
)''')
flow.Define("NSoft5New2","SoftActivityJetNjets5-Sum(FootprintSAJet_pt>5)")
flow.Define("NSoft2New2","SoftActivityJetNjets2-Sum(FootprintSAJet_pt>2)")
flow.Define("NSoft10New2","SoftActivityJetNjets10-Sum(FootprintSAJet_pt>10)")
flow.Define("FootHT","Sum(FootprintSAJet_pt)")
flow.Define("SAHT","SoftActivityJetHT-Sum(FootprintSAJet_pt)")
flow.Define("SAHT5","SoftActivityJetHT5-Sum(FootprintSAJet_pt*(FootprintSAJet_pt>5.f))")



flow.Define("qq","QJet0_p4+QJet1_p4")
flow.Define("Mqq","qq.M()")
flow.Define("MqqGenJet","(QJet0_genJetIdx>=0&&QJet1_genJetIdx>=0&&QJet0_genJetIdx<nGenJet&&QJet1_genJetIdx<nGenJet)?(At(GenJet_p4,QJet0_genJetIdx)+At(GenJet_p4,QJet1_genJetIdx)).M():-99")
flow.Define("qq_pt","qq.Pt()")
flow.Define("qqDeltaEta","abs(QJet0_eta-QJet1_eta)")
flow.Define("qqDeltaPhi","abs(ROOT::Math::VectorUtil::DeltaPhi(QJet0_p4,QJet1_p4))")

#QQ vs ll kinematic
flow.Define("ll_ystar","Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())/2.f")
flow.Define("ll_zstar"," abs( ll_ystar/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )) ")
flow.Define("ll_zstar_log"," log(ll_zstar) ")
flow.Define("ll_zstarbug_log","log(abs( (  Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())  )/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )))")
flow.Define("DeltaEtaQQSum","abs(QJet0_eta) +  abs(QJet1_eta)")
flow.Define("PhiZQ1","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet0_p4))")
flow.Define("PhiZQ2","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet1_p4))")
flow.Define("EtaHQ1","abs(Higgs.Eta() - QJet0_eta)")
flow.Define("EtaHQ2","abs(Higgs.Eta() - QJet1_eta)")
flow.Define("DeltaRelQQ","(QJet0_p4+QJet1_p4).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt())")
flow.Define("Rpt","(QJet0_p4+QJet1_p4+ Higgs).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt() + Higgs.Pt())")
flow.Define("mmjj","Higgs+qq")
flow.Define("theta2","Higgs.Vect().Dot(QJet1_p4.Vect())/QJet1_p4.Vect().R()/Higgs.Vect().R()")
flow.ObjectAt("LeadMuon","SelectedMuon","0",requires=["twoMuons"])
flow.ObjectAt("SubMuon","SelectedMuon","1",requires=["twoMuons"])
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("Higgs_m_uncalib","HiggsUncalib.M()")
flow.Define("Mqq_log","log(Mqq)")
flow.Define("mmjj_pt","mmjj.Pt()")
flow.Define("mmjj_pt_log","log(mmjj_pt)")
flow.Define("mmjj_pz","mmjj.Pz()")
flow.Define("mmjj_pz_logabs","log(abs(mmjj_pz))")
flow.Define("MaxJetAbsEta","std::max(std::abs(QJet0_eta), std::abs(QJet1_eta))")


flow.DefaultConfig(higgsMassWindowWidth=10,mQQcut=250,nominalHMass=125.03,btagCut=0.8)
flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Selection("MassWindowZ","abs(Higgs.M()-91)<15")
flow.Selection("VBFRegion","Mqq > mQQcut && QJet0_pt_touse> 35 && QJet1_pt_touse > 25")
flow.Selection("PreSel","VBFRegion && twoOppositeSignMuons && Max(SelectedJet_btagDeepB) < btagCut && LeadMuon_pt > 30 && SubMuon_pt > 20 && HLT_IsoMu24 && abs(SubMuon_eta) <2.4 && abs(LeadMuon_eta) < 2.4",requires=["VBFRegion","twoOppositeSignMuons"])
flow.Selection("SideBand","Higgs_m < 150 && Higgs_m > 110 && ! MassWindow && VBFRegion &&  qqDeltaEta > 2.5",requires=["VBFRegion"])
flow.Selection("SignalRegion","VBFRegion && MassWindow && Mqq > 400 &&  qqDeltaEta > 2.5", requires=["VBFRegion","MassWindow","PreSel"])
flow.Selection("ZRegion","VBFRegion && MassWindowZ  && qqDeltaEta > 2.5", requires=["VBFRegion","MassWindowZ","PreSel"])
flow.Selection("ZRegionNadya","VBFRegion && MassWindowZ && QJet0_pt_touse> 50 && QJet1_pt_touse > 30 ", requires=["VBFRegion","MassWindowZ","PreSel"])
flow.Selection("TwoJetsTwoMu","twoJets && twoOppositeSignMuons", requires=["twoJets","twoOppositeSignMuons"])

flow.Define("SBClassifier","mva.eval(__slot,{Higgs_m,Mqq_log,mmjj_pt_log,qqDeltaEta,float(NSoft5),ll_zstarbug_log,Higgs_pt,theta2,mmjj_pz_logabs,MaxJetAbsEta})") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])
flow.Define("SBClassifierNoMass","mva.eval(__slot,{125.,Mqq_log,mmjj_pt_log,qqDeltaEta,float(NSoft5),ll_zstarbug_log,Higgs_pt,theta2,mmjj_pz_logabs,MaxJetAbsEta})") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])
#flow.Define("SBClassifierZ","mva.eval({125.,log(Mqq),mmjj.Pt(),qqDeltaEta,float(NSoft5),ll_zstar,Higgs.Pt(),theta2,mmjj.Pz(),std::max(std::abs(QJet0_eta), std::abs(QJet1_eta))})") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])
flow.Define("BDTAtan","atanh((SBClassifier+1.)/2.)")
flow.Define("BDTAtanNoMass","atanh((SBClassifierNoMass+1.)/2.)")
flow.Selection("BDT0p8","BDTAtan>0.8")
flow.Selection("BDT1p0","BDTAtan>1.0")
flow.Selection("BDT1p1","BDTAtan>1.1")
flow.Selection("BDT1p2","BDTAtan>1.2")

