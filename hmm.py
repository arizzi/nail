from nail import *
import ROOT
from histobinning import binningrules

#muons="Muon"
#jets="Jet"


def buildFlow(flow):
    flow.Define("Muon_m","0*Muon_pfRelIso04_all+0.1056f")
    flow.Define("Muon_p4","vector_map_t<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >        >(Muon_pt , Muon_eta, Muon_phi, Muon_m)")
    flow.Define("Muon_iso","Muon_pfRelIso04_all")
    flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < 0.25 && Muon_mediumId && Muon_pt > 20. && abs(Muon_eta) < 2.4") 
    flow.Selection("twoUnpreselMuons","nMuon>=2")
    flow.Selection("twoMuons","nSelectedMuon==2") 
    flow.Distinct("MuMu","SelectedMuon")
    flow.Define("OppositeSignMuMu","Nonzero(MuMu0_charge != MuMu1_charge)",requires=["twoMuons"])
    flow.Selection("twoOppositeSignMuons","OppositeSignMuMu.size() > 0")
    flow.TakePair("Mu","SelectedMuon","MuMu","At(OppositeSignMuMu,0,-200)",requires=["twoOppositeSignMuons"])
    flow.Define("Higgs","Mu0_p4+Mu1_p4")
    flow.Define("Higgs_m","Higgs.M()")
    flow.Define("Jet_p4","vector_map_t<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >        >(Jet_pt , Jet_eta, Jet_phi, Jet_mass)")

    flow.Match("Jet","Muon") #associate Jets to muon

    flow.SubCollection("SelectedJet","Jet",'''
    Jet_pt > 25 && Jet_btagDeepB >= 0 && /* ( Jet_pt > 50 
    || (  Jet_puId  > 0 ))
     &&   Jet_jetId  > 0  && */ abs(Jet_eta) < 4.7  &&  
    (Jet_MuonIdx==-1 || Jet_MuonDr>0.4 || TakeDef(Muon_iso,Jet_MuonIdx,100) > 0.25 || abs(TakeDef(Muon_pt,Jet_MuonIdx,0)) < 20 || abs(TakeDef(Muon_mediumId,Jet_MuonIdx,0) == 0 )) 
     
    ''')
    flow.Selection("twoJets","nSelectedJet>=2")
    flow.Distinct("JetPair","SelectedJet")
    flow.Define("SortedSelectedJetIndices","Argsort(-SelectedJet_pt)")
    flow.ObjectAt("QJet0","SelectedJet",'At(SortedSelectedJetIndices,0)',requires=["twoJets"])
    flow.ObjectAt("QJet1","SelectedJet",'At(SortedSelectedJetIndices,1)',requires=["twoJets"])
    


    flow.Define("qq","QJet0_p4+QJet1_p4")
    flow.Define("Mqq","qq.M()")
    flow.Define("qq_pt","qq.Pt()")
    flow.Define("qqDeltaEta","abs(QJet0_eta-QJet1_eta)")
    flow.Define("qqDeltaPhi","abs(ROOT::Math::VectorUtil::DeltaPhi(QJet0_p4,QJet1_p4))")

    #QQ vs ll kinematic
    flow.Define("ll_ystar","Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())/2.f")
    flow.Define("ll_zstar"," abs( ll_ystar/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )) ")
    flow.Define("ll_zstar_log"," log(ll_zstar) ")
    flow.Define("ll_zstarbug_log","log(abs( (  Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())  )/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )))")
    flow.Define("DeltaEtaQQSum","abs(QJet0_eta) +  abs(QJet1_eta)")
    flow.Define("PhiHQ1","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet0_p4))")
    flow.Define("PhiHQ2","abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet1_p4))")
    flow.Define("EtaHQ1","abs(Higgs.Eta() - QJet0_eta)")
    flow.Define("EtaHQ2","abs(Higgs.Eta() - QJet1_eta)")
    flow.Define("DeltaRelQQ","(QJet0_p4+QJet1_p4).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt())")
    flow.Define("Rpt","(QJet0_p4+QJet1_p4+ Higgs).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt() + Higgs.Pt())")
    flow.Define("mmjj","Higgs+qq")
    flow.Define("theta2","Higgs.Vect().Dot(QJet1_p4.Vect())/QJet1_p4.Vect().R()/Higgs.Vect().R()")
    flow.Define("SortedSelectedMuonIndices","Argsort(-SelectedMuon_pt)")
    flow.ObjectAt("LeadMuon","SelectedMuon","SortedSelectedMuonIndices[0]",requires=["twoMuons"])
    flow.ObjectAt("SubMuon","SelectedMuon","SortedSelectedMuonIndices[1]",requires=["twoMuons"])
    flow.Define("Higgs_pt","Higgs.Pt()")
    flow.Define("Higgs_rapidity","Higgs.Rapidity()")

    flow.Define("Higgs_eta","Higgs.Eta()")
    flow.Define("Higgs_m_uncalib","HiggsUncalib.M()")
    flow.Define("Mqq_log","log(Mqq)")
    flow.Define("Mqq_over400_log","log(Mqq/400)")
    flow.Define("mmjj_pt","mmjj.Pt()")
    flow.Define("mmjj_pt_log","log(mmjj_pt)")
    flow.Define("mmjj_pz","mmjj.Pz()")
    flow.Define("mmjj_pz_logabs","log(abs(mmjj_pz))")
    flow.Define("MaxJetAbsEta","std::max(std::abs(QJet0_eta), std::abs(QJet1_eta))")
    flow.Define("minEtaHQ","std::min(abs(EtaHQ1),(EtaHQ2))")
    flow.Define("minPhiHQ","std::min(abs(PhiHQ1),abs(PhiHQ2))")

    flow.DefaultConfig(higgsMassWindowWidth=10,mQQcut=400,nominalHMass=125.03) #,btagCut=0.8)
    flow.Define("btagCut","0.4184f")
    flow.Define("btagCutL","0.1241f")
    #adding for sync
    flow.Define("nbtagged","int(Nonzero(SelectedJet_btagDeepB > btagCut && abs(SelectedJet_eta)< 2.5).size())")
    flow.Define("nbtaggedL","int(Nonzero(SelectedJet_btagDeepB > btagCutL && abs(SelectedJet_eta)< 2.5).size())")
    flow.Define("nelectrons","int(0)")
#    flow.Define("nelectrons","int(Nonzero(Electron_pt > 20 && abs(Electron_eta) < 2.5 && Electron_mvaFall17V2Iso_WP90 ).size())")


    flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
    flow.Selection("MassWindowZ","abs(Higgs.M()-91)<15")
    flow.Selection("VBFRegion","Mqq > mQQcut && QJet0_pt> 35 && QJet1_pt > 25")
    flow.Selection("PreSel","nelectrons==0 && nbtaggedL < 2 && VBFRegion && twoOppositeSignMuons && nbtagged < 1 && LeadMuon_pt > 26 && SubMuon_pt > 20  && abs(SubMuon_eta) <2.4 && abs(LeadMuon_eta) < 2.4",requires=["VBFRegion","twoOppositeSignMuons"])
    flow.Selection("PreSelSimple","VBFRegion && twoOppositeSignMuons && LeadMuon_pt > 26 && SubMuon_pt > 20  && abs(SubMuon_eta) <2.4 && abs(LeadMuon_eta) < 2.4",requires=["VBFRegion","twoOppositeSignMuons"])
    flow.Selection("SideBand","Higgs_m < 150 && Higgs_m > 110 && ! MassWindow && VBFRegion &&  qqDeltaEta > 2.5",requires=["VBFRegion","PreSel"])
    flow.Selection("SignalRegion","VBFRegion && MassWindow &&  qqDeltaEta > 2.5", requires=["VBFRegion","MassWindow","PreSel"])
    flow.Selection("ZRegion","VBFRegion && MassWindowZ  && qqDeltaEta > 2.5", requires=["VBFRegion","MassWindowZ","PreSel"])
    return flow


flow=SampleProcessing("Simple Test","/gpfs/ddn/users/fvaselli/muonData/synth_nanoaod.root")
flow=buildFlow(flow)

flow1=SampleProcessing("Simple Test","/scratch/arizzi/0088F3A1-0457-AB4D-836B-AC3022A0E34F.root")
flow1=buildFlow(flow1)
#flow=SampleProcessing("Simple Test"," /store/mc/RunIISummer20UL18NanoAODv2/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/230000/0088F3A1-0457-AB4D-836B-AC3022A0E34F.root")
flow.binningRules = binningrules
flow1.binningRules = binningrules

nthreads=50
histos={
   "SignalRegion" :["LeadMuon_pt","Higgs_m","ll_zstar_log","Rpt"],
   "ZRegion" :["LeadMuon_pt","Higgs_m","ll_zstar_log","Rpt"],
    "PreSel" :["LeadMuon_pt","ll_zstar_log","Rpt","Higgs_pt","Higgs_eta","minEtaHQ","LeadMuon_pt","SubMuon_pt","Higgs_m"],
    "": ["SelectedMuon_pt"],
    "VBFRegion": ["SelectedMuon_pt"],
    "twoOppositeSignMuons": ["SelectedMuon_pt","LeadMuon_pt","SubMuon_pt"],
    "PreSelSimple": ["SelectedMuon_pt","LeadMuon_pt","SubMuon_pt"]
}
#targets=["nMuon","Muon_charge"]
targets=["SelectedMuon","Mqq_log","Rpt","qqDeltaEta","ll_zstar_log","minEtaHQ","Higgs_pt","Higgs_eta","Mqq","QJet0_pt","QJet1_pt","QJet0_eta","QJet1_eta","QJet0_phi","QJet1_phi","QJet0_qgl","QJet1_qgl","Higgs_m"]
print flow.Describe(targets+["PreSel"])
#argets=["nGenJet","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_GenMuonDr","GenJet_GenMuonIdx","nGenMuon","GenMuon_pt","GenMuon_eta","GenMuon_phi","GenJet_goodGenMuonIdx"]
processor = flow.CreateProcessor("eventProcessor", targets , histos, [], "", nthreads)
processor1 = flow1.CreateProcessor("eventProcessor1", targets , histos, [], "", nthreads)

#rdf = ROOT.RDataFrame("Events", "/gpfs/ddn/users/fvaselli/muonData/synth_nanoaod.root") 
#rdf1 = ROOT.RDataFrame("Events", "/scratch/arizzi/0088F3A1-0457-AB4D-836B-AC3022A0E34F.root") 

for sam in ["DY","HMM","EWK"] :  
    rdf = ROOT.RDataFrame("Events", "/gpfs/ddn/users/fvaselli/muonData/synth_%s_nanoaod.root"%sam) 
    rdf1 = ROOT.RDataFrame("Events", "/gpfs/ddn/users/fvaselli/nanoaods/%s1.root"%sam) 
    result=processor(rdf)
    result1=processor1(rdf1)

    cc=[]
    for (o,s) in zip(result1.histos, result.histos) :
       cc.append(ROOT.TCanvas())
       o.Draw()
       s.SetLineColor(2)
       s.Draw("sames")
       cc[-1].SaveAs("/afs/pi.infn.it/user/arizzi/public_html/flashsim/%s--%s.png"%(sam,o.GetName()))

    processed_rdf=result.rdf.find("").second
    processed_rdf1=result1.rdf.find("").second

    processed_rdf.Snapshot("Events","out%s.root"%sam,targets)
    processed_rdf1.Snapshot("Events","out1%s.root"%sam,targets)

