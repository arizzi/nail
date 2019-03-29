from nail import *
import ROOT
import sys

#f=ROOT.TFile.Open("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root")
#f=ROOT.TFile.Open("/dev/shm/VBF_HToMuMu_nano2016.root")
#e=f.Get("Events")
#allbranches=[(x.GetName(),x.GetListOfLeaves()[0].GetTypeName()) for x in e.GetListOfBranches()]
#df=ROOT.RDataFrame("Events","/dev/shm/VBF_HToMuMu_nano2016.root")
#dftypes={x[0]:df.GetColumnType(x[0]) for x in allbranches}


#print >>sys.stderr, "ROOT loaded"
#flow=SampleProcessing("",["Muon_pt","Muon_eta","Muon_phi","Muon_tightId","Muon_looseId","Jet_pt","Muon_iso","Jet_muonIdx1","Jet_eta","Jet_phi","Jet_mass"])

#flow=SampleProcessing("VBF Hmumu Analysis","/dev/shm/VBF_HToMuMu_nano2016.root")
#flow=SampleProcessing("VBF Hmumu Analysis","root://xrootd-cms.infn.it:1194//store/mc/RunIIFall17NanoAODv4/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/00000/62A8C0F2-EB7B-2D45-B7E2-35DA0D54A774.root")
flow=SampleProcessing("VBF Hmumu Analysis","/scratch/mandorli/Hmumu/samplePerAndrea/GluGlu_HToMuMu_skim_nano2016.root")

#Higgs to mumu reconstruction
flow.Selection("hasHiggs","Sum(GenPart_pdgId == 25) > 0")
flow.Define("GenHiggs_idx","Nonzero(GenPart_pdgId == 25)", requires=["hasHiggs"])
flow.SubCollection("QParton","GenPart",sel="GenPart_genPartIdxMother==Take(GenPart_genPartIdxMother,GenHiggs_idx)[0] && GenPart_pdgId!= 25")
flow.Define("QParton_p4","@p4v(QParton)")
flow.Distinct("QQ","QParton")
flow.Selection("twoQ","nQParton>=2")
flow.Define("QQ_p4","QQ0_p4+QQ1_p4",requires=["twoQ"])
flow.Define("QQ_mass","MemberMap(QQ_p4,M())")
flow.Define("HighestGenQQMass","QQ_mass[Argmax(QQ_mass)]")

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
flow.TakePair("Mu","SelectedMuon","MuMu","OppositeSignMuMu[0]",requires=["twoOppositeSignMuons"])
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
flow.MatchDeltaR("Jet","Lepton")

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
flow.MatchDeltaR("SelectedJet","SoftActivityJet") #associate signal jets
flow.MatchDeltaR("SelectedMuon","SoftActivityJet") #associate signal muons
flow.Define("NSoft2",'''SoftActivityJetNjets2-Sum( 
(	   (SoftActivityJet_SelectedJetDr<0.2 && ( SoftActivityJet_SelectedJetIdx == QJet0 ||  SoftActivityJet_SelectedJetIdx == QJet1)) ||
	   (SoftActivityJet_SelectedMuonDr<0.2 && ( SoftActivityJet_SelectedJetIdx == Mu0 || SoftActivityJet_SelectedJetIdx == Mu1) ) || 
	   (SoftActivityJet_eta > std::max(QJet0_eta, QJet1_eta) || SoftActivityJet_eta < std::min(QJet0_eta, QJet1_eta))
)&&
	   SoftActivityJet_pt > 2. 
)''')


flow.Define("qq","QJet0_p4+QJet1_p4")
flow.Define("Mqq","qq.M()")
flow.Define("MqqGenJet","(QJet0_genJetIdx>=0&&QJet1_genJetIdx>=0&&QJet0_genJetIdx<nGenJet&&QJet1_genJetIdx<nGenJet)?(GenJet_p4[QJet0_genJetIdx]+GenJet_p4[QJet1_genJetIdx]).M():-99")
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
flow.Selection("SideBand","! MassWindow")
flow.Selection("VBFRegion","Mqq > mQQcut && QJet0_pt > 35")
flow.Selection("SignalRegion","VBFRegion && MassWindow")

#flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion") 
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("SBClassifier","Higgs_pt+Higgs_m+Mqq+Rpt+DeltaRelQQ+NSoft2") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])


#define some event weights
flow.Define("SelectedJet_btagWeight","vector_map(btagWeight,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
flow.Define("btagEventWeight","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
flow.CentralWeight("genWeight")
flow.CentralWeight("btagEventWeight")
flow.ObjectAt("LeadMuon","SelectedMuon","0")
flow.ObjectAt("SubMuon","SelectedMuon","1")
flow.Define("muEffWeight","effMu2016(LeadMuon_pt,LeadMuon_eta)*effMu2016(SubMuon_pt,SubMuon_eta)")
flow.CentralWeight("muEffWeight",["twoMuons"])


#Systematic weights
flow.VariationWeightArray("LHEScaleWeight",9,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic
#this is not obvious as N replicas can change... think about it
#flow.AddVariationWeightArray("LHEPdfWeight",30,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic


#create btag systematics
#this should be simplified
flow.Define("SelectedJet_btagWeight_up","vector_map(btagWeightUp,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
#flow.Define("btagEventWeightUp","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
flow.Systematic("BTagUp","SelectedJet_btagWeight","SelectedJet_btagWeight_up")

flow.createVariationBranch("BTagUp",["btagEventWeight"])
flow.VariationWeight("btagEventWeight__syst__BTagUP","btagEventWeight")

#or x in  flow.validCols :
#   if x[:len("defaultWeight")]=="defaultWeight" :
#if x!="defaultWeight":
#	flow.AddVariationWeight(x,filt=lambda hname,wname : "__syst__" not in hname,nodefault=True)

#Define Systematic variations
flow.Define("Muon_pt_scaleUp","Muon_pt*1.01f") #this should be protected against systematic variations
flow.Define("Muon_pt_scaleDown","Muon_pt*0.97f")
flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

for i in range(5):
  flow.Define("Jet_pt_JEC%s"%i,"Jet_pt+%s/100.f"%i)
  flow.Systematic("JEC%s"%i,"Jet_pt","Jet_pt_JEC%s"%i) #name, target, replacement

for i in range(5):
   print >> sys.stderr, "JEC",i
   flow.createVariationBranch("JEC%s"%i,["SBClassifier"])
  
flow.createVariationBranch("MuScaleUp",["SBClassifier"])
flow.createVariationBranch("MuScaleDown",["SBClassifier"])

print >> sys.stderr, "Number of known columns", len(flow.validCols)
#flow.printRDFCpp(["GenQQ_mass","QJet_indices","QJet0","QJet1","Rpt","SBClassifier","qqDeltaEta","MqqGenJet"])
for x in ["QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","HighestGenQQMass","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"] :
	flow.Histo(x)



flow.printRDFCpp(["SideBand","nSoftActivityJet","SoftActivityJet_pt","SoftActivityJet_eta","SoftActivityJet_phi","SoftActivityJet_SelectedJetDr","SoftActivityJet_SelectedJetIdx","SoftActivityJet_SelectedMuonDr","SoftActivityJet_SelectedMuonIdx","VBFRegion","QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","HighestGenQQMass","QJet0","QJet1","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"],debug=False,outname=sys.argv[1],selections=["VBFRegion","SideBand"])


#flow.printRDF(list(flow.allNodesTo("SBClassifier")))


