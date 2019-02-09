from nail import *
import ROOT
f=ROOT.TFile.Open("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root")
e=f.Get("Events")
allbranches=[(x.GetName(),x.GetListOfLeaves()[0].GetTypeName()) for x in e.GetListOfBranches()]

flow=SampleProcessing("",allbranches)
#flow=SampleProcessing("",["Muon_pt","Muon_eta","Muon_phi","Muon_tightId","Muon_looseId","Jet_pt","Muon_iso","Jet_muonIdx1","Jet_eta","Jet_phi","Jet_mass"])
print "Start"
flow.DefaultConfig(muIsoCut=0.13,muIdCut=3,muPtCut=25) #cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible

#Higgs to mumu reconstruction
flow.Define("Muon_mass","0.106+0*Muon_pt") #ensure same lenght of Muon_pt
flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Filter("twoOppositeSignMuons","nSelectedMuon>=2 && SelectedMuon_charge[0]*SelectedMuon_charge[1] < 0")
#flow.Define("SelectedMuon_p4","@p4(SelectedMuon)")
flow.Define("Higgs","@p4(SelectedMuon)[0]+@p4(SelectedMuon)[1]",requires=["twoOppositeSignMuons"]) 
#flow.Combinations("JetPairs","Jets","Jets",filter="[](auto a,auto b){a != b}")
#flow.Reduce("HighestPtPair","MaxArg(first(p4)+second(p4))")

#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet","Jet_pt > jetPtCut && (Jet_muonIdx1 == -1 || Take(Muon_iso,Jet_muonIdx1) > muIsoCut || Take(Muon_id,Jet_muonIdx1) > 0)")
flow.Filter("twoJets","nSelectedJet>=2")
flow.Define("Qjet1","@p4(SelectedJet)[0]",requires=["twoJets"])
flow.Define("Qjet2","@p4(SelectedJet)[1]",requires=["twoJets"])
flow.Define("qq","Qjet1+Qjet2")
flow.Define("Mqq","qq.M()")
flow.Define("qq_pt","qq.Pt()")
flow.Define("qqDeltaEta","TMath::Abs(Qjet1.Eta()-Qjet2.Eta())")
flow.Define("qqDeltaPhi","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(Qjet1,Qjet2))")

#QQ vs ll kinematic
flow.Define("ll_ystar","Higgs.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity())")
flow.Define("ll_zstar"," TMath::Abs( ll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ")
flow.Define("DeltaEtaQQSum","TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta())")
flow.Define("PhiZQ1","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,Qjet1))")
flow.Define("PhiZQ2","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,Qjet2))")
flow.Define("EtaHQ1","TMath::Abs(Higgs.Eta() - Qjet1.Eta())")
flow.Define("EtaHQ2","TMath::Abs(Higgs.Eta() - Qjet2.Eta())")
flow.Define("DeltaRelQQ","(Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt())")
flow.Define("Rpt","(Qjet1+Qjet2+ Higgs).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Higgs.Pt())")

flow.DefaultConfig(higgsMassWindowWidth=15,mQQcut=400,nominalHMass=125.03)
flow.Filter("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Filter("SideBand","! MassWindow")
flow.Filter("VBFRegion","Mqq > mQQcut")
flow.Filter("SignalRegion","VBFRegion && MassWindow")

#flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion") 
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("SBClassifier","Higgs_pt+Higgs_m+Mqq+Rpt+DeltaRelQQ",inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])


#Define Systematic variations
flow.Define("Muon_pt_scaleUp","Muon_pt*1.01")
flow.Define("Muon_pt_scaleDown","Muon_pt*0.97")
flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

flow.createSystematicBranch("MuScaleUp","SBClassifier")
flow.createSystematicBranch("MuScaleDown","SBClassifier")

#flow.printRDF(["Higgs_m","SBClassifier"])
flow.printRDF(["Higgs_m","SBClassifier","SBClassifier__syst__MuScaleUp"])

exit(1)
import os
class Analaysis:
    def __init__(self,flow,samples,store="analysis_store"):
	self.flow=flow
	self.samples=samples
	self.store=store
	os.mkdir(store)
	self.workers={}
	
    def clearActions():
        pass    

    def estimateTime():
	pass

    def run():
	pass

    def Draw(self,observables,samples,sampleref,normalization="area",systematics=[]) :
	pass

    def ComparisonChi2(self,observables,samples,sampleref,normalization="area",systematics=[]) :
	flow.fetchHistos(observables,samples+[sampleref],normalization,systematics)

    def fetchHistos(self,observables,samples,systematics=[]):	
	for s in samples :
	    w=self.getWorker(s)


from samples import *
ana=Analysis(flow,samples)


ana.ComparisonChi2(flow.inputs["SBClassifier"],["b","s"],"d" ,norm="area")
ana.ComparisonChi2(["SBClassifier"],["s+b","b"],"d", norm="xsec",systematics=[".*"])

#interactive?
#ana=Analysis(flow,samples)
#ana.stateLoad()
#ana.



#print "list of vars to update with MuScaleUp"
if 0:
  import networkx as nx
  nodes=set(flow.findAffectedNodesForSystematicOnTarget("MuScaleUp","SBClassifier")+["Muon_pt"])
  G = nx.DiGraph()
  for k in flow.inputs :
   if k in nodes:
    G.add_node(k)
  for k in flow.inputs :
    for i in flow.inputs[k] :
      if k in nodes and i in nodes :
       G.add_edge(k,i)

  import matplotlib.pyplot as plt
  plt.subplot(121)
  nx.draw_networkx(G, with_labels=True, font_weight='bold')
  plt.show()


#flow.printRDF()


import networkx as nx
if 0:
  print "digraph {"
  for k in flow.inputs :
    for i in flow.inputs[k] :
       print k,"->",i,";"
  print "}"

if 0 :
  G = nx.DiGraph()
  for k in flow.inputs :
    G.add_node(k)
  for k in flow.inputs :
    for i in flow.inputs[k] :
       G.add_edge(k,i)
 
  import matplotlib.pyplot as plt
  plt.subplot(121)
  nx.draw_networkx(G, with_labels=True, font_weight='bold')
  plt.show()

