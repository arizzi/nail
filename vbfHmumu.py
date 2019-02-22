from nail import *
import ROOT
import sys

#f=ROOT.TFile.Open("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root")
f=ROOT.TFile.Open("/dev/shm/VBF_HToMuMu_nano2016.root")
e=f.Get("Events")
allbranches=[(x.GetName(),x.GetListOfLeaves()[0].GetTypeName()) for x in e.GetListOfBranches()]
df=ROOT.RDataFrame("Events","/dev/shm/VBF_HToMuMu_nano2016.root")
dftypes={x[0]:df.GetColumnType(x[0]) for x in allbranches}


print >>sys.stderr, "ROOT loaded"
flow=SampleProcessing("",allbranches,dftypes)
#flow=SampleProcessing("",["Muon_pt","Muon_eta","Muon_phi","Muon_tightId","Muon_looseId","Jet_pt","Muon_iso","Jet_muonIdx1","Jet_eta","Jet_phi","Jet_mass"])
flow.DefaultConfig(muIsoCut=0.13,muIdCut=3,muPtCut=25) #cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible

#Higgs to mumu reconstruction
flow.Define("GenHiggs_idx","Nonzero(GenPart_pdgId == 25)")
flow.SubCollection("QParton","GenPart",sel="GenPart_genPartIdxMother==Take(GenPart_genPartIdxMother,GenHiggs_idx)[0] && GenPart_pdgId!= 25")
flow.Define("QParton_p4","@p4v(QParton)")
flow.Distinct("QQ","QParton")
flow.Selection("twoQ","nQParton>=2")
flow.Define("QQ_p4","QQ0_p4+QQ1_p4",requires=["twoQ"])
flow.Define("GenQQ_mass","MemberMap(QQ_p4,M())")
flow.Define("HighestGenQQMass","QQ_mass[Argmax(QQ_mass)]")

flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Define("SelectedMuon_p4","@p4v(SelectedMuon)")
flow.Selection("twoMuons","nSelectedMuon>=2")
flow.Distinct("MuMu","SelectedMuon")
flow.Define("OppositeSignMuMu","Nonzero(MuMu0_charge != MuMu1_charge)",requires=["twoMuons"])
flow.Selection("twoOppositeSignMuons","OppositeSignMuMu.size() > 0")
flow.TakePair("Mu","SelectedMuon","MuMu","OppositeSignMuMu[0]",requires=["twoOppositeSignMuons"])
flow.Define("Higgs","Mu0_p4+Mu1_p4")

#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet","Jet_pt > jetPtCut && (Jet_muonIdx1 == -1 || Take(Muon_iso,Jet_muonIdx1) > muIsoCut || Take(Muon_id,Jet_muonIdx1) > 0)")
flow.Selection("twoJets","nSelectedJet>=2")
flow.Define("GenJet_p4","@p4v(GenJet)")
flow.Define("SelectedJet_p4","@p4v(SelectedJet)")
flow.Distinct("JetPair","SelectedJet")
flow.TakePair("QJet","SelectedJet","JetPair","Argmax(MemberMap((JetPair0_p4+JetPair1_p4),M() ))",requires=["twoJets"])
#flow.ObjectAt("QJet0","SelectedJet","0",requires=["twoJets"])
#flow.ObjectAt("QJet1","SelectedJet","1",requires=["twoJets"])

flow.Define("qq","QJet0_p4+QJet1_p4")
flow.Define("Mqq","qq.M()")
flow.Define("MqqGenJet","(QJet0_genJetIdx>=0&&QJet1_genJetIdx>=0&&QJet0_genJetIdx<nGenJet&&QJet1_genJetIdx<nGenJet)?(GenJet_p4[QJet0_genJetIdx]+GenJet_p4[QJet1_genJetIdx]).M():-99")
flow.Define("qq_pt","qq.Pt()")
flow.Define("qqDeltaEta","abs(QJet0_eta-QJet1_eta)")
flow.Define("qqDeltaPhi","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(QJet0_p4,QJet1_p4))")

#QQ vs ll kinematic
flow.Define("ll_ystar","Higgs.Rapidity() - (QJet0_p4.Rapidity() + QJet1_p4.Rapidity())")
flow.Define("ll_zstar"," TMath::Abs( ll_ystar/ (QJet0_p4.Rapidity()-QJet1_p4.Rapidity() )) ")
flow.Define("DeltaEtaQQSum","TMath::Abs(QJet0_p4.Eta()) +  TMath::Abs(QJet1_p4.Eta())")
flow.Define("PhiZQ1","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet0_p4))")
flow.Define("PhiZQ2","TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(Higgs,QJet1_p4))")
flow.Define("EtaHQ1","TMath::Abs(Higgs.Eta() - QJet0_p4.Eta())")
flow.Define("EtaHQ2","TMath::Abs(Higgs.Eta() - QJet1_p4.Eta())")
flow.Define("DeltaRelQQ","(QJet0_p4+QJet1_p4).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt())")
flow.Define("Rpt","(QJet0_p4+QJet1_p4+ Higgs).Pt()/( QJet0_p4.Pt()+QJet1_p4.Pt() + Higgs.Pt())")

flow.DefaultConfig(higgsMassWindowWidth=15,mQQcut=400,nominalHMass=125.03)
flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Selection("SideBand","! MassWindow")
flow.Selection("VBFRegion","Mqq > mQQcut")
flow.Selection("SignalRegion","VBFRegion && MassWindow")

#flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion") 
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("SBClassifier","Higgs_pt+Higgs_m+Mqq+Rpt+DeltaRelQQ") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])


#Define Systematic variations
flow.Define("Muon_pt_scaleUp","Muon_pt*1.01") #this should be protected against systematic variations
flow.Define("Muon_pt_scaleDown","Muon_pt*0.97")
flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

#for i in range(1):
#  flow.Define("Jet_pt_JEC%s"%i,"Jet_pt+%s/100."%i)
#  flow.Systematic("JEC%s"%i,"Jet_pt","Jet_pt_JEC%s"%i) #name, target, replacement

#flow.createSystematicBranch("MuScaleUp","SBClassifier")
#for i in range(1):
#   print >> sys.stderr, "JEC",i
#   flow.createSystematicBranch("JEC%s"%i,"SBClassifier")
  
flow.createSystematicBranch("MuScaleDown","SBClassifier")

print '''
#include <Math/VectorUtil.h>

template <typename type>
auto Argmax(const type & v){
 return ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(v))[0];
}

template <typename type, typename Vec,typename... OtherVecs>
auto vector_map_t(const Vec & v,  const OtherVecs &... args) {
  ROOT::VecOps::RVec<type> res(v.size());
  for(size_t i=0;i<v.size(); i++) res[i]=type(v[i],args[i]...);
  return res;
}

template <typename func, typename Vec,typename... OtherVecs>
auto vector_map(func f, const Vec & v,  const OtherVecs &... args) {
  ROOT::VecOps::RVec<decltype(f(std::declval<typename Vec::value_type>(),std::declval<typename OtherVecs::value_type>()...))> res(v.size());
  for(size_t i=0;i<v.size(); i++) res[i]=f(v[i],args[i]...);
  return res;
}


/*template <typename func, typename Vec,typename... OtherVecs>
auto matrix_map(shape, int axis,func f, const Vec & v,  const OtherVecs &... args) {

}*/

auto pt(const ROOT::Math::PtEtaPhiMVector &i){
 return i.pt();
}

auto mass(const ROOT::Math::PtEtaPhiMVector &i){
 return i.M();
}

#define MemberMap(vector,member) Map(vector,[](auto x){return x.member;})
#define RETURN return 
''' 
print >> sys.stderr, "Number of known columns", len(flow.validCols)
#flow.printRDFCpp(["GenQQ_mass","QJet_indices","QJet0","QJet1","Rpt","SBClassifier","qqDeltaEta","MqqGenJet"])
flow.printRDFCpp(["GenQQ_mass","QJet0","QJet1","Rpt","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"])


#flow.printRDF(list(flow.allNodesTo("SBClassifier")))


