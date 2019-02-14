from nail import *
import ROOT

f=ROOT.TFile.Open("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root")
e=f.Get("Events")
allbranches=[(x.GetName(),x.GetListOfLeaves()[0].GetTypeName()) for x in e.GetListOfBranches()]

print "ROOT loaded"
flow=SampleProcessing("",allbranches)
#flow=SampleProcessing("",["Muon_pt","Muon_eta","Muon_phi","Muon_tightId","Muon_looseId","Jet_pt","Muon_iso","Jet_muonIdx1","Jet_eta","Jet_phi","Jet_mass"])
flow.DefaultConfig(muIsoCut=0.13,muIdCut=3,muPtCut=25) #cuts value should not be hardcoded below but rather being declared here so that scans and optimizations are possible

#Higgs to mumu reconstruction
flow.Define("GenHiggs_idx","Nonzero(GenPart_pdgId == 25)")
flow.SubCollection("QQParton","GenPart",sel="GenPart_genPartIdxMother==Take(GenPart_genPartIdxMother,GenHiggs_idx)[0] && GenPart_pdgId!= 25")
#flow.Define("QQParton_p4","vector_map_t<ROOT::Math::PtEtaPhiMVector>(QQParton_pt, QQParton_eta, QQParton_phi,QQParton_mass)")
flow.Define("QQParton_p4","@p4v(QQParton)")

flow.Distinct("QQ","QQParton")
flow.Define("QQ_p4","QQ0_p4+QQ1_p4")
flow.Define("QQ_mass","Map(QQ_p4,mass)")
flow.Define("HighestQQMass","QQ_mass[Reverse(Argsort(QQ_mass))[0]]")


#flow.Define("QQPair_p4","QQPair0_p4+QQPair1_p4")
#flow.Define("SelectedQQPair","Reverse(Argsort(QQPair_pt))[0]"
#flow.Define("QQPairAll","Combinations(Nonzero(QQParton),Nonzero(QQParton))")
#flow.Define("QQPair","Nonzero(QQPairAll[0] > QQPairAll[1])")
#flow.Define("QQPair0","Take(QQPairAll[0],QQPair)")
#flow.Define("QQPair1","Take(QQPairAll[1],QQPair)")
#flow.Define("QQPair1_pt","Take(QQParton_pt,QQPair0)")
#flow.Define("QQPair2_pt","Take(QQParton_pt,QQPair1)")
#flow.Define("QQPair_p4","(Take(QQParton_p4,QQPair0)+Take(QQParton_p4,QQPair1))")
#flow.Define("QQPair_mass","Map(QQPair_p4,mass)")
#flow.Define("QQPair_pt","Map(QQPair_p4,pt)")
#flow.Define("QQPairHighestPt_mass","QQPair_mass[Reverse(Argsort(QQPair_pt))[0]]") 



flow.Define("Muon_mass","0.106+0*Muon_pt") #ensure same lenght of Muon_pt
flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Filter("twoOppositeSignMuons","nSelectedMuon>=2 && SelectedMuon_charge[0]*SelectedMuon_charge[1] < 0")
flow.Define("SelectedMuon_p4","@p4v(SelectedMuon)")
flow.Define("Higgs","@p4(SelectedMuon)[0]+@p4(SelectedMuon)[1]",requires=["twoOppositeSignMuons"]) 

#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet","Jet_pt > jetPtCut && (Jet_muonIdx1 == -1 || Take(Muon_iso,Jet_muonIdx1) > muIsoCut || Take(Muon_id,Jet_muonIdx1) > 0)")
flow.Filter("twoJets","nSelectedJet>=2")

flow.Define("GenRecoJets","Combinations(SelectedJet,GenJet_pt)")
flow.Define("GenJet_p4","vector_map_t<ROOT::Math::PtEtaPhiMVector>(GenJet_pt, GenJet_eta, GenJet_phi,GenJet_mass)")
flow.Define("SelectedJet_p4","vector_map_t<ROOT::Math::PtEtaPhiMVector>(SelectedJet_pt, SelectedJet_eta, SelectedJet_phi,SelectedJet_mass)")
flow.Define("GenRecoJets_dr","vector_map(ROOT:Math::DeltaR,SelectedJet_p4[GenRecoJets[0]],GenJet_p4[GenRecoJets[1]])")

#flow.Define("JetPairs","Combinations(SelectedJet,SelectedJet)")
#flow.SubCollection("SelectedJet0","SelectedJet","SelectedJet[JetPairs[0]]")
#flow.SubCollection("SelectedJet1","SelectedJet","SelectedJet[JetPairs[1]]")
#flow.Define("PairPt","(@p4(SelectedJet0)+@p4(SelectedJet1)).Pt()")
#flow.Define("SelectedPair","ArgSort(PairPt)[0]")
#flow.Reduce("HighestPtPair","MaxArg(first(p4)+second(p4))")
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
#flow.printRDF(["Higgs_m","SBClassifier","SBClassifier__syst__MuScaleUp","QQPartons_eta","nQQPartons"])
print '''
#include <Math/VectorUtil.h>



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

auto pt(const ROOT::Math::PtEtaPhiMVector &i){
 return i.pt();
}

auto mass(const ROOT::Math::PtEtaPhiMVector &i){
 return i.M();
}

 
''' 
flow.printRDF(["HighestQQMass","nQQParton","QQ_mass"])


