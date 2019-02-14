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
flow.SubCollection("QParton","GenPart",sel="GenPart_genPartIdxMother==Take(GenPart_genPartIdxMother,GenHiggs_idx)[0] && GenPart_pdgId!= 25")
flow.Define("QParton_p4","@p4v(QParton)")
flow.Distinct("QQ","QParton")
flow.Define("QQ_p4","QQ0_p4+QQ1_p4")
flow.Define("GenQQ_mass","Map(QQ_p4,[] (const auto & x){RETURN x.M();})")
flow.Define("HighestGenQQMass","QQ_mass[Argmax(QQ_mass)]")

flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("SelectedMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Distinct("MuMu","Muon")
flow.Define("OppositeSignMuMu","MuMu0_charge != MuMu1_charge")
flow.Selection("twoOppositeSignMuons","Sum(OppositeSignMuMu) > 0")
flow.Define("Mu0","OppositeSignMuMu[0][0]", requires=["twoOppositeSignMuons"])
flow.Define("Mu1","OppositeSignMuMu[1][0]", requires=["twoOppositeSignMuons"])
flow.Define("SelectedMuon_p4","@p4v(SelectedMuon)")
flow.Define("Higgs","@p4(SelectedMuon)[Mu0]+@p4(SelectedMuon)[Mu1]")

#VBF Jets kinematics
flow.DefaultConfig(jetPtCut=25)
flow.SubCollection("SelectedJet","Jet","Jet_pt > jetPtCut && (Jet_muonIdx1 == -1 || Take(Muon_iso,Jet_muonIdx1) > muIsoCut || Take(Muon_id,Jet_muonIdx1) > 0)")
flow.Selection("twoJets","nSelectedJet>=2")

flow.Define("GenJet_p4","@p4v(GenJet)")
flow.Define("SelectedJet_p4","@p4v(SelectedJet)")
#flow.Define("GenRecoJets","Combinations(SelectedJet,GenJet_pt)")
#flow.Match("SelectedJet","GenJet","ROOT:Math::DeltaR(SelectedJet_p4,GenJet_p4)",threshold=0.4,needwrapper=True,unique=False,bidirectional=True)
#flow.Define("GenRecoJets_dr","vector_map(ROOT:Math::DeltaR,SelectedJet_p4[GenRecoJets[0]],GenJet_p4[GenRecoJets[1]])")

flow.Distinct("JetPair","SelectedJet")
flow.Define("JetPair_p4","JetPair0_p4+JetPair1_p4",requires=["twoJets"])
flow.Define("JetPair_mass","Map(JetPair_p4,[](auto x){RETURN x.M()'})")
flow.Define("HighestPairMass","Argmax(JetPair_mass)")
flow.Define("HighestPairMass1","Argmax(Map( (JetPair0_p4+JetPair1_p4) ,[](auto x){RETURN x.M()'}))")
flow.Define("Jet0","JetPair[0][HighestPairMass]")
flow.Define("Jet1","JetPair[1][HighestPairMass]")
flow.Define("Jet0_p4","@p4(SelectedJet)[Jet0]")
flow.Define("Jet1_p4","@p4(SelectedJet)[Jet1]")


#flow.DefineObject("QJet0","SelectedJet","Jet0")
#flow.DefineObject("QJet1","SelectedJet","Jet1")
flow.Define("Qjet1","@p4(SelectedJet)[Jet0]")
flow.Define("Qjet2","@p4(SelectedJet)[Jet1]")
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
flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Selection("SideBand","! MassWindow")
flow.Selection("VBFRegion","Mqq > mQQcut")
flow.Selection("SignalRegion","VBFRegion && MassWindow")

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

print '''
#include <Math/VectorUtil.h>

template <typename type>
auto Argmax(const v & v){
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

#define RETURN return 
''' 
flow.printRDF(["GenQQ_mass","HighestPairMass1","HighestPairMass"])

#flow.printRDF(["Higgs_m","SBClassifier","SBClassifier__syst__MuScaleUp"])


