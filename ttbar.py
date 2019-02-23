from nail import *
import ROOT
import sys


flow=SampleProcessing("ttbar","/dev/shm/VBF_HToMuMu_nano2016.root")
# Toy Analysis for dileptonic and semileptonic ttbar
# * Loose Lepton selection (requires pt>20, an “Id Loose” flag, relative isolation < 0.25)
#    * Both electrons and muons
#Muons
flow.DefaultConfig(muIsoCut=0.13,muIdCut=3,muPtCut=25) 
flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("LooseMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Define("LooseMuon_p4","@p4v(LooseMuon)")

#Electrons
flow.DefaultConfig(eIsoCut=0.13,eIdCut=3,ePtCut=25) 
flow.Define("Electron_id","Electron_tightId*3+Electron_mediumId") 
flow.Define("Electron_iso","Electron_miniPFRelIso_all")
flow.SubCollection("LooseEle","Electron",sel="Electron_iso < eIsoCut && Electron_id > eIdCut && Electron_pt > ePtCut") 
flow.Define("LooseEle_p4","@p4v(LooseEle)")

#Lepton collection
flow.Define("Electron_p4","@p4v(Electron)")
flow.Define("SelectedMuon_pid","(SelectedMuon_pt*0)+13")
flow.Define("Electron_pid","(Electron_pt*0)+11")
flow.MergeCollections("Lepton",["SelectedMuon","Electron"])

# * Jet select/cleaning against loose leptons , jet pt > 25 , jet id
flow.DefaultConfig(jetPtCut=25,jetIdCut=0,jetPUIdCut=0)
flow.SubCollection("CleanJet","Jet",'''
 Jet_pt > jetPtCut &&
 Jet_jetId > jetIdCut &&
 Jet_puId > jetPUIdCut &&
 (Jet_muonIdx1 == -1 || Take(Muon_iso,Jet_muonIdx1) > muIsoCut || Take(Muon_id,Jet_muonIdx1) > muIdCut)
 (Jet_muonIdx2 == -1 || Take(Muon_iso,Jet_muonIdx2) > muIsoCut || Take(Muon_id,Jet_muonIdx2) > muIdCut)
 (Jet_electronIdx1 == -1 || Take(Electron_iso,Jet_electronIdx1) > eIsoCut || Take(Electron_id,Jet_electronIdx1) > eIdCut)
 (Jet_electronIdx2 == -1 || Take(Electron_iso,Jet_electronIdx2) > eIsoCut || Take(Electron_id,Jet_electronIdx2) > eIdCut)
''')
flow.Define("CleanedJet_btag","Jet_btagCSVV2") #this allows to later just use "btag" and to study changes of btag algorithm
flow.Define("CleanJet_p4","@p4v(CleanJet)")
flow.Define("CleanJet_pt","MemberMap(CleanJet_p4,pt())")
# * Compute event variables using selected cleaned jets
flow.Define("JetSum_p4","std::accumulate(CleanJet_p4.begin(), CleanJet_p4.end(), ROOT::Math::PtEtaPhiMVector())") #Feature request=> Sum with stardvalue
flow.Define("LeptonSum_p4","std::accumulate(Lepton_p4.begin(), Lepton_p4.end(), ROOT::Math::PtEtaPhiMVector())") #Feature request=> Sum with stardvalue
#    * HT, HT+loose leptons, MHT
flow.Define("HT","Sum(CleanJet_pt)")
flow.Define("MHTL","-(JetSum_p4+LeptonSum_p4).pt()")
flow.Define("MHT","-JetSum_p4.pt()")

# * Find opposite sign, same flavour pairs
flow.Selection("twoLeptons","nLepton>=2")
flow.Distinct("LPair","Lepton")
flow.Define("isOSSF","LPair0_charge != LPair1_charge && LPair0_pid == LPair1_pid",requires=["twoLeptons"])
flow.Selection("hasOSSF","Sum(isOSSF) > 0")
#    * Closest to Z
flow.TakePair("ZLep","Lepton","LPair","Argmax(-abs(MemberMap((LPair0_p4+LPair1_p4),M() )-91.2)*isOSSF)",requires=["hasOSSF"])
flow.Define("Z","ZLep0_p4+ZLep1_p4")
flow.Define("Z_mass","Z.M()")
#    * Highest pt pair
#    * Highest mass pair

# * Find opposite sign, opposite flavour pairs
flow.Define("isOSOF","LPair0_charge != LPair1_charge && LPair0_pid != LPair1_pid",requires=["twoLeptons"])
flow.Selection("hasOSOF","Sum(isOSOF)")
#    * Highest pt pair
#    * Highest mass pair

# * Find jj on W mass
flow.Distinct("JJ","CleanedJet")
flow.Selection("oneJets","nCleanedJet>=1")
flow.Selection("twoJets","nCleanedJet>=2")
#    * Closest to W mass
flow.TakePair("JJBestWMass","CleanedJet","JJ","Argmax(-abs(MemberMap((JJPair0_p4+JJPair1_p4),M() )-80.4))",requires=["twoJets"])
#    * Minimizing (Mjj-Mw)**2/(20)**2 + j1_btag+j2_btag (proxy for an actual likelihood) 
flow.TakePair("JJBestLike","CleanedJet","JJ","Argmax(-abs(MemberMap((JJPair0_p4+JJPair1_p4),M() )-80.4)/20.-JJPair0_btag-JJPair1_btag)",requires=["twoJets"])

# * Define B tagged jets
#    * Tight tagged
flow.DefaultConfig(tightOp=0.9,mediumOp=0.5)
flow.SubCollection("TightBTagged","CleanedJet","CleanedJet_btag > tightOp")
#    * Medium tagged
flow.SubCollection("MediumBTagged","CleanedJet","CleanedJet_btag > mediumOp")
#    * Leading, subleading btagged jet
flow.Define("BTagRank","Argsort(-CleanedJet_btag)")
flow.ObjectAt("LeadingBJet","CleanedJet","CleanedJet_btag[0]",requires=["oneJet"])
flow.ObjectAt("SubLeadingBJet","CleanedJet","CleanedJet_btag[1]",requires=["twoJets"])
# * All triplets jjB 
#    * Minimizing (Mjj-Mw)**2/(20)**2 + j1_btag+j2_btag -j3_btag
#    * Minimizing mass diff to top nominal

# * Tight leptons for trigger eff and signal/control region definition
#    * pt>30 (fiducial region)
#    * Id tight, rel iso < 0.15
flow.DefaultConfig(tightPtCut=30, tightIdCut=2,tightIso=0.15)
flow.SubCollection("TightLepton","Lepton","Lepton_iso < tightIso && Lepton_id > tightIdCut && Lepton_pt > tightPtCut")
# * Trigger eff weight
#    * From leading lepton pt (no trigger obj matching)
flow.ObjectAt("LeadingLepton","TightLepton","-Argsort(TightLepton_pt)")

#TODO: how to handle different weights in different regions?
flow.AddDefaultWeight("SingleLeptonEfficiency","efficiency(LeadingLepton_pt,LeadingLepton_eta,LeadingLepton_pid)")


# * Regions to be defined:
#    * Signal Dile: 2 loose leptons, opposite charge
flow.Define("isOS","LPair0_charge != LPair1_charge",requires=["twoLeptons"])
flow.Selection("hasOS","Sum(isOS)")
flow.Selection("DileptonRegion","hasOS") #this is noop
#    * Signal Semi: 1 tight lepton, no second loose lepton
flow.Selection("SemileptonRegion","nTightLepton==1 && nLepton == 1"
#    * DY control: 2 tight, closest to Z within 10 GeV of nominal Z 
flow.DefaultConfig(ZmassWindow=10)
flow.Selection("DYControl","nTightLepton==2 && abs(Z_mass-91.2) < ZmassWindow")
#    * W+jet: 1 tight, no second loose, leading btag < medium OP
flow.Selection("DYControl","nTightLepton==1 && nLepton ==1 & nMediumBTagged==0")
#    * Some other top selection (e.g. good candidate for mass, loose/tight lept)
#       * ???

# * Systematics:
#    * 50 JEC variations: no need for actual variations let’s scan from -25% to +25% with the 50 variations (steps of 1%)
#    * JER: use TRand smearing
#    * Mu Scale: up/down by 3%
#    * Ele Scale: up/down by 3%
#    * Btag weights x 10:
#       * Compute eventBtagWeights for 10 different values of the individual jetweights (from 0.9 to 1.1 in steps of 0.02)
#    * Lept eff weight
#       * ???
# * Histo to event weight
#    * Trigger eff
#    * Lepton eff
# * Weight arrays
#    * 100 PDF
#    * 8 qcdscale
# * Good runs/lumi selection
# * Long “OR” Trigger bit selection (expanding from wildcards?)
# * Output results (in multiple regions, with systematics):
#    * Cutflow: all selections individually (if allowed by dependencies), sorted as in declaration down to signal region
#    * Pt,eta for jets,leptons (tight and loose),W,Top cands(triples)
#    * MHT, HT, MET, rho, tkmet, nPVs
#    * Top and W mass in regions where it makes sense
#    * Z(dilepton) mass in regions where it makes sense



# Define B tagged jets
# Tight tagged
# Medium tagged
# Leading, subleading btagged jet
# All triplets jjB 
# Minimizing (Mjj-Mw)**2/(20)**2 + j1_btag+j2_btag -j3_btag
# Minimizing mass diff to top nominal
# Tight leptons for trigger eff and signal/control region definition
# pt>30 (fiducial region)
# Id tight, rel iso < 0.15
# Trigger eff weight
# From leading lepton pt (no trigger obj matching)
# Regions to be defined:
# Signal Dile: 2 loose leptons, opposite charge
# Signal Semi: 1 tight lepton, no second loose lepton, opposite charge
# DY control: 2 tight, closest to Z within 10 GeV of nominal Z 
# W+jet: 1 tight, no second loose, leading btag < medium OP
# Some other top selection (e.g. good candidate for mass, loose/tight lept)
# ???
# Systematics:
# 50 JEC variations: no need for actual variations let’s scan from -25% to +25% with the 50 variations (steps of 1%)
# JER: use TRand smearing
# Mu Scale: up/down by 3%
# Ele Scale: up/down by 3%
# Btag weights x 10:
# Compute eventBtagWeights for 10 different values of the individual jetweights (from 0.9 to 1.1 in steps of 0.02)
# Lept eff weight
# ???
# Histo to event weight
# Trigger eff
# Lepton eff
# Weight arrays
# 100 PDF
# 8 qcdscale
# Good runs/lumi selection
# Long “OR” Trigger bit selection (expanding from wildcards?)
# Output results (in multiple regions, with systematics):
# Cutflow: all selections individually (if allowed by dependencies), sorted as in declaration down to signal region
# Pt,eta for jets,leptons (tight and loose),W,Top cands(triples)
# MHT, HT, MET, rho, tkmet, nPVs
# Top and W mass in regions where it makes sense
# Z(dilepton) mass in regions where it makes sense
# N btags, N leptons


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

flow.DefaultConfig(higgsMassWindowWidth=15,mQQcut=400,nominalHMass=125.03)
flow.Selection("MassWindow","abs(Higgs.M()-nominalHMass)<higgsMassWindowWidth")
flow.Selection("SideBand","! MassWindow")
flow.Selection("VBFRegion","Mqq > mQQcut")
flow.Selection("SignalRegion","VBFRegion && MassWindow")

#flow.Trainable("SBClassifier","evalMVA",["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"],splitMode="TripleMVA",requires="VBFRegion") 
flow.Define("Higgs_pt","Higgs.Pt()")
flow.Define("Higgs_m","Higgs.M()")
flow.Define("SBClassifier","Higgs_pt+Higgs_m+Mqq+Rpt+DeltaRelQQ") #,inputs=["Higgs_pt","Higgs_m","Mqq","Rpt","DeltaRelQQ"])

#define some event weights
flow.Define("SelectedJet_btagWeight","vector_map(btagWeight,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
flow.Define("btagEventWeight","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
flow.AddDefaultWeight("genWeight")
flow.AddDefaultWeight("btagEventWeight")


#Systematic weights
flow.AddWeightArray("LHEScaleWeight",9,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic
#this is not obvious as N replicas can change... think about it
#flow.AddWeightArray("LHEPdfWeight",30,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic


#create btag systematics
#this should be simplified
flow.Define("SelectedJet_btagWeight_up","vector_map(btagWeightUp,SelectedJet_btagCSVV2,SelectedJet_pt,SelectedJet_eta)")
flow.Define("btagEventWeightUp","std::accumulate(SelectedJet_btagWeight.begin(),SelectedJet_btagWeight.end(),1, std::multiplies<double>())")
flow.Systematic("BTagUp","SelectedJet_btagWeight","SelectedJet_btagWeight_up")
flow.createVariationBranch("BTagUp","defaultWeight")
for x in  flow.validCols :
    if x[:len("defaultWeight")]=="defaultWeight" :
	if x!="defaultWeight":
		flow.AddWeight(x,filt=lambda hname,wname : "__syst__" not in hname,nodefault=True)

#Define Systematic variations
flow.Define("Muon_pt_scaleUp","Muon_pt*1.01") #this should be protected against systematic variations
flow.Define("Muon_pt_scaleDown","Muon_pt*0.97")
flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

for i in range(10):
  flow.Define("Jet_pt_JEC%s"%i,"Jet_pt+%s/100."%i)
  flow.Systematic("JEC%s"%i,"Jet_pt","Jet_pt_JEC%s"%i) #name, target, replacement

for i in range(10):
   print >> sys.stderr, "JEC",i
   flow.createVariationBranch("JEC%s"%i,"SBClassifier")
  
flow.createVariationBranch("MuScaleUp","SBClassifier")
flow.createVariationBranch("MuScaleDown","SBClassifier")


print '''
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>

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

float btagWeight(float csv,float pt,float eta){
 return 1.01;
}
float btagWeightUp(float csv,float pt,float eta){
 return 1.03;
}

float efficiency(float pt,float eta,int pid){
 if(pid==11) return 0.99;
 if(pid==13) return 0.92;
}

#define MemberMap(vector,member) Map(vector,[](auto x){return x.member;})
#define RETURN return 
''' 
print >> sys.stderr, "Number of known columns", len(flow.validCols)
#flow.printRDFCpp(["GenQQ_mass","QJet_indices","QJet0","QJet1","Rpt","SBClassifier","qqDeltaEta","MqqGenJet"])
for x in ["QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","HighestGenQQMass","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"] :
	flow.Histo(x)



flow.printRDFCpp(["defaultWeight","QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","HighestGenQQMass","QJet0","QJet1","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"]+flow.weights.keys(),debug=False)


#flow.printRDF(list(flow.allNodesTo("SBClassifier")))


