from nail import *
import ROOT
import sys


flow=SampleProcessing("ttbar","/dev/shm/VBF_HToMuMu_nano2016.root")
# Toy Analysis for dileptonic and semileptonic ttbar
# * Loose Lepton selection (requires pt>20, a  Loosflag, relative isolation < 0.25)
#    * Both electrons and muons
#Muons
flow.DefaultConfig(muIsoCut=0.13,muIdCut=3,muPtCut=25) 
flow.Define("Muon_id","Muon_tightId*3+Muon_mediumId") 
flow.Define("Muon_iso","Muon_miniPFRelIso_all")
flow.SubCollection("LooseMuon","Muon",sel="Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut") 
flow.Define("LooseMuon_p4","@p4v(LooseMuon)")

#Electrons
flow.Define("Electron_p4","@p4v(Electron)")
flow.DefaultConfig(eIsoCut=0.13,eIdCut=3,ePtCut=25) 
flow.Define("Electron_id","Electron_cutBased") 
flow.Define("Electron_iso","Electron_miniPFRelIso_all")
flow.SubCollection("LooseEle","Electron",sel="Electron_iso < eIsoCut && Electron_id > eIdCut && Electron_pt > ePtCut") 
flow.Define("LooseEle_p4","@p4v(LooseEle)")

#Lepton collection
flow.Define("LooseMuon_pid","(LooseMuon_pt*0)+13")
flow.Define("LooseEle_pid","(LooseEle_pt*0)+11")
flow.MergeCollections("Lepton",["LooseMuon","LooseEle"])

#Match jets to selected loose Leptons
flow.Define("Jet_p4","@p4v(Jet)")
flow.MatchDeltaR("Jet","Lepton",embed=(["pt","pid"],[]))

# * Jet select/cleaning against loose leptons , jet pt > 25 , jet id
flow.DefaultConfig(jetPtCut=25,jetIdCut=0,jetPUIdCut=0)
flow.SubCollection("CleanJet","Jet",'''
 Jet_pt > jetPtCut &&
 Jet_jetId > jetIdCut &&
 Jet_puId > jetPUIdCut &&
 (Jet_LeptonIdx==-1 || Jet_LeptonDr > 0.3)
''')
flow.Define("CleanJet_btag","Jet_btagCSVV2") #this allows to later just use "btag" and to study changes of btag algorithm
flow.Define("CleanJet_p4","@p4v(CleanJet)")



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
flow.TakePair("LepHPt","Lepton","LPair","Argmax((MemberMap((LPair0_p4+LPair1_p4),pt() ))*isOSSF)",requires=["hasOSSF"])
flow.Define("OSSFHPt","LepHpt0_p4+LepHPt1_p4")
flow.Define("OSSFHPt_mass","OSSFHPt.M()")

#    * Highest mass pair
flow.TakePair("LepHM","Lepton","LPair","Argmax(MemberMap((LPair0_p4+LPair1_p4),M() )*isOSSF)",requires=["hasOSSF"])
flow.Define("OSSFHM","LepHM0_p4+LepHM1_p4")
flow.Define("OSSFHM_mass","OSSFHM.M()")

# * Find opposite sign, opposite flavour pairs
flow.Define("isOSOF","LPair0_charge != LPair1_charge && LPair0_pid != LPair1_pid",requires=["twoLeptons"])
flow.Selection("hasOSOF","Sum(isOSOF)")
#    * Highest pt pair
flow.TakePair("LepOFHPt","Lepton","LPair","Argmax((MemberMap((LPair0_p4+LPair1_p4),pt() ))*isOSOF)",requires=["hasOSOF"])
flow.Define("OSOFHPt","LepOFHpt0_p4+LepOFHPt1_p4")
flow.Define("OSOFHPt_mass","OSOFHPt.M()")

#    * Highest mass pair
flow.TakePair("LepOFHM","Lepton","LPair","Argmax(MemberMap((LPair0_p4+LPair1_p4),M() )*isOSOF)",requires=["hasOSOF"])
flow.Define("OSOFHM","LepOFHM0_p4+LepOFHM1_p4")
flow.Define("OSOFHM_mass","OSOFHM.M()")


# * Find jj on W mass
flow.Distinct("JJ","CleanJet")
flow.Selection("oneJets","nCleanJet>=1")
flow.Selection("twoJets","nCleanJet>=2")
#    * Closest to W mass
flow.TakePair("JJBestWMass","CleanJet","JJ","Argmax(-abs(MemberMap((JJ0_p4+JJ1_p4),M() )-80.4))",requires=["twoJets"])
flow.Define("JJBestWMass_mass","(JJBestWMass0_p4+JJBestWMass1_p4).M()")
#    * Minimizing (Mjj-Mw)**2/(20)**2 + j1_btag+j2_btag (proxy for an actual likelihood) 
flow.TakePair("JJBestLike","CleanJet","JJ","Argmax(-(pow(MemberMap((JJ0_p4+JJ1_p4),M() )-80.4,2)/400.)-JJ0_btag-JJ1_btag)",requires=["twoJets"])
flow.Define("JJBestLike_mass","(JJBestLike0_p4+JJBestLike1_p4).M()")
# * Define B tagged jets
#    * Tight tagged
flow.DefaultConfig(tightOp=0.9,mediumOp=0.5)
flow.SubCollection("TightBTagged","CleanJet","CleanJet_btag > tightOp")
#    * Medium tagged
flow.SubCollection("MediumBTagged","CleanJet","CleanJet_btag > mediumOp")
#    * Leading, subleading btagged jet
flow.Define("BTagRank","Argsort(-CleanJet_btag)")
flow.ObjectAt("LeadingBJet","CleanJet","CleanJet_btag[0]",requires=["oneJet"])
flow.ObjectAt("SubLeadingBJet","CleanJet","CleanJet_btag[1]",requires=["twoJets"])
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
#TODO: defaultWeight now has a dependency!!! arghhhhhh
#flow.Define("SingleLeptonEfficiency","efficiency(LeadingLepton_pt,LeadingLepton_eta,LeadingLepton_pid)")
#flow.AddDefaultWeight("SingleLeptonEfficiency")


# * Regions to be defined:
#    * Signal Dile: 2 loose leptons, opposite charge
flow.Define("isOS","LPair0_charge != LPair1_charge",requires=["twoLeptons"])
flow.Selection("hasOS","Sum(isOS)")
flow.Selection("DileptonRegion","hasOS") #this is noop
#    * Signal Semi: 1 tight lepton, no second loose lepton
flow.Selection("SemileptonRegion","nTightLepton==1 && nLepton == 1")
#    * DY control: 2 tight, closest to Z within 10 GeV of nominal Z 
flow.DefaultConfig(ZmassWindow=10)
flow.Selection("DYControl","nTightLepton==2 && abs(Z_mass-91.2) < ZmassWindow")
#    * W+jet: 1 tight, no second loose, leading btag < medium OP
flow.Selection("WControl","nTightLepton==1 && nLepton ==1 & nMediumBTagged==0")
#    * Some other top selection (e.g. good candidate for mass, loose/tight lept)
#       * ???

# * Systematics:
#    * 50 JEC variations: no need for actual variations let's scan from -25% to +25% with the 50 variations (steps of 1%)
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
# * Long OR Trigger bit selection (expanding from wildcards?)


#define some event weights
flow.Define("CleanJet_btagWeight","vector_map(btagWeight,CleanJet_btag,CleanJet_pt,CleanJet_eta)")
flow.Define("btagEventWeight","std::accumulate(CleanJet_btagWeight.begin(),CleanJet_btagWeight.end(),1, std::multiplies<double>())")
flow.AddDefaultWeight("genWeight")
flow.AddDefaultWeight("btagEventWeight")


#Systematic weights
flow.AddWeightArray("LHEScaleWeight",9,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic
#this is not obvious as N replicas can change... think about it
flow.AddWeightArray("LHEPdfWeight",100,filt=lambda hname,wname : "__syst__" not in hname ) #systematic variations are 1D, let's avoid systematics of systematic


#create btag systematics
#this should be simplified
flow.Define("CleanJet_btagWeight_up","vector_map(btagWeightUp,CleanJet_btag,CleanJet_pt,CleanJet_eta)")
flow.Define("btagEventWeightUp","std::accumulate(CleanJet_btagWeight.begin(),CleanJet_btagWeight.end(),1, std::multiplies<double>())")
flow.Systematic("BTagUp","CleanJet_btagWeight","CleanJet_btagWeight_up")
flow.createVariationBranch("BTagUp",["defaultWeight"])

#takes all syst variations of defaultWeight
for x in  flow.validCols :
    if x[:len("defaultWeight")]=="defaultWeight" :
	if x!="defaultWeight":
		flow.AddWeight(x,filt=lambda hname,wname : "__syst__" not in hname,nodefault=True)

#Define Systematic variations
flow.Define("Muon_pt_scaleUp","Muon_pt*1.01") #this should be protected against systematic variations
flow.Define("Muon_pt_scaleDown","Muon_pt*0.97")
flow.Systematic("MuScaleDown","Muon_pt","Muon_pt_scaleDown") #name, target, replacement
flow.Systematic("MuScaleUp","Muon_pt","Muon_pt_scaleUp") #name, target, replacement

for i in range(50):
  flow.Define("Jet_pt_JEC%s"%i,"Jet_pt*(1.+(%s-24.5)/100.)"%i)
  flow.Systematic("JEC%s"%i,"Jet_pt","Jet_pt_JEC%s"%i) #name, target, replacement

# * Output results (in multiple regions, with systematics):
#    * Cutflow: all selections individually (if allowed by dependencies), sorted as in declaration down to signal region
#    * Pt,eta for jets,leptons (tight and loose),W,Top cands(triples)
#    * MHT, HT, MET, rho, tkmet, nPVs
#    * Top and W mass in regions where it makes sense
#    * Z(dilepton) mass in regions where it makes sense
colsToPlot=["nJet","nLepton","MHT","HT","Z_mass","OSOFHM_mass","OSSFHM_mass","JJBestLike_mass","JJBestWMass_mass"]
for attr in ["_pt","_eta"]:
  for coll in ["LooseMuon","TightLepton","LooseEle","CleanJet"]:
   colsToPlot.append(coll+attr)


for i in range(50):
   print >> sys.stderr, "JEC",i
   flow.createVariationBranch("JEC%s"%i,colsToPlot)
flow.createVariationBranch("MuScaleUp",colsToPlot)
flow.createVariationBranch("MuScaleDown",colsToPlot)



colsToPlotWithSyst=[]
for c in colsToPlot :
    colsToPlotWithSyst+=[x for x in flow.validCols if x[:len(c)]==c]


print '''
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
#include "helpers.h"
#define MemberMap(vector,member) Map(vector,[](auto x){return x.member;})
#define P4DELTAR ROOT::Math::VectorUtil::DeltaR<ROOT::Math::PtEtaPhiMVector,ROOT::Math::PtEtaPhiMVector> 

''' 
print >> sys.stderr, "Number of known columns", len(flow.validCols)
#flow.printRDFCpp(["GenQQ_mass","QJet_indices","QJet0","QJet1","Rpt","SBClassifier","qqDeltaEta","MqqGenJet"])
for x in colsToPlotWithSyst :
	flow.Histo(x)



flow.printRDFCpp(colsToPlotWithSyst+flow.weights.keys()+["Jet_LeptonDr","Lepton_JetDr","Jet_LeptonIdx","Jet_pt","Jet_eta","Jet_phi","Lepton_eta","Lepton_phi","Lepton_JetIdx","Lepton_jetIdx"],debug=False)
#["defaultWeight","QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","HighestGenQQMass","QJet0","QJet1","qqDeltaEta","MqqGenJet"]+[x for x in flow.validCols if x[:len("SBClassifier")]=="SBClassifier"]+flow.inputs["SBClassifier"]+flow.weights.keys(),debug=False)


#flow.printRDF(list(flow.allNodesTo("SBClassifier")))


