
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
#include "helpers.h"
#define MemberMap(vector,member) Map(vector,[](auto x){return x.member;})
#define P4DELTAR ROOT::Math::VectorUtil::DeltaR<ROOT::Math::PtEtaPhiMVector,ROOT::Math::PtEtaPhiMVector> 


auto func__muPtCut() { 
 return 25; }
using type__muPtCut = ROOT::TypeTraits::CallableTraits<decltype(func__muPtCut)>::ret_type;
auto func__muIdCut() { 
 return 3; }
using type__muIdCut = ROOT::TypeTraits::CallableTraits<decltype(func__muIdCut)>::ret_type;
auto func__muIsoCut() { 
 return 0.13; }
using type__muIsoCut = ROOT::TypeTraits::CallableTraits<decltype(func__muIsoCut)>::ret_type;
auto func__Muon_id(const ROOT::VecOps::RVec<Bool_t> & Muon_mediumId, const ROOT::VecOps::RVec<Bool_t> & Muon_tightId) { 
 return Muon_tightId*3+Muon_mediumId; }
using type__Muon_id = ROOT::TypeTraits::CallableTraits<decltype(func__Muon_id)>::ret_type;
auto func__Muon_iso(const ROOT::VecOps::RVec<Float_t> & Muon_miniPFRelIso_all) { 
 return Muon_miniPFRelIso_all; }
using type__Muon_iso = ROOT::TypeTraits::CallableTraits<decltype(func__Muon_iso)>::ret_type;
auto func__LooseMuon(const type__muPtCut & muPtCut, const type__Muon_id & Muon_id, const type__Muon_iso & Muon_iso, const ROOT::VecOps::RVec<Float_t> & Muon_pt, const type__muIdCut & muIdCut, const type__muIsoCut & muIsoCut) { 
 return Muon_iso < muIsoCut && Muon_id > muIdCut && Muon_pt > muPtCut; }
using type__LooseMuon = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon)>::ret_type;
auto func__LooseMuon_eta(const ROOT::VecOps::RVec<Float_t> & Muon_eta, const type__LooseMuon & LooseMuon) { 
 return Muon_eta[LooseMuon]; }
using type__LooseMuon_eta = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_eta)>::ret_type;
auto func__LooseMuon_mass(const ROOT::VecOps::RVec<Float_t> & Muon_mass, const type__LooseMuon & LooseMuon) { 
 return Muon_mass[LooseMuon]; }
using type__LooseMuon_mass = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_mass)>::ret_type;
auto func__LooseMuon_phi(const ROOT::VecOps::RVec<Float_t> & Muon_phi, const type__LooseMuon & LooseMuon) { 
 return Muon_phi[LooseMuon]; }
using type__LooseMuon_phi = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_phi)>::ret_type;
auto func__LooseMuon_pt(const type__LooseMuon & LooseMuon, const ROOT::VecOps::RVec<Float_t> & Muon_pt) { 
 return Muon_pt[LooseMuon]; }
using type__LooseMuon_pt = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_pt)>::ret_type;
auto func__LooseMuon_jetIdx(const ROOT::VecOps::RVec<Int_t> & Muon_jetIdx, const type__LooseMuon & LooseMuon) { 
 return Muon_jetIdx[LooseMuon]; }
using type__LooseMuon_jetIdx = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_jetIdx)>::ret_type;
auto func__nLooseMuon(const type__LooseMuon & LooseMuon) { 
 return Sum(LooseMuon); }
using type__nLooseMuon = ROOT::TypeTraits::CallableTraits<decltype(func__nLooseMuon)>::ret_type;
auto func__LooseMuon_p4(const type__LooseMuon_mass & LooseMuon_mass, const type__LooseMuon_pt & LooseMuon_pt, const type__LooseMuon_phi & LooseMuon_phi, const type__LooseMuon_eta & LooseMuon_eta) { 
 return vector_map_t<ROOT::Math::PtEtaPhiMVector>(LooseMuon_pt , LooseMuon_eta, LooseMuon_phi, LooseMuon_mass); }
using type__LooseMuon_p4 = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_p4)>::ret_type;
auto func__Electron_p4(const ROOT::VecOps::RVec<Float_t> & Electron_mass, const ROOT::VecOps::RVec<Float_t> & Electron_pt, const ROOT::VecOps::RVec<Float_t> & Electron_phi, const ROOT::VecOps::RVec<Float_t> & Electron_eta) { 
 return vector_map_t<ROOT::Math::PtEtaPhiMVector>(Electron_pt , Electron_eta, Electron_phi, Electron_mass); }
using type__Electron_p4 = ROOT::TypeTraits::CallableTraits<decltype(func__Electron_p4)>::ret_type;
auto func__eIdCut() { 
 return 3; }
using type__eIdCut = ROOT::TypeTraits::CallableTraits<decltype(func__eIdCut)>::ret_type;
auto func__eIsoCut() { 
 return 0.13; }
using type__eIsoCut = ROOT::TypeTraits::CallableTraits<decltype(func__eIsoCut)>::ret_type;
auto func__ePtCut() { 
 return 25; }
using type__ePtCut = ROOT::TypeTraits::CallableTraits<decltype(func__ePtCut)>::ret_type;
auto func__Electron_id(const ROOT::VecOps::RVec<Int_t> & Electron_cutBased) { 
 return Electron_cutBased; }
using type__Electron_id = ROOT::TypeTraits::CallableTraits<decltype(func__Electron_id)>::ret_type;
auto func__Electron_iso(const ROOT::VecOps::RVec<Float_t> & Electron_miniPFRelIso_all) { 
 return Electron_miniPFRelIso_all; }
using type__Electron_iso = ROOT::TypeTraits::CallableTraits<decltype(func__Electron_iso)>::ret_type;
auto func__LooseEle(const type__Electron_id & Electron_id, const type__ePtCut & ePtCut, const type__eIdCut & eIdCut, const type__eIsoCut & eIsoCut, const type__Electron_iso & Electron_iso, const ROOT::VecOps::RVec<Float_t> & Electron_pt) { 
 return Electron_iso < eIsoCut && Electron_id > eIdCut && Electron_pt > ePtCut; }
using type__LooseEle = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle)>::ret_type;
auto func__LooseEle_eta(const type__LooseEle & LooseEle, const ROOT::VecOps::RVec<Float_t> & Electron_eta) { 
 return Electron_eta[LooseEle]; }
using type__LooseEle_eta = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_eta)>::ret_type;
auto func__LooseEle_phi(const ROOT::VecOps::RVec<Float_t> & Electron_phi, const type__LooseEle & LooseEle) { 
 return Electron_phi[LooseEle]; }
using type__LooseEle_phi = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_phi)>::ret_type;
auto func__LooseEle_pt(const ROOT::VecOps::RVec<Float_t> & Electron_pt, const type__LooseEle & LooseEle) { 
 return Electron_pt[LooseEle]; }
using type__LooseEle_pt = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_pt)>::ret_type;
auto func__LooseEle_jetIdx(const ROOT::VecOps::RVec<Int_t> & Electron_jetIdx, const type__LooseEle & LooseEle) { 
 return Electron_jetIdx[LooseEle]; }
using type__LooseEle_jetIdx = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_jetIdx)>::ret_type;
auto func__LooseEle_p4(const type__Electron_p4 & Electron_p4, const type__LooseEle & LooseEle) { 
 return Electron_p4[LooseEle]; }
using type__LooseEle_p4 = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_p4)>::ret_type;
auto func__nLooseEle(const type__LooseEle & LooseEle) { 
 return Sum(LooseEle); }
using type__nLooseEle = ROOT::TypeTraits::CallableTraits<decltype(func__nLooseEle)>::ret_type;
auto func__LooseMuon_pid(const type__LooseMuon_pt & LooseMuon_pt) { 
 return (LooseMuon_pt*0)+13; }
using type__LooseMuon_pid = ROOT::TypeTraits::CallableTraits<decltype(func__LooseMuon_pid)>::ret_type;
auto func__LooseEle_pid(const type__LooseEle_pt & LooseEle_pt) { 
 return (LooseEle_pt*0)+11; }
using type__LooseEle_pid = ROOT::TypeTraits::CallableTraits<decltype(func__LooseEle_pid)>::ret_type;
auto func__nLepton(const type__nLooseEle & nLooseEle, const type__nLooseMuon & nLooseMuon) { 
 return nLooseMuon+nLooseEle; }
using type__nLepton = ROOT::TypeTraits::CallableTraits<decltype(func__nLepton)>::ret_type;
auto func__Lepton_pid(const type__LooseMuon_pid & LooseMuon_pid, const type__LooseEle_pid & LooseEle_pid) { 
 return Concat(LooseMuon_pid,LooseEle_pid); }
using type__Lepton_pid = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_pid)>::ret_type;
auto func__Lepton_phi(const type__LooseEle_phi & LooseEle_phi, const type__LooseMuon_phi & LooseMuon_phi) { 
 return Concat(LooseMuon_phi,LooseEle_phi); }
using type__Lepton_phi = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_phi)>::ret_type;
auto func__Lepton_p4(const type__LooseEle_p4 & LooseEle_p4, const type__LooseMuon_p4 & LooseMuon_p4) { 
 return Concat(LooseMuon_p4,LooseEle_p4); }
using type__Lepton_p4 = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_p4)>::ret_type;
auto func__Lepton_eta(const type__LooseEle_eta & LooseEle_eta, const type__LooseMuon_eta & LooseMuon_eta) { 
 return Concat(LooseMuon_eta,LooseEle_eta); }
using type__Lepton_eta = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_eta)>::ret_type;
auto func__Lepton_jetIdx(const type__LooseMuon_jetIdx & LooseMuon_jetIdx, const type__LooseEle_jetIdx & LooseEle_jetIdx) { 
 return Concat(LooseMuon_jetIdx,LooseEle_jetIdx); }
using type__Lepton_jetIdx = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_jetIdx)>::ret_type;
auto func__Jet_p4(const ROOT::VecOps::RVec<Float_t> & Jet_mass, const ROOT::VecOps::RVec<Float_t> & Jet_pt, const ROOT::VecOps::RVec<Float_t> & Jet_phi, const ROOT::VecOps::RVec<Float_t> & Jet_eta) { 
 return vector_map_t<ROOT::Math::PtEtaPhiMVector>(Jet_pt , Jet_eta, Jet_phi, Jet_mass); }
using type__Jet_p4 = ROOT::TypeTraits::CallableTraits<decltype(func__Jet_p4)>::ret_type;
auto func__JetLeptonPair(const type__nLepton & nLepton, const UInt_t & nJet) { 
 return Combinations(nJet,nLepton); }
using type__JetLeptonPair = ROOT::TypeTraits::CallableTraits<decltype(func__JetLeptonPair)>::ret_type;
auto func__JetLeptonPair_dr(const type__Lepton_p4 & Lepton_p4, const type__Jet_p4 & Jet_p4, const type__JetLeptonPair & JetLeptonPair) { 
 return vector_map(P4DELTAR,Take(Jet_p4,JetLeptonPair[0]),Take(Lepton_p4,JetLeptonPair[1])); }
using type__JetLeptonPair_dr = ROOT::TypeTraits::CallableTraits<decltype(func__JetLeptonPair_dr)>::ret_type;
auto func__Jet_LeptonDr(const type__nLepton & nLepton, const type__JetLeptonPair_dr & JetLeptonPair_dr, const UInt_t & nJet) { 
 return matrix_map(nJet,nLepton,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):-99;},JetLeptonPair_dr); }
using type__Jet_LeptonDr = ROOT::TypeTraits::CallableTraits<decltype(func__Jet_LeptonDr)>::ret_type;
auto func__Lepton_JetDr(const type__nLepton & nLepton, const type__JetLeptonPair_dr & JetLeptonPair_dr, const UInt_t & nJet) { 
 return matrix_map(nJet,nLepton,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):-99;},JetLeptonPair_dr); }
using type__Lepton_JetDr = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_JetDr)>::ret_type;
auto func__Jet_LeptonIdx(const type__nLepton & nLepton, const type__JetLeptonPair_dr & JetLeptonPair_dr, const UInt_t & nJet) { 
 return matrix_map(nJet,nLepton,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):-1;},JetLeptonPair_dr); }
using type__Jet_LeptonIdx = ROOT::TypeTraits::CallableTraits<decltype(func__Jet_LeptonIdx)>::ret_type;
auto func__Lepton_JetIdx(const type__nLepton & nLepton, const type__JetLeptonPair_dr & JetLeptonPair_dr, const UInt_t & nJet) { 
 return matrix_map(nJet,nLepton,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):-1;},JetLeptonPair_dr); }
using type__Lepton_JetIdx = ROOT::TypeTraits::CallableTraits<decltype(func__Lepton_JetIdx)>::ret_type;
auto func__jetPtCut() { 
 return 25; }
using type__jetPtCut = ROOT::TypeTraits::CallableTraits<decltype(func__jetPtCut)>::ret_type;
auto func__jetIdCut() { 
 return 0; }
using type__jetIdCut = ROOT::TypeTraits::CallableTraits<decltype(func__jetIdCut)>::ret_type;
auto func__jetPUIdCut() { 
 return 0; }
using type__jetPUIdCut = ROOT::TypeTraits::CallableTraits<decltype(func__jetPUIdCut)>::ret_type;
auto func__CleanJet(const type__jetPUIdCut & jetPUIdCut, const type__Jet_LeptonIdx & Jet_LeptonIdx, const type__jetPtCut & jetPtCut, const ROOT::VecOps::RVec<Int_t> & Jet_puId, const ROOT::VecOps::RVec<Int_t> & Jet_jetId, const ROOT::VecOps::RVec<Float_t> & Jet_pt, const type__Jet_LeptonDr & Jet_LeptonDr, const type__jetIdCut & jetIdCut) { 
 return 
 Jet_pt > jetPtCut &&
 Jet_jetId > jetIdCut &&
 Jet_puId > jetPUIdCut &&
 (Jet_LeptonIdx==-1 || Jet_LeptonDr > 0.3)
; }
using type__CleanJet = ROOT::TypeTraits::CallableTraits<decltype(func__CleanJet)>::ret_type;
auto func__CleanJet_eta(const type__CleanJet & CleanJet, const ROOT::VecOps::RVec<Float_t> & Jet_eta) { 
 return Jet_eta[CleanJet]; }
using type__CleanJet_eta = ROOT::TypeTraits::CallableTraits<decltype(func__CleanJet_eta)>::ret_type;
auto func__CleanJet_pt(const type__CleanJet & CleanJet, const ROOT::VecOps::RVec<Float_t> & Jet_pt) { 
 return Jet_pt[CleanJet]; }
using type__CleanJet_pt = ROOT::TypeTraits::CallableTraits<decltype(func__CleanJet_pt)>::ret_type;
auto func__CleanJet_btag(const ROOT::VecOps::RVec<Float_t> & Jet_btagCSVV2) { 
 return Jet_btagCSVV2; }
using type__CleanJet_btag = ROOT::TypeTraits::CallableTraits<decltype(func__CleanJet_btag)>::ret_type;
auto func__CleanJet_btagWeight(const type__CleanJet_pt & CleanJet_pt, const type__CleanJet_eta & CleanJet_eta, const type__CleanJet_btag & CleanJet_btag) { 
 return vector_map(btagWeight,CleanJet_btag,CleanJet_pt,CleanJet_eta); }
using type__CleanJet_btagWeight = ROOT::TypeTraits::CallableTraits<decltype(func__CleanJet_btagWeight)>::ret_type;
auto func__btagEventWeight(const type__CleanJet_btagWeight & CleanJet_btagWeight) { 
 return std::accumulate(CleanJet_btagWeight.begin(),CleanJet_btagWeight.end(),1, std::multiplies<double>()); }
using type__btagEventWeight = ROOT::TypeTraits::CallableTraits<decltype(func__btagEventWeight)>::ret_type;
auto func__defaultWeight(const type__btagEventWeight & btagEventWeight, const Float_t & genWeight) { 
 return genWeight*btagEventWeight; }
using type__defaultWeight = ROOT::TypeTraits::CallableTraits<decltype(func__defaultWeight)>::ret_type;
auto func__LHEScaleWeight0(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[0])*defaultWeight; }
using type__LHEScaleWeight0 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight0)>::ret_type;
auto func__LHEScaleWeight1(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[1])*defaultWeight; }
using type__LHEScaleWeight1 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight1)>::ret_type;
auto func__LHEScaleWeight2(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[2])*defaultWeight; }
using type__LHEScaleWeight2 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight2)>::ret_type;
auto func__LHEScaleWeight3(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[3])*defaultWeight; }
using type__LHEScaleWeight3 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight3)>::ret_type;
auto func__LHEScaleWeight4(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[4])*defaultWeight; }
using type__LHEScaleWeight4 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight4)>::ret_type;
auto func__LHEScaleWeight5(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[5])*defaultWeight; }
using type__LHEScaleWeight5 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight5)>::ret_type;
auto func__LHEScaleWeight6(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[6])*defaultWeight; }
using type__LHEScaleWeight6 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight6)>::ret_type;
auto func__LHEScaleWeight7(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[7])*defaultWeight; }
using type__LHEScaleWeight7 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight7)>::ret_type;
auto func__LHEScaleWeight8(const ROOT::VecOps::RVec<Float_t> & LHEScaleWeight, const type__defaultWeight & defaultWeight) { 
 return (LHEScaleWeight[8])*defaultWeight; }
using type__LHEScaleWeight8 = ROOT::TypeTraits::CallableTraits<decltype(func__LHEScaleWeight8)>::ret_type;

int main(int argc, char** argv)
{
   auto n_cores = 0;
   if (argc > 1)
      n_cores = std::atoi(argv[1]);
   if (n_cores > 0)
      ROOT::EnableImplicitMT(n_cores);


ROOT::RDataFrame rdf("Events","/dev/shm/VBF_HToMuMu_nano2016.root");
auto toplevel =
rdf.Define("muPtCut",func__muPtCut,{})
.Define("muIdCut",func__muIdCut,{})
.Define("muIsoCut",func__muIsoCut,{})
.Define("Muon_id",func__Muon_id,{"Muon_mediumId","Muon_tightId"})
.Define("Muon_iso",func__Muon_iso,{"Muon_miniPFRelIso_all"})
.Define("LooseMuon",func__LooseMuon,{"muPtCut","Muon_id","Muon_iso","Muon_pt","muIdCut","muIsoCut"})
.Define("LooseMuon_eta",func__LooseMuon_eta,{"Muon_eta","LooseMuon"})
.Define("LooseMuon_mass",func__LooseMuon_mass,{"Muon_mass","LooseMuon"})
.Define("LooseMuon_phi",func__LooseMuon_phi,{"Muon_phi","LooseMuon"})
.Define("LooseMuon_pt",func__LooseMuon_pt,{"LooseMuon","Muon_pt"})
.Define("LooseMuon_jetIdx",func__LooseMuon_jetIdx,{"Muon_jetIdx","LooseMuon"})
.Define("nLooseMuon",func__nLooseMuon,{"LooseMuon"})
.Define("LooseMuon_p4",func__LooseMuon_p4,{"LooseMuon_mass","LooseMuon_pt","LooseMuon_phi","LooseMuon_eta"})
.Define("Electron_p4",func__Electron_p4,{"Electron_mass","Electron_pt","Electron_phi","Electron_eta"})
.Define("eIdCut",func__eIdCut,{})
.Define("eIsoCut",func__eIsoCut,{})
.Define("ePtCut",func__ePtCut,{})
.Define("Electron_id",func__Electron_id,{"Electron_cutBased"})
.Define("Electron_iso",func__Electron_iso,{"Electron_miniPFRelIso_all"})
.Define("LooseEle",func__LooseEle,{"Electron_id","ePtCut","eIdCut","eIsoCut","Electron_iso","Electron_pt"})
.Define("LooseEle_eta",func__LooseEle_eta,{"LooseEle","Electron_eta"})
.Define("LooseEle_phi",func__LooseEle_phi,{"Electron_phi","LooseEle"})
.Define("LooseEle_pt",func__LooseEle_pt,{"Electron_pt","LooseEle"})
.Define("LooseEle_jetIdx",func__LooseEle_jetIdx,{"Electron_jetIdx","LooseEle"})
.Define("LooseEle_p4",func__LooseEle_p4,{"Electron_p4","LooseEle"})
.Define("nLooseEle",func__nLooseEle,{"LooseEle"})
.Define("LooseMuon_pid",func__LooseMuon_pid,{"LooseMuon_pt"})
.Define("LooseEle_pid",func__LooseEle_pid,{"LooseEle_pt"})
.Define("nLepton",func__nLepton,{"nLooseEle","nLooseMuon"})
.Define("Lepton_pid",func__Lepton_pid,{"LooseMuon_pid","LooseEle_pid"})
.Define("Lepton_phi",func__Lepton_phi,{"LooseEle_phi","LooseMuon_phi"})
.Define("Lepton_p4",func__Lepton_p4,{"LooseEle_p4","LooseMuon_p4"})
.Define("Lepton_eta",func__Lepton_eta,{"LooseEle_eta","LooseMuon_eta"})
.Define("Lepton_jetIdx",func__Lepton_jetIdx,{"LooseMuon_jetIdx","LooseEle_jetIdx"})
.Define("Jet_p4",func__Jet_p4,{"Jet_mass","Jet_pt","Jet_phi","Jet_eta"})
.Define("JetLeptonPair",func__JetLeptonPair,{"nLepton","nJet"})
.Define("JetLeptonPair_dr",func__JetLeptonPair_dr,{"Lepton_p4","Jet_p4","JetLeptonPair"})
.Define("Jet_LeptonDr",func__Jet_LeptonDr,{"nLepton","JetLeptonPair_dr","nJet"})
.Define("Lepton_JetDr",func__Lepton_JetDr,{"nLepton","JetLeptonPair_dr","nJet"})
.Define("Jet_LeptonIdx",func__Jet_LeptonIdx,{"nLepton","JetLeptonPair_dr","nJet"})
.Define("Lepton_JetIdx",func__Lepton_JetIdx,{"nLepton","JetLeptonPair_dr","nJet"})
.Define("jetPtCut",func__jetPtCut,{})
.Define("jetIdCut",func__jetIdCut,{})
.Define("jetPUIdCut",func__jetPUIdCut,{})
.Define("CleanJet",func__CleanJet,{"jetPUIdCut","Jet_LeptonIdx","jetPtCut","Jet_puId","Jet_jetId","Jet_pt","Jet_LeptonDr","jetIdCut"})
.Define("CleanJet_eta",func__CleanJet_eta,{"CleanJet","Jet_eta"})
.Define("CleanJet_pt",func__CleanJet_pt,{"CleanJet","Jet_pt"})
.Define("CleanJet_btag",func__CleanJet_btag,{"Jet_btagCSVV2"})
.Define("CleanJet_btagWeight",func__CleanJet_btagWeight,{"CleanJet_pt","CleanJet_eta","CleanJet_btag"})
.Define("btagEventWeight",func__btagEventWeight,{"CleanJet_btagWeight"})
.Define("defaultWeight",func__defaultWeight,{"btagEventWeight","genWeight"})
.Define("LHEScaleWeight0",func__LHEScaleWeight0,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight1",func__LHEScaleWeight1,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight2",func__LHEScaleWeight2,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight3",func__LHEScaleWeight3,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight4",func__LHEScaleWeight4,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight5",func__LHEScaleWeight5,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight6",func__LHEScaleWeight6,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight7",func__LHEScaleWeight7,{"LHEScaleWeight","defaultWeight"})
.Define("LHEScaleWeight8",func__LHEScaleWeight8,{"LHEScaleWeight","defaultWeight"})
;
auto nJet=toplevel.Histo1D("nJet","defaultWeight");
auto nJet__weight__LHEScaleWeight8=toplevel.Histo1D("nJet","LHEScaleWeight8");
auto nJet__weight__LHEScaleWeight0=toplevel.Histo1D("nJet","LHEScaleWeight0");
auto nJet__weight__LHEScaleWeight1=toplevel.Histo1D("nJet","LHEScaleWeight1");
auto nJet__weight__LHEScaleWeight2=toplevel.Histo1D("nJet","LHEScaleWeight2");
auto nJet__weight__LHEScaleWeight3=toplevel.Histo1D("nJet","LHEScaleWeight3");
auto nJet__weight__LHEScaleWeight4=toplevel.Histo1D("nJet","LHEScaleWeight4");
auto nJet__weight__LHEScaleWeight5=toplevel.Histo1D("nJet","LHEScaleWeight5");
auto nJet__weight__LHEScaleWeight6=toplevel.Histo1D("nJet","LHEScaleWeight6");
auto nJet__weight__LHEScaleWeight7=toplevel.Histo1D("nJet","LHEScaleWeight7");
auto nLepton=toplevel.Histo1D("nLepton","defaultWeight");
auto nLepton__weight__LHEScaleWeight8=toplevel.Histo1D("nLepton","LHEScaleWeight8");
auto nLepton__weight__LHEScaleWeight0=toplevel.Histo1D("nLepton","LHEScaleWeight0");
auto nLepton__weight__LHEScaleWeight1=toplevel.Histo1D("nLepton","LHEScaleWeight1");
auto nLepton__weight__LHEScaleWeight2=toplevel.Histo1D("nLepton","LHEScaleWeight2");
auto nLepton__weight__LHEScaleWeight3=toplevel.Histo1D("nLepton","LHEScaleWeight3");
auto nLepton__weight__LHEScaleWeight4=toplevel.Histo1D("nLepton","LHEScaleWeight4");
auto nLepton__weight__LHEScaleWeight5=toplevel.Histo1D("nLepton","LHEScaleWeight5");
auto nLepton__weight__LHEScaleWeight6=toplevel.Histo1D("nLepton","LHEScaleWeight6");
auto nLepton__weight__LHEScaleWeight7=toplevel.Histo1D("nLepton","LHEScaleWeight7");


   auto tr=toplevel.Snapshot("ot", "outputFile.root", {"nJet","nLepton","Jet_LeptonDr","Lepton_JetDr","Jet_LeptonIdx","Jet_pt","Jet_eta","Jet_phi","Lepton_eta","Lepton_phi","Lepton_JetIdx","Lepton_jetIdx"});
   

   return 0;
}

