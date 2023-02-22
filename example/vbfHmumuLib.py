from nail.nail import *
import ROOT
import sys

from eventprocessing import flow

#define some event weights
from weights import *
addDefaultWeights(flow)
addMuEffWeight(flow)

from systematics import *
#addLheScale(flow)
addBtag(flow)
addMuScale(flow)
addCompleteJecs(flow)

from histograms import histosPerSelection

systematics=flow.variations #take all systematic variations
histosWithSystematics=flow.createSystematicBranches(systematics,histosPerSelection)

print("The following histograms will be created in the following regions")
for sel in  histosWithSystematics:
	print(sel,":",histosWithSystematics[sel])


print("Number of known columns", len(flow.validCols), file=sys.stderr)


#snap=["SideBand","nSoftActivityJet","SoftActivityJet_pt","SoftActivityJet_eta","SoftActivityJet_phi","SoftActivityJet_SelectedJetDr","SoftActivityJet_SelectedJetIdx","SoftActivityJet_SelectedMuonDr","SoftActivityJet_SelectedMuonIdx","VBFRegion","QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","QJet0","QJet1","qqDeltaEta","MqqGenJet"]

snap=[] #x[0] for x in flow.originalCols ]
snaplist=["QJet0_pt","QJet1_pt","QJet0_eta","QJet1_eta","Mqq","Higgs_pt","twoJets","twoOppositeSignMuons","PreSel","VBFRegion","MassWindow","SignalRegion","qqDeltaEta","event","HLT_IsoMu24","QJet0_pt_nom","QJet1_pt_nom","QJet0_puId","QJet1_puId","SBClassifier","Higgs_m","Mqq_log","mmjj_pt_log","NSoft5","ll_zstar","theta2","mmjj_pz_logabs","MaxJetAbsEta","ll_zstar_log"]

from histobinning import binningrules
flow.binningRules = binningrules
flow.printRDFCpp(snaplist,debug=False,outname="tmp.C",selections=histosWithSystematics,snap=snap,snapsel="PreSel",lib=True)
flow.printRDFCpp(snaplist,debug=False,outname="tmp1.C",selections=histosWithSystematics,snap=snap,snapsel="PreSel",lib=True)

#compile and process
import os
if sys.argv[1] != "nobuild" :
  os.system("rm eventProcessor.so")
  os.system("g++ -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o eventProcessor.so --shared -lTMVA -I..")

ROOT.gROOT.ProcessLine('''
ROOT::EnableImplicitMT(20);
''')
ROOT.gInterpreter.Declare('''
#include "ext.h"
''')
ROOT.gInterpreter.Declare('''
template <typename T>
class NodeCaster {
   public:
   static ROOT::RDF::RNode Cast(T rdf)
   {
      return ROOT::RDF::RNode(rdf);
   }
};
''')
ROOT.gSystem.Load("eventProcessor.so")

def CastToRNode(node):
   return ROOT.NodeCaster(node.__cppname__).Cast(node)


#from samples2016Zuncomp import samples
#from tttests import samples
#from vbfdy import samples
from samples2016Zshm import samples

#from samples2016Z import samples
#from samples import samples
import psutil
def f(ar):
#f,s,i=ar
     p = psutil.Process()
     print("Affinity",p.cpu_affinity())
     p.cpu_affinity( list(range(psutil.cpu_count())))
     print("After set Affinity",p.cpu_affinity())
     s,f=ar
     print(f)
     ROOT.gROOT.ProcessLine('''
     ROOT::EnableImplicitMT(20);
     ''')
     vf=ROOT.vector("string")()
     list(map(lambda x : vf.push_back(x), f))
     for x in vf:
	print(x)
     rdf=ROOT.RDataFrame("Events",vf)
     if rdf :
       try:
	 if s=="data" :
	   rdf=rdf.Define("isMC","false")
	   rdf=rdf.Define("genWeight","1.0f")
	   rdf=rdf.Define("puWeight","1.0f")
	   rdf=rdf.Define("btagWeight_CSVV2","1.0f")
	   rdf=rdf.Define("Jet_pt_nom","Jet_pt")
	   print("Defined")
	 else :
	   rdf=rdf.Define("isMC","true")
	 if "filter" in samples[s] :
	   rdf=rdf.Filter(samples[s]["filter"])

         ou=ROOT.processRDF(CastToRNode(rdf))
         fff=ROOT.TFile.Open("out/%sHistos.root"%(s),"recreate")
         #fff=ROOT.TFile.Open("out/%s%sHistos.root"%(s,i),"recreate")
#	 h=ou.Histo1D('Jet_pt')
#	 h.Write()
         # inferred.
         #snaplist=["QJet0_pt","QJet1_pt","QJet0_eta","QJet1_eta","Mqq","Higgs_pt","twoJets","twoOppositeSignMuons","PreSel","VBFRegion","MassWindow","SignalRegion"]
         branchList = ROOT.vector('string')()
	 list(map(lambda x : branchList.push_back(x), snaplist))

  #       ou.rdf.Filter("twoMuons","twoMuons").Filter("twoOppositeSignMuons","twoOppositeSignMuons").Filter("twoJets","twoJets").Filter("MassWindow","MassWindow").Filter("VBFRegion","VBFRegion").Filter("PreSel","PreSel").Filter("SignalRegion","SignalRegion").Snapshot("Events","out/%sSnapshot.root"%(s),branchList)
  #       ou.rdf.Filter("event==24331988").Snapshot("Events","out/%sEventPick.root"%(s),branchList)
         print(ou.histos.size())
         for h in ou.histos :
#	    print h
 	    h.Write()
         fff.Write()
         fff.Close()
	 return 0
       except Exception as e: 
	 print(e)
	 print("FAIL",f)
	 return 1
     else :
	print("Null file",f)

#     return  os.system("./eventProcessor %s %s out/%s%s "%(4,f,s,i))  

from multiprocessing import Pool
runpool = Pool(20)

print(list(samples.keys()))
sams=list(samples.keys())

#sams=["DY2J","TTlep"]
#toproc=[(x,y,i) for y in sams for i,x in enumerate(samples[y]["files"])]
toproc=[ (s,samples[s]["files"]) for s in sams ]
print(list(zip(runpool.map(f, toproc ),sams)))
