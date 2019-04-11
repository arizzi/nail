from nail import *
import ROOT
import sys

from eventprocessing import flow

#define some event weights
from weights import *
addBtagWeight(flow)
addMuEffWeight(flow)

from systematics import *
addLheScale(flow)
addBtag(flow)
addMuScale(flow)
addCompleteJecs(flow)

from histograms import histosPerSelection

systematics=flow.variations #take all systematic variations
histosWithSystematics=flow.createSystematicBranches(systematics,histosPerSelection)

print "The following histograms will be created in the following regions"
for sel in  histosWithSystematics:
	print sel,":",histosWithSystematics[sel]


print >> sys.stderr, "Number of known columns", len(flow.validCols)


snap=["SideBand","nSoftActivityJet","SoftActivityJet_pt","SoftActivityJet_eta","SoftActivityJet_phi","SoftActivityJet_SelectedJetDr","SoftActivityJet_SelectedJetIdx","SoftActivityJet_SelectedMuonDr","SoftActivityJet_SelectedMuonIdx","VBFRegion","QJet0_pt","QJet0_eta","QJet0_btagCSVV2","QJet1_pt","QJet1_eta","QJet1_btagCSVV2","Mu0_pt","Mu0_eta","Mu1_pt","Mu1_eta","QJet0","QJet1","qqDeltaEta","MqqGenJet"]

flow.printRDFCpp(snap,debug=False,outname="tmp.C",selections=histosWithSystematics,snap=snap,snapsel="TwoJetsTwoMu")

#compile and process
import os
os.system("rm eventProcessor")
os.system("g++ -fPIC -Wall -O3 tmp.C $(root-config --libs --cflags)  -o eventProcessor")

from samples import samples
for s in samples :
  i=0
  for f in samples[s]["files"]:
     os.system("./eventProcessor %s %s out/%s%s "%(48,f,s,i))
     i+=1

