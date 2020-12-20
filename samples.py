from nail import Sample

s1 = Sample("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root", "Events")
s2 = Sample("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/GluGlu_HToMuMu_nano2016.root", "Events")
signal = s1+s2
dy = Sample("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160_VBFFilter-madgraphMLM_nano2016.root", "Events")
