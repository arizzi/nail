from nail import Sample

s1=Sample("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root","Events")
s2=Sample("/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/GluGlu_HToMuMu_nano2016.root","Events")
signal=s1+s2


