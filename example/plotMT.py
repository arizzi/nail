import ROOT
#from samples2016 import samples
from samples2017 import samples #as samples2
#samples.update(samples2)
from models2017H import *
from math import *
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
import re

totev={}
totevCount={}
totevSkim={}
hnForSys={}
def findSyst(hn,sy,f) :
    if hn in hnForSys and sy in hnForSys[hn]:
	return hnForSys[hn][sy]
    if hn not in hnForSys :
	hnForSys[hn]={}
    allh=list([x.GetName() for x in f.GetListOfKeys()])
    h1=hn+"__syst__"+sy    
    h2=re.sub("___","__syst__"+sy+"___",hn)    
    h3=re.sub("___","__syst__"+sy+"___",hn)+"__syst__"+sy    
    print("Syst options",h1,h2,h3)
    if h1 in allh:
	 hnForSys[hn]=h1
	 return h1
    if h2 in allh:
	 hnForSys[hn]=h2
	 return h2
    if h3 in allh:
	 hnForSys[hn]=h3
	 return h3
    print("none matching")
    return ""

def totevents(s):
    if s not in totev:
       totev[s]=1e-9
       totevCount[s]=1e-9
       totevSkim[s]=1e-9
       for fn in samples[s]["files"]:
	  f=ROOT.TFile.Open(fn)
	  run=f.Get("Runs")
	  totevSkim[s]+=f.Get("Events").GetEntries()
	  if run :
  	    hw=ROOT.TH1F("hw","", 5,0,5)
	    run.Project("hw","1","genEventSumw")
	    totev[s]+=hw.GetSumOfWeights()
	    run.Project("hw","1","genEventCount")
	    totevCount[s]+=hw.GetSumOfWeights()
	    print(totev[s])
#    print "returning",totev[s], "for",s
    return totev[s]



f={}
for group in signal :
    for s in signal[group] :
        f[s]=ROOT.TFile.Open("out/%sHistos.root"%s)
for group in background :
    for b in background[group] :
        f[b]=ROOT.TFile.Open("out/%sHistos.root"%b)
for group in data :
    for d in data[group] :
        f[d]=ROOT.TFile.Open("out/%sHistos.root"%d)

histoNames=[x.GetName() for x in f["data"].GetListOfKeys() ]
canvas={}
datastack={}
datasum={}
histos={}
histosum={}
histosSig={}
histoSigsum={}

datasumSyst={}
histosumSyst={}
histoSigsumSyst={}

#i=1
ROOT.gStyle.SetOptStat(0)
def makeplot(hn):
 if "__syst__" not in hn and "LHE" not in hn :
   print("Making histo",hn)
   histos[hn]=ROOT.THStack(hn,hn) 
   histosSig[hn]=ROOT.THStack(hn,hn) 
   datastack[hn]=ROOT.THStack(hn,hn) 
     
   canvas[hn]=ROOT.TCanvas("canvas_"+hn,"",900,750) 
   canvas[hn].Divide(1,2)
   canvas[hn].GetPad(2).SetPad(0,0,1,0.25) 
   canvas[hn].GetPad(1).SetPad(0,0.25,1,1) 
   #for gr in sorted(background,key=lambda g:):
   lumitot=0
   for gr in data:
     for d in data[gr]:
      lumitot+=samples[d]["lumi"]	
      if f[d] :
        h=f[d].Get(hn)
	if h:
   	   h.SetMarkerStyle(21)
	   datastack[hn].Add(h)
	   if hn not in datasum :
		datasum[hn]=h.Clone()
		datasumSyst[hn]={}
   		for sy in systematicsToPlot :
                  datasumSyst[hn][sy]=h.Clone()
	   else :
		datasum[hn].Add(h)	
   		for sy in systematicsToPlot :
		  hs=f[d].Get(findSyst(hn,sy,f[d]))
		  if hs:
                    datasumSyst[hn][sy].Add(hs)
		  else :
                    datasumSyst[hn][sy].Add(h)
#   print "Lumi tot", lumitot

   for gr in backgroundSorted:
     for b in background[gr]: 
      nevents=totevents(b)
      if f[b] :
	h=f[b].Get(hn)
	if h :
	   h.Scale(samples[b]["xsec"]/nevents*lumitot)
#	   print "adding", b, "to", hn 
	   h.SetFillColor(fillcolor[gr])
	   h.SetMarkerColor(markercolor[gr])
	   h.SetLineColor(linecolor[gr])
#	   i+=1
	   if hn not in histosum :
		histosum[hn]=h.Clone()
		histosumSyst[hn]={}
   		for sy in systematicsToPlot :
      		  hs=f[b].Get(findSyst(hn,sy,f[b]))
	          if hs:
        	     hs.Scale(samples[b]["xsec"]/nevents*lumitot)
	             histosumSyst[hn][sy]=hs.Clone()
		  else :
                     histosumSyst[hn][sy]=h.Clone()
	   else :
		histosum[hn].Add(h)	
   		for sy in systematicsToPlot :
		  hs=f[b].Get(findSyst(hn,sy,f[d]))
		  if hs:
	   	    hs.Scale(samples[b]["xsec"]/nevents*lumitot)
                    histosumSyst[hn][sy].Add(hs)
		  else :
		    print("using unchanged histo")
                    histosumSyst[hn][sy].Add(h)
	   histos[hn].Add(h)

   for gr in signal:
     for b in signal[gr]:
      nevents=totevents(b)
      if f[b] :
        h=f[b].Get(hn)
        if h :
           h.Scale(samples[b]["xsec"]/nevents*lumitot)
#           print "adding", b, "to", hn
           h.SetFillColor(fillcolor[gr])
           h.SetMarkerColor(markercolor[gr])
           h.SetLineColor(linecolor[gr])
 #          i+=1
           if hn not in histoSigsum :
                histoSigsum[hn]=h.Clone()
		histoSigsumSyst[hn]={}
   		for sy in systematicsToPlot :
                  hs=f[b].Get(findSyst(hn,sy,f[b]))
                  if hs:
                    hs.Scale(samples[b]["xsec"]/nevents*lumitot)
                    histoSigsumSyst[hn][sy]=hs.Clone()
                  else :
	            histoSigsumSyst[hn][sy]=h.Clone()

           else :
                histoSigsum[hn].Add(h)
   		for sy in systematicsToPlot :
		  hs=f[b].Get(findSyst(hn,sy,f[b]))
		  if hs:
           	    hs.Scale(samples[b]["xsec"]/nevents*lumitot)
                    histoSigsumSyst[hn][sy].Add(hs)
		  else :
		    print("Not found",hn+"__syst__"+sy)
                    histoSigsumSyst[hn][sy].Add(h)
           histosSig[hn].Add(h)
           histos[hn].Add(h)
     
   firstBlind=100000
   lastBlind=-1
   for i in range(histosSig[hn].GetStack().Last().GetNbinsX()) :
	if histosSig[hn].GetStack().Last().GetBinContent(i) > 0.1*sqrt(abs(histos[hn].GetStack().Last().GetBinContent(i))) and histosSig[hn].GetStack().Last().GetBinContent(i)/(0.1+abs(histos[hn].GetStack().Last().GetBinContent(i))) > 0.05 :
		print("to blind",hn,i,abs(histos[hn].GetStack().Last().GetBinContent(i)), histosSig[hn].GetStack().Last().GetBinContent(i))	
	        if i < firstBlind:
		    firstBlind=i
                lastBlind=i
   for i in range(firstBlind,lastBlind) :
       datastack[hn].GetStack().Last().SetBinContent(i,0)
       datasum[hn].SetBinContent(i,0)
       print("blinded",i,hn)
   canvas[hn].cd(1)
   datastack[hn].Draw("E P")
   datastack[hn].GetXaxis().SetTitle(hn)
   histos[hn].Draw("hist same")
   datastack[hn].Draw("E P same")
   datastack[hn].GetHistogram().SetMarkerStyle(20)

   canvas[hn].Update()
   ratio=datasum[hn].Clone()
   ratio.Add(histosum[hn],-1.) 
   ratio.Divide(histosum[hn])
   ratio.SetMarkerStyle(20)
   canvas[hn].cd(2)
   canvas[hn].cd(2)
   ratio.SetLabelSize(datastack[hn].GetHistogram().GetLabelSize()*3)
   ratio.GetYaxis().SetLabelSize(datastack[hn].GetHistogram().GetLabelSize()*3)
   ratio.Draw()
   ratio.SetAxisRange(-0.5,0.5,"Y")
   ratio.GetYaxis().SetNdivisions(5)
   ratiosy=[]
   for j,sy in enumerate(systematicsToPlot):
       ratiosy.append(histosumSyst[hn][sy].Clone())
       ratiosy[-1].Add(histosum[hn],-1.)
       ratiosy[-1].Divide(histosum[hn])
       ratiosy[-1].SetLineColor(1+j)
       ratiosy[-1].SetLineStyle(j)
       ratiosy[-1].SetFillStyle(0)
       ratiosy[-1].Draw("same hist")
       print("Heu",hn,sy,histosumSyst[hn][sy].Integral(),histosum[hn].Integral(),lumitot,ratiosy[-1])
#   systematics=[x for x in histoNames if x[:hn.find("___")]==hn[:hn.find("___")] and "__syst__" in x]
#   print "available systematics",hn,systematics
#  for s in systematics:
#  getsum up , down
   canvas[hn].GetPad(2).SetGridy()
   canvas[hn].SaveAs("%s.png"%hn)	   
   #canvas[hn].SaveAs("%s.root"%hn)	   
   canvas[hn].GetPad(1).SetLogy(True)
   canvas[hn].SaveAs("%s_log.png"%hn)	   


his=[x for x in histoNames if "__syst__" not in x]
print(his[0])
makeplot(his[0]) #do once for caching normalizations

if True:
 from multiprocessing import Pool
 runpool = Pool(3)
 #toproc=[(x,y,i) for y in sams for i,x in enumerate(samples[y]["files"])]
 runpool.map(makeplot, his[1:])
else :
 for x in his[1:] :
    makeplot(x)

tot=0
for s in totevCount:
  tot+=totevSkim[s]

print(tot, "input events") 
