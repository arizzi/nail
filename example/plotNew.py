import ROOT
#from samples2016 import samples
from samples2018 import samples #as samples2
#samples.update(samples2)
from models2018Z import *
from labelDict import *  
year=2018             
lumi = "%2.1f fb^{-1}" 
from math import *
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
import re

totev={}
totevCount={}
totevSkim={}
hnForSys={}
def makeText (x, y, someText, font) :
    tex = ROOT.TLatex(x,y,someText);
    tex.SetNDC();
    tex.SetTextAlign(35);
    tex.SetTextFont(font);
    tex.SetTextSize(0.045);
    tex.SetLineWidth(2);
    return tex


    
def setHistoStyle (h, gr) :
    h.SetFillColor(fillcolor[gr])
    h.SetTitle("")
    h.SetLineColor(linecolor[gr])
    h.SetFillStyle(1001) #NEW
    h.SetLineStyle(1) #NEW    
    
    
def makeRatioMCplot(h) :
    hMC = h.Clone()
    hMC.SetLineWidth(1)
    for n in range(hMC.GetNbinsX()) :
       # hMC.SetBinError(n,  hMC.GetBinError(n)/hMC.GetBinContent(n) if hMC.GetBinContent(n)>0 else 0 )
        e = hMC.GetBinError(n)/hMC.GetBinContent(n) if hMC.GetBinContent(n)>0 else 0     
        hMC.SetBinError(n,  e if e<0.5 else 0.5 )                                        
        hMC.SetBinContent(n, 0.)
    return hMC
    
def setStyle(h, isRatio=False) :
    h.SetTitle("")
    w = 0.055 * (2. if isRatio else 1.)
    h.GetYaxis().SetLabelSize(w)
    h.GetXaxis().SetLabelSize(w)
    h.GetYaxis().SetTitleSize(w)
    h.GetXaxis().SetTitleSize(w)
    if isRatio : 
        h.GetYaxis().SetTitle("Data/MC - 1")
        h.GetYaxis().SetTitleOffset(0.5)
#        h.GetXaxis().SetTitle(str(h.GetName()).split("___")[0])
	xKey = str(h.GetName()).split("___")[0]                                              
	h.GetXaxis().SetTitle(labelVariable[xKey] if xKey in labelVariable.keys() else xKey) 
    else :
        binWidht = str(h.GetBinWidth(1))[:4]
        if binWidht.endswith(".") : binWidht = binWidht[:3]
        h.GetXaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle("Entries/"+binWidht)
        h.GetXaxis().SetLabelSize(0)
        h.GetXaxis().SetTitleSize(0)


def findSyst(hn,sy,f) :
  #  print hnForSys.keys()
    if hn in hnForSys and sy in hnForSys[hn]:
 #       print hn,sy,hnForSys[hn]
	return hnForSys[hn][sy]
    if hn not in hnForSys :
	hnForSys[hn]={}
    allh=list([x.GetName() for x in f.GetListOfKeys()])
    h1=hn+"__syst__"+sy    
    h2=re.sub("___","__syst__"+sy+"___",hn)    
    h3=re.sub("___","__syst__"+sy+"___",hn)+"__syst__"+sy    
#    print "Syst options",h1,h2,h3
    if h1 in allh:
	 hnForSys[hn][sy]=h1
	 return h1
    if h2 in allh:
	 hnForSys[hn][sy]=h2
	 return h2
    if h3 in allh:
	 hnForSys[hn][sy]=h3
	 return h3
    print "none matching"
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
	    print totev[s]
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

histoNames=list(set([x.GetName() for y in f.keys() for x in f[y].GetListOfKeys() ]))

canvas={}
datastack={}
datasum={}
histos={}
histosum={}
histosSig={}
histoSigsum={}
histoTH = {} 

datasumSyst={}
histosumSyst={}
histoSigsumSyst={}
histosSignal={}

#i=1
ROOT.gStyle.SetOptStat(0)


def makeplot(hn):
 myLegend= ROOT.TLegend(0.85, 0.4, 1, 0.9, "") 
 myLegend.SetFillColor(0);                          
 myLegend.SetBorderSize(0);                         
 myLegend.SetTextFont(42);                          
 myLegend.SetTextSize(0.025);                       
 if "__syst__" not in hn and "LHE" not in hn :
   #print "Making histo",hn
   histos[hn]=ROOT.THStack(hn,hn) 
   histosSig[hn]=ROOT.THStack(hn,hn) 
   datastack[hn]=ROOT.THStack(hn,hn) 

   canvas[hn]=ROOT.TCanvas("canvas_"+hn,"",900,750)       
   canvas[hn].SetRightMargin(.05);                        
   canvas[hn].Divide(1,2)                                 
   canvas[hn].GetPad(2).SetPad(0.0,0.,0.85,0.3)           
   canvas[hn].GetPad(1).SetPad(0.0,0.25,0.85,1.) 

   ROOT.gStyle.SetPadLeftMargin(0.15)                     
   canvas[hn].GetPad(2).SetBottomMargin(0.3)              
   canvas[hn].GetPad(2).SetTopMargin(0.05)                
                                                          

#   canvas[hn]=ROOT.TCanvas("canvas_"+hn,"",900,750) 
#   canvas[hn].Divide(1,2)
#   canvas[hn].GetPad(2).SetPad(0,0,1,0.25) 
#   canvas[hn].GetPad(1).SetPad(0,0.25,1,1) 
   #for gr in sorted(background,key=lambda g:):
   lumitot=0
   for gr in data:
     for d in data[gr]:
      lumitot+=samples[d]["lumi"]	
      if f[d] :
        h=f[d].Get(hn)
	if h:
   	   h.SetMarkerStyle(10)
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
        else:
	   print "Cannot open",d,hn
	   exit(1)
     myLegend.AddEntry(h,"data","P")
#   print "Lumi tot", lumitot

   for gr in backgroundSorted:
     for b in background[gr]: 
      nevents=totevents(b)
      if f[b] :
	h=f[b].Get(hn)
	if h :
	   h.Scale(samples[b]["xsec"]/nevents*lumitot)
	   setHistoStyle (h, gr) 
#	   dprint "adding", b, "to", hn 
#	   i+=1
	   if hn not in histosum :
		histosum[hn]=h.Clone()
	        histoTH[hn]=h.Clone()   
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
		histoTH[hn].Add(h)
   		for sy in systematicsToPlot :
		  hs=f[b].Get(findSyst(hn,sy,f[d]))
		  if hs:
	   	    hs.Scale(samples[b]["xsec"]/nevents*lumitot)
                    histosumSyst[hn][sy].Add(hs)
		  else :
		    #print "using unchanged histo"
                    histosumSyst[hn][sy].Add(h)
	   histos[hn].Add(h)
        else:
	   print "Cannot open",b,hn
	   exit(1)
     myLegend.AddEntry(h,gr,"f")

   histosSignal[hn]={} 
   for gr in signal:
     for b in signal[gr]:
      nevents=totevents(b)
      if f[b] :
        h=f[b].Get(hn)
        if h :
           h.Scale(samples[b]["xsec"]/nevents*lumitot)
#           print "adding", b, "to", hn
	   setHistoStyle (h, gr) 
 #          i+=1
           histosSignal[hn][b] = h.Clone()  	
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
		    #print "Not found",hn+"__syst__"+sy
                    histoSigsumSyst[hn][sy].Add(h)
           histosSig[hn].Add(h)
           histos[hn].Add(h)
           histoTH[hn].Add(h)               
        else:
	   print "Cannot open",b,hn
	   exit(1)
     myLegend.AddEntry(h,gr,"f")            
                                            
   for gr in signal:                        
     for b in signal[gr]:                   
        h=histosSignal[hn][b]               
        h.SetLineColor(linecolor[gr])       
        h.SetFillStyle(0)                   
        h.SetLineWidth(3)                   
        h.SetLineStyle(2)                   
        h.Scale(20.)                        
        myLegend.AddEntry(h,gr+" x20","l")       
   firstBlind=100000
   lastBlind=-1
   for i in range(histosSig[hn].GetStack().Last().GetNbinsX()) :
	if histosSig[hn].GetStack().Last().GetBinContent(i) > 0.1*sqrt(abs(histos[hn].GetStack().Last().GetBinContent(i))) and histosSig[hn].GetStack().Last().GetBinContent(i)/(0.1+abs(histos[hn].GetStack().Last().GetBinContent(i))) > 0.05 :
		print "to blind",hn,i,abs(histos[hn].GetStack().Last().GetBinContent(i)), histosSig[hn].GetStack().Last().GetBinContent(i)	
	        if i < firstBlind:
		    firstBlind=i
                lastBlind=i
   for i in range(firstBlind,lastBlind) :
       datastack[hn].GetStack().Last().SetBinContent(i,0)
       datasum[hn].SetBinContent(i,0)
       print "blinded",i,hn
   myLegend.Draw() #NEW  
   canvas[hn].cd(1)
   histos[hn].SetTitle("") 
   datastack[hn].Draw("E P")
   datastack[hn].GetXaxis().SetTitle(hn)
   histos[hn].Draw("hist same")
#  histos[hn].Draw("hist")                                                               
   histoTH[hn].SetLineWidth(0)                                                           
   histoTH[hn].SetFillColor(ROOT.kBlack);                                                
   histoTH[hn].SetFillStyle(3004);                                                       
   #if "log" in hn : histoTH[hn].GetYaxis().SetRangeUser(0.1,histoTH[hn].GetMaximum()*5) 
   #else : histoTH[hn].GetYaxis().SetRangeUser(0.,histoTH[hn].GetMaximum()*1.2)          
   #histoTH[hn].GetYaxis().SetRangeUser(0.1,histoTH[hn].GetMaximum()*2)                  
   setStyle(histos[hn].GetHistogram())
   canvas[hn].Update()                
   histoTH[hn].Draw("same E2")        


   datastack[hn].Draw("E P same")
   for gr in signal:                                                     
     for b in signal[gr]: histosSignal[hn][b].Draw("hist same")          
                                                                         
   t0 = makeText(0.75,0.80,labelRegion[hn.split("___")[1]] if hn.split("___")[1] in labelRegion.keys() else hn.split("___")[1], 61)  
   t1 = makeText(0.15,0.95,"CMS", 61)                                                           
   t2 = makeText(0.30,0.95,str(year), 42)                                                       
   t3 = makeText(0.92,0.95,lumi%(lumitot/1000.)+" (13 TeV)", 42)                                
   t0.Draw()                                                             
   t1.Draw()                                                             
   t2.Draw()                                                             
   t3.Draw()                                     
   datastack[hn].GetHistogram().SetMarkerStyle(10)                       

   canvas[hn].Update()
   ratio=datasum[hn].Clone()
   ratio.Add(histosum[hn],-1.) 
   ratio.Divide(histosum[hn])
   ratio.SetMarkerStyle(10)

   canvas[hn].cd(2)
   setStyle(ratio, isRatio=True)

#   ratio.SetLabelSize(datastack[hn].GetHistogram().GetLabelSize()*3)
#   ratio.GetYaxis().SetLabelSize(datastack[hn].GetHistogram().GetLabelSize()*3)
   ratio.Draw()
   ratioError = makeRatioMCplot(histoTH[hn])  
   ratioError.Draw("same E2")                 

   ratio.SetAxisRange(-0.5,0.5,"Y")
   ratio.GetYaxis().SetNdivisions(5)
   ratiosy=[]
   for j,sy in enumerate(systematicsToPlot):
       ratiosy.append(histosumSyst[hn][sy].Clone())
       ratiosy[-1].Add(histosum[hn],-1.)
       ratiosy[-1].Divide(histosum[hn])
       ratiosy[-1].SetLineColor(1+j)
       #ratiosy[-1].SetLineStyle(j)
       ratiosy[-1].SetFillStyle(0)
       ratiosy[-1].Draw("same hist")
       print "Heu",hn,sy,histosumSyst[hn][sy].Integral(),histosum[hn].Integral(),lumitot,ratiosy[-1]
#   systematics=[x for x in histoNames if x[:hn.find("___")]==hn[:hn.find("___")] and "__syst__" in x]
#   print "available systematics",hn,systematics
#  for s in systematics:
#  getsum up , down
   canvas[hn].GetPad(2).SetGridy()
   canvas[hn].SaveAs("figure/2018/Z/%s.png"%hn)	   
   #canvas[hn].SaveAs("%s.root"%hn)	   
   canvas[hn].GetPad(1).SetLogy(True)
   canvas[hn].SaveAs("figure/2018/Z/%s_log.png"%hn)	   


his=[x for x in histoNames if "__syst__" not in x]
print his[0]
makeplot(his[0]) #do once for caching normalizations

print "Preload"
for ff in f:
   for h in histoNames :
     f[ff].Get(h)
print "Preload-done"

if True:
 from multiprocessing import Pool
 runpool = Pool(20)
 #toproc=[(x,y,i) for y in sams for i,x in enumerate(samples[y]["files"])]
 runpool.map(makeplot, his[1:])
else :
 for x in his[1:] :
    makeplot(x)

tot=0
for s in totevCount:
  tot+=totevSkim[s]

print tot, "input events" 
