import clang.cindex
import re
class AggregatedSample:
    def __init__(self,*args):
      self.samples=list(args)

class Sample:
    def __init__(self,filename,treename):
	self.filename=filename
	self.treename=treename
    def __add__(self,other):
	return AggregatedSamples(self,other)


class SampleProcessing:

    def __init__(self,name,cols):
	self.name=name
        self.obs={} 
        self.filters={} 
        self.code={} 
        self.inputs={} 
        self.selections={} 
	self.syst={}
	self.conf={}
	self.regexps=[]
	self.validCols=cols
	for c in cols:
	    self.inputs[c]=[]
	    self.selections[c]=[]
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_]+)\]","ROOT::Math::PtEtaPhiMVector(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","ROOT::Math::PtEtaPhiMVector(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))

    def AddCodeRegex(self,regexp):
	self.regexps.append(regexp)

    def DefaultConfig(self,**kwargs):
	self.conf.update(kwargs)
	for k in kwargs.keys() :
	    self.Define(k,"%s"%(kwargs[k]))

    def SubCollection(self,name,existing,sel):
	self.Define(name,sel)
	l=len(existing)
	additionalCols= [ (name+c[l:],c) for c in self.validCols  if c[0:l+1]==existing+"_" ]
	for (ac,oc) in additionalCols :	
	   self.Define(ac,"%s[%s]"%(oc,name))

    def Define(self,name,code,inputs=[],requires=[]):
	if name not in self.validCols :
  	    self.validCols.append(name)
	    pcode=self.preprocess(code)
            self.obs[name]={}
            self.code[name]=pcode
            self.inputs[name]=self.findCols(pcode)+inputs
            self.selections[name]=list(set(requires+[y for x in self.inputs[name] if x in self.selections for y in self.selections[x]]))
	else :
	    print "Attempt to redefine column", name," => noop"

    def Filter(self,name,code,inputs=[]) :
        if name not in self.validCols :
            self.validCols.append(name)
	    pcode=self.preprocess(code)
	    self.filters[name]={}
            self.code[name]=pcode
            self.inputs[name]=self.findCols(pcode)
            self.selections[name]=list(set([y for x in self.inputs[name] if x in self.selections for y in self.selections[x]]))

	else :
	    print "Attempt to redefine column", name," => noop"

    def Systematic(self,name,original,modified): 
	self.syst[name]={}
        self.syst[name]["original"]=original
        self.syst[name]["modified"]=modified

	

    def findCols(self,code) :
	idx = clang.cindex.Index.create()
	tu = idx.parse('tmp.cpp', args=['-std=c++11'], unsaved_files=[('tmp.cpp', code)],  options=0)
	identifiers=set()
	for t in tu.get_tokens(extent=tu.cursor.extent):
	   if t.kind==clang.cindex.TokenKind.IDENTIFIER :
	     if t.spelling in self.validCols:
	            identifiers.add(t.spelling)

	ret=[]
	regBound="[^a-zA-Z0-9_]"
	for c in identifiers:
	    reg=regBound+c+regBound
	    if re.search(reg," "+code+" ") :
		ret.append(c)
	#print "In: #######\n",code
	#print "####\nFound",ret
	return ret

    def preprocess(self,code) :
        for s,r in self.regexps :
	    code=re.sub(s,r,code)
	return code

    def printRDF(self):
	print 'ROOT::RDataFrame toplevel_("Events","/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root");'
#print "ROOT::RDataFrame d(treeName, fileName);"
	definedSels=["toplevel"]
	last={"toplevel":"toplevel"}
	for c in self.validCols:
	    if(len(self.selections[c])==0):
		rdf="toplevel";
	    else :
   		rdf="_".join(sorted(self.selections[c]))
	    if rdf not in definedSels:
		print "\nauto node_%s=toplevel"%rdf
		for f in sorted(self.selections[c]) :
		    print '.Filter("%s","%s")'%(f,self.code[f])
		print ";"
	    	definedSels.append(rdf)
		last[rdf]=rdf
	    if c in self.obs:
		print 'auto node_%s=node_%s.Define("%s","%s");'%(c,last[rdf],c,self.code[c])
		last[rdf]=c
	    if c in self.filters:
		print '\nauto node_%s=node_%s.Filter("%s","%s");'%(c,last[rdf],c,self.code[c])
		last[c]=c
	    	definedSels.append(c)



    def baseInputs(self,x) :
        if len(self.inputs[x]) == 0 :
          return [x]
        else :
          ret=[]
          for i in self.inputs[x] :
             ret.extend(self.baseInputs(i))
          return ret
    
    def allNodesTo(self,x) :   
          ret=[x]
          for i in self.inputs[x] :
             ret.extend(self.allNodesTo(i))
          for i in self.selections[x] :
             ret.extend(self.allNodesTo(i))
          return list(set(ret))

    def allNodesFrom(self,x) :   
          children=[n for n in self.inputs.keys() if x in self.inputs[n]+self.selections[n]]
	  ret=children
	  for c in children :
             ret.extend(self.allNodesFrom(c))
          return list(set(ret))

    def findAffectedNodesForSystematicOnTarget(self,name,target):
	 return [x for x in self.allNodesFrom(self.syst[name]["original"]) if x in self.allNodesTo(target)]
    
    def createSystematicBranch(self,name,target):
	 affected=(self.findAffectedNodesForSystematicOnTarget(name,target))
	 for x in affected:
	     x_syst=x+"__syst__"+name
	     if x in self.obs:
		 selections=[]
		 for s in self.selections[x] :
		     if s in affected :	
			selections.append(s+"__syst__"+name)
		     else :
			selections.append(s)
	         self.Define(x_syst,self.code[x],requires=selections) 
	     if x in self.filters:
	         self.Filter(x_syst,self.code[x]) 
                 
	     #FIXME: not clear what to do here, the code should be the same just bound to different inputs...so:
	     #  - clone the code?
             #  - provide translation table?
	     # should also adjust the "requires" and detect filter vs define	



    




class AnalysisYields:
    def __init__(self):
	self.hisots={}
	self.counters={}
	self.tuples={}
	self.fillers={}


    def Histo(self,name,what,binHint=None):
	pass

class Interperations:
    def __init__(self):
	self.fits={}
	self.tables={}
  
    def Compare(self) :
	pass

    def Fit():
	pass
