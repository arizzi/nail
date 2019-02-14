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
	self.histos={}
	self.regexps=[]
	self.validCols=[x[0] for x in cols]
	self.inputTypes={x[0]:x[1] for x in cols}
#	print self.inputTypes
	for c,t in cols:
	    self.inputs[c]=[]
	    self.selections[c]=[]
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_]+)\]","makeP4(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","makeP4(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
#	print 'TLorentzVector makeP4(float pt,float eta,float phi,float m) { TLorentzVector r; r.SetPtEtaPhiM(pt,eta,phi,m); return r;}'

	print "gSystem->Load(\"libGenVector.so\")"
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_\[\]]+)\]","ROOT::Math::PtEtaPhiMVector(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","ROOT::Math::PtEtaPhiMVector(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.AddCodeRegex(("@p4v\(([a-zA-Z0-9_]+)\)","vector_map_t<ROOT::Math::PtEtaPhiMVector>(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
    def AddCodeRegex(self,regexp):
	self.regexps.append(regexp)

    def DefaultConfig(self,**kwargs):
	self.conf.update(kwargs)
	for k in kwargs.keys() :
	    self.Define(k,"%s"%(kwargs[k]))

    def SubCollection(self,name,existing,sel=""):
	if sel!="":
		self.Define(name,sel)
	else :
		sel=name
	l=len(existing)
	additionalCols= [ (name+c[l:],c) for c in self.validCols  if c[0:l+1]==existing+"_" ]
	for (ac,oc) in additionalCols :	
	   if oc in self.inputTypes and self.inputTypes[oc] =='Bool_t' :
	       self.Define(ac,"(1*%s)[%s]"%(oc,name)) #FIX RDF BUG
	   else:
	       self.Define(ac,"%s[%s]"%(oc,name))
	self.Define("n%s"%name,"Sum(%s)"%(name))

    def SubCollectionFromIndices(self,name,existing,sel=""):
        if sel!="":
                self.Define(name,sel)
        else :
                sel=name
        l=len(existing)
        additionalCols= [ (name+c[l:],c) for c in self.validCols  if c[0:l+1]==existing+"_" ]
        for (ac,oc) in additionalCols :
            self.Define(ac,"Take(%s,%s)"%(oc,name))
        self.Define("n%s"%name,"%s.size()"%(name))

    def Distinct(self,name,collection,selection="") :
	if collection not in self.validCols :
	   if "n%s"%collection in self.validCols:
		self.Define(collection,"RVec<bool>(n%s,true)"%collection)
	   else :
		print "Cannot find collection",collection
		return
	self.Define("%s_allpairs"%name,"Combinations(Nonzero(%s),Nonzero(%s))"%(collection,collection))
	self.Define(name,"%s_allpairs[0] > %s_allpairs[1]"%(name,name))
	self.Define("%s0"%name,"%s_allpairs[0][%s]"%(name,name))
	self.Define("%s1"%name,"%s_allpairs[1][%s]"%(name,name))
	self.SubCollectionFromIndices("%s0"%name,collection)
	self.SubCollectionFromIndices("%s1"%name,collection)

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

    def Selection(self,name,code,inputs=[]) :
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

    def Histo(self,name,binHint=None):
	self.histos[name]={}

	

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

    def printRDF(self,to,optimizeFilters=False):
	#print 'ROOT::RDataFrame rdf("Events","/gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root");'
	print 'ROOT::RDataFrame rdf("Events","/dev/shm/VBF_HToMuMu_nano2016.root");'
	print "auto toplevel ="
	rdf="rdf"
	if optimizeFilters :
	    allfilters=set([s for t in to for s in self.selections[t]])
	    allfiltersinputs=set([i for t in to for s in self.selections[t] for i in self.allNodesTo(s)])
            filterstring="||".join( ["(%s)"%("&&".join(self.selections[t])) for t in to])
      	    orderedColumns=[x for x in self.validCols if x in allfiltersinputs] + [x for x in self.validCols if x not in allfiltersinputs]
	    commonfilters=set(self.selections[to[0]]) if len(to) > 0 else set()
	    for s in to[1:]:
	        commonfilters.intersection_update(self.selections[s])
	
	#or t in to :
	#  sels=self.selections[t]
	#  filters.add('&&'.join(self.selections[t])
	
	toprint=set([x for t in to for x in self.allNodesTo(t)])
	cols=self.validCols if not optimizeFilters else orderedColumns
	for c in cols :
           if c in toprint:
	    if c in self.obs or c in self.filters :
	        print '%s.Define("%s","%s")'%(rdf,c,self.code[c])
		if optimizeFilters :
		    if c in allfilters :
			allfilters.remove(c)
		    if len(allfilters) == 0 and filterstring!="":
			print '%s.Filter("%s")'%(rdf,filterstring)
			filterstring=""
		    if c in commonfilters :
		        print '%s.Filter("%s")'%(rdf,c)
		
		rdf=""
	print ";"
        for t in to :
	    if len(self.selections[t]) > 0 :
	  	    print 'auto %s=toplevel.Filter("%s").Histo1D("%s");'%(t,'&&'.join(self.selections[t]),t)
	    else:
	  	    print 'auto %s=toplevel.Histo1D("%s");'%(t,t)

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
	 affected.sort(key=lambda x: self.validCols.index(x)) #keep original sorting
         replacementTable=[(x,x+"__syst__"+name) for x in affected]
         for x,x_syst in replacementTable:

             ncode=" "+self.code[x]+" "  #FIXME: we should avoid duplicating the code
             for y,y_syst in replacementTable:
	         regBound="([^a-zA-Z0-9_])"
                 reg=regBound+y+regBound
                 ncode=re.sub(reg,"\\1"+y_syst+"\\2",ncode)
             if x in self.obs:
                 selections=[]
                 for s in self.selections[x] :
                     if s in affected :
                        selections.append(s+"__syst__"+name)
                     else :
                        selections.append(s)
                 self.Define(x_syst,ncode,requires=selections)
             if x in self.filters:
                 self.Selection(x_syst,ncode)



    




class AnalysisYields:
    def __init__(self):
	self.hisots={}
	self.counters={}
	self.tuples={}
	self.fillers={}



    

class Interperations:
    def __init__(self):
	self.fits={}
	self.tables={}
  
    def Compare(self) :
	pass

    def Fit():
	pass
