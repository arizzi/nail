import clang.cindex
import re
import sys
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

    def __init__(self,name,cols,dftypes):
	self.name=name
        self.obs={} 
        self.filters={} 
        self.code={} 
        self.inputs={} 
        self.selections={} 
	self.variations={}
	self.conf={}
	self.histos={}
	self.weights={}
	self.defweights=[]
	self.regexps=[]
	self.validCols=[x[0] for x in cols]
	self.inputTypes={x[0]:x[1] for x in cols}
	self.dftypes=dftypes
	print >> sys.stderr, "Start"
#	print self.inputTypes
	for c,t in cols:
	    self.inputs[c]=[]
	    self.selections[c]=[]
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_]+)\]","makeP4(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","makeP4(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
#	print 'TLorentzVector makeP4(float pt,float eta,float phi,float m) { TLorentzVector r; r.SetPtEtaPhiM(pt,eta,phi,m); return r;}'

#	print "gSystem->Load(\"libGenVector.so\")"
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_\[\]]+)\]","ROOT::Math::PtEtaPhiMVector(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","ROOT::Math::PtEtaPhiMVector(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.AddCodeRegex(("@p4v\(([a-zA-Z0-9_]+)\)","vector_map_t<ROOT::Math::PtEtaPhiMVector>(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.nodesto={}
	self.Define("defaultWeight","1.")

    def AddCodeRegex(self,regexp):
	self.regexps.append(regexp)

    def DefaultConfig(self,**kwargs):
	self.conf.update(kwargs)
	for k in kwargs.keys() :
	    self.Define(k,"%s"%(kwargs[k]))

    def SubCollection(self,name,existing,sel="",requires=[],singleton=False):
	if sel!="":
		self.Define(name,sel,requires=requires)
	else :
		sel=name
	l=len(existing)
	additionalCols= [ (name+c[l:],c) for c in self.validCols  if c[0:l+1]==existing+"_" ]
	for (ac,oc) in additionalCols :	
	   if oc in self.inputTypes and self.inputTypes[oc] =='Bool_t' :
	       self.Define(ac,"(1*%s)[%s]"%(oc,name),requires=requires) #FIX RDF BUG
	   else:
	       self.Define(ac,"%s[%s]"%(oc,name),requires=requires)
	if not singleton:
		self.Define("n%s"%name,"Sum(%s)"%(name),requires=requires)

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

    def ObjectAt(self,name,existing,index="",requires=[]):
	self.SubCollection(name,existing,index,requires,True)

    def Distinct(self,name,collection,selection="") :
	if collection not in self.validCols :
	   if "n%s"%collection in self.validCols:
		self.Define(collection,"ROOT::VecOps::RVec<unsigned int>(n%s,true)"%collection)
	   else :
		print "Cannot find collection",collection
		return
	self.Define("%s_allpairs"%name,"Combinations(Nonzero(%s),Nonzero(%s))"%(collection,collection))
	self.Define(name,"%s_allpairs[0] > %s_allpairs[1]"%(name,name))
	self.Define("%s0"%name,"%s_allpairs[0][%s]"%(name,name))
	self.Define("%s1"%name,"%s_allpairs[1][%s]"%(name,name))
	self.SubCollectionFromIndices("%s0"%name,collection)
	self.SubCollectionFromIndices("%s1"%name,collection)

    def TakePair(self,name,existing,pairs,index,requires=[]):
	self.Define("%s_indices"%(name),index,requires=requires)
	for i in [0,1]:
	    self.ObjectAt("%s%s"%(name,i), existing,"int(%s%s[%s_indices])"%(pairs,i,name),requires=requires)


    def Define(self,name,code,inputs=[],requires=[]):
	print >> sys.stderr, name
	if name not in self.validCols or name=="defaultWeight":
	    if name == "defaultWeight" and name in self.validCols:
		self.validCols.remove(name)
	    self.validCols.append(name)
	    pcode=self.preprocess(code)
            self.obs[name]={}
            self.code[name]=pcode
            self.inputs[name]=list(set(self.findCols(pcode)+inputs))
            self.selections[name]=list(set(requires+[y for x in self.inputs[name] if x in self.selections for y in self.selections[x]]))
	else :
	    print "Attempt to redefine column", name," => noop"

    def Selection(self,name,code,inputs=[]) :
        if name not in self.validCols :
            self.validCols.append(name)
	    pcode=self.preprocess(code)
	    self.filters[name]={}
            self.code[name]=pcode
            self.inputs[name]=list(set(self.findCols(pcode)))
            self.selections[name]=list(set([y for x in self.inputs[name] if x in self.selections for y in self.selections[x]]))

	else :
	    print "Attempt to redefine column", name," => noop"

    def Systematic(self,name,original,modified, exceptions=[]):
	self.Variation(name,original,modified, exceptions)

    def Variation(self,name,original,modified,exceptions=[]):
	self.variations[name]={}
        self.variations[name]["original"]=original
        self.variations[name]["modified"]=modified
	self.variations[name]["exceptions"]=exceptions

    def AddDefaultWeight(self,name):
	if not name in self.defweights:
	    self.defweights.append(name)
	    #update defaultWeight string
	    self.Define("defaultWeight","*".join(self.defweights))

    def AddWeight(self,name,code="",filt=lambda x,y: True,nodefault=False):
	if code:
            if nodefault :
		self.Define(name,"%s"%code)
	    else:
		self.Define(name,"(%s)*defaultWeight"%code)
	self.weights[name]={}
	self.weights[name]["filter"]=filt

    def AddWeightArray(self,name,nentries=1,filt=lambda x,y: True):
	for i in range(nentries):
		self.AddWeight("%s%s"%(name,i),"%s[%s]"%(name,i),filt=filt)

    def Histo(self,name,binHint=None):
	self.histos[name]={}
	self.histos["bin"]=binHint

    

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

    def printRDFCpp(self,to,optimizeFilters=False):

        toprint=set([x for t in to for x in self.allNodesTo(t)])
        cols=self.validCols if not optimizeFilters else orderedColumns
        for c in cols :
           if c in toprint:
            if c in self.obs or c in self.filters :
		inputs=""
		for i in self.inputs[c]:
		   if inputs!="":
		       inputs+=", "
		   inputs+="const %s & %s"%(self.dftypes[i],i)
		#print "auto func__%s = [](%s) { return %s; };" %(c,inputs,self.code[c])
		print "auto func__%s(%s) { return %s; }" %(c,inputs,self.code[c])
		print "using type__%s = ROOT::TypeTraits::CallableTraits<decltype(func__%s)>::ret_type;"%(c,c)
		self.dftypes[c]="type__%s"%(c)

	print '''
int main(int argc, char** argv)
{
   auto n_cores = 0;
   if (argc > 1)
      n_cores = std::atoi(argv[1]);
   if (n_cores > 0)
      ROOT::EnableImplicitMT(n_cores);

'''
        print 'ROOT::RDataFrame rdf("Events","/dev/shm/VBF_HToMuMu_nano2016.root");'
        rdf="rdf"

        print "auto toplevel ="
        for c in cols :
           if c in toprint:
            if c in self.obs or c in self.filters :
                print '%s.Define("%s",func__%s,{%s})'%(rdf,c,c,",".join([ '"%s"'%x for x in self.inputs[c]]))
                rdf=""
        print ";"

        for t in to :
	    if t in self.histos:
                if len(self.selections[t]) > 0 :
                    print 'auto %s_neededselection=toplevel.Filter("%s");'%(t,'&&'.join(self.selections[t]))
		    rdf="%s_neededselection"%t
                else:
		    rdf="toplevel"
                print 'auto %s=%s.Histo1D("%s","defaultWeight");'%(t,rdf,t)
		for w in self.weights :
		    if self.weights[w]["filter"](t,w):
	                    print 'auto %s__weight__%s=%s.Histo1D("%s","%s");'%(t,w,rdf,t,w)

	print '''
   TStopwatch s;
   SBClassifier.OnPartialResult(0, [&s](TH1D &) {
      std::cout << "starting event loop" << std::endl;
      s.Start();
      return true;
   });
   std::cout << "now jitting+event loop..." << std::endl;
   *SBClassifier;
   s.Stop();
   std::cout << "elapsed time: " << s.RealTime() << "s" << std::endl;
   return 0;
}
'''




    def baseInputs(self,x) :
        if len(self.inputs[x]) == 0 :
          return [x]
        else :
          ret=[]
          for i in self.inputs[x] :
             ret.extend(self.baseInputs(i))
          return ret
    
    def allNodesTo(self,x) :   
	  if x in self.nodesto :
		return self.nodesto[x]
	  print >> sys.stderr, "Nodes to ",x
          ret=set([x])
          for i in self.inputs[x] :
             ret.update(self.allNodesTo(i))
          for i in self.selections[x] :
             ret.update(self.allNodesTo(i))
	  #this is cached because it cannot change
	  self.nodesto[x]=ret;
          return ret
    def allNodesFromWithWhiteList(self,x,wl,cache={}) :
          #cache is local because it changes with more defines
          if x in cache :
             return cache[x]
          print >> sys.stderr, "Nodes from ",x
          children=set([n for n in self.inputs.keys() if ((x in self.inputs[n]+self.selections[n]) and (n in wl))])
          ret=set(children)
          for c in children :
             ret.update(self.allNodesFromWithWhiteList(c,wl,cache))
          cache[x]=ret
          return ret


    def allNodesFrom(self,x,cache={}) :   
	  #cache is local because it changes with more defines
	  if x in cache :
	     return cache[x]
	  print >> sys.stderr, "Nodes from ",x
          children=set([n for n in self.inputs.keys() if x in self.inputs[n]+self.selections[n]])
	  ret=set(children)
	  for c in children :
             ret.update(self.allNodesFrom(c,cache))
	  cache[x]=ret
          return ret

    def findAffectedNodesForVariationOnTarget(self,name,target):
	 nodesTo= [x for x in self.allNodesTo(target) if x not in self.variations[name]["exceptions"] ]
	 return [x for x in self.allNodesFromWithWhiteList(self.variations[name]["original"],nodesTo) if x in nodesTo]
   
    def createVariationBranch(self,name,target):
	 print >> sys.stderr, "Find affected"
         affected=(self.findAffectedNodesForVariationOnTarget(name,target))
	 print >> sys.stderr, "Found affected"
	 affected.sort(key=lambda x: self.validCols.index(x)) #keep original sorting
         replacementTable=[(x,x+"__syst__"+name) for x in affected]
         for x,x_syst in replacementTable:

             ncode=" "+self.code[x]+" "  #FIXME: we should avoid duplicating the code
             for y,y_syst in replacementTable+[(self.variations[name]["original"],self.variations[name]["modified"])]:
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
