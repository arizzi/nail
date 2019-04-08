import ROOT
import clang.cindex
import re
import sys
import copy
from clang.cindex import CursorKind
from clang.cindex import Index
from clang.cindex import TypeKind
headerstring='''
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
#include "helpers.h"
#define MemberMap(vector,member) Map(vector,[](auto x){return x.member;})
#define P4DELTAR ROOT::Math::VectorUtil::DeltaR<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>,ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>> 
//ROOT::Math::PtEtaPhiMVector,ROOT::Math::PtEtaPhiMVector> 
'''


ROOT.gInterpreter.Declare(headerstring)
def returnType(code,function):
   myFuncCode = code
   ROOT.gInterpreter.Declare(myFuncCode)
   getTypeNameF = "auto %s_typestring=ROOT::Internal::RDF::TypeID2TypeName(typeid(ROOT::TypeTraits::CallableTraits<decltype(%s)>::ret_type));"%(function,function)
   ROOT.gInterpreter.ProcessLine(getTypeNameF)
   return getattr(ROOT,"%s_typestring"%function)

whatever='''  idx = clang.cindex.Index.create()
  tu = idx.parse('tmp.cpp', args=['-std=c++14'], unsaved_files=[('tmp.cpp', headerstring+code)],  options=0)
  for cursor in tu.cursor.get_children():
    # Ignore AST elements not from the main source file (e.g.
    # from included files).
    if not cursor.location.file or cursor.location.file.name != 'tmp.cpp':
        continue

    # Ignore AST elements not a function declaration.
    if cursor.kind != CursorKind.FUNCTION_DECL:
        continue

    # Obtain the return Type for this function.
    result_type = cursor.type.get_result()
    print "function", cursor.spelling,function
    if cursor.spelling == function:
	return result_type.kind.spelling
  return "void"
'''

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

    def __init__(self,name,filename) :
        f=ROOT.TFile.Open(filename)
	e=f.Get("Events")
	allbranches=[(x.GetName(),x.GetListOfLeaves()[0].GetTypeName()) for x in e.GetListOfBranches()]
	df=ROOT.RDataFrame("Events",filename)
	dftypes={x[0]:df.GetColumnType(x[0]) for x in allbranches}
        self.init(name,allbranches,dftypes) 	
	self.defFN=filename
	self.dupcode=False
    def init(self,name,cols,dftypes):
	self.name=name
        self.lazyParse=True
	self.cols={}
        self.obs={} 
        self.selections={} 
        self.code={} 
        self.inputs={} 
        self.explicitInputs={} 
        self.originals={} 
        self.requirements={} 
        self.explicitRequirements={} 
	self.variations={}
	self.conf={}
	self.histos={}
	self.variationWeights={} #variation weights
	self.centralWeights={"":[]} #default weights for each selection
	self.regexps=[]
	self.originalCols=cols
	self.validCols=[x[0] for x in cols]
	self.inputTypes={x[0]:x[1] for x in cols}
	self.dftypes=dftypes
	print >> sys.stderr, "Start"
#	print self.inputTypes
	for c,t in cols:
	    self.inputs[c]=[]
	    self.requirements[c]=[]
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_]+)\]","makeP4(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
#	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","makeP4(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
#	print 'TLorentzVector makeP4(float pt,float eta,float phi,float m) { TLorentzVector r; r.SetPtEtaPhiM(pt,eta,phi,m); return r;}'

#	print "gSystem->Load(\"libGenVector.so\")"
	#self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_\[\]]+)\]","ROOT::Math::PtEtaPhiMVector(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
	#self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","ROOT::Math::PtEtaPhiMVector(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	#self.AddCodeRegex(("@p4v\(([a-zA-Z0-9_]+)\)","vector_map_t<ROOT::Math::PtEtaPhiMVector>(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)\[([a-zA-Z0-9_\[\]]+)\]","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >>(\\1_pt[\\2] , \\1_eta[\\2], \\1_phi[\\2], \\1_mass[\\2])"))
	self.AddCodeRegex(("@p4\(([a-zA-Z0-9_]+)\)","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >>(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.AddCodeRegex(("@p4v\(([a-zA-Z0-9_]+)\)","vector_map_t<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >	>(\\1_pt , \\1_eta, \\1_phi, \\1_mass)"))
	self.nodesto={}
#	self.Define("defaultWeight","1.")

    def Inputs(self,colName):
	if colName not in self.inputs:
#	    print "parse! ", colName
            self.inputs[colName]=sorted(set(self.findCols(self.code[colName],colName)+self.explicitInputs[colName]))
        return self.inputs[colName]

    def Requirements(self,colName):
	if colName not in self.requirements:
            self.requirements[colName]=self.sortedUniqueColumns(self.explicitRequirements[colName]+[y for x in self.Inputs(colName) for y in self.Requirements(x)])
        return self.requirements[colName]

	
    def AddCodeRegex(self,regexp):
	self.regexps.append(regexp)

    def DefaultConfig(self,**kwargs):
	self.conf.update(kwargs)
	for k in kwargs.keys() :
	    self.Define(k,"%s"%(kwargs[k]))

    def colsForObject(self,name):
	l=len(name)
	return [ c[l+1:] for c in self.validCols  if c[0:l+1]==name+"_" ]

    def MergeCollections(self,name,collections,requires=[]):
 	commonCols=set(self.colsForObject(collections[0]))
	#print commonCols
	for coll in collections[1:] :
	    commonCols=commonCols.intersection(self.colsForObject(coll))
	self.Define("n"+name,"+".join(["n"+coll for coll in collections]))
	for ac in commonCols :
	    self.Define(name+"_"+ac,"Concat(%s)"%(",".join([c+"_"+ac for c in collections]) ))

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

    def SubCollectionFromIndices(self,name,existing,sel="",requires=[]):
        if sel!="":
                self.Define(name,sel)
        else :
                sel=name
        l=len(existing)
        additionalCols= [ (name+c[l:],c) for c in self.validCols  if c[0:l+1]==existing+"_" ]
        for (ac,oc) in additionalCols :
            self.Define(ac,"Take(%s,%s)"%(oc,name),requires=requires)
        self.Define("n%s"%name,"int(%s.size())"%(name),requires=requires) #FIXME: cast to unsigned due to ROOT problem with uint leaf

    def ObjectAt(self,name,existing,index="",requires=[]):
	self.SubCollection(name,existing,index,requires,True)

    def Distinct(self,name,collection,selection="",requires=[]) :
	if collection not in self.validCols :
	   if "n%s"%collection in self.validCols:
		self.Define(collection,"ROOT::VecOps::RVec<unsigned int>(n%s,true)"%collection,requires=requires)
	   else :
		print "Cannot find collection",collection
		return
	self.Define("%s_allpairs"%name,"Combinations(Nonzero(%s),Nonzero(%s))"%(collection,collection),requires=requires)
	self.Define(name,"%s_allpairs[0] < %s_allpairs[1]"%(name,name),requires=requires)
	self.Define("%s0"%name,"%s_allpairs[0][%s]"%(name,name),requires=requires)
	self.Define("%s1"%name,"%s_allpairs[1][%s]"%(name,name),requires=requires)
	self.SubCollectionFromIndices("%s0"%name,collection,requires=requires)
	self.SubCollectionFromIndices("%s1"%name,collection,requires=requires)

    def TakePair(self,name,existing,pairs,index,requires=[]):
	self.Define("%s_indices"%(name),index,requires=requires)
	for i in [0,1]:
	    self.ObjectAt("%s%s"%(name,i), existing,"int(%s%s[%s_indices])"%(pairs,i,name),requires=requires)

#    def DefineWithWildCards(self,name,function,inputs,requires=[])
	# need for example for HLT bits
	#efineWithWildCars("HLTMuon","multiOR"
#	print "Not implemented"
#	pass

    def Define(self,name,code,inputs=[],requires=[],original=""):
#	print >> sys.stderr, name
	if name not in self.validCols :
#	if name not in self.validCols or name[:len(defaultWeight)]=="defaultWeight":
#	    if name[:len(defaultWeight)] == "defaultWeight" and name in self.validCols:
#		self.validCols.remove(name)
	    self.validCols.append(name)
	    pcode=self.preprocess(code)
            self.obs[name]={}
            self.originals[name]=original
            self.code[name]=pcode
	    self.explicitInputs[name]=inputs
	    self.explicitRequirements[name]=requires
	    if not self.lazyParse:
                self.inputs[name]=sorted(list(set(self.findCols(pcode,name)+inputs)))
                self.requirements[name]=list(set(requires+[y for x in self.inputs[name] if x in self.requirements for y in self.requirements[x]]))

	else :
	    print "Attempt to redefine column", name," => noop"

    def Selection(self,name,code,inputs=[],requires=[],original="") :
        if name not in self.validCols :
            self.validCols.append(name)
	    pcode=self.preprocess(code)
	    self.selections[name]={}
            self.originals[name]=original
            self.code[name]=pcode
            self.explicitInputs[name]=inputs
            self.explicitRequirements[name]=requires
            if not self.lazyParse:
                self.inputs[name]=sorted(list(set(self.findCols(pcode,name)+inputs)))
                self.requirements[name]=list(set(requires+[y for x in self.inputs[name] if x in self.requirements for y in self.requirements[x]]))

	else :
	    print "Attempt to redefine column", name," => noop"

    def Match(self,name1,name2,metric="P4DELTAR",metricName="Dr",embed=([],[]),defIdx=-1,defVal=-99):
        name=name1+name2+"Pair"
        self.Define(name,"Combinations(n%s,n%s)"%(name1,name2))
        self.Define("%s_%s"%(name,metricName),"vector_map(%s,Take(%s_p4,%s[0]),Take(%s_p4,%s[1]))"%(metric,name1,name,name2,name))
        self.Define("%s_%s%s"%(name1,name2,metricName),"matrix_map(n%s,n%s,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):%s;},%s_%s)"%(name1,name2,defVal,name,metricName))
        self.Define("%s_%s%s"%(name2,name1,metricName),"matrix_map(n%s,n%s,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):%s;},%s_%s)"%(name1,name2,defVal,name,metricName))
        self.Define("%s_%sIdx"%(name1,name2),"matrix_map(n%s,n%s,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):%s;},%s_%s)"%(name1,name2,defIdx,name,metricName))
        self.Define("%s_%sIdx"%(name2,name1),"matrix_map(n%s,n%s,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):%s;},%s_%s)"%(name1,name2,defIdx,name,metricName))

    def MatchDeltaR(self,name1,name2,embed=([],[]),defIdx=-1,defVal=-99):
        name=name1+name2+"Pair"
        self.Define(name,"Combinations(n%s,n%s)"%(name1,name2))
        self.Define("%s_dr"%name,"vector_map(P4DELTAR,Take(%s_p4,%s[0]),Take(%s_p4,%s[1]))"%(name1,name,name2,name))
        self.Define("%s_%sDr"%(name1,name2),"matrix_map(n%s,n%s,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):%s;},%s_dr)"%(name1,name2,defVal,name))
        self.Define("%s_%sDr"%(name2,name1),"matrix_map(n%s,n%s,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?(-Max(-v)):%s;},%s_dr)"%(name1,name2,defVal,name))
        self.Define("%s_%sIdx"%(name1,name2),"matrix_map(n%s,n%s,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):%s;},%s_dr)"%(name1,name2,defIdx,name))
        self.Define("%s_%sIdx"%(name2,name1),"matrix_map(n%s,n%s,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?Argmax(-v):%s;},%s_dr)"%(name1,name2,defIdx,name))
#FIXME: embedding is broken
#	for attr in embed[0] :
#            self.Define("%s_%s%s"%(name1,name2,attr.capitalize()),"matrix_map(n%s,n%s,1,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?%s_%s[Argmax(-v)]:%s;},%s_dr)"%(name1,name2,name2,attr,defVal,name))
#	for attr in embed[1] :
#            self.Define("%s_%s%s"%(name2,name1,attr.capitalize()),"matrix_map(n%s,n%s,0,[](const ROOT::VecOps::RVec<float> & v) {return v.size()>0?%s_%s[Argmax(-v)]:%s;},%s_dr)"%(name2,name1,name1,attr,defVal,name))

    def Systematic(self,name,original,modified, exceptions=[]):
	self.Variation(name,original,modified, exceptions)

    def Variation(self,name,original,modified,exceptions=[]):
	self.variations[name]={}
        self.variations[name]["original"]=original
        self.variations[name]["modified"]=modified
	self.variations[name]["exceptions"]=exceptions

#    def addInput(self,node,inp):
 #       if inp not in self.inputs[node]:
#	    self.inputs[node]=sorted(list(set(self.inputs[node]+[inp])))

    def CentralWeight(self,name,selections=[""]):
	for s in selections :
	   missing=[x for x in self.Requirements(name) if x not in (self.Requirements(s)+[s] if s != "" else []) ]
	   if len(missing) > 0 :
		print "Cannot add weight",name,"on selection",s,"because",name,"requires the following additional selections"
		print missing
		exit(1)
	   if s not in self.centralWeights:
		self.centralWeights[s]=[]
	   self.centralWeights[s].append(name)
#	   if s!="": #FIXME: seems a bug to me
#	       self.addInput(s,name)

    def VariationWeight(self,name,replacing="",filt=lambda sname,hname,wname : "__syst__" not in hname and "__syst__" not in sname):
        self.variationWeights[name]={}
        self.variationWeights[name]["replacing"]=replacing
        self.variationWeights[name]["filter"]=filt
	     
    def VariationWeightArray(self,name,nentries=1,replacing="",filt=lambda sname,hname,wname : "__syst__" not in hname and "__syst__" not in sname):
        for i in range(nentries):
		self.Define("%s%s"%(name,i),"%s[%s]"%(name,i))
                self.VariationWeight("%s%s"%(name,i),replacing,filt)

    def Histo(self,name,binHint=None):
	self.histos[name]={}
	self.histos["bin"]=binHint

    def recursiveGetWeights(self,sel):
#	print "Weights for",sel
#	if "__syst__" in sel :
#	    return self.recursiveGetWeights(sel[:sel.find("__syst__")])
	res=set(self.centralWeights[""])
	if sel in self.centralWeights :
	    res.update(self.centralWeights[sel])
	for dep in self.Requirements(sel) :
	    res.update(self.recursiveGetWeights(dep))
#	print res
	return res

    def sortedUniqueColumns(self,cols):
	return sorted(list(set(cols)),key=lambda c:self.validCols.index(c))
	
    def selSetName(self,sels):
	sels=self.sortedUniqueColumns(sels)
	return "_".join(sels)

    def replaceWeightWithVariation(self,selname,variation,weights,hist="" ):
	v=self.variationWeights[variation]
	if v["replacing"] == "" and v["filter"](selname,variation,hist) :
	   return list(weights)+[variation]
	res=[]
	for w in weights:
	   if w==v["replacing"] and v["filter"](selname,variation,hist) :
		res.append(variation)
	   else :
		res.append(w)
	return res


    def defineWeights(self,selectionsSets):
	res=[]
	for name,selections in selectionsSets.iteritems():
#    name=selSetName(selections)	   
	    print "define weight",name
	    weights=set()
	    if name=="" and "" in self.centralWeights:
	        weights.update(self.centralWeights[""])
	    for s in selections:
	        weights.update(self.recursiveGetWeights(s))
	    for variation in ["Central"]+self.variationWeights.keys():
	        replacedWeights=weights
		same=False
		if variation != "Central":
		    replacedWeights=self.replaceWeightWithVariation(name,variation,weights)
		  #  same=(set(replacedWeights)==weights)
	        if variation == "Central" or "__syst__" not in name :	
 		  if not replacedWeights:
	                self.Define("%sWeight__%s"%(name,variation),"1.")
		  else:	
	                self.Define("%sWeight__%s"%(name,variation),"*".join(replacedWeights))
		  res.append("%sWeight__%s"%(name,variation))
#		print "%sWeight__%s"%(name,variation)
	return res
   


    def findCols(self,code,colname) :
	idx = clang.cindex.Index.create()
	tu = idx.parse('tmp.cpp', args=['-std=c++11'], unsaved_files=[('tmp.cpp', code)],  options=0)
	identifiers=set()
        boundary=self.validCols.index(colname) if colname in self.validCols else -1
	for t in tu.get_tokens(extent=tu.cursor.extent):
	   if t.kind==clang.cindex.TokenKind.IDENTIFIER :
	     if t.spelling in self.validCols[:boundary]:
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

    def printRDFCpp(self,to,debug=False,outname="out.C",selections={},snap=[],snapsel=""):
	histos=[x for y in selections for x in selections[y]]
	to=list(set(to+histos+selections.keys()))
	sels={self.selSetName(self.Requirements(x)):self.sortedUniqueColumns(self.Requirements(x)) for x in to}
#	print "Update with",{x:self.sortedUniqueColumns(self.Requirements(x)+[x]) for x in selections}
	sels.update({x:self.sortedUniqueColumns(self.Requirements(x)+[x]) for x in selections})

##AR###	selections.append("")

#	print "Sels are",sels
	weights=self.defineWeights(sels)
#	print "Weights to print", weights
	f=open(outname,"w")
	ftxt=open(outname[:-2]+"-data.txt","w")
	f.write(headerstring)
        toprint=set(selections.keys()+[x for t in to+weights for x in self.allNodesTo([t])])
#	print "toprint:",toprint
        cols=self.validCols # if not optimizeFilters else orderedColumns
        for c in cols :
           if c in toprint:
#	    print '.... doing..',c
            if c in self.obs or c in self.selections :
	      if self.code[c] != "" :
		inputs=""
		for i in self.Inputs(c):
		   if inputs!="":
		       inputs+=", "
		   if i in self.dftypes :
			   inputs+="const %s & %s"%(self.dftypes[i],i)
		   else:
			   print i, self.dftypes, c
			   inputs+="const %s & %s"%(self.dftypes[i[:i.rfind("__syst__")]],i)
		#print "auto func__%s = [](%s) { return %s; };" %(c,inputs,self.code[c])
		debugcode=""
		if debug :
			debugcode='std::cout << "%s" << std::endl;\n'%c
		cppcode="auto func__%s(%s) { %s return %s; }\n" %(c,inputs,debugcode,self.code[c])
		f.write(cppcode)
		self.dftypes[c]=returnType(cppcode,"func__%s"%c) #"type__%s"%(c)
		#f.write("using type__%s = ROOT::TypeTraits::CallableTraits<decltype(func__%s)>::ret_type;\n"%(c,c))
		#self.dftypes[c]="type__%s"%(c)

	f.write('''
int main(int argc, char** argv)
{
   auto n_cores = 0;
   if (argc > 1)
      n_cores = std::atoi(argv[1]);
   if (n_cores > 0)
      ROOT::EnableImplicitMT(n_cores);
   std::string fname="%s";
   if(argc > 2)
     fname=argv[2];
   std::string out="out";
   if(argc > 3)
     out=argv[3];


'''%self.defFN)
        f.write('ROOT::RDataFrame rdf("Events",fname.c_str());\n')
        rdf="rdf"
	if debug:
	   rdf="rdf.Range(1000)"
	i=0
        f.write("auto rdf0 =")
        for c in cols :
           if c in toprint:
            if c in self.obs or c in self.selections :
		if self.code[c] != "" :
	               f.write('%s.Define("%s",func__%s,{%s})\n'%(rdf,c,c,",".join([ '"%s"'%x for x in self.Inputs(c)])))
		else:
	               f.write('%s.Define("%s",func__%s,{%s})\n'%(rdf,c,self.originals[c],",".join([ '"%s"'%x for x in self.Inputs(c)])))
		rdf="rdf%s"%i
		i+=1
#		if debug :
#		   f.write('std::cout << "%s" << std::endl;\n'%rdf);
#		f.write("auto rdf%s ="%i)
                rdf=""
        f.write(rdf+";\n")
        rdf="rdf0"
 	f.write("auto toplevel=%s;\n"%rdf)
	f.write("std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;\n")
	#rdflast=rdf
	rdflast="rdf0"
	selsprinted=[]
	selname=""
	f.write("{")
	for s in selections.keys() :
	    for t in selections[s] :
#        for t in to :
#	    if t in self.histos:
#	      for s in selections.keys():
	          if s=="":
		    if len(self.Requirements(t)) > 0 :
		      selname=self.selSetName(self.Requirements(t))
		      sellist=self.Requirements(t)
		    else :
		      selname=""  
		  else:
           	      missing=[x for x in self.Requirements(t) if x not in self.Requirements(s)+[s]  ]
		      if missing :
			 print "Cannot make ",t," in ",s," because the following selections are needed",missing
		 	 continue
		      selname=s
		      sellist=self.Requirements(s)+[s]
#		      print "making ",t, "in",selname
 
                  if selname!="" :
#                    print 'auto %s_neededselection=%s.Filter("%s");'%(t,rdf,'").Filter("'.join(self.requirements[t]))
		    if selname not in selsprinted :
#			print "Printing",selname
			f.write("}")
 		        f.write('auto selection_%s=%s.Filter("%s","%s")'%(selname,rdflast,sellist[0],sellist[0]))
 	 	        for ss in sellist[1:]:
                            f.write('.Filter("%s","%s")'%(ss,ss))
		        f.write(";\n")
			selsprinted.append(selname)
			f.write("{")
			f.write('loadHistograms(%s,"%s",histos);\n'%("selection_%s"%selname,"selection_%s"%selname))
		    rdf="selection_%s"%selname
                  else:
		    rdf=rdflast


		  ftxt.write('%s,%s,%s,2000,0,2000,%s,%sWeight__Central\n'%(rdf,t,t,t,selname))
                  f.write('histos.emplace_back(%s.Histo1D({"%s--%s", "%s {%s}", 2000, 0, 2000},"%s","%sWeight__Central"));\n'%(rdf,t,s,t,s,t,selname))
                  #f.write('histos.emplace_back(%s.Histo1D({"%s%s", "%s {%s}", 500, 0, 0},"%s","%sWeight__Central"));\n'%(rdf,t,s,t,s,t,selname))
		  for w in self.variationWeights :
		    if self.variationWeights[w]["filter"](selname,t,w):
			    ww="%sWeight__%s"%(selname,w)
	                    ftxt.write('%s,%s__weight__%s,%s,1000,0,100,%s,%s\n'%(rdf,t,w,t,t,ww))
	                    f.write('histos.emplace_back(%s.Histo1D({"%s%s__weight__%s", "%s {%s}", 1000, 0, 100},"%s","%s"));\n'%(rdf,t,s,w,t,s,t,ww))

	
	if snapsel=="":
		rdf=rdflast
	else :
		rdf="selection_%s"%snapsel

	f.write("}")
	f.write('auto snap=%s.Snapshot("ot", (out+"Snap.root").c_str(),{%s})\n;'%(rdf,",".join(['"%s"'%x for x in snap])))

	f.write('''

//   auto tr=selection_twoMuons_twoOppositeSignMuons_twoJets.Snapshot("ot", (out+"Snap.root").c_str(), {"nJet","SBClassifier","NSoft2","event","VBFRegion","nSoftActivityJet","SoftActivityJet_pt","SoftActivityJet_eta","SoftActivityJet_phi","SoftActivityJet_SelectedJetDr","SoftActivityJet_SelectedJetIdx","SoftActivityJet_SelectedMuonDr","SoftActivityJet_SelectedMuonIdx"});
//   selection_twoMuons_twoOppositeSignMuons_twoJets.Report()->Print();



   //   auto tr=SBClassifier_neededselection.Snapshot("ot", "outputFile.root", {"nJet","SBClassifier"});

   //auto tr=toplevel.Snapshot("ot", (out+"File.root").c_str(), {"nJet","nLepton","Jet_LeptonDr","Lepton_JetDr","Jet_LeptonIdx","Jet_pt","Jet_eta","Jet_phi","Lepton_eta","Lepton_phi","Lepton_JetIdx","Lepton_jetIdx"});
   
/*TStopwatch s;
   SBClassifier.OnPartialResult(0, [&s](TH1D &) {
      std::cout << "starting event loop" << std::endl;
      s.Start();
      return true;
   });
   std::cout << "now jitting+event loop..." << std::endl;
   *SBClassifier;
   s.Stop();
   std::cout << "elapsed time: " << s.RealTime() << "s" << std::endl;*/
   auto fff=TFile::Open((out+"Histos.root").c_str(),"recreate");
   for(auto h : histos) h->Write();
   fff->Write();  
   fff->Close();
   return 0;
}
''')

#        for t in to :
#	   if t in self.histos:
#	  	print '%s->Write();'%(t)
#		for w in self.weights :
#		    if self.weights[w]["filter"](t,w):
#	                    print '%s__weight__%s->Write();'%(t,w)




    def baseInputs(self,x) :
        if len(self.Inputs(x)) == 0 :
          return [x]
        else :
          ret=[]
          for i in self.Inputs(x) :
             ret.extend(self.baseInputs(i))
          return ret
    
    def allNodesTo(self,nodes) :   
	  ret=set()
#	  print "Looking for nodes",nodes
	  for x in nodes:
#	      print "x is ",x
#	      print "     inputs ", self.inputs[x]
#	      print "     selections ", self.requirements[x]
#	      print "     weights ", (self.centralWeights[x] if x in self.centralWeights else [])
	      if x in self.nodesto :
		  toset=self.nodesto[x]
	      else:
#	          print >> sys.stderr, "Nodes to ",x
                  self.nodesto[x]=set([x]) 
                  self.nodesto[x].update(self.allNodesTo(self.Inputs(x)))
#		  print x,self.requirements[x]
                  self.nodesto[x].update(self.allNodesTo(self.Requirements(x)))
                  self.nodesto[x].update(self.allNodesTo((self.centralWeights[x] if x in self.centralWeights else [])))
	          #this is cached because it cannot change
	          toset=self.nodesto[x]
              ret.update(toset)
          return ret
    def allNodesFromWithWhiteList(self,nodes,wl,cache=None) :
	  if cache is None:
		cache={}
          #cache is local because it changes with more defines
	  ret=set()
	  for x in nodes :
              if x in cache :
                  childrenset=cache[x]
#	          print >> sys.stderr, "Nodes from ",x
	      else:
#                 children=set([n for n in self.validCols if ((x in self.Inputs(n)+self.Requirements(n)+(self.centralWeights[n] if n in self.centralWeights else [])) and (n in wl))])
                  children=set([n for n in wl if ((x in self.Inputs(n)+self.Requirements(n)+(self.centralWeights[n] if n in self.centralWeights else [])) )])
#		  print "children of ",x,children
                  cache[x]=set(children)
                  cache[x].update(self.allNodesFromWithWhiteList(children,wl,cache))
                  childrenset=cache[x]
	      ret.update(childrenset)
          return ret


    def allNodesFrom(self,nodes,cache=None) :   
	  if cache is None:
		cache={}
	  #cache is local because it changes with more defines
	  ret=set()
	  for x in nodes :
	      if x in cache :
	          childrenset=cache[x]
	      else :
#	          print >> sys.stderr, "Nodes from ",x
                  children=set([n for n in self.validCols if x in self.Inputs(n)+self.Requirements(n)+(self.centralWeights[n] if n in self.centralWeights else [])])
	          childrenset=set(children)
                  childrenset.update(self.allNodesFrom(children,cache))
	          cache[x]=childrenset
	      ret.update(childrenset)
          return ret

    def findAffectedNodesForVariationOnTargets(self,name,targets):
	 nodesTo=set(targets)
 	 nodesTo.update([x for x in self.allNodesTo(targets) if x not in self.variations[name]["exceptions"] ])
         #print  >> sys.stderr, "VarNodesTo:",nodesTo
	 return [x for x in self.allNodesFromWithWhiteList([self.variations[name]["original"]],nodesTo) if x in nodesTo]

   
    def createVariationBranch(self,name,target):
#	 print >> sys.stderr, "Find affected"
         affected=(self.findAffectedNodesForVariationOnTargets(name,target))
#	 print >> sys.stderr, "Found affected\n", affected
	 affected.sort(key=lambda x: self.validCols.index(x)) #keep original sorting
	 res=[y+"__syst__"+name for y in affected if y in target] 
         replacementTable=[(x,x+"__syst__"+name) for x in affected]
         for x,x_syst in replacementTable:
	     ncode=""
	     if self.dupcode :
               ncode=" "+self.code[x]+" "  #FIXME: we should avoid duplicating the code
               for y,y_syst in replacementTable+[(self.variations[name]["original"],self.variations[name]["modified"])]:
	         regBound="([^a-zA-Z0-9_])"
                 reg=regBound+y+regBound
                 ncode=re.sub(reg,"\\1"+y_syst+"\\2",ncode)
	     	
	     originalinputs=self.Inputs(x)
	     repdict={c[0]:c[1] for c in replacementTable+[(self.variations[name]["original"],self.variations[name]["modified"])]}
	     replacedinputs=[(c if c not in repdict else repdict[c]) for c in originalinputs ]
             if x in self.obs:
                 selections=[]
                 for s in self.Requirements(x) :
                     if s in affected :
                        selections.append(s+"__syst__"+name)
                     else :
                        selections.append(s)

 		 if x+"__syst__"+name not in self.validCols:
 		   if self.dupcode:
	                 self.Define(x_syst,ncode,requires=selections)
		   else :
	                 self.Define(x_syst,ncode,requires=selections,original=x,inputs=replacedinputs)
	
             if x in self.selections:
                 selections=[]
                 for s in self.Requirements(x) :
                     if s in affected :
                        selections.append(s+"__syst__"+name)
                     else :
                        selections.append(s)


 		 if x+"__syst__"+name not in self.validCols:
	             if self.dupcode :
	                 self.Selection(x_syst,ncode,requires=selections)
			 #FIXME: weights
		     else:
	                 self.Selection(x_syst,ncode,original=x,inputs=replacedinputs,requires=selections)
			 if x in self.centralWeights :
	 	     		 replacedweights=[(c if c not in repdict else repdict[c]) for c in self.centralWeights[x] ]
        	         	 self.centralWeights[x_syst]=replacedweights

#	 print "Recomputing for",res
         return res

    def createSystematicBranches(self,systematics,selWithHistos):
	res=copy.deepcopy(selWithHistos)
	for syst in systematics :
	    print "creating systematic",syst
 	    for sel in selWithHistos :
		selWsyst=self.createVariationBranch(syst,[sel])
	   	histos=self.createVariationBranch(syst,selWithHistos[sel])
		selKey=sel
		if selWsyst:
		   selKey = selWsyst[0]
		if selKey not in res:
		    res[selKey]=[x for x in selWithHistos[sel] if x+"__syst__"+syst not in histos] #copy.deepcopy(res[sel])
		res[selKey]+=histos
        return res
    




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
