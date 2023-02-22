# Natual Analysis Implementation Language

NAIL is an analysis language that should allow to define in an abstract way a data analsysis of a typical HEP experiment such as CMS or ATLAS.
NAIL assumes an input data model for the event to process (see below) and allow to specify the event by event processing actions in a declarative form.
Analysis variations for optimizations and systematics do not need to be explicitelly coded but are automatically derived from the event processing computational graph.
Currently ROOT's RDataFrame is used as backend for a concrete implementation of the event processing as it allows parallelization and lazy evaluation.

The final goal is to be able to query for a result and be able to automatically trigger the needed calculation with any optimization happening automatically.

The idea is that the user can employ NAIL not only to describe the analysis at a given step but to develop the analysis, querying for intermediate results, updating the event processing strategy, adding systematics or exploring optimizations, all the way to finally statistically interepted data.

## Current status
A prototype for the Event processing exists. As it is the more delicate part it was developed first. We anticipate that the implementation of datasets handling, stat interpretation etc.. is a simpler problem that needs less R&D


## Input Data model
The event data used is of 4 different types
- scalars (one number per event)
- vectors (multiple entries per event, with variable number of entries)
- objects (sets of a fixed number of scalars per events representing different properties of the same object)
- collections (a variable number of objects per event)

Each property of an object in a collection can be seen as a simple vector, but all the vectors representing individual properties have the same length.
Each property of singleton object is just a scalar.

I.e. the data of a single event can be seen as split on an fixed size horizontal axis, with some "per object" or "per collection" grouping, and on a variable size vertical axis that is not only different event by event but also different for different objects.

The NANOAOD naming convention is assumed so that length of vectors and collections are named _nNameOfCollection_ and properties of an object or a collection are named *NameOfCollection_property*

![datamodel](image.png)



## Event processing actions

##### DefaultConfig
DefaultConfig command allows to specify in one line multiple fixed values columns such as cut values that can be later variated in systematic shifts or optimization studies.

##### Define
Define a variable of any type given some code snippet to compute it. It can reference any know variable (either a defined one or one known to be available in the input)
Requirements can be specified, a _requirement_ is a selection that must be fullfilled in order for the code to make sense (e.g. Muon_pt[1] make sense only if there are at least 2 muons)

##### Selection
An event selection defined with an expression that can be evaluated as True or False

##### SubCollection
While each column can be considered independently some sets of columns are properties of the same analysis object and hence can be considered as aligned (i.e. same index in the column correspond to the same object). The SubCollection action allows to thin the collection applying a selection on it. SubCollection generates a new variable for each input property of the original collection.

##### Distinct
Create a list of pairs of the given input collection and new colulmns with all object properties for the element 0 and 1 of the pairs.

##### TakePair
Extract an individual pair from a list of pairs creating new scalar columns for the properties of element 0 and 1 of the pair

##### ObjectAt
Extract an individual object from a collection creating new scalar columns for the properties of the object

##### MergeCollections
Merge two collections having all properties shared by the two input collections. If the user want to put dummy values for non shared properties it must do so with a prior define (e.g. merging Muon and Electron the user can _Define_ Electron_nMuonStations=-99 to have this variable appearing also in the merged collection)

##### Match
Matches two collections via DeltaR or a user specified metric in a non-unique way. For each element of collection A the closest element of collection B is computed and viceversa. Delta R and index of the matched object are defined. 

##### CentralWeight
Add a weight that has to be used when making nominal distribution (central value). CentralWeight can be called several times, the product of the weights will be the resulting nominal weight. The weight can be applied only to some phase space regions by specifying the list of selections the weight belong to (a selection inherit the product of the weight from selections that are listed as requirements)

##### VariationWeight
Add a systematic/variation weight. Those weights are applied in addition to central weight to evaluate systematic shifts. By default distribution made with those weight are not computed for distribution obtained with systematic variation of some input (i.e. we do not evaluate the effect of applying multiple systematics in one go)

##### VariationWeightArray
Same as VariationWeight but taking an array of weights in input and creating a different weight for each element. The number of elements to consider must be specified

### Inner looping operations and special keywords
Several operations can be defined to (inner) loop on an collection (e.g. the Muons).
In addition some special operations starting with "@" are implemented to simplify common operations
Currently *nail* uses ROOT RDataFrame as a backend hence some of the syntax is taken directly from it.
For example logical operations can be defined across columns and result in a vector of 1 and 0, e.g.
"Muon_pt > 20 && abs(Muon_eta) < 2" is a vector of 1 and 0 corresponding to muons passing/failing such selection.
A selection like that can be directly used as Muon_phi[Muon_pt > 20 && abs(Muon_eta) < 2] returning a vector with only the phi of the muons passing the selection.
it should be noted that as of today Muon_phi[2] also works and return phi of the 3rd muon, while in order to obtain the muons and indices {1,3,5} the special function "Take" should be used instead of operator[] (feature request to ROOT was made to improve this https://sft.its.cern.ch/jira/browse/ROOT-10071)

###### Special keywords: @p4, @p4v
@p4 and @p4v expand a collection name to get its pt,eta,phi,mass property and build a 4 vector. E.g. one can do
Define("Muon_p4","@p4v(Muon)") or Define("Muon0_p4","@p4(Muon[0])") or Define("Higgs_p4","@p4(Higgs)")

###### vector_map(function, collection_properties...) 
calls *function* on each element of the collection passing the listed properties as argument.
 This is now available in the new Map function of VecOps


A variant vector_map_t<Type> exists to call a constructor of Type for each entry in the collection instead of calling the *function*.
  

###### MemberMap (preprocessor macro)
*MemberMap* allows to easily call a member of class A for each entry in a vector<A>
  

###### operations defined in ROOT::VecOps
other operations are listed here https://root.cern.ch/doc/master/namespaceROOT_1_1VecOps.html 
and can be applied to inner vectors.
For example Argsort,Max,Min,Argmax,Argmin, DeltaPhi, DeltaR, InvariantMass ...





## Yield and Queries
The output of the processing consists of histograms and ntuples. Histograms can be specified as a dictionary where for each selection a list of histograms is given. Ntuples are specified with a selection and a list of columns to store.

## Histogram binning
Rules should be specified for binning and range using regexps


# Testing recipe

```
pip3 install libclang clang
git clone git@github.com:arizzi/nail.git
cd nail/
git checkout py3
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh 
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
python3 simple.py 
```
