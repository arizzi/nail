#ifndef NAIL_EXT_h_
#define NAIL_EXT_h_
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <utility>
using RNode = ROOT::RDF::RNode;
struct Result {
 Result(RNode  rdf_): rdf(rdf_){}
 RNode rdf;
 ROOT::RDF::RResultPtr<TH1D> histo;
 std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
};
Result processRDF(RNode rdf);
#endif


