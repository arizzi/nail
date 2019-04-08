#include <TFile.h>
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
#include <TH1F.h>
template <typename type>
auto Argmax(const type & v){
 return ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(v))[0];
}

template <typename type>
auto Max(const type & v){
 return v[ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(v))[0]];
}

template <typename type, typename Vec,typename... OtherVecs>
auto vector_map_t(const Vec & v,  const OtherVecs &... args) {
  ROOT::VecOps::RVec<type> res(v.size());
  for(size_t i=0;i<v.size(); i++) res[i]=type(v[i],args[i]...);
  return res;
}

template <typename func, typename Vec,typename... OtherVecs>
auto vector_map(func f, const Vec & v,  const OtherVecs &... args) {
  ROOT::VecOps::RVec<decltype(f(std::declval<typename Vec::value_type>(),std::declval<typename OtherVecs::value_type>()...))> res(v.size());
  for(size_t i=0;i<v.size(); i++) res[i]=f(v[i],args[i]...);
  return res;
}

template <typename func, typename Vec>
auto matrix_map(size_t xsize, size_t ysize, size_t axis, func f, const Vec & v) {
  ROOT::VecOps::RVec<decltype(f(std::declval<Vec>()))> res(axis==1?xsize:ysize );
  for(size_t i=0;i<res.size(); i++){
	Vec part(axis==0?xsize:ysize);
  	for(size_t j=0;j<part.size(); j++) {
	   part[j]=v[axis==1?(i*ysize+j):(i+j*ysize)];

	}
	res[i]=f(part);
  }
  return res;
}



/*template <typename func, typename Vec,typename... OtherVecs>
auto matrix_map(shape, int axis,func f, const Vec & v,  const OtherVecs &... args) {

}*/

auto pt(const ROOT::Math::PtEtaPhiMVector &i){
 return i.pt();
}

auto mass(const ROOT::Math::PtEtaPhiMVector &i){
 return i.M();
}

float btagWeight(float csv,float pt,float eta){
 return 1.01;
}
float btagWeightUp(float csv,float pt,float eta){
 return 1.03;
}

float efficiency(float pt,float eta,int pid){
 if(pid==11) return 0.99;
 if(pid==13) return 0.92;
}
template <typename T, typename T2>
ROOT::VecOps::RVec<T> Concat(const ROOT::VecOps::RVec<T> & v1,  const ROOT::VecOps::RVec<T2> & v2){
ROOT::VecOps::RVec<T> v;
for(auto i:v1) {v.push_back(i);}
for(auto i:v2) {v.push_back(i);}
return v;
} 



template <typename size_type,typename size_type2>
ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_type>> Combinations(size_type size1, size_type2 size2) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_type>> r(2);
    r[0].resize(size1*size2);
    r[1].resize(size1*size2);
    size_type c = 0;
    for(size_type i=0; i<size1; i++) {
       for(size_type j=0; j<size2; j++) {
          r[0][c] = i;
          r[1][c] = j;
          c++;
       }
    }
    return r;
 }


void loadHistograms(const ROOT::RDF::RNode &rdf, const std::string & name, std::vector<ROOT::RDF::RResultPtr<TH1D>>  & histos ){

}

#include "hmmtools/hmm_code.h"
