#ifndef HELPERS_H
#define HELPERS_H
#include <TFile.h>
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
#include <TH1F.h>
static int verbosecount = 100;

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

template <typename type>
auto At(const type &v, size_t i,  typename type::value_type def){
// if(i>=v.size() || i < 0) {std::cout << "ERROR out of boundaries" << i << " vs " <<  long(i) <<  v.size() << std::endl; return def; }
 if(i>=v.size() || i < 0) { return def; }
 return v[i];
}

template <typename type>
auto At(const type &v, size_t i){
 typedef typename type::value_type rettype;
 if(i>=v.size() || i < 0) {if(verbosecount-->0) {std::cout << "ERROR out of boundaries" << i  << " " <<  long (i) << " vs " <<  v.size() <<  std::endl;}  return  rettype();}
 return v[i];
}


template <typename type,typename masktype>
auto At(const type &v, const ROOT::VecOps::RVec<masktype> &m){
 if(v.size() != m.size()) {if(verbosecount-->0){std::cout << "ERROR mismatch mask length" << std::endl;  return  v[ROOT::VecOps::RVec<masktype>(v.size())];}}
 return v[m];
}

template <typename type,typename indextype>
auto TakeDef(const type &v, const ROOT::VecOps::RVec<indextype> &m,  typename type::value_type def){
 type ret;
 for(auto i : m) {
   if(i>=0 and i < v.size()) ret.push_back(v[i]); else ret.push_back(def);
 }
 return ret;
}
template <typename type,typename indextype>
auto TakeDef(const type &v, const ROOT::VecOps::RVec<indextype> &m,  const type &defs){
 type ret;
 size_t ii=0;
 for(auto i : m) {
   if(i>=0 and i < v.size()) ret.push_back(v[i]); else ret.push_back(defs[ii]);
   ii++;
 }
 return ret;
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
 return 1.0;
}
float btagWeightUp(float csv,float pt,float eta){
 return 1.0;
}

//float efficiency(float pt,float eta,int pid){
// if(pid==11) return 0.99;
// if(pid==13) return 0.92;
//}
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


template <typename size_type,typename size_type2,typename size_type3>
ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_type>> Combinations(size_type size1, size_type2 size2,size_type3 size3) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_type>> r(3);
    r[0].resize(size1*size2*size3);
    r[1].resize(size1*size2*size3);
    r[2].resize(size1*size2*size3);
    size_type c = 0;
    for(size_type i=0; i<size1; i++) {
       for(size_type j=0; j<size2; j++) {
          for(size_type k=0; k<size3; k++) {
             r[0][c] = i;
             r[1][c] = j;
             r[2][c] = k;
             c++;
	  }
       }
    }
    return r;
 }


ROOT::VecOps::RVec<size_t> Range(size_t n){
  ROOT::VecOps::RVec<size_t>  res;
  for(size_t i=0;i<n;i++) res.push_back(i);
  return res;
}

void loadHistograms(const ROOT::RDF::RNode &rdf, const std::string & name, std::vector<ROOT::RDF::RResultPtr<TH1D>>  & histos ){

}
template <typename T>
T SumDef(const ROOT::VecOps::RVec<T> &v)
{
   return std::accumulate(v.begin(), v.end(), T());
}
#endif
