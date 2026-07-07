// types_nf_rdf.cxx
// I put neighbor finder and rdf related data structs here.
module;

#include <vector>
#include <algorithm>
#include <numeric>
#include <span>

import logger;

export module types:nf_rdf;

export {
  struct NeighList{
      std::vector<int> n;
      std::vector<std::vector<int>> ids;

      NeighList(){}

      NeighList(int natom, int cap){
          n.resize(natom, 0);
          ids.resize(natom);
          for (std::vector<int>& kk:ids){
              kk.resize(cap, 0);
          }
      }

      NeighList(const NeighList& other) : n(other.n), ids(other.ids){}

      // Addition assignment operator (optional but recommended)
      NeighList& operator+=(const NeighList& other) {
          if (n.size() != other.n.size()) {
              pklog.fatal("NeighList sizes don't match for addition");
          }

          for (size_t i = 0; i < n.size(); ++i) {
              n[i] += other.n[i];
              if (n[i] > ids.at(i).size()){
                  pklog.fatal("NeighList capacity not large enough!");
              }
              for (size_t j=0; j<other.n.at(i); ++j){
                  ids[i].at(n[i]-other.n[i]+j) = other.ids[i].at(j);
              }
          }

          return *this;
      }

      NeighList operator+(const NeighList& other) const {
          NeighList result(*this);   // copy (uses your existing copy ctor)
          result += other;           // reuse the correct logic
          return result;
      }
  };

  // bin is the mesh when doing neighbor building.
  using BinList = NeighList;
  
  // Couting struct for rdf
  struct RDFCount{  // store the type-wise count in each bin.
      int ntype, nbin;
      //std::vector<std::vector<std::vector<int>>> count;
      std::vector<int> count;

      RDFCount() : ntype(0), nbin(0) {}

      RDFCount(int n, int cap) : ntype(n), nbin(cap){
          count.resize(n*n*cap, 0);
      }
      
      auto view(int i, int j) -> std::span<int> {
          size_t start = i*ntype*nbin + j*nbin;
          return std::span<int>(count.data() + start, nbin);
      }

      auto at(int i, int j, int k) -> int& {
          return count.at(i*ntype*nbin + j*nbin + k);
      }

      RDFCount(const RDFCount& other) : ntype(other.ntype), nbin(other.nbin), count(other.count){}
      
      RDFCount& operator+= (const RDFCount& other){
        if (count.empty() || other.count.empty()){
          pklog.fatal("Empty count");
        }

        if(nbin != other.nbin || ntype != other.ntype){
          pklog.fatal("Adding different size RDF data!");
        }

        std::transform(count.begin(), count.end(),
                      other.count.begin(), count.begin(), std::plus<int>());
        return *this;
      }
      
  };

  struct HistPos{
      int hist_num;
      std::vector<float> lbound, ubound, center;

      HistPos(float range, int num) : hist_num(num) {
          const float binsize = range/ static_cast<float>(hist_num);

          lbound.resize(hist_num);
          ubound.resize(hist_num);
          center.resize(hist_num);

          std::iota(lbound.begin(), lbound.end(), 0.0f);
          std::transform(lbound.begin(), lbound.end(), lbound.begin(),
                        [binsize](float v) { return v * binsize; });

          std::transform(lbound.begin(), lbound.end(), ubound.begin(),
                        [binsize](float v) { return v + binsize; });

          std::transform(lbound.begin(), lbound.end(), center.begin(),
                        [binsize](float v) { return v + binsize * 0.5f; });
      }
  };
  
  struct RDF{
    int hist_num, ntype;
    std::vector<float> center;
    std::vector<std::vector<std::vector<float>>> rdf;

    RDF(){}
    
    RDF(float r, int nbin, int ntype): hist_num(nbin) {
      rdf.resize(ntype);

      for(std::vector<std::vector<float>>& kk : rdf){
          kk.resize(ntype);
          for(std::vector<float>& pp : kk){
              pp.resize(nbin);
          }
      }
    }
  };
}
