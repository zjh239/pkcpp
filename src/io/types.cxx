// types.cxx
module;
#include <vector>

import logger;
export module types;

export{
using mreal = float;

struct Atom {
    int type;
    mreal x, y, z;
};

struct Box{
    mreal xmin, ymin, zmin;
    mreal xmax, ymax, zmax;
    mreal lx, ly, lz;
};

struct Bin{
    int n;
    std::vector<int> ids;
};

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

    // Addition operator
    NeighList operator+(const NeighList& other) const {
        if (n.size() != other.n.size()) {
            pklog.fatal("NeighList/BinList sizes don't match for addition");
        }

        NeighList result;
        result.n.resize(n.size());
        result.ids.resize(ids.size());

        for (size_t i = 0; i < n.size(); ++i) {
            // Sum the n values
            result.n[i] = n[i] + other.n[i];
            if (result.n[i] > ids.at(i).size()){
                pklog.fatal("NeighList/BinList capacity not large enough!");
            }

            // Combine ids vectors (append other.ids[i] to ids[i])
            result.ids[i].reserve(ids[i].size());
            result.ids[i] = ids[i];  // Copy first
            result.ids[i].insert(result.ids[i].end(),
                                  other.ids[i].begin(),
                                  other.ids[i].end());
            
            result.ids[i].resize(ids.at(i).size(), 0);
        }

        return result;
    }

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
};
using BinList = NeighList;

struct DeltaList{};

}
