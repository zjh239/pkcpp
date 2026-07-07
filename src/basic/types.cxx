// types.cxx
module;
#include <vector>

import logger;
export module types;

export import :para;
export import :nf_rdf;

export{


using mreal = float;

//struct Atom {
//    int type;
//    mreal x, y, z;
//};

struct AtomList {
    std::vector<mreal> x, y, z;
    std::vector<int> ptype;
    size_t ntype;

    AtomList(){}

    AtomList(int n){
        x.resize(n, 0.0);
        y.resize(n, 0.0);
        z.resize(n, 0.0);
        ptype.resize(n, 0);
    }
};

struct Box{
    mreal xmin, ymin, zmin;
    mreal xmax, ymax, zmax;
    mreal lx, ly, lz;
};

//struct Bin{
//    int n;
//    std::vector<int> ids;
//};


struct DeltaList{};

}
