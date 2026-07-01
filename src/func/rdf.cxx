// rdf.cxx
module;

#include <vector>
#include <thread>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

import binning;
import types;
import logger;

export module rdf;
export{

struct RawRDF{
    int ntype, nb;
    std::vector<std::vector<std::vector<int>>> count;

    RawRDF() : ntype(0), nb(0) {}

    RawRDF(int n, int cap) : ntype(n), nb(cap){
        count.resize(ntype);
        for(std::vector<std::vector<int>>& kk : count){
            kk.resize(ntype);
            for(std::vector<int>& pp : kk){
                pp.resize(cap);
            }
        }
    }
};

struct Histogram{
    int cap;
    std::vector<float> lbound, ubound, center;

    Histogram(float cutoff, int capacity) : cap(capacity) {
        const float binsize = cutoff / static_cast<float>(cap);

        lbound.resize(cap);
        ubound.resize(cap);
        center.resize(cap);

        std::iota(lbound.begin(), lbound.end(), 0.0f);
        std::transform(lbound.begin(), lbound.end(), lbound.begin(),
                      [binsize](float v) { return v * binsize; });

        std::transform(lbound.begin(), lbound.end(), ubound.begin(),
                      [binsize](float v) { return v + binsize; });

        std::transform(lbound.begin(), lbound.end(), center.begin(),
                      [binsize](float v) { return v + binsize * 0.5f; });
    }
};

void compute_dist(const BinList& bins,
                  const std::vector<std::array<mreal, 3>> xyz0,
                  int start, int end,
                  const std::array<int, 3> bin_num,
                  mreal cutoff,
                  const std::vector<bool> is_border,
                  const std::vector<int>& ptype,
                  RawRDF& rdf_raw){

    const mreal r = cutoff*cutoff;
    const Histogram his(cutoff, rdf_raw.nb);

    for (int i=start;i<end;++i){
        if(!is_border.at(i)){
        int z_bin = i / ((bin_num[0]+2)*(bin_num[1]+2));
        int y_bin = (i % ((bin_num[0]+2)*(bin_num[1]+2))) / (bin_num[0]+2);
        int x_bin = i % (bin_num[0]+2);
        // std::cout<<x_bin<<" "<<y_bin<<" "<<z_bin<<std::endl;
        for (int atom = 0; atom < bins.n[i]; ++atom){
            int id =  bins.ids[atom].at(i);
            for (int p = -1; p<=1; ++p){
                for (int q = -1; q<=1; ++q){
                    for (int o = -1; o<=1; ++o){
                        int checked_index = x_bin+p
                                        + (y_bin+q)*(bin_num[0]+2)
                                        + (z_bin+o)*(bin_num[0]+2)*(bin_num[1]+2);
                        int checked_n = bins.n[checked_index];
                        for (int atom2 = 0; atom2 < checked_n; ++atom2){
                            int checkid = bins.ids[atom2].at(checked_index);
                            if (checkid < id){
                                mreal d = (xyz0[id][0]-xyz0[checkid][0])*(xyz0[id][0]-xyz0[checkid][0])
                                        + (xyz0[id][1]-xyz0[checkid][1])*(xyz0[id][1]-xyz0[checkid][1])
                                        + (xyz0[id][2]-xyz0[checkid][2])*(xyz0[id][2]-xyz0[checkid][2]);
                                if (d < r) {
                                    int type_1 = ptype.at(id);
                                    int type_2 = ptype.at(checkid);

                                    for(int k=0; k<400;++k){
                                        if(k < his.ubound.at(k)){
                                            rdf_raw.count[type_1][type_2][k] += 1;
                                            rdf_raw.count[type_2][type_1][k] += 1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    }

};

void merge_histogram(){};

// Divide the box into meshes to accelarate the neighbor list construction.
void compute_rdf(std::vector<mreal> cutoff, std::vector<Atom>& atoms, Box mbox, std::vector<bool>& pbc){

    pklog.warn("Starting radial distribution function analysis;");

    const int num_threads = 4;

    const int his_num = 400;
    const int ntype = 3;

    // Get the number of bins on each dimension.
    std::array<int,3> bin_num;
    const mreal r = cutoff.at(0);

    bin_num[0] = std::floor(mbox.lx/r);
    bin_num[1] = std::floor(mbox.ly/r);
    bin_num[2] = std::floor(mbox.lz/r);

    std::cout<<"Bin number on each dimension: "<<bin_num[0]+2<<" * "<<bin_num[1]+2<<" * "<<bin_num[2]+2<<std::endl;

    const int bin_num_1d = (bin_num[0]+2)*(bin_num[1]+2)*(bin_num[2]+2);
    int cap = std::ceil(r*r*r);

    std::vector<BinList> multi_bins(num_threads);
    for(BinList& aa:multi_bins){
        aa = BinList(bin_num_1d, cap);
    }

    const int natom = atoms.size();

    std::vector<std::array<mreal, 3>> realxyz(natom);

    for(int i=0;i<natom;++i){
    realxyz[i][0] = atoms[i].x - mbox.xmin;
    realxyz[i][1] = atoms[i].y - mbox.ymin;
    realxyz[i][2] = atoms[i].z - mbox.zmin;
    }

    int per_cpu_atom = natom / num_threads;

    std::vector<std::thread> threads(num_threads);
    //
    for (int t = 0; t < num_threads; ++t) {

        int start = t * per_cpu_atom;
        int end = (t == num_threads - 1) ? natom : start + per_cpu_atom;
        std::cout<<"This is thread "<<t<<std::endl;
        threads[t] = std::thread(sort_atom_to_bin,
                                 std::cref(realxyz),
                                 start, end,
                                 bin_num,
                                 r,
                                 std::ref(multi_bins.at(t)));
    }

    for (int t = 0; t < num_threads; ++t){
        threads[t].join();
    }

    BinList bins = merge_bins(multi_bins);

    multi_bins.clear();
    multi_bins.shrink_to_fit();

    std::vector<bool> is_border = set_bins_pbc(bins, pbc, bin_num);

    // Histogram count raw data for each pair type.
    std::vector<RawRDF> data_3d;
    data_3d.resize(num_threads);
    for(auto& kk: data_3d){
        kk = RawRDF(3, 400);
    }

    // Distribute raw rdf histogram to thread and construct rdf.
    int per_cpu_bin = his_num/ num_threads;
    std::vector<int> ptype(natom);
    std::transform(atoms.begin(), atoms.end(), ptype.begin(),
                      [](const Atom& a) { return a.type; });

    for (int t = 0; t < num_threads; ++t) {
        int start = t * per_cpu_bin;
        int end = (t == num_threads - 1) ? bin_num_1d : start + per_cpu_bin;
        std::cout<<"This is thread "<<t<<std::endl;
        threads[t] = std::thread(compute_dist,
                                 std::cref(bins),
                                 std::cref(realxyz),
                                 start, end,
                                 std::cref(bin_num),
                                 r,
                                 std::cref(is_border),
                                 std::cref(ptype),
                                 std::ref(data_3d.at(t)));
    }

    for (int t = 0; t < num_threads; ++t){
        threads[t].join();
    }

};

}
