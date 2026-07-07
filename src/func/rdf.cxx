// rdf.cxx
module;

#include <vector>
#include <thread>
#include <iostream>
#include <cmath>

import binning;
import types;
import logger;

export module rdf;
export{

// like neighbor finder, construct a pair-wise counting histogram
void sort_to_rdf_hist(const Box& mbox,
                      const BinList& bins,
                      const AtomList& xyz0,
                      const int start, const int end,
                      const std::array<int, 3> bin_num,
                      const mreal cutoff,
                      const std::vector<bool>& is_border,
                      RDFCount& rdf_raw){
    pklog.debug("are you okay?");
    const int nbin = rdf_raw.nbin;
    HistPos his(cutoff, nbin);
    
    // i: bin index
    // j: first atom index
    // h: second atom index
    // m: adjacent bin index
    for (int i=start;i<end;++i){
        if(!is_border.at(i)){
            int z_bin = i / ((bin_num[0]+2)*(bin_num[1]+2));
            int y_bin = (i % ((bin_num[0]+2)*(bin_num[1]+2))) / (bin_num[0]+2);
            int x_bin = i % (bin_num[0]+2);
            
            pklog.debug("Analyzing bin {} {} {};", x_bin, y_bin, z_bin);

            for (int j = 0; j < bins.n[i]; ++j){
                int id =  bins.ids[i].at(j);
                for (int p = -1; p<=1; ++p){
                    for (int q = -1; q<=1; ++q){
                        for (int o = -1; o<=1; ++o){
                            int m = x_bin+p
                                            + (y_bin+q)*(bin_num[0]+2)
                                            + (z_bin+o)*(bin_num[0]+2)*(bin_num[1]+2);
                            int checked_n = bins.n[m];
                            for (int h = 0; h < checked_n; ++h){
                                int checkid = bins.ids[h].at(m);
                                if (checkid < id){
                                    mreal dx = xyz0.x[id] - xyz0.x[checkid];
                                    mreal dy = xyz0.y[id] - xyz0.y[checkid];
                                    mreal dz = xyz0.z[id] - xyz0.z[checkid];
                                    mreal d = std::hypot(dx, dy, dz);

                                    if (d < cutoff) {
                                        int type_1 = xyz0.ptype.at(id);
                                        int type_2 = xyz0.ptype.at(checkid);

                                        for(int k=0; k<nbin;++k){
                                            if(k < his.ubound.at(k)){
                                                //rdf_raw.count[type_1][type_2][k] += 1;
                                                rdf_raw.at(type_1, type_2, k) += 1;
                                                //rdf_raw.count[type_2][type_1][k] += 1;
                                                rdf_raw.at(type_2, type_1, k) += 1;
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

auto merge_histogram(std::vector<RDFCount>& multi_count) -> RDFCount {
    RDFCount result(multi_count.at(0));

    int n = multi_count.size();
    for (int i=1;i<n;++i){
        result += multi_count[i];
    }

    multi_count.clear();
    multi_count.shrink_to_fit();

    for (int i=0; i<result.nbin; ++i){
        int j = result.at(0,0,i);
        pklog.info("RDF at: {};", j);
    }

    return result;
};

// Divide the box into meshes to accelarate the neighbor list construction.
auto compute_rdf(const std::vector<mreal> cutoff,
                 const AtomList& atoms,
                 const Box mbox,
                 const std::vector<bool>& pbc,
                 const int hist_num = 400) -> RDFCount{

    pklog.warn("Starting radial distribution function analysis;");
    const int num_threads = thread_number;

    // Get the number of bins on each dimension.
    std::array<int,3> bin_num;
    const mreal r = cutoff.at(0);

    bin_num[0] = std::floor(mbox.lx/r);
    bin_num[1] = std::floor(mbox.ly/r);
    bin_num[2] = std::floor(mbox.lz/r);

    pklog.info("Bin number on each dimension: {} x {} x {};", bin_num[0]+2, bin_num[1]+2, bin_num[2]+2);

    const int bin_num_1d = (bin_num[0]+2)*(bin_num[1]+2)*(bin_num[2]+2);
    pklog.info("Total bin number: {};", bin_num_1d);
    int cap = std::ceil(r*r*r);

    std::vector<BinList> multi_bins(num_threads);
    for(BinList& aa:multi_bins){
        aa = BinList(bin_num_1d, cap);
    }

    const int natom = atoms.x.size();

    AtomList xyz_real(natom);
    
    xyz_real.ntype = atoms.ntype;
    for(int i=0;i<natom;++i){
        xyz_real.x[i] = atoms.x[i] - mbox.xmin;
        xyz_real.y[i] = atoms.y[i] - mbox.ymin;
        xyz_real.z[i] = atoms.z[i] - mbox.zmin;
    }
    xyz_real.ptype = atoms.ptype;

    int per_cpu_atom = natom / num_threads;

    std::vector<std::thread> threads(num_threads);
    //
    for (int t = 0; t < num_threads; ++t) {

        int start = t * per_cpu_atom;
        int end = (t == num_threads - 1) ? natom : start + per_cpu_atom;
        pklog.debug("Thread {} will handle xyz from {} to {};", t, start, end-1);
        threads[t] = std::thread(sort_atom_to_bin,
                                 std::cref(xyz_real),
                                 start, end,
                                 bin_num,
                                 r,
                                 std::ref(multi_bins.at(t)));
    }

    for (int t = 0; t < num_threads; ++t){
        threads[t].join();
    }

    BinList bins = merge_bins(multi_bins);

    std::vector<bool> is_border = set_bins_pbc(bins, pbc, bin_num);
    int tc = std::count(is_border.begin(), is_border.end(), true);
    pklog.info("{} border bins modified;", tc);

    // multi thread data
    std::vector<RDFCount> multi_count;
    multi_count.reserve(num_threads);
    multi_count.emplace_back(atoms.ntype, hist_num);
    
    int per_cpu_bin = bin_num_1d/ num_threads;

    for (int t = 0; t < num_threads; ++t) {
        int start = t * per_cpu_bin;
        int end = (t == num_threads - 1) ? bin_num_1d : start + per_cpu_bin;
        pklog.debug("Thread {} is handling bins from {} to {};", t, start, end);
        threads[t] = std::thread(sort_to_rdf_hist, std::cref(mbox),
                                 std::cref(bins),
                                 std::cref(xyz_real),
                                 start, end,
                                 std::cref(bin_num),
                                 r,
                                 std::cref(is_border),
                                 std::ref(multi_count.at(t)));
    }

    for (int t = 0; t < num_threads; ++t){
        threads[t].join();
    }
    
    RDF result(r, hist_num, atoms.ntype);

    return merge_histogram(multi_count);

};

}
