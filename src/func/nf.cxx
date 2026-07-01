// nf.cxx
// Neighbor list construction.
module;

#include <print>

#include <vector>
#include <thread>
#include <iostream>
#include <cmath>
#include <algorithm>

import types;
import binning;
import logger;

export module nf;

export{

size_t nt = 2; // set thread number from xmake in future

// Print neighbor list.
void print_bin(Bin& kk){
    int n = kk.n;
    std::cout<<"(num): "<<n<<" (ids):";
    for(int c=0;c<n;++c){
        std::cout<<" "<<kk.ids.at(c);
    }
    std::cout<<std::endl;
};

void find_neighbors(const BinList& bins,
                    const std::vector<std::array<mreal, 3>> xyz0,
                    const int start, const int end,
                    const std::array<int, 3> bin_num,
                    mreal cutoff,
                    const std::vector<bool> is_border,
                    bool cal_delta,
                    NeighList& neigh_list,
                    std::vector<std::vector<std::array<float, 3>>>& delta){

    const mreal r = cutoff*cutoff;

    int cap=std::ceil(cutoff*cutoff*cutoff);

    if (cal_delta){
        delta.resize(neigh_list.n.size());
        for(std::vector<std::array<float, 3>>& kk : delta){
        kk.resize(cap);
    }
    }

    for (int i=start;i<end;++i){
        if(!is_border.at(i)){
            int x_bin = i % (bin_num[0]+2);
            int y_bin = (i % ((bin_num[0]+2)*(bin_num[1]+2))) / (bin_num[0]+2);
            int z_bin = i / ((bin_num[0]+2)*(bin_num[1]+2));

        for (int atom = 0; atom < bins.n.at(i); ++atom){
            int id =  bins.ids.at(i).at(atom);

            pklog.trace("{} {} {} {} {} {} {}", id, x_bin, y_bin, z_bin, xyz0.at(id)[0], xyz0.at(id)[1], xyz0.at(id)[2]);

            for (int p = -1; p<=1; ++p){
                for (int q = -1; q<=1; ++q){
                    for (int o = -1; o<=1; ++o){
                        int checked_index = x_bin+p
                                        + (y_bin+q)*(bin_num[0]+2)
                                        + (z_bin+o)*(bin_num[0]+2)*(bin_num[1]+2);
                        int checked_n = bins.n.at(checked_index);
                        for (int atom2 = 0; atom2 < checked_n; ++atom2){
                            int checkid = bins.ids.at(checked_index).at(atom2);
                            if (checkid < id){
                                mreal d = (xyz0[id][0]-xyz0[checkid][0])*(xyz0[id][0]-xyz0[checkid][0])
                                        + (xyz0[id][1]-xyz0[checkid][1])*(xyz0[id][1]-xyz0[checkid][1])
                                        + (xyz0[id][2]-xyz0[checkid][2])*(xyz0[id][2]-xyz0[checkid][2]);
                                if (d < r) {
                                    neigh_list.ids.at(checkid).at(neigh_list.n.at(checkid)) = id;
                                    neigh_list.ids.at(id).at(neigh_list.n.at(id)) = checkid;

                                    if(cal_delta){
                                        delta[checkid].at(neigh_list.n.at(checkid))[0] = xyz0[id][0]-xyz0[checkid][0];
                                        delta[checkid].at(neigh_list.n.at(checkid))[1] = xyz0[id][1]-xyz0[checkid][1];
                                        delta[checkid].at(neigh_list.n.at(checkid))[2] = xyz0[id][2]-xyz0[checkid][2];

                                        delta[id][neigh_list.n.at(id)][0] = xyz0[checkid][0]-xyz0[id][0];
                                        delta[id][neigh_list.n.at(id)][1] = xyz0[checkid][1]-xyz0[id][1];
                                        delta[id][neigh_list.n.at(id)][2] = xyz0[checkid][2]-xyz0[id][2];
                                    }
                                    neigh_list.n.at(checkid) += 1;
                                    neigh_list.n.at(id) += 1;
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

auto merge_neighbor(std::vector<NeighList>& mn) -> NeighList{

    NeighList result(mn.at(0));

    int n = mn.size();
    for (int i=1;i<n;++i){
        result+=mn[i];
    }

    mn.clear();
    mn.shrink_to_fit();

    auto [at, bt] = std::minmax_element(result.n.begin(), result.n.end());
    int a = *at;
    int b = *bt;

    for (int i=a;i<=b;++i){
        int j = std::count(result.n.begin(), result.n.end(), i);
        pklog.info("{} atoms with CN = {};", j, i);
    }

    return result;
};

std::vector<std::vector<std::array<float, 3>>> a(0); // placeholder

// Divide the box into meshes to accelarate the neighbor list construction.
auto build_neighbors(const std::vector<mreal> cutoff,
                     const std::vector<Atom>& atoms,
                     const Box mbox,
                     const std::vector<bool> pbc,
                     // NeighList& neighbor_list,
                     std::vector<std::vector<std::array<float, 3>>> delta=a) -> NeighList{

    pklog.warn("Starting neighbor list constructing");
    const int num_threads = 1;

    // Get the number of bins on each dimension.
    std::array<int,3> bin_num;
    const mreal r = cutoff.at(0);
    pklog.error("Remember to use pair-wise cutoff in future!");

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
        pklog.debug("Thread {} will handle xyz from {} to {};", t, start, end-1);
        threads.at(t) = std::thread(sort_atom_to_bin,
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

    std::vector<bool> is_border = set_bins_pbc(bins, pbc, bin_num);
    int tc = std::count(is_border.begin(), is_border.end(), true);
    pklog.info("{} border bins modified;", tc);

    std::vector<NeighList> multi_nlist(num_threads);
    for (NeighList& kk: multi_nlist){
        kk = NeighList(natom, cap);
    }

    bool calculate_delta = delta.size();
    pklog.error("Delta calculation indicator is {}, check in future!", calculate_delta);

    // Distribute bins to thread and construct neighbor list.
    int per_cpu_bin = bin_num_1d / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        int start = t * per_cpu_bin;
        int end = (t == num_threads - 1) ? bin_num_1d : start + per_cpu_bin;
        pklog.debug("Thread {} will handle bins from {} to {};", t, start, end-1);
        threads[t] = std::thread(find_neighbors,
                                 std::cref(bins),
                                 std::cref(realxyz),
                                 start, end,
                                 std::cref(bin_num),
                                 r,
                                 std::cref(is_border),
                                 calculate_delta,
                                 std::ref(multi_nlist.at(t)),
                                 std::ref(delta));
    }

    for (int t = 0; t < num_threads; ++t){
        threads[t].join();
    }

    return merge_neighbor(multi_nlist);


};


}
