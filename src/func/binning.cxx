// binning.cxx
module;

#include <vector>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
    #include <omp.h>
#endif

import types;
import logger;
export module binning;

export{

// Put atoms to the bin.
void sort_atom_to_bin(std::vector<std::array<mreal, 3>> xyz0, int start, int end, std::array<int,3> binnum, mreal r, BinList& thread_bins){

    #pragma omp parallel for
    for (int id=start; id < end; ++id){

        int thread_id = omp_get_thread_num();

        int xbin = std::ceil(xyz0[id][0]/r);
        int ybin = std::ceil(xyz0[id][1]/r);
        int zbin = std::ceil(xyz0[id][2]/r);

        if (xbin > binnum[0]) xbin = binnum[0];
        if (ybin > binnum[1]) ybin = binnum[1];
        if (zbin > binnum[2]) zbin = binnum[2];

        if (xbin == 0) xbin = 1;
        if (ybin == 0) ybin = 1;
        if (zbin == 0) zbin = 1;

        int binidex = xbin + ybin*(binnum[0]+2) + zbin*(binnum[0]+2)*(binnum[1]+2);
        int k = thread_bins.n.at(binidex);
        thread_bins.ids.at(binidex).at(k) = id;
        thread_bins.n.at(binidex) += 1;
    };
};

// Merge the bin data
auto merge_bins(std::vector<BinList>& tmp_bins) -> BinList{
    BinList bins(tmp_bins.at(0));
    int n = tmp_bins.size();
    for (int i=1;i<n;++i){
        bins+=tmp_bins[i];
    }

    tmp_bins.clear();
    tmp_bins.shrink_to_fit();

    return bins;
};

// modify the periodic images and return border index.
auto set_bins_pbc(BinList& bins, const std::vector<bool>& pbc, const std::array<int,3>& bin_num)
    -> std::vector<bool>{

    std::vector<bool> is_border(bins.n.size(), 0);
    pklog.error("Now only using the first PBC, check in future!");

    for (int xbin = 0; xbin <= bin_num[0]+1; ++xbin){
        for (int ybin = 0; ybin <= bin_num[1]+1; ++ybin){
            for (int zbin = 0; zbin <= bin_num[2]+1; ++zbin){
                int x_pbc=0, y_pbc=0, z_pbc=0;
                int this_index = xbin
                                + ybin*(bin_num[0]+2)
                                + zbin*(bin_num[0]+2)*(bin_num[1]+2);
                if (xbin == 0)           {x_pbc = 1; is_border.at(this_index)=1;}
                if (xbin == bin_num[0]+1){x_pbc = -1; is_border.at(this_index)=1;}
                if (ybin == 0)           {y_pbc = 1; is_border.at(this_index)=1;}
                if (ybin == bin_num[1]+1){y_pbc = -1; is_border.at(this_index)=1;}
                if (zbin == 0)           {z_pbc = 1; is_border.at(this_index)=1;}
                if (zbin == bin_num[2]+1){z_pbc = -1; is_border.at(this_index)=1;}

                if (is_border.at(this_index)){
                    if(!pbc[0] && x_pbc != 0) x_pbc = 0;
                    if(!pbc[1] && y_pbc != 0) y_pbc = 0;
                    if(!pbc[2] && z_pbc != 0) z_pbc = 0;

                    int t_xbin = xbin + x_pbc*bin_num[0];
                    int t_ybin = ybin + y_pbc*bin_num[1];
                    int t_zbin = zbin + z_pbc*bin_num[2];
                    pklog.trace("{}, {}, {} bin copying from {}, {}, {};", xbin, ybin, zbin, t_xbin, t_ybin, t_zbin);
                    int target = xbin + x_pbc*bin_num[0]
                                + (ybin + y_pbc*bin_num[1])*(bin_num[0]+2)
                                + (zbin + z_pbc*bin_num[2])*(bin_num[0]+2)*(bin_num[1]+2);
                    bins.n.at(this_index) = bins.n.at(target);
                    bins.ids.at(this_index) = bins.ids.at(target);
                }
            }
        }
    }

    auto [at, bt] = std::minmax_element(bins.n.begin(), bins.n.end());
    int a = *at;
    int b = *bt;

    for (int i=a;i<=b;++i){
        int j = std::count(bins.n.begin(), bins.n.end(), i);
        pklog.info("{} bins with atom number = {};", j, i);
    }

    return is_border;
};

}
