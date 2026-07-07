// init.cxx
// the full command is
// g++ -std=c++20 -fmodules-ts init.cpp init.cxx main.cpp
module;

#include <iostream>
#include <string>
#include <vector>

import parser;
import flags;
import read_file;
import nf;
import rdf;
import types;
import logger;

export module init;

export {
void init_logger(){
    pklog.enable_console(true);
    // pklog.enable_file("game.log");
    pklog.set_min_level(LogLevel::DEBUG);

};

void initialize(int argc, char** argv){

    set_thread(1);  // empty for auto-detecting
    pklog.warn("Thread number set to {};", thread_number);
    // Set a default thread count, e.g., use hardware threads, but fallback to 2.
    // int thread_count = hardware_threads != 0 ? hardware_threads : 2;

    init_logger();

    get_input_options(argc, argv);
    save_config();

    Box box;
    AtomList atoms = read_xyz_file(file_name, box);

    if (run_nf) {
        DeltaList delta_list;
        NeighList neighbor_list = build_neighbors(cutoffs, atoms, box, pbcs);
    }
    if (run_rdf) {
        //NeighList neighbor_list = build_neighbors(rdf_cut, atoms, box, pbcs);
        RDFCount rdf_data = compute_rdf(rdf_cut, atoms, box, pbcs);
    }
    // std::cout<<atoms[1].x;
    // run_analysis();
};

}

// void run_analysis(float** &xyz){
//     std::cout<<"deal with xyz"<<std::endl;
//     get_xyz_from_file(d_name+"/"+f_name);
//
//     if (run_nfd){
//         neighbor_finder_vector();
//     }else if(run_nf){
//         neighbor_finder();
//     };
//     if(run_poly){polyhedral_analysis();};
//     if(run_rs){ring_analysis();};
//     if(run_nc){neighbor_change();};
// };
//
//
//
//
// export void neighbor_finder(){
//
// };
//
// };
//
// export struct bins{
//
// int atomNum;
// int atomList[30];
//
// };



