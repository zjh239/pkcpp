// parser.cxx
module;
#include <print>
#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>

#include <toml++/toml.h>

import flags;
import logger;
import types;

export module parser;

export {
auto vars = toml::table{};

std::string path;
std::string file_name;
std::string dir_name;
std::string* fnames;

int fnumber;
int frame_interval;
std::vector<bool> pbcs;
std::vector<mreal> cutoffs;
std::vector<mreal> rdf_cut;
bool dynamic = 0;

void save_config(){
    std::ofstream file("config.toml");
    file << vars;
};

void show_input_err(char* in){
    pklog.error("Unknown variable '{}'!", in);
    std::terminate();
};

int show_help_msg(){
    pklog.info("Example usage:");
    pklog.info("  poam -f abc.xyz -p 1 -r 2.32 -c pt -o ntpl       (for single-file analysis)");
    pklog.info("  poam -d xyzfiles 10 -p 1 -r 2.32 -c pt -o ntpl   (for multi-files analysis)");
    pklog.info("  -f [filename(char)]   Load data which is a single file, incompatible with '-d' option.");
    pklog.info("  -d [dir(char)]        Load files in the directory, incompatible with '-f' option.");
    pklog.info("      [interval(int)]     ");
    pklog.info("  -r [cutoff(float)]    The cutoff value.");
    pklog.info("  -p [pbc(bool)]        Apply periodic boundary condition or not.");
    return 0;
};

int show_version_msg(){
    pklog.info("PoAM - Polyhedral Analysis Module (Version {})", 0.1);
    pklog.info("  Please kindly cite:");
    pklog.info("  ---");
    pklog.info("  Room temperature plasticity in amorphous SiO2 and amorphous Al2O3 : A computational and topological study. Zhang, J., Frankberg, E. J., Kalikka, J. & Kuronen, A.  Acta Mater. 259, 119223. https://doi.org/10.1016/j.actamat.2023.119223 (2023).");
    pklog.info("  ---");
    pklog.info("  Bug report: zjh239@foxmail.com");
    return 0;
};

auto get_cutoffs(const std::string& a) -> std::vector<mreal>{

    std::vector<mreal> result;
    std::stringstream ss(a);
    std::string tmp;
    mreal item;

    while (std::getline(ss, tmp, ',')) {
        item = std::stof(tmp);
        result.push_back(item);
    }
    return result;
}

auto get_pbcs(const std::string& a) -> std::vector<bool> {
    std::vector<bool> result;
    std::stringstream ss(a);
    std::string tmp;
    bool item;
    while (std::getline(ss, tmp, ',')) {
        item = std::stoi(tmp);
        if (item != 1 && item != 0){
            pklog.fatal("PBC parameter must be 1 or 0!");
        }
        result.push_back(item);
    }
    if(result.size() != 3){
        if(result.size() != 1){
            pklog.fatal("Only 1 or 3 PBC parameters should be given!");
        }
        result.push_back(result.at(0));
        result.push_back(result.at(0));
    }
    pklog.info("Periodic boundary condition on each dimension is {}, {}, {};",result.at(0),result.at(1),result.at(2));
    return result;
}

void get_input_options(int argc, char** argv){

    if(argc==1){
        show_help_msg();
        show_version_msg();
    }

    for(int c=1;c<argc;++c){

        std::string cli(argv[c]);
        // std::cout<<cli<<std::endl;

        if(cli=="-f"){
            path = "";
            file_name = argv[c+1];

            vars.insert("path", "");
            vars.insert("file", argv[c+1]);

            pklog.warn("Filename to be analyzed: {}", file_name);
            ++c;
        }
        else if(cli=="-d"){
            dir_name = argv[c+1];
            vars.insert("dir", argv[c+1]);
            ++c;
        }
        else if(cli=="-p"){
            std::string cli(argv[c+1]);
            // vars.insert("pbc", toml::array{pbcs})
            pbcs = get_pbcs(cli);
            ++c;
        }
        else if(cli=="-rdf"){
            run_rdf = 1;
            std::string cli(argv[c+1]);
            rdf_cut = get_cutoffs(cli);
            ++c;
        }
        else if(cli=="-nf"){
            run_nf = 1;
            std::string cli(argv[c+1]);
            cutoffs = get_cutoffs(cli);
            // pklog.warn("Run neighbor construction with cutoff: {}", cutoffs[0]);
            ++c;
        }
        else if(cli=="-h" || cli=="--help"){
            show_help_msg();
            ++c;
        }
        else if(cli== "-v" || cli=="--version"){
            show_version_msg();
            ++c;
        }
        else if(cli== "-i" || cli=="--input"){
            // parser input config toml file.
            std::string cli(argv[c+1]);
            vars = toml::parse_file(cli);
            std::cout<<vars["basic"]<<std::endl;
            ++c;
        }
        else{
            show_input_err(argv[c]);
            ++c;
        };

    };
};

}

// void Flow::runAnalysis(){
//     Methods pipe;
//
//     // pipe.getCoord();
//     if(flag_nf){
//         // pipe.nf();
//         cout<<"Run neighbor finder .."<<endl;
//     };
//     if(flag_nfd){
//         // pipe.nf();
//         // neighborFinderDirect();
//         cout<<"Run neighbor finder with vectors.."<<endl;
//     };
//     if(flag_poly){
//         // pipe.nf();
//         // polyhedralAnalysis();
//         cout<<"Run polyhedral analysis .."<<endl;
//     };
// }
//
// void Methods::nf(mreal** xyz0, mreal lx, mreal ly, mreal lz, mreal rCut, bins** cellbins)
// {
//     vector<vector<int>> v2;
//
//     int natom;
//     mreal xyz[natom];
//
//
//
//     vector<int> vtr;
//     //初始化容器
//     for (int i = 0; i < 10; ++i)
//     {
//         vtr.push_back(i);
//     }
//     //利用迭代器遍历容器
//
//     cout << "方式1：";
//     for (auto it = vtr.begin(); it != vtr.end(); ++it)
//     {
//         cout << *it << " ";
//     }
//     cout << "\n方式1：";
//     for (auto it = begin(vtr); it != end(vtr); ++it)
//     {
//         cout << *it << " ";
//     }
//     cout << endl;
// }

