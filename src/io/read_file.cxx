module;
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

import types;
import logger;
export module read_file;

export{

auto read_types(const std::vector<std::string>& typechar) -> std::vector<int>{
    // Create dict of types
    std::vector<std::string> type_name;
    int type_num=1;
    int natom=typechar.size();
    std::vector<int> partype(natom, 0);

    type_name.push_back(typechar[0]);
    for(int k=0;k<natom;++k){
        for(int i=0;i<type_num;++i){
            // std::cout<<typechar[k]<<" "<<type_name[i]<<std::endl;
            if(typechar[k]==type_name[i]){
                // Assign integer type.
                partype[k] = i+1;
                break;
            }else if(i==type_num-1 && typechar[k]!=type_name[i]){
                // Create new type.
                type_name.push_back(typechar[k]);
                type_num++;
            }
        }
    }

    for(size_t k=1;k<=type_num;++k){
        int number=std::count(partype.begin(), partype.end(), k);
        pklog.info("{} atoms with atom type {};", number, k);
    }

    return partype;
};

void read_xyz_file(std::string filename, Box& box, std::vector<Atom>& atoms){

    pklog.warn("Reading data file;");
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    int atom_num;
    std::string line;

    file >> atom_num;
    pklog.info("{} atoms read from file;", atom_num);
    std::getline(file, line);   // read the \n symbol.
    std::getline(file, line);   // read comment line.

    atoms.resize(atom_num);   // allocate space.

    std::vector<std::string> typechar(atom_num);

    // read atom xyz.
    for (int i = 0; i < atom_num; ++i) {
        file >> typechar[i] >> atoms[i].x >> atoms[i].y >> atoms[i].z;

        if(i==0){
            box.xmin = atoms[i].x;
            box.ymin = atoms[i].y;
            box.zmin = atoms[i].z;

            box.xmax = atoms[i].x;
            box.ymax = atoms[i].y;
            box.zmax = atoms[i].z;
        }
        box.xmin = std::min(box.xmin, atoms[i].x);
        box.ymin = std::min(box.ymin, atoms[i].y);
        box.zmin = std::min(box.zmin, atoms[i].z);

        box.xmax = std::max(box.xmax, atoms[i].x);
        box.ymax = std::max(box.ymax, atoms[i].y);
        box.zmax = std::max(box.zmax, atoms[i].z);
    }

    file.close();

    // read box info.
    pklog.info("Box size on dimension X: {}, {};", box.xmin, box.xmax);
    pklog.info("Box size on dimension Y: {}, {};", box.ymin, box.ymax);
    pklog.info("Box size on dimension Z: {}, {};", box.zmin, box.zmax);

    box.lx = box.xmax - box.xmin;
    box.ly = box.ymax - box.ymin;
    box.lz = box.zmax - box.zmin;
    pklog.info("Simulation box dimension: {}, {}, {};", box.lx, box.ly, box.lz);

    // read types.
    std::vector<int> ptype = read_types(typechar);

    // check if pair wise cutoff is given with correct number and report the value for each pair.

    return;
};
};
