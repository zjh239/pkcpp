module;
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

import types;
import logger;
export module read_file;

export{

auto read_types(const std::vector<std::string>& typechar, int atoms_ntype) -> std::vector<int>{
    // Create dict of types
    std::vector<std::string> type_name;
    int type_num=1;
    int natom=typechar.size();
    std::vector<int> partype(natom, 0);

    type_name.push_back(typechar[0]);
    for(int k = 0; k < natom; ++k){
        for(int i = 0; i < type_num; ++i){
            // std::cout<<typechar[k]<<" "<<type_name[i]<<std::endl;
            if(typechar[k] == type_name[i]){
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
    atoms_ntype = type_num;
    pklog.info("{} types of atoms read from file;", atoms_ntype);
    for(size_t k=1; k<= type_num;++k){
        int number=std::count(partype.begin(), partype.end(), k);
        pklog.info("{} atoms with atom type {};", number, k);
    }

    return partype;
};

auto read_xyz_file(const std::string filename, Box& box) -> AtomList{

    pklog.warn("Reading data file;");
    std::ifstream file(filename);
    if (!file.is_open()) {
        pklog.fatal("Unable to open file {}, stop.", filename);
    }

    int atom_num;
    std::string line;

    file >> atom_num;
    pklog.info("{} atoms read from file;", atom_num);
    std::getline(file, line);   // read the \n symbol.
    std::getline(file, line);   // read comment line.

    AtomList atoms(atom_num);

    std::vector<std::string> typechar(atom_num);

    // read atom xyz.
    for (int i = 0; i < atom_num; ++i) {
        file >> typechar[i] >> atoms.x[i] >> atoms.y[i] >> atoms.z[i];

        if(i==0){
            box.xmin = atoms.x[i];
            box.ymin = atoms.y[i];
            box.zmin = atoms.z[i];

            box.xmax = atoms.x[i];
            box.ymax = atoms.y[i];
            box.zmax = atoms.z[i];
        }
        box.xmin = std::min(box.xmin, atoms.x[i]);
        box.ymin = std::min(box.ymin, atoms.y[i]);
        box.zmin = std::min(box.zmin, atoms.z[i]);

        box.xmax = std::max(box.xmax, atoms.x[i]);
        box.ymax = std::max(box.ymax, atoms.y[i]);
        box.zmax = std::max(box.zmax, atoms.z[i]);
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
    atoms.ptype = read_types(typechar, atoms.ntype);

    // check if pair wise cutoff is given with correct number and report the value for each pair.

    return atoms;
};
};
