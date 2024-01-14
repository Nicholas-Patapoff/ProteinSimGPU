#include <iostream>
#include <fstream>
#include <string>
#include "Env.h"

Environment::Environment(std::string pdb){
    std::fstream pdb_file;
    pdb_file.open(pdb, std::ios::in);
    if(pdb_file.is_open()){
        std::string temp;
        while(getline(pdb_file, temp)){
            //std::cout << temp << std::endl;
            if(temp.substr(0,4) == "ATOM"){
                append_name(temp.substr(13, 4));
                append_residue(temp.substr(18, 2));
                append_coords(temp.substr(31,8),temp.substr(39, 8),temp.substr(47, 8));
                //std::cout << "atom!" << std::endl;            
                }
            else if(temp.substr(0, 6) == "HETATM"){
                append_Hname(temp.substr(13, 4));
                append_Hresidue(temp.substr(18, 2));
                append_Hcoords(temp.substr(31,8),temp.substr(39, 8),temp.substr(47, 8));
                //class name.HETATM.append(string)
                //std::cout << "hetatm!" << std::endl;
            }
        }
    } else{
        printf("did not open\n");
    }
    
    
}

void Environment::append_name(std::string a){
    Aatom_name.push_back(a);
}
void Environment::append_Hname(std::string a){
    Hatom_name.push_back(a);
}
void Environment::append_residue(std::string a){
    Aresidue.push_back(a);
}
void Environment::append_Hresidue(std::string a){
    Hresidue.push_back(a);
}
void Environment::append_coords(std::string a, std::string b, std::string c){
    double x = std::stof(a);
    double y = std::stof(b);
    double z = std::stof(c);
    //std::cout << x << " " << y << " " << z << std::endl;
    Acoords.push_back(x);
    Acoords.push_back(y);
    Acoords.push_back(z);
}

void Environment::append_Hcoords(std::string a, std::string b, std::string c){
    double x = std::stof(a);
    double y = std::stof(b);
    double z = std::stof(c);
    Hcoords.push_back(x);
    Hcoords.push_back(y);
    Hcoords.push_back(z);
}
