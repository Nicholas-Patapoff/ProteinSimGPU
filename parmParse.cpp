#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <variant>
#include "parm.h"


parm::parm(std::string parm7file){
    std::fstream nfile;
    nfile.open(parm7file, std::ios::in);

    if(nfile.is_open()){

        std::string temp;
        std::string current_flag;
        std::string current_format;

        while(getline(nfile, temp)){
            
            if(temp.substr(0,5) == "%FLAG"){ // this is good but problematic for files which have comments between the flags/formats and thier identifiers
                std::string format = temp;
                std::istringstream splitstring(temp);
                splitstring >> current_flag; // saves as "%flag"
                splitstring >> current_flag; // saves as "flag_type"

                getline(nfile, temp);
                temp = temp.substr(7, 10); // makes a substring after %format and then separates by whitespace- incase of comments- likely not a perfect solution
                std::istringstream splitstring2(temp);
                splitstring2 >> current_format; //this saves as "(format_type)"

                if(current_format == "(20a4)" || current_format == "(1a80)"){
                    values.emplace(current_flag, std::vector<T>());  

                } else if(current_format == "(5E16.8)"){
                    values.emplace(current_flag, std::vector<T>());

                } else if(current_format == "(10I8)" || current_format == "(1I8)"){
                    values.emplace(current_flag, std::vector<T>());

                }
            } 
            
            else if(current_format == "(20a4)" && temp != ""){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;

                for(int i = 0; i < temp.size() - 1; i+= 4){
                    vec.push_back(temp.substr(i, 4));
                }
            }

            else if(current_format == "(5E16.8)" && temp != ""){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;

                for(int i = 0; i < temp.size() - 1; i+= 16){
                    vec.push_back(std::stof(temp.substr(i, 16)));
                }
            }

            else if(current_format == "(10I8)" && temp != ""){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;

                for(int i = 0; i < temp.size() - 1; i += 8){
                    vec.push_back(std::stoi(temp.substr(i, 8)));
                }
            }

            else if(current_format == "(1a80)" && temp != ""){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;
                vec.push_back(temp);
            }

            else if(current_format == "(1I8)" && temp != ""){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;
                vec.push_back(std::stoi(temp)); 

            }
        } 
    }
}