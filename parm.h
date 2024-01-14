#ifndef parm_H
#define parm_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <variant>
#include <unordered_map>

class parm{
private:

public: 
using T = std::variant<int, double, std::string>;
std::unordered_map<std::string, std::vector<T> > values;
parm(std::string pdb);
};


#endif