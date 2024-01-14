#include <iostream>
#include <fstream> 
#include <string>
#include <variant>
#include <unordered_map>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
#include "GPU.h"
#include <time.h>
#include "GPU_Data.h"






int main()
{

Environment coords("Data/chignolin_ext.pdb");
parm data("Data/chig_ext.parm7"); 
time_t start, end;
cpu_simulation small_probe(coords, data, 1, "test_prot");
start = time(NULL);
small_probe.update_coord(0.01, 10000, 1000);
end = time(NULL);
std::cout<< "Completed CPU in " << end - start << " seconds!" << std::endl;


gpu_sim(10000, 1000, 0.01, ready_data(coords, data));



return 0;
}





