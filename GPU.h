#include "GPU_Data.h"
 
#ifndef GPUs
#define GPUs

void gpu_sim(int frames, int exp_num, float step_size, Gdata GPU_ready);

void exports(float* data, int rows, int size);

#endif