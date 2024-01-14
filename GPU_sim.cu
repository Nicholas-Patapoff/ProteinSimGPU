#include "GPU_Data.h"
#include <cuda_runtime.h>
#include "GPU.h"
#include <cooperative_groups.h>
#include <stdio.h>
#include <iomanip>

//put constants onto GPU 
//set up GPU so that program runs entirely on gpu until 

//globals
std::fstream temp_file;


__global__ void force_calculation(float step_size, int Natoms, float* Mass,
 float* atomcoord, float* atomvel, float* atomforces,
 int* BWoutH, int BWoutH_size, int* BIH, int BIH_size, float* BForceC, float* BEQV, 
int* AWoutH, int AWoutH_size, int* AIH, int AIH_size, float* AForceC, float* AEQV, int* DWoutH,
int DWoutH_size,int* DincH,int DincH_size,float* DForceC,float* DPeriod,float* SCEE_SF,
float* SCNB_SF,float* DPhase);
__global__ void wait_completion();
__global__ void setZero(float* arr, int size);
__global__ void Vhalf_step(float* atomvel, float* atomforces, float* Mass, float step_size, int max);
__global__ void Vcoord_update(float* atomcoord, float* atomvel, float step_size, int max );
__global__ void harmonics(float step_size, int* atom_set, float* BForceC, float* BEQV, float* Acoords, float* forces, int max);
__global__ void angle_forces(float step_size, int* atom_set, float* AForceC, float* AEQB, float* Acoords, float* forces, int max);
__global__ void dih_forces(float step_size, int* atom_set, float* DForceC, float* DPeriod, float* SCNB_SF, float* DPhase, float* Acoords, float* forces, int max);
__global__ void export_init(float* buffer, float* exported, int size, int* flag, int flag_value);
__global__ void stopper(int* prog_flag);

__device__ void disp_vect(float* d, float* Acoords , int atom1, int atom2 );
__device__ void unit_vector(float* mag, float* d, float* unitv);
__device__ void magnitude(float* object, float* mag);
__device__ void cross (float* vect1, float* vect2, float* cprod );
__device__ void theta_from_dot(int* atom1, int* atom2, int* atom3, float* theta, float* Acoords);
__device__ void dot(float* disp1, float* disp2, float* val);
__device__ void DHrotTheta_from_dot(float* disp1, float* disp2, float* mag_disp1, float* mag_disp2, float* theta);
__device__ void vect_add(float* v1, float* v2, float* product, int size);
__device__ void resize(float* vect, float scale, int vect_length);

__global__ void sim_limited(int frames, float step_size, int Natoms, float* Mass,
 float* atomcoord, float* atomvel, float* atomforces,
 int* BWoutH, int BWoutH_size, int* BIH, int BIH_size, float* BForceC, float* BEQV, 
int* AWoutH, int AWoutH_size, int* AIH, int AIH_size, float* AForceC, float* AEQV, int* DWoutH,
int DWoutH_size,int* DincH,int DincH_size,float* DForceC,float* DPeriod,float* SCEE_SF,
float* SCNB_SF,float* DPhase, float* exp_coords,int* exp_counter, int export_num, int rows, float* coordbuffer, int* row_progress, int* prog_flag){



cudaError_t err;

//Vhalf_step init
int V_threads = 64;
int V_blocks = 1 + (Natoms / V_threads);
int current_row = 0;


cudaStream_t wait;

cudaStreamCreateWithFlags(&wait, cudaStreamNonBlocking);
err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("CUDA Error1: %s\n", cudaGetErrorString(err));
}
for(int i = 0; i < frames; i++){

//*prog_flag = Natoms;
Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, Mass, step_size, Natoms);
 err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("CUDA Error1: %s\n", cudaGetErrorString(err));
}

Vcoord_update<<<V_blocks, V_threads, 0, wait>>>(atomcoord, atomvel, step_size, Natoms);
err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("CUDA Error2: %s\n", cudaGetErrorString(err));
}


setZero<<< V_blocks, V_threads, 0, wait>>>(atomforces, Natoms);
err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("CUDA Error3: %s\n", cudaGetErrorString(err));
}


force_calculation<<<1, 1, 0, wait>>>(step_size, Natoms, Mass , atomcoord, atomvel, atomforces,  BWoutH, BWoutH_size, BIH, BIH_size,  BForceC, BEQV, 
AWoutH, AWoutH_size, AIH, AIH_size, AForceC, AEQV, DWoutH, DWoutH_size, DincH, DincH_size, DForceC, DPeriod, SCEE_SF, SCNB_SF, DPhase);

err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("CUDA Error4: %s\n", cudaGetErrorString(err));
}





Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, Mass, step_size, Natoms);
if (err != cudaSuccess) {
    printf("CUDA Error5: %s\n", cudaGetErrorString(err));
}

if(i % 10 == 0){ // need to make this adapt so that max queued items is less than items queued in this loop
    *prog_flag = 1;
   stopper<<<1, 1, 0, wait>>>(prog_flag);
while(*prog_flag != 0){
;
}
 
}



if( i % export_num == 0){


cudaMemcpyAsync(&coordbuffer[current_row * Natoms * 3], atomcoord, Natoms * 3 * sizeof(float), cudaMemcpyDeviceToDevice, wait);
current_row ++;
    if(current_row == rows){
        export_init<<<1, 1, 0, wait>>>(coordbuffer, exp_coords, Natoms * 3 * rows * sizeof(float), exp_counter, 1); 
        current_row = 0;
    }

}

}

if(current_row > 0){
    *row_progress = current_row;
    int temp = Natoms * 3 * current_row * sizeof(float);
        export_init<<<1, 1, 0, wait>>>(coordbuffer, exp_coords, temp, exp_counter, 2); 

} 


}


__global__ void export_init(float* buffer, float* exported, int size, int* flag, int flag_value){

    
    cudaStream_t wait;
    int t = flag_value;
    cudaStreamCreateWithFlags(&wait, cudaStreamNonBlocking);
    cudaMemcpyAsync(exported, buffer, size, cudaMemcpyDeviceToDevice, wait);
    *flag = flag_value;
    cudaStreamDestroy(wait);
   
    
    
    

}

__global__ void stopper(int* prog_flag) {
    
        atomicSub(prog_flag, 1);
}

__global__ void setZero(float* array, int max) {
    int ID = threadIdx.x + blockDim.x * blockIdx.x;
    if(ID < max){
        ID = ID * 3; 
        for(int atom = ID; atom < ID + 3; atom++){
                    array[atom] = 0.0f;

    }
        //atomicSub(prog_flag, 1);
    }
}


__global__ void force_calculation(float step_size, int Natoms, float* Mass,
 float* atomcoord, float* atomvel, float* atomforces,
 int* BWoutH, int BWoutH_size, int* BIH, int BIH_size, float* BForceC, float* BEQV, 
int* AWoutH, int AWoutH_size, int* AIH, int AIH_size, float* AForceC, float* AEQV, int* DWoutH,
int DWoutH_size,int* DincH,int DincH_size,float* DForceC,float* DPeriod,float* SCEE_SF,
float* SCNB_SF,float* DPhase){

int default_threads = 32;

int BWH_blocks = (BWoutH_size/3 + default_threads - 1) / default_threads;
int BIH_blocks = (BIH_size/3 + default_threads - 1) / default_threads;
int AWH_blocks = (AWoutH_size/4 + default_threads - 1) / default_threads;
int AIH_blocks = (AIH_size/4 + default_threads - 1) / default_threads;
int DWH_blocks = (DWoutH_size / 5 + default_threads - 1) / default_threads;
int DIH_blocks = (DincH_size / 5 + default_threads - 1) / default_threads;



dih_forces<<<DWH_blocks, default_threads, 0, cudaStreamFireAndForget>>> (step_size, DWoutH, DForceC, DPeriod, SCNB_SF, DPhase, atomcoord, atomforces, DWoutH_size / 5 );
dih_forces<<<DIH_blocks, default_threads, 0, cudaStreamFireAndForget>>> (step_size, DincH, DForceC, DPeriod, SCNB_SF, DPhase, atomcoord, atomforces, DincH_size / 5);
angle_forces<<<AWH_blocks, default_threads, 0, cudaStreamFireAndForget>>> (step_size, AWoutH, AForceC, AEQV, atomcoord, atomforces, AWoutH_size / 4);
angle_forces<<<AIH_blocks, default_threads, 0, cudaStreamFireAndForget>>> (step_size, AIH, AForceC, AEQV, atomcoord, atomforces, AIH_size / 4);
harmonics<<<BWH_blocks, default_threads, 0, cudaStreamFireAndForget>>>(step_size, BWoutH, BForceC, BEQV, atomcoord, atomforces, BWoutH_size / 3);
harmonics<<<BIH_blocks, default_threads, 0, cudaStreamFireAndForget>>>(step_size, BIH, BForceC, BEQV, atomcoord, atomforces, BIH_size / 3);





}

__global__ void Vhalf_step(float* atomvel, float* atomforces, float* Mass, float step_size, int max){
    int ID = threadIdx.x + blockDim.x * blockIdx.x;
    if(ID < max){
        ID = ID * 3; 
        for(int atom = ID; atom < ID + 3; atom++){

        atomvel[atom] = atomvel[atom] + (atomforces[atom] * step_size / (2 * (Mass[atom/3])));

        }
        //atomicSub(prog_flag, 1);
    }
    

}
__global__ void Vcoord_update(float* atomcoord, float* atomvel, float step_size, int max){
    int ID = threadIdx.x + blockDim.x * blockIdx.x;
    if(ID < max){
        ID = ID * 3;
    for(int atom = ID; atom < ID + 3; atom++){
    atomcoord[atom] = atomcoord[atom] + atomvel[atom] * step_size;
    }
        //atomicSub(prog_flag, 1);
}

}

__global__ void harmonics(float step_size, int* atom_set, float* BForceC, float* BEQV, float* Acoords, float* forces, int max){

int arrayloc = (threadIdx.x + blockDim.x * blockIdx.x);

if(arrayloc < max){
    //printf("%d ", arrayloc);
arrayloc = arrayloc * 3;

int atom1 = atom_set[arrayloc] / 3;
int atom2 = atom_set[arrayloc + 1] / 3;
int constIndex = atom_set[arrayloc + 2] - 1; 

float dispmag,  dispvectr[3], unitvect[3];
disp_vect(dispvectr, Acoords, atom1, atom2);


magnitude(dispvectr, &dispmag);
unit_vector(&dispmag, dispvectr, unitvect);

for(int i = 0; i < 3; i++){

    float force = -0.5 * (-BForceC[constIndex] * (dispvectr[i] - BEQV[constIndex] * unitvect[i]));
    atomicAdd(&forces[atom2 * 3 + i], force);
    atomicAdd(&forces[atom1 * 3 + i], -force);

    //printf("%f \n", force);
}
        //atomicSub(prog_flag, 1);
}

}

__global__ void angle_forces(float step_size, int* atom_set, float* AForceC, float* AEQB, float* Acoords, float* forces, int max){

int arrayloc = (threadIdx.x + blockDim.x * blockIdx.x);

if(arrayloc < max){
    
    arrayloc = arrayloc * 4;
    int atom1 = atom_set[arrayloc] / 3 ;
    int atom2 = atom_set[arrayloc + 1] / 3;
    int atom3 = atom_set[arrayloc + 2] / 3;
    int constIndex = atom_set[arrayloc + 3] - 1;
     //printf( "%d ", constIndex);
    float disp1[3], disp2[3], dispac[3];
    disp_vect(disp1, Acoords, atom1, atom2);
    disp_vect(disp2, Acoords, atom3, atom2);
    disp_vect(dispac, Acoords, atom1, atom3);

    //printf("%f %f %f", disp1[0], disp1[1], disp1[2]);

    float orth_abc[3], inp_ba[3], inp_bc[3];
    cross(disp1, disp2, orth_abc);
    cross(disp1, orth_abc, inp_ba);
    cross(disp2, orth_abc, inp_bc);

    float magba, magbc, mag_inp_ba, mag_inp_bc;
    magnitude(disp1, &magba);
    magnitude(disp2, &magbc);
    magnitude(inp_ba, &mag_inp_ba);
    magnitude(inp_bc, &mag_inp_bc);

    //printf("%f %f %f %f ", magba, magbc, mag_inp_ba, mag_inp_bc);

    float unit_inp_ba[3], unit_inp_bc[3];
    unit_vector(&mag_inp_ba, inp_ba, unit_inp_ba);
    unit_vector(&mag_inp_bc, inp_bc, unit_inp_bc);

    float theta; 
    theta_from_dot(&atom1, &atom2, &atom3, &theta, Acoords);
            
    float force_ba, force_bc;
    for(int i = 0; i < 3; i++){
        force_ba =  AForceC[constIndex] * (AEQB[constIndex] - theta)/magba * unit_inp_ba[i];
        
        force_bc =  -AForceC[constIndex] * (AEQB[constIndex] - theta)/magbc * unit_inp_bc[i];
        

        atomicAdd(&forces[atom1 * 3 + i],force_ba);
        atomicAdd(&forces[atom3 * 3 + i],force_bc);
        atomicAdd(&forces[atom2 * 3 + i], -force_ba - force_bc);
        
    }    
        //atomicSub(prog_flag, 1);
}

}



__global__ void dih_forces(float step_size, int* atom_set, float* DForceC, float* DPeriod, float* SCNB_SF, float* DPhase, float* Acoords, float* forces, int max){
int index = threadIdx.x + blockDim.x * blockIdx.x;
if(index < max){
index = index * 5; 
int atom1 = atom_set[index] / 3;
int atom2 = atom_set[index + 1] / 3;
int atom3 = atom_set[index + 2] / 3;
int atom4 = atom_set[index + 3] / 3;
int constIndex = atom_set[index + 4] - 1;




float dispba[3], dispbc[3], dispcb[3], dispcd[3];
disp_vect(dispba, Acoords, atom1, atom2);
disp_vect(dispbc, Acoords, abs(atom3), atom2);
disp_vect(dispcb, Acoords, atom2, abs(atom3));
disp_vect(dispcd, Acoords, abs(atom4), abs(atom3));

float magba, magbc, magcd, baDOTbc = 0, cdDOTbc = 0, magorthabc, magorthbcd;
magnitude(dispbc, &magbc);
magnitude(dispba, &magba);
magnitude(dispcd, &magcd);

float R[3], Rnorm[3], dispbcnorm[3], S[3], Snorm[3]; // (r(ij) - (r(ij) dot r(kj)))
unit_vector(&magbc, dispbc, dispbcnorm);
dot(dispba, dispbcnorm, &baDOTbc);
dot(dispbc, dispcd, &cdDOTbc);

for(int i = 0; i < 3; i++){
    R[i] = (dispba[i] - (baDOTbc * dispbcnorm[i]));
    S[i] = (dispcd[i] - (cdDOTbc * dispbcnorm[i]));
}

float magS, magR; 
magnitude(S, &magS);
magnitude(R, &magR);

unit_vector(&magS, S, Snorm);
unit_vector(&magR, R, Rnorm);

float Rnorm_DOT_Snorm = 0, DH_theta;
dot(Snorm, Rnorm, &Rnorm_DOT_Snorm);

DH_theta = acos(Rnorm_DOT_Snorm);


float ab_theta, cd_theta;
DHrotTheta_from_dot(dispba, dispcb, &magba, &magbc, &ab_theta);
DHrotTheta_from_dot(dispbc, dispcd, &magbc, &magcd, &cd_theta);


float orthaabc[3],orthabcd[3];
cross(dispba, dispbc, orthaabc);
cross(dispcd, dispcb, orthabcd);
magnitude(orthaabc, &magorthabc);
magnitude(orthabcd, &magorthbcd);

float northabc[3], northbcd[3], Fa[3], Fd[3];
unit_vector(&magorthabc, orthaabc, northabc);
unit_vector(&magorthbcd, orthabcd, northbcd);

for(int i = 0; i < 3; i++){

    Fa[i] = ( 0.5 * DPeriod[constIndex] * DForceC[constIndex] * sin(DPeriod[constIndex] * DH_theta + DPhase[constIndex])/(magba * sin(ab_theta)));
    Fd[i] = ( 0.5 * DPeriod[constIndex] * DForceC[constIndex] * sin(DPeriod[constIndex] * DH_theta + DPhase[constIndex])/(magcd * sin(cd_theta)));
    Fa[i] *= northabc[i];
    Fd[i] *= northbcd[i];

}

float dispoc[3], tc[3], tb[3], ocXFd[3], cdXFd[3], baXFa[3], Fc[3];
disp_vect(dispoc, Acoords, abs(atom3), atom2);
resize(dispoc, 0.5, 3);
resize(dispcd, 0.5, 3);
resize(dispba, 0.5, 3);

cross(dispoc, Fd, ocXFd);
cross(dispcd, Fd, cdXFd);
cross(dispba, Fa, baXFa);

for(int i = 0; i < 3; i++){
    tc[i] = 0;
    tb[i] = 0;
}

vect_add(ocXFd, cdXFd, tc, 3);
vect_add(tc, baXFa, tc, 3);

resize(tc, -1, 3);


float ocmag;
magnitude(dispoc, &ocmag);
ocmag = 1/(ocmag * ocmag);
resize(tc, ocmag, 3);

cross(tc, dispoc, Fc);

vect_add(Fa, Fd, tb, 3);
vect_add(tb, Fc, tb, 3);

resize(tb, -1, 3);


//SUM OF TORQUES
resize(dispba, 2, 3);
resize(dispcd, 2, 3);


float dispob[3], dispoa[3], dispod[3], dispoc2[3];     
for(int i = 0; i < 3; i++){ // can make more efficient with a proper init
    dispoa[i] = 0;
    dispod[i] = 0;
}
disp_vect(dispoc2, Acoords, abs(atom3), atom2);
disp_vect(dispob, Acoords, atom2, abs(atom3));
resize(dispoc2, 0.5, 3);
resize(dispob, 0.5, 3);
vect_add(dispba, dispob, dispoa, 3);
vect_add(dispoc2, dispcd, dispod, 3);

float torqoa[3], torqob[3], torqoc[3], torqod[3];

cross(dispoa, Fa, torqoa);
cross(dispob, tb, torqob);
cross(dispoc2, Fc, torqoc);
cross(dispod, Fd, torqod);

double total = 0;
for(int i = 0; i < 3; i++){
total += torqoa[i] + torqob[i] + torqoc[i] + torqod[i];
}
//SUM OF TROQUES


if(ab_theta){
    for(int i = 0; i < 3; i++){
        atomicAdd(&forces[atom1 * 3 + i], Fa[i]);
        atomicAdd(&forces[atom2 * 3 + i], tb[i]);
        atomicAdd(&forces[abs(atom3) * 3 + i], Fc[i]);
        atomicAdd(&forces[abs(atom4) * 3 + i], Fd[i]);
        }
    }

        //atomicSub(prog_flag, 1);

}


}


__device__ void vect_add(float* v1, float* v2, float* product, int size){

    for(int i = 0; i < size; i++){
        product[i] = v1[i] + v2[i]; 
    }

}

__device__ void resize(float* vect, float scale, int vect_length){

for(int i = 0; i < vect_length; i++){
        vect[i] *= scale;
    }

}
__device__ void DHrotTheta_from_dot(float* disp1, float* disp2, float* mag_disp1, float* mag_disp2, float* theta){
    float dotted = 0;
    dot(disp1, disp2, &dotted);
    double cosine_value = dotted / ((*mag_disp1) * (*mag_disp2));
    cosine_value = max((double)-1.0f, min((double)1.0f, cosine_value));
    *theta = acos(cosine_value);
}




__global__ void wait_completion(){ //method is used after force_calulation as a way to prevent race conditions for the final velocity algorithm
    ;
}
__device__ void disp_vect(float* d, float* Acoords , int atom1, int atom2 ){

for(int i = 0; i < 3; i++){
        d[i] = Acoords[atom1 * 3 + i] - Acoords[atom2 * 3 + i];
    
    //printf("%d, %d , coords %f %f ",atom1, atom2, Acoords[atom1 * 3 + i], Acoords[atom2 * 3 + i] );
}
}

__device__ void unit_vector(float* mag, float* d, float* unitv){
    
    for(int i = 0; i < 3; i++){
        unitv[i] = (d[i]/(*mag));
    }
}

__device__ void magnitude(float* object, float* mag){
    double temp = 0;
    for(int i = 0; i < 3; i++){
        temp+= object[i] * object[i];
    }
    temp = std::sqrt(temp);
    *mag = temp;
}

__device__ void cross (float* vect1, float* vect2, float* cprod ){
    cprod[0] = (vect1[1] * vect2[2] - vect1[2] * vect2[1]);
    cprod[1] = (-(vect1[0] * vect2[2] - vect1[2] * vect2[0]));
    cprod[2] = (vect1[0] * vect2[1] - vect1[1] * vect2[0]);
}

__device__ void theta_from_dot(int* atom1, int* atom2, int* atom3, float* theta, float* Acoords){
    float disp1[3], disp2[3];
    disp_vect(disp1, Acoords, *atom1, *atom2);
    disp_vect(disp2, Acoords, *atom3, *atom2);
    float magba;
    magnitude(disp1, &magba);
    float magbc;
    magnitude(disp2, &magbc);
    float dotac = 0;
    dot(disp1, disp2, &dotac);

    float cosine_value = dotac/(magba*magbc);
    cosine_value = max((double)-1.0f, min((double)1.0f, cosine_value));
    *theta = acos(cosine_value);
}

__device__ void dot(float* disp1, float* disp2, float* val){
    for(int i = 0; i < 3; i++){
        *val += disp1[i] * disp2[i];
    }
}

void gpu_sim(int frames, int exp_num, float step_size, Gdata final){
temp_file.open("/home/nich/Documents/GitHub/ProteinSim/gpuOUT.crd", std::ios::out);
temp_file << std::endl;


//mass
float* mass;

cudaMalloc((void**)&mass, sizeof(float) * final.MASS_size);

cudaMemcpy(mass, final.MASS, sizeof(float) * final.MASS_size, cudaMemcpyHostToDevice);

//atom forces, velocities, location
float* atomcoord; 
float* atomvel;
float* atomforces;
int* row_progress;

cudaMalloc((void**)&atomcoord, sizeof(float) * final.atom_coords_size);  //float
cudaMalloc((void**)&atomvel, sizeof(float) * final.atom_coords_size);    //float
cudaMalloc((void**)&atomforces, sizeof(float) * final.atom_coords_size); //float
cudaMallocManaged((void**)&row_progress, sizeof(int)); 

cudaMemcpy(atomcoord, final.atom_coords, sizeof(float) * final.atom_coords_size, cudaMemcpyHostToDevice);
cudaMemset(atomvel, 0, sizeof(float) * final.atom_coords_size);
cudaMemset(atomforces, 0, sizeof(float) * final.atom_coords_size);
cudaMemset(row_progress, 0, sizeof(int));




//bond force constants
int* BWoutH;
int BWoutH_size = final.BONDS_WITHOUT_HYDROGEN_size;
int* BIH;
int BIH_size = final.BONDS_INC_HYDROGEN_size;
float* BForceC;
float* BEQV;
  
cudaMalloc((void**)&BWoutH, sizeof(int) * final.BONDS_WITHOUT_HYDROGEN_size);  //int
cudaMalloc((void**)&BIH, sizeof(int) * final.BONDS_INC_HYDROGEN_size); //int
cudaMalloc((void**)&BForceC, sizeof(float) * final.BOND_FORCE_CONSTANT_size);    //float
cudaMalloc((void**)&BEQV, sizeof(float) * final.BOND_EQUIL_VALUE_size);  //float

cudaMemcpy(BWoutH, final.BONDS_WITHOUT_HYDROGEN, sizeof(int) * final.BONDS_WITHOUT_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(BIH, final.BONDS_INC_HYDROGEN, sizeof(int) * final.BONDS_INC_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(BForceC, final.BOND_FORCE_CONSTANT, sizeof(float) * final.BOND_FORCE_CONSTANT_size, cudaMemcpyHostToDevice);
cudaMemcpy(BEQV, final.BOND_EQUIL_VALUE, sizeof(float) * final.BOND_EQUIL_VALUE_size, cudaMemcpyHostToDevice);

//angle force constants
int* AWoutH;
int AWoutH_size = final.ANGLES_WITHOUT_HYDROGEN_size; 
int* AIH;
int AIH_size = final.ANGLES_INC_HYDROGEN_size;
float* AForceC;
float* AEQV;

cudaMalloc((void**)&AWoutH, sizeof(int) * final.ANGLES_WITHOUT_HYDROGEN_size); //int
cudaMalloc((void**)&AIH, sizeof(int) * final.ANGLES_INC_HYDROGEN_size); //int
cudaMalloc((void**)&AForceC, sizeof(float) * final.ANGLE_FORCE_CONSTANT_size);   //float
cudaMalloc((void**)&AEQV, sizeof(float) * final.ANGLE_EQUIL_VALUE_size); //float

cudaMemcpy(AWoutH, final.ANGLES_WITHOUT_HYDROGEN, sizeof(int) * final.ANGLES_WITHOUT_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(AIH, final.ANGLES_INC_HYDROGEN, sizeof(int) * final.ANGLES_INC_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(AForceC, final.ANGLE_FORCE_CONSTANT, sizeof(float) * final.ANGLE_FORCE_CONSTANT_size, cudaMemcpyHostToDevice);
cudaMemcpy(AEQV, final.ANGLE_EQUIL_VALUE, sizeof(float) * final.ANGLE_EQUIL_VALUE_size, cudaMemcpyHostToDevice);

//DIH force constants

int* DWoutH;
int DWoutH_size = final.DIHEDRALS_WITHOUT_HYDROGEN_size;
int* DincH;
int DincH_size = final.DIHEDRALS_INC_HYDROGEN_size;
float* DForceC;
float* DPeriod;
float* SCEE_SF;
float* SCNB_SF;
float* DPhase;

cudaMalloc((void**)&DWoutH, sizeof(int) * final.DIHEDRALS_WITHOUT_HYDROGEN_size); //int
cudaMalloc((void**)&DincH, sizeof(int) * final.DIHEDRALS_INC_HYDROGEN_size); //int
cudaMalloc((void**)&DForceC, sizeof(float) * final.DIHEDRAL_FORCE_CONSTANT_size);   //float
cudaMalloc((void**)&DPeriod, sizeof(float) * final.DIHEDRAL_PERIODICITY_size); //float
cudaMalloc((void**)&SCEE_SF, sizeof(int) * final.SCEE_SCALE_FACTOR_size); //int
cudaMalloc((void**)&SCNB_SF, sizeof(float) * final.SCNB_SCALE_FACTOR_size);   //float
cudaMalloc((void**)&DPhase, sizeof(float) * final.DIHEDRAL_PHASE_size); //float

cudaMemcpy(DWoutH, final.DIHEDRALS_WITHOUT_HYDROGEN, sizeof(int) * final.DIHEDRALS_WITHOUT_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(DincH, final.DIHEDRALS_INC_HYDROGEN, sizeof(int) * final.DIHEDRALS_INC_HYDROGEN_size, cudaMemcpyHostToDevice);
cudaMemcpy(DForceC, final.DIHEDRAL_FORCE_CONSTANT, sizeof(float) * final.DIHEDRAL_FORCE_CONSTANT_size, cudaMemcpyHostToDevice);
cudaMemcpy(DPeriod, final.DIHEDRAL_PERIODICITY, sizeof(float) * final.DIHEDRAL_PERIODICITY_size, cudaMemcpyHostToDevice);
cudaMemcpy(SCEE_SF, final.SCEE_SCALE_FACTOR, sizeof(int) * final.SCEE_SCALE_FACTOR_size, cudaMemcpyHostToDevice);
cudaMemcpy(SCNB_SF, final.SCNB_SCALE_FACTOR, sizeof(int) * final.SCNB_SCALE_FACTOR_size, cudaMemcpyHostToDevice);
cudaMemcpy(DPhase, final.DIHEDRAL_PHASE, sizeof(float) * final.DIHEDRAL_PHASE_size, cudaMemcpyHostToDevice);




//LJ forces
int* ATI;
int* NBPIndex;
float* LJA;
float* LJB; 

cudaMalloc((void**)&ATI, sizeof(int) * final.ATOM_TYPE_INDEX_size); //int
cudaMalloc((void**)&NBPIndex, sizeof(int) * final.NONBONDED_PARM_INDEX_size); //int
cudaMalloc((void**)&LJA, sizeof(float) * final.LENNARD_JONES_ACOEF_size);   //float
cudaMalloc((void**)&LJB, sizeof(float) * final.LENNARD_JONES_BCOEF_size); //float

cudaMemcpy(ATI, final.ATOM_TYPE_INDEX, sizeof(int) * final.ATOM_TYPE_INDEX_size, cudaMemcpyHostToDevice);
cudaMemcpy(NBPIndex, final.NONBONDED_PARM_INDEX, sizeof(int) * final.NONBONDED_PARM_INDEX_size, cudaMemcpyHostToDevice);
cudaMemcpy(LJA, final.LENNARD_JONES_ACOEF, sizeof(float) * final.LENNARD_JONES_ACOEF_size, cudaMemcpyHostToDevice);
cudaMemcpy(LJB, final.LENNARD_JONES_BCOEF, sizeof(float) * final.LENNARD_JONES_BCOEF_size, cudaMemcpyHostToDevice);

//multiframe buffer size determination
int rows;
if(sizeof(float) * final.atom_coords_size < 10000000){
    rows =   10000000 / (sizeof(float) * final.atom_coords_size);
} else{
    rows = 1;
}




//return and managed memory
int* exp_counter;
float* exp_coords;
float* buffer;
cudaMallocManaged(&exp_coords,rows * final.atom_coords_size * sizeof(float));
cudaMallocManaged(&exp_counter, sizeof(int));
cudaMalloc((void**)&buffer, sizeof(float) * final.atom_coords_size * rows);
*exp_counter = 0;


int atom_tot = final.atom_coords_size / 3;



struct timespec start, finish;
double elapsed;










//ALL CPU MANAGED

cudaStream_t export_hold;
cudaStreamCreateWithFlags(&export_hold, cudaStreamNonBlocking);

cudaStream_t stream;
cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);

// Get the start time





int* prog_flag;
cudaMalloc((void**) &prog_flag, sizeof(int));

//Vhalf_step init
int V_threads = 64;
int V_blocks = (atom_tot + V_threads - 1) / V_threads;
int current_row = 0;
cudaError_t err;

cudaStream_t wait;
cudaEvent_t event;
    
    cudaStreamCreateWithFlags(&wait, cudaStreamNonBlocking);
    cudaEventCreateWithFlags(&event, cudaEventDefault);
clock_gettime(CLOCK_MONOTONIC, &start);

for(int i = 0; i < frames; i++){

Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, mass, step_size, atom_tot);


Vcoord_update<<<V_blocks, V_threads, 0, wait>>>(atomcoord, atomvel, step_size, atom_tot);

//reset frame forces

setZero<<< V_blocks, V_threads, 0, wait>>>(atomforces, atom_tot);

//find new forces

force_calculation<<<1, 1, 0, wait>>>(step_size, atom_tot, mass , atomcoord, atomvel, atomforces,  BWoutH, BWoutH_size, BIH, BIH_size,  BForceC, BEQV, 
AWoutH, AWoutH_size, AIH, AIH_size, AForceC, AEQV, DWoutH, DWoutH_size, DincH, DincH_size, DForceC, DPeriod, SCEE_SF, SCNB_SF, DPhase );


Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, mass, step_size, 1);
err = cudaGetLastError();

cudaEventRecord(event, wait);
cudaStreamWaitEvent(wait, event);
int t = 0;
cudaMemcpy(prog_flag, &t, sizeof(int), cudaMemcpyHostToDevice);

if( i % atom_tot == 0){

cudaMemcpyAsync(&exp_coords[current_row * atom_tot * 3],atomcoord, atom_tot * 3 * sizeof(float), cudaMemcpyDeviceToHost, wait);
cudaEventRecord(event, wait);
cudaStreamWaitEvent(wait, event);

current_row ++;
    if(current_row == rows){
        exports(exp_coords, rows, atom_tot * 3 * rows * sizeof(float));
        current_row = 0;
    }

}

}

if(current_row > 0){

        exports(exp_coords, current_row, atom_tot * 3 * current_row * sizeof(float));

} 

clock_gettime(CLOCK_MONOTONIC, &finish);

// Calculate the elapsed time in nanoseconds
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  
    printf("Elapsed time GPU1: %.9f seconds\n", elapsed);








// GPU MANAGED 

clock_gettime(CLOCK_MONOTONIC, &start);

sim_limited<<<1, 1>>>(frames, step_size,  atom_tot,  mass,
 atomcoord,  atomvel,  atomforces,
 BWoutH,  BWoutH_size,  BIH,  BIH_size,  BForceC,  BEQV, 
 AWoutH, AWoutH_size,  AIH,  AIH_size,  AForceC,  AEQV,  DWoutH,
 DWoutH_size, DincH,DincH_size, DForceC, DPeriod, SCEE_SF,
 SCNB_SF, DPhase,  exp_coords, exp_counter,  exp_num, rows, buffer, row_progress, prog_flag);



while(*exp_counter != 2){
if(*exp_counter == 1){
    cudaDeviceSynchronize();
    exports(exp_coords, rows, atom_tot * 3 * rows * sizeof(float));
    *exp_counter = 0;
}
}

if(*exp_counter == 2){
    cudaDeviceSynchronize();
    exports(exp_coords, *row_progress, atom_tot * 3 * *row_progress * sizeof(float));
}

clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time GPU2: %.9f seconds\n", elapsed);
















//GPU GRAPH

bool graphCreated = false;
cudaGraph_t graph;
cudaGraphExec_t instance;

clock_gettime(CLOCK_MONOTONIC, &start);
for(int i = 0; i < frames; i++){
if(!graphCreated){
    cudaStreamBeginCapture(wait, cudaStreamCaptureModeGlobal);
Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, mass, step_size, atom_tot);

Vcoord_update<<<V_blocks, V_threads, 0, wait>>>(atomcoord, atomvel, step_size, atom_tot);

//reset frame forces


setZero<<< V_blocks, V_threads, 0, wait>>>(atomforces, atom_tot);
err = cudaGetLastError();


//find new forces

force_calculation<<<1, 1, 0, wait>>>(step_size, atom_tot, mass , atomcoord, atomvel, atomforces,  BWoutH, BWoutH_size, BIH, BIH_size,  BForceC, BEQV, 
AWoutH, AWoutH_size, AIH, AIH_size, AForceC, AEQV, DWoutH, DWoutH_size, DincH, DincH_size, DForceC, DPeriod, SCEE_SF, SCNB_SF, DPhase );


Vhalf_step<<<V_blocks, V_threads, 0, wait>>>(atomvel, atomforces, mass, step_size, 1);


cudaEventRecord(event, wait);
cudaStreamWaitEvent(wait, event);
int t = 0;
cudaMemcpy(prog_flag, &t, sizeof(int), cudaMemcpyHostToDevice);

if( i % atom_tot == 0){



cudaMemcpyAsync(&exp_coords[current_row * atom_tot * 3],atomcoord, atom_tot * 3 * sizeof(float), cudaMemcpyDeviceToHost, wait);


}

cudaStreamEndCapture(wait, &graph);
cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);
graphCreated = true;

   

}

cudaGraphLaunch(instance, wait);
cudaStreamSynchronize(wait);
if( i % atom_tot == 0){



cudaMemcpyAsync(&exp_coords[current_row * atom_tot * 3],atomcoord, atom_tot * 3 * sizeof(float), cudaMemcpyDeviceToHost, wait);
cudaEventRecord(event, wait);
cudaStreamWaitEvent(wait, event);
current_row ++;
if(current_row == rows){
        exports(exp_coords, rows, atom_tot * 3 * rows * sizeof(float));
        current_row = 0;
    }
}
 
}

if(current_row > 0){

exports(exp_coords, current_row, atom_tot * 3 * current_row * sizeof(float));

} 
clock_gettime(CLOCK_MONOTONIC, &finish);

elapsed = (finish.tv_sec - start.tv_sec);
elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  
printf("Elapsed time GPU3: %.9f seconds\n", elapsed);











cudaFree(atomcoord);
cudaFree(atomvel);
cudaFree(atomforces);
cudaFree(BWoutH);
cudaFree(BIH);
cudaFree(BForceC);
cudaFree(BEQV);
cudaFree(AWoutH);
cudaFree(AIH);
cudaFree(AForceC);
cudaFree(AEQV);
cudaFree(mass);
cudaFree(DWoutH);
cudaFree(DincH);
cudaFree(DForceC);
cudaFree(DPeriod);
cudaFree(DPhase);
cudaFree(SCEE_SF);
cudaFree(SCNB_SF);
cudaFree(ATI);
cudaFree(LJA);
cudaFree(LJB);
cudaFree(NBPIndex);
cudaFree(exp_coords);
cudaFree(exp_counter);
cudaFree(row_progress);







}



void exports(float* data, int rows, int size){

    int row_length = size / (rows * sizeof(float));
    for(int y = 0; y < rows; y++){
    if(!temp_file.is_open()){
        std::cout << "failed open external file" << std::endl;
    }
    for(int i = 0; i < size / (rows * sizeof(float)); i++){
        if(i % 10 == 0 && i != 0){
            temp_file << std::endl;
        }
        temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (data[i + (y * row_length)]);
        if(i == row_length - 1){
            temp_file << std::endl;
            for(int z = 1; z < 4; z++){
                temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (1.000);
            }
            
        }
        
    }
    temp_file << std::endl;
    }
}

