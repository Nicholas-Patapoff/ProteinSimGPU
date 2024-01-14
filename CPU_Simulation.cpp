#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
#include "Vector_Math.h"
#include <math.h>
#include <thread>
#include <iomanip>
#include <libgen.h>
#include <eigen3/Eigen/Dense>


#define RATTLE 1


cpu_simulation::cpu_simulation(Environment& temp1, parm& temp2, double step, std::string export_name){

    coord = std::make_unique<Environment>(temp1);
    top = std::make_unique<parm>(temp2);


    velocities = std::vector<double>(coord->Acoords.size(), 0);
    forces = std::vector<double>(coord->Acoords.size(), 0);


    temp_file.open( export_name +".crd", std::ios::out);
    temp_file << std::endl;
    }


void cpu_simulation::update_coord(double step_size, int frames, int export_step){
    //initialize atom properties
    std::unordered_map<std::string, std::vector<T> >::iterator ms = top->values.find("MASS");
        std::vector<T>& Mass = ms->second;

    
    //initialize bond force constants 
    std::unordered_map<std::string, std::vector<T> >::iterator bwh = top->values.find("BONDS_WITHOUT_HYDROGEN");
        std::vector<T>& BWoutH = bwh->second;
    
    std::unordered_map<std::string, std::vector<T> >::iterator bih = top->values.find("BONDS_INC_HYDROGEN");
        std::vector<T>& BIH = bih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bfc = top->values.find("BOND_FORCE_CONSTANT");
        std::vector<T>& BForceC = bfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bev = top->values.find("BOND_EQUIL_VALUE");
        std::vector<T>& BEQV = bev->second;


    //initialize angle force constants
    std::unordered_map<std::string, std::vector<T> >::iterator awh = top->values.find("ANGLES_WITHOUT_HYDROGEN");
        std::vector<T>& AWoutH = awh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aih = top->values.find("ANGLES_INC_HYDROGEN");
        std::vector<T>& AIH = aih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator afc = top->values.find("ANGLE_FORCE_CONSTANT");
        std::vector<T>& AForceC = afc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aev = top->values.find("ANGLE_EQUIL_VALUE");
        std::vector<T>& AEQV = aev->second;


    //initialize DIH force constants
    std::unordered_map<std::string, std::vector<T> >::iterator dih = top->values.find("DIHEDRALS_INC_HYDROGEN");
        std::vector<T>& DincH = dih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dwh = top->values.find("DIHEDRALS_WITHOUT_HYDROGEN");
        std::vector<T>& DWoutH = dwh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dfc = top->values.find("DIHEDRAL_FORCE_CONSTANT");
        std::vector<T>& DForceC = dfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dp = top->values.find("DIHEDRAL_PERIODICITY");
        std::vector<T>& DPeriod = dp->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dph = top->values.find("DIHEDRAL_PHASE");
        std::vector<T>& DPhase = dph->second;

    std::unordered_map<std::string, std::vector<T> >::iterator scee = top->values.find("SCEE_SCALE_FACTOR");
        std::vector<T>& SCEE_SF = scee->second;      

    std::unordered_map<std::string, std::vector<T> >::iterator scnb = top->values.find("SCNB_SCALE_FACTOR");
        std::vector<T>& SCNB_SF = scnb->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator lja = top->values.find("LENNARD_JONES_ACOEF");
        std::vector<T>& LJAC = lja->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator ljb = top->values.find("LENNARD_JONES_BCOEF");
        std::vector<T>& LJBC = ljb->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator ati = top->values.find("ATOM_TYPE_INDEX");
        std::vector<T>& ATI = ati->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator nbpi = top->values.find("NONBONDED_PARM_INDEX");
        std::vector<T>& NBPIndex = nbpi->second;
  

    //IMPORT EXCLUDED ATOMS 
    std::unordered_map<std::string, std::vector<T> >::iterator nea = top->values.find("NUMBER_EXCLUDED_ATOMS");
        std::vector<T>& NEA = nea->second;
    std::unordered_map<std::string, std::vector<T> >::iterator eal = top->values.find("EXCLUDED_ATOMS_LIST");
        std::vector<T>& EAL = eal->second;
    //IMPORT EXCLUDED ATOMS


    // INIT EXCLUSION LIST
    std::vector<std::vector<int> > excluded;
    find_excluded(NEA.size(), excluded, NEA, EAL);
    // INIT EXCLUSION LIST


    //INIT BOX PARAMETERS
    std::vector<double> COM;
    //find center of mass
    center_of_mass(Mass, coord->Acoords, COM);
    
    //measure atomic distances from center
    //farthest is found and R become atom + 15 anstrums 
    double insphereR, seg_length;
    farthest_from_center(coord->Acoords, COM, insphereR);
    insphereR += 15; //check units to assure that untis are in Angstroms
    seg_length = 2 * insphereR /(sqrtf64(3 * (5 + sqrt(5))));

    //generate face vectors 

    std::vector<Eigen::Vector3d> dodeca_faces = generateFaceVectors(insphereR);
    //INIT BOX PARAMETERS
    
    
    //INIT BOXED VECTORS
    //each vector will cover a 3d space and sort based on Z axis
    


    
    int reorder_atoms_freq;

    if(step_size >= 0.01){
        reorder_atoms_freq = 5;
    } else if (step_size < 0.01 && step_size >= 0.001){
        reorder_atoms_freq = 25;
    } else if(step_size < 0.001){
        reorder_atoms_freq = 50;
    }



for(int i = 0; i <= frames; i++){
    //refresh atoms in box vectors every reorder_atoms_freq
    if(i % reorder_atoms_freq == 0 ){
        reorder_atoms_freq;
    }
        VerletAlg(step_size,  BWoutH,  BIH,  BForceC, BEQV,  AWoutH,
  AIH, AForceC, AEQV,  DincH,  DWoutH,  DForceC, 
 DPeriod ,  DPhase,  SCEE_SF,  SCNB_SF, LJAC,  LJBC, ATI,
 NBPIndex, excluded);
        
        if(i % export_step == 0) {
            std::cout<< (double) i / frames * 100 << "% complete!" << std::endl;
            double total_f = 0;
            exports();
        }
    }




    
}


std::vector<Eigen::Vector3d> cpu_simulation::generateFaceVectors(double inradius) {
    // Define a vector to store the face vectors
    std::vector<Eigen::Vector3d> faceVectors;

    //define dih_angle and rotation angle
    double dih_angle = acosf64(-1.0/ sqrtf64(5)); 
    double rot_angle = M_PI_2f64 * 2/5; //72 deg

    Eigen::Vector3d startingFaceVector(0, 0, inradius);
    faceVectors.push_back(startingFaceVector);

        // Eigen's rotation matrices for rotations around the Z-axis and around an edge
    Eigen::Matrix3d rotationZ;
    rotationZ = Eigen::AngleAxisd(rot_angle, Eigen::Vector3d::UnitZ());

    Eigen::Matrix3d rotationEdge;
    rotationEdge = Eigen::AngleAxisd(dih_angle, Eigen::Vector3d::UnitX());

        // Generate face vectors for the 5 surrounding faces (the belt around the dodecahedron)
        for(int i = 0; i < 5; i++){
            startingFaceVector = rotationZ * startingFaceVector;
            faceVectors.push_back(startingFaceVector);
        }
    // generates top-most and bottom-most face vectors
    faceVectors.push_back(rotationEdge * startingFaceVector);
    faceVectors.push_back(rotationEdge.inverse() * startingFaceVector);


    for (int i = 1; i <= 5; ++i) {
        Eigen::Vector3d vector = faceVectors[i];
        Eigen::Matrix3d rotationAroundFaceVector;
        rotationAroundFaceVector = Eigen::AngleAxisd(dih_angle, Eigen::Vector3d::UnitZ().cross(vector));
        faceVectors.push_back(rotationAroundFaceVector * vector);
    }


    return faceVectors;
}


void cpu_simulation::center_of_mass(std::vector<T>& Mass,std::vector<double>& Acoords,std::vector<double>& center_of_mass ){
std::vector<double> added_coordsxmass, added_masses;
    added_coordsxmass = std::vector<double>(3, 0);
    added_masses.push_back(0);
    for(int i = 0; i < Acoords.size(); i += 3){
        added_coordsxmass[0] += Acoords[i] * std::get<double>(Mass[i/3]);
        added_coordsxmass[1] += Acoords[i + 1] * std::get<double>(Mass[i/3]);
        added_coordsxmass[2] += Acoords[i + 2] * std::get<double>(Mass[i/3]);
        added_masses[0] += std::get<double>(Mass[i/3]);
    }



    center_of_mass = std::vector<double>(3,0);
    center_of_mass[0] = added_coordsxmass[0]/added_masses[0];
    center_of_mass[1] = added_coordsxmass[1]/added_masses[0];
    center_of_mass[2] = added_coordsxmass[2]/added_masses[0];
}

void cpu_simulation::farthest_from_center(std::vector<double>& Acoords, std::vector<double> COM, double& max_dist){
std::vector<double> atom = std::vector<double>(3,0), current_disp;
            atom[0] = Acoords[0];
            atom[1] = Acoords[1];
            atom[2] = Acoords[2];
            for(int i = 0; i < 3; i++){
                current_disp.push_back(atom[i] - COM[i]);
            }
            magnitude(current_disp, max_dist);
    for(int i = 0; i < Acoords.size(); i += 3){
        double current;
            atom[0] = Acoords[i];
            atom[1] = Acoords[i + 1];
            atom[2] = Acoords[i + 2];
            for(int i = 0; i < 3; i++){
                current_disp.push_back(atom[i] - COM[i]);
            }
            magnitude(current_disp, current);
            if(current > max_dist){
                max_dist = current;
            }
    }

}

void cpu_simulation::find_excluded(int atoms, std::vector<std::vector<int> >& excluded, std::vector<T>& NEA, std::vector<T>& EAL){
        int count = 0;
    
    for(int i = 0; i < NEA.size(); i++){
        std::vector<int> temp;

        for(int y = count; y < count + std::get<int> (NEA[i]); y++){
            temp.push_back(std::get<int>(EAL[y]));

    }
        count += std::get<int> (NEA[i]);
        excluded.push_back(temp);
    }
}

void cpu_simulation::force_additions(std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex, std::vector<std::vector<int> >& excluded){ //all forces on protien for bond/FF interactions


    //bond forces
    
    for(int i = 0; i < BWoutH.size() - 1; i+=3){ //BWoutH.size()
        spring_force( std::get<int>(BWoutH[i]) / 3,std::get<int>(BWoutH[i + 1]) / 3, std::get<double>(BForceC[std::get<int>(BWoutH[i + 2]) - 1]), std::get<double>(BEQV[std::get<int>(BWoutH[i + 2]) - 1]) );
    
    }
    
    for(int i = 0; i < BIH.size(); i+=3){ // this is a bond which include hydrogen. I want to stiffen interaction so as to reduce computational complexity. I can set up a method for selecting a model and integrate it into there
        spring_force( std::get<int>(BIH[i]) / 3,std::get<int>(BIH[i + 1]) / 3, std::get<double>(BForceC[std::get<int>(BIH[i + 2]) - 1]), std::get<double>(BEQV[std::get<int>(BIH[i + 2]) - 1]) );
    }


    //angle forces
    for(int i = 0; i < AWoutH.size(); i+=4){
        angle_force( std::get<int>(AWoutH[i]) / 3,std::get<int>(AWoutH[i + 1]) / 3, std::get<int>(AWoutH[i + 2]) / 3, std::get<double>(AForceC[std::get<int>(AWoutH[i + 3]) - 1]), std::get<double>(AEQV[std::get<int>(AWoutH[i + 3]) - 1]));

    }

    for(int i = 0; i < AIH.size() - 1; i+=4){ 
        angle_force( std::get<int>(AIH[i]) / 3,std::get<int>(AIH[i + 1]) / 3, std::get<int>(AIH[i + 2]) / 3, std::get<double>(AForceC[std::get<int>(AIH[i + 3]) - 1]), std::get<double>(AEQV[std::get<int>(AIH[i + 3]) - 1]));

    }

    //Dihedral forces
    
    
    
    for(int i = 0; i < DWoutH.size(); i+=5){
        if((std::get<int>(DWoutH[i + 2]) / 3) >= 0){
            dihedral_force( std::get<int>(DWoutH[i]) / 3, std::get<int>(DWoutH[i + 1]) / 3,
             (std::get<int>(DWoutH[i + 2]) / 3), (std::get<int>(DWoutH[i + 3]) / 3), 
            std::get<double>(DForceC[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(DPeriod[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(SCEE_SF[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(SCNB_SF[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(DPhase[std::get<int>(DWoutH[i + 4]) - 1]));
            double LJA, LJB;
            int atom1, atom4;
            atom1 = std::get<int>(DWoutH[i]) / 3;
            atom4 = std::get<int>(DWoutH[i + 3]) / 3;
            int temp, ati1, ati2;
            ati1 = std::get<int>(ATI[atom1]);
            ati2 = std::get<int>(ATI[abs(atom4)]);

            temp = std::get<int>(NBPIndex[sqrt(NBPIndex.size()) * (ati1 - 1) + ati2 - 1]) - 1;

            LJB = std::get<double>(LJBC[temp]);
            LJA = std::get<double>(LJAC[temp]);
            DH_LJF(atom1, abs(atom4), std::get<double>(SCNB_SF[std::get<int>(DWoutH[i + 4]) - 1]), LJA, LJB);
        } else if((std::get<int>(DWoutH[i + 3]) / 3) < 0){
            improper_dih_force(std::get<int>(DWoutH[i]) / 3, std::get<int>(DWoutH[i + 1]) / 3,
             (std::get<int>(DWoutH[i + 2]) / 3), (std::get<int>(DWoutH[i + 3]) / 3), 
            std::get<double>(DForceC[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(DPeriod[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(SCEE_SF[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(SCNB_SF[std::get<int>(DWoutH[i + 4]) - 1]), 
            std::get<double>(DPhase[std::get<int>(DWoutH[i + 4]) - 1]));
        }


    }    

   
    for(int i = 0; i < DincH.size(); i+=5){

        
        if((std::get<int>(DincH[i + 2]) / 3) >= 0){
            dihedral_force( std::get<int>(DincH[i]) / 3,std::get<int>(DincH[i + 1]) / 3, 
            (std::get<int>(DincH[i + 2]) / 3), (std::get<int>(DincH[i + 3]) / 3), 
            std::get<double>(DForceC[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(DPeriod[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(SCEE_SF[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(SCNB_SF[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(DPhase[std::get<int>(DincH[i + 4]) - 1]));
        
            double LJA, LJB;
            int atom1, atom4;
            atom1 = abs(std::get<int>(DincH[i]) / 3);
            atom4 = abs(std::get<int>(DincH[i + 3]) / 3);

            int ati1, ati2;
            int temp;
            ati1 = std::get<int>(ATI[atom1]);
            ati2 = std::get<int>(ATI[abs(atom4)]);

            temp = std::get<int>(NBPIndex[sqrt(NBPIndex.size()) * (ati1 - 1) + ati2 - 1]) - 1 ;
            
            LJA = std::get<double>(LJAC[temp]);
            LJB = std::get<double>(LJBC[temp]);


            
 
            DH_LJF(atom1, abs(atom4), std::get<double>(SCNB_SF[std::get<int>(DincH[i + 4]) - 1]), LJA, LJB);
        } else if((std::get<int>(DincH[i + 3]) / 3) < 0){
            improper_dih_force(std::get<int>(DincH[i]) / 3,std::get<int>(DincH[i + 1]) / 3, 
            (std::get<int>(DincH[i + 2]) / 3), (std::get<int>(DincH[i + 3]) / 3), 
            std::get<double>(DForceC[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(DPeriod[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(SCEE_SF[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(SCNB_SF[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<double>(DPhase[std::get<int>(DincH[i + 4]) - 1]));
        }
    }
    
    //LJ/hbond forces
    /*
    
    
    
    for(int i = 0; i < ATI.size() - 1; i++){
        for(int y = i + 1; y < ATI.size(); y++){
            int check = 0;
            for(int exc = 0; exc < excluded[i].size(); exc++){
                if(excluded[i][exc] == y + 1){
                    check = 1;
                }
                }
            if(check == 1){
                ;
            } else {
            double A, B;
            int a1 = i; 
            int a2 = y;
            a1 = std::get<int>(ATI[a1]);
            a2 = std::get<int>(ATI[a2]);
            int temp = std::get<int>(NBPIndex[(sqrt(NBPIndex.size()) * (a1 - 1) + a2) - 1]) - 1;
            if(temp >= 0){
            A = std::get<double>(LJAC[temp]); 
            B = std::get<double>(LJBC[temp]); 
            LennardJ_force(i, y, A, B);
            }
        
           
            }

        }
    }
    */
    
    

    
}


void cpu_simulation::spring_force(int atom1, int atom2, double k, double eq){
    std::vector<double> disp;
    double dmag;
    std::vector<double> dunit_vect;

    displacement_vect(disp,coord->Acoords, atom1, atom2);
    magnitude(disp, dmag);
    unit_vector(dmag, disp, dunit_vect);

    for(int i = 0; i < disp.size(); i ++){
        double force = -0.5 * (-k * (disp[i] - eq * dunit_vect[i]));
        forces[atom2 * 3 + i] += force;
        forces[atom1 * 3 + i] -= force;

    }

}



void cpu_simulation::angle_force(int atom1, int atom2, int atom3, double k, double eq){ //this needs to be optimized so as not to make disp and mags two times
    //first find the change of potential with respect the change of angle : 2k(thetaEQ - current_theta)
    //then find the change of potential with respect to the change of bond length
    // we must then make a nomalized vector which is in plane to atoms 1,2,3 and points in a direction orth to vector 1-2, and then 2-3
    // from there we can apply forces via F12 = - du/dr dot orth_vect12 which: -2k(thetaEQ - current_theta)/mag(vector1-2) dot orth_vect12.
    //F23 is the same 
    //then take the total of both as a negative ( - f12 - f23) and then we obtain force on atom 2. since -FA + -FC = FB;
    
    std::vector<double> disp1, disp2, dispac;
    displacement_vect(disp1,coord->Acoords, atom1, atom2);
    displacement_vect(disp2,coord->Acoords, atom3, atom2);
    displacement_vect(dispac,coord->Acoords, atom1, atom3);
    std::vector<double> orth_abc, inp_ba, inp_bc;
    cross(disp1, disp2, orth_abc);
    cross(disp1, orth_abc, inp_ba);
    cross(disp2, orth_abc, inp_bc);

    double magba, magbc, mag_inp_ba, mag_inp_bc;
    magnitude(disp1, magba);
    magnitude(disp2, magbc);
    magnitude(inp_ba, mag_inp_ba);
    magnitude(inp_bc, mag_inp_bc);

    std::vector<double> unit_inp_ba, unit_inp_bc;
    unit_vector(mag_inp_ba, inp_ba, unit_inp_ba);
    unit_vector(mag_inp_bc, inp_bc, unit_inp_bc);

    double theta; 
    theta_from_dot(atom1, atom2, atom3, theta);
    //std::cout << theta << " ";
    double force_ba, force_bc;
    for(int i = 0; i < unit_inp_ba.size(); i++){
        force_ba =  k * (eq - theta)/magba * unit_inp_ba[i];
        force_bc =  -k * (eq - theta)/magbc * unit_inp_bc[i];
        forces[atom1 * 3 + i] += force_ba;
        forces[atom3 * 3 + i] += force_bc;
        forces[atom2 * 3 + i] += -force_ba - force_bc;
        
    }    

}

void cpu_simulation::improper_dih_force(int atom1, int atom2, int atom3, int atom4, double k,double period, double sceef, double scnbf, double phase){
std::vector<double> dispba, dispbc, dispcd;
atom4 = abs(atom4);
displacement_vect(dispba,coord->Acoords, abs(atom3), atom1);
displacement_vect(dispbc,coord->Acoords, abs(atom3), atom2);
displacement_vect(dispcd,coord->Acoords, abs(atom3), atom4);

double magba, magbc, magcd, baDOTbc = 0, cdDOTbc = 0, magorthacb, magorthbcd, magorthdca;
magnitude(dispbc, magbc);
magnitude(dispba, magba);
magnitude(dispcd, magcd);

std::vector<double> orthaacb,orthabcd, orthadca;
cross(dispba, dispbc, orthaacb);
cross(dispbc, dispcd, orthabcd);
cross(dispcd, dispba, orthadca);
magnitude(orthaacb, magorthacb);
magnitude(orthabcd, magorthbcd);
magnitude(orthadca, magorthdca);

std::vector<double> Northaacb,Northabcd, Northadca;
unit_vector(magorthacb, orthaacb, Northaacb);
unit_vector(magorthbcd, orthabcd, Northabcd);
unit_vector(magorthdca, orthadca, Northadca);

double dot1 = 0, dot2 = 0, dot3 = 0;
dot(dispcd, Northaacb, dot1);
dot(dispba, Northabcd, dot2);
dot(dispbc, Northadca, dot3);

double theta1, theta2, theta3;
theta1 = acosf64(dot1/magcd);
theta2 = acosf64(dot2/magba);
theta3 = acosf64(dot3/magbc);

std::vector<double> Fa, Fb, Fc, Fd;
for(int i = 0; i < 3; i++){ //0.5 * period * k * sin(period*DH_theta + phase)
    Fa.push_back(0.5 * period * k * sin(period* (theta1 - 1.5708) + phase) * Northaacb[i]);
    Fb.push_back(0.5 * period * k * sin(period* (theta2 - 1.5708) + phase) * Northabcd[i]);
    Fd.push_back(0.5 * period * k * sin(period* (theta3 - 1.5708) + phase) * Northadca[i]);
    Fc.push_back(-Fa[i] - Fb[i] - Fd[i]);

}

    for(int i = 0; i < 3; i++){
        forces[atom1 * 3 + i] += Fa[i];
        forces[atom2 * 3 + i] += Fb[i];
        forces[abs(atom3) * 3 + i] += Fc[i];
        forces[abs(atom4) * 3 + i] += Fd[i];
        }




}

void cpu_simulation::dihedral_force(int atom1, int atom2, int atom3, int atom4, double k, double period, double sceef, double scnbf, double phase){
//torsion potential is defined as U= 0.5[A1(1 + cos(theta)) + A2(1 - cos2theta)) + A3(1 + cost(3theta)) + A4]
//the partial derivative of tors w/ respect to pos of atom a (pa) is: dU/dpa = du/dtheta * dtheta/dpa
//dU/dpa = +/- 0.5(-A1sin(theta) + 2A2sing(2theta) - 3A3sin(3theta))
//dtheta/dpa = 1/(mag(end_vector) * sin(theta))
//steps to find forces on end atoms: 1.create unit vectors whose direction is orthogonal abc and bcd 2. find the angle between the planes 
//3. calculate the forces for atoms a and d usign (0.5/mag(ab)sin(theta)) * (1A1sin(theta) - 2A2sin(2theta) + 3A3(sin(theta)) * norm vector orth to plane. 

//1. finding DIH angle without cross products where theta = arccose(Rnorm dot Snorm)
std::vector<double> dispba, dispbc, dispcb, dispcd;
displacement_vect(dispba,coord->Acoords, atom1, atom2);
displacement_vect(dispbc,coord->Acoords, abs(atom3), atom2);
displacement_vect(dispcb,coord->Acoords, atom2, abs(atom3));
displacement_vect(dispcd,coord->Acoords, abs(atom4), abs(atom3));

double magba, magbc, magcd, baDOTbc = 0, cdDOTbc = 0, magorthabc, magorthbcd;
magnitude(dispbc, magbc);
magnitude(dispba, magba);
magnitude(dispcd, magcd);

std::vector<double> R, Rnorm, dispbcnorm, S, Snorm; // (r(ij) - (r(ij) dot r(kj)))
unit_vector(magbc, dispbc, dispbcnorm);
dot(dispba, dispbcnorm, baDOTbc);
dot(dispbc, dispcd, cdDOTbc);

for(int i = 0; i < 3; i++){
    R.push_back(dispba[i] - (baDOTbc * dispbcnorm[i]));
    S.push_back(dispcd[i] - (cdDOTbc * dispbcnorm[i]));
}

double magS, magR; 
magnitude(S, magS);
magnitude(R, magR);

unit_vector(magS, S, Snorm);
unit_vector(magR, R, Rnorm);


double Rnorm_DOT_Snorm = 0, DH_theta;
dot(Snorm, Rnorm, Rnorm_DOT_Snorm);

DH_theta = acosf64(Rnorm_DOT_Snorm);


double ab_theta, cd_theta;
DHrotTheta_from_dot(dispba, dispcb, magba, magbc, ab_theta);
DHrotTheta_from_dot(dispbc, dispcd, magbc, magcd, cd_theta);


std::vector<double> orthaabc,orthabcd;
cross(dispba, dispbc, orthaabc);
cross(dispcd, dispcb, orthabcd);
magnitude(orthaabc, magorthabc);
magnitude(orthabcd, magorthbcd);

std::vector<double> northabc, northbcd, Fa, Fd;
unit_vector(magorthabc, orthaabc, northabc);
unit_vector(magorthbcd, orthabcd, northbcd);

for(int i = 0; i < 3; i++){

    Fa.push_back( 0.5 * period * k * sin(period*DH_theta + phase)/(magba * sin(ab_theta)));
    Fd.push_back( 0.5 * period * k * sin(period*DH_theta + phase)/(magcd * sin(cd_theta)));
    Fa[i] *= northabc[i];
    Fd[i] *= northbcd[i];

}




// can I reduce the amount of cross products below? 
std::vector<double> dispoc, tc, tb, ocXFd, cdXFd, baXFa, Fc;
displacement_vect(dispoc, coord->Acoords, abs(atom3), atom2);
resize(dispoc, 0.5);
resize(dispcd, 0.5);
resize(dispba, 0.5);

cross(dispoc, Fd, ocXFd);
cross(dispcd, Fd, cdXFd);
cross(dispba, Fa, baXFa);

for(int i = 0; i < 3; i++){
    tc.push_back(0);
    tb.push_back(0);
}

vect_add(ocXFd, cdXFd, tc);
vect_add(tc, baXFa, tc);

resize(tc, -1);


double ocmag;
magnitude(dispoc, ocmag);
ocmag = 1/(ocmag * ocmag);
resize(tc, ocmag);

cross(tc, dispoc, Fc);

vect_add(Fa, Fd, tb);
vect_add(tb, Fc, tb);

resize(tb, -1);


//SUM OF TORQUES
resize(dispba, 2);
resize(dispcd, 2);


std::vector<double> dispob, dispoa, dispod, dispoc2;     
for(int i = 0; i < 3; i++){
    dispoa.push_back(0);
    dispod.push_back(0);
}
displacement_vect(dispoc2,coord->Acoords, abs(atom3), atom2);
displacement_vect(dispob, coord->Acoords, atom2, abs(atom3));
resize(dispoc2, 0.5);
resize(dispob, 0.5);
vect_add(dispba, dispob, dispoa);
vect_add(dispoc2, dispcd, dispod);

std::vector<double> torqoa, torqob, torqoc, torqod;

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
        forces[atom1 * 3 + i] += Fa[i];
        forces[atom2 * 3 + i] += tb[i];
        forces[abs(atom3) * 3 + i] += Fc[i];
        forces[abs(atom4) * 3 + i] += Fd[i];
        }
    }

}

void cpu_simulation::LennardJ_force(int atom1,int atom2, double LJA, double LJB){
std::vector<double> dispad, normad;
    displacement_vect(dispad, coord->Acoords, atom2, atom1);
    double magad;
    magnitude(dispad, magad);
    unit_vector(magad, dispad, normad);

    double Fa = 6/magad * (2 * LJA/(pow(magad, 12)) - LJB/pow(magad, 6));
for(int i = 0; i < 3; i++){
        forces[atom1 * 3 + i] -= Fa * normad[i];
        forces[atom2 * 3 + i] += Fa * normad[i];
    }

}

void cpu_simulation::DH_LJF(int atom1, int atom4, double SCNBF, double LJA, double LJB){

    std::vector<double> dispad, normad;
    displacement_vect(dispad, coord->Acoords, atom4, atom1);
    double magad;
    magnitude(dispad, magad);
    unit_vector(magad, dispad, normad);

    double Fa = 6/magad * (2 * LJA/(pow(magad, 12)) - LJB/pow(magad, 6))/SCNBF;

    


    if(SCNBF != 0){
        for(int i = 0; i < 3; i++){
        forces[atom1 * 3 + i] -=Fa * normad[i];
        forces[atom4 * 3 + i] +=Fa * normad[i];
    }
    }


}







void cpu_simulation::DHrotTheta_from_dot(std::vector<double>& disp1, std::vector<double>& disp2, double& mag_disp1, double& mag_disp2, double& theta){
    double dotted = 0;
    dot(disp1, disp2, dotted);
    double cosine_value = dotted / (mag_disp1 * mag_disp2);
    cosine_value = std::max((double)-1.0f, std::min((double)1.0f, cosine_value));
    theta = acos(cosine_value);
    
}


void cpu_simulation::DHtheta_from_dot(std::vector<double>& nplane1, std::vector<double>& nplane2, double np1mag, double np2mag, double& theta){
    double dot12 = 0;
    dot(nplane1, nplane2, dot12);
    double cosine_value = dot12 / (np1mag * np2mag);
    cosine_value = std::max((double)-1.0f, std::min((double)1.0f, cosine_value));
    theta = acos(cosine_value);
}



void cpu_simulation::theta_from_dot(int& atom1, int& atom2, int& atom3, double& theta){
    std::vector<double> disp1, disp2;
    displacement_vect(disp1, coord->Acoords,atom1, atom2);
    displacement_vect(disp2, coord->Acoords, atom3, atom2);
    double magba;
    magnitude(disp1, magba);
    double magbc;
    magnitude(disp2, magbc);
    double dotac = 0;
    dot(disp1, disp2, dotac);
    std::vector<double> ba_unitvect, bc_unitvect;

    double cosine_value = dotac/(magba*magbc);
    cosine_value = std::max((double)-1.0f, std::min((double)1.0f, cosine_value));
    theta = acos(cosine_value);
}



void cpu_simulation::VerletAlg(double& step_size, std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex, std::vector<std::vector<int> >& excluded){
    std::unordered_map<std::string, std::vector<T> >::iterator ms = top->values.find("MASS");
        std::vector<T>& Mass = ms->second;

for(int atom = 0; atom < velocities.size(); atom++){
    velocities[atom] = velocities[atom] + (forces[atom] * step_size/(2 * std::get<double>(Mass[atom/3])));
}
for(int atom = 0; atom < coord->Acoords.size(); atom++){
    coord->Acoords[atom] = coord->Acoords[atom] + velocities[atom]*step_size;
}

forces.assign(forces.size(), 0);

force_additions(BWoutH,  BIH,  BForceC, BEQV,  AWoutH,
 AIH, AForceC, AEQV,  DincH,  DWoutH,  DForceC, 
 DPeriod ,  DPhase,  SCEE_SF,  SCNB_SF, LJAC,  LJBC, ATI,
 NBPIndex, excluded);



for(int atom = 0; atom < velocities.size(); atom++){
   velocities[atom] = velocities[atom] + forces[atom]*step_size / (2 * std::get<double>(Mass[atom/3]));
}


}

void cpu_simulation::VerletAlgRattle(double& step_size, std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex, std::vector<std::vector<int> >& excluded){
    std::unordered_map<std::string, std::vector<T> >::iterator ms = top->values.find("MASS");
        std::vector<T>& Mass = ms->second;

std::vector<double> predicted_vel = std::vector<double>(0,velocities.size());
std::vector<double> predicted_loc = std::vector<double>(0,coord->Acoords.size());;

for(int atom = 0; atom < velocities.size(); atom++){
    velocities[atom] = (velocities[atom] + (forces[atom] * step_size/(2 * std::get<double>(Mass[atom/3]))));
}
for(int atom = 0; atom < coord->Acoords.size(); atom++){
    coord->Acoords[atom] = coord->Acoords[atom] + velocities[atom]*step_size;
}

forces.assign(forces.size(), 0);

force_additions(BWoutH,  BIH,  BForceC, BEQV,  AWoutH,
 AIH, AForceC, AEQV,  DincH,  DWoutH,  DForceC, 
 DPeriod ,  DPhase,  SCEE_SF,  SCNB_SF, LJAC,  LJBC, ATI,
 NBPIndex, excluded);



for(int atom = 0; atom < velocities.size(); atom++){
   velocities[atom] = velocities[atom] + forces[atom]*step_size / (2 * std::get<double>(Mass[atom/3]));
}

if(RATTLE == 1){
    
    for(int atom = 0; atom < coord->Acoords.size(); atom++){
    predicted_loc[atom] = coord->Acoords[atom] + velocities[atom] * step_size;
}




}
}




void cpu_simulation::Rattle_positions(std::vector<double>& predicted, const std::vector<double>& original, double step_size, const std::vector<T>& Mass) {
    // Your RATTLE position correction algorithm here



}








void cpu_simulation::exports(){
    if(!temp_file.is_open()){
        std::cout << "failed open external file" << std::endl;
    }
    for(int i = 0; i < coord->Acoords.size(); i++){
        if(i % 10 == 0 && i != 0){
            temp_file << std::endl;
        }
        temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (coord->Acoords[i]);
        if(i == coord->Acoords.size() - 1){
            temp_file << std::endl;
            for(int z = 1; z < 4; z++){
                temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (1.000);
            }
            
        }
        
    }
    temp_file << std::endl;

}



