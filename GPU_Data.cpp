#include "GPU_Data.h"
#include "Env.h"
#include "parm.h"
#include <vector>

Gextract::Gextract(Environment& coord, parm& top){

 
std::unordered_map<std::string, std::vector<T> >::iterator ms = top.values.find("MASS");
        std::vector<T>& Mass = ms->second;

    
    //initialize bond force constants 
    std::unordered_map<std::string, std::vector<T> >::iterator bwh = top.values.find("BONDS_WITHOUT_HYDROGEN");
        std::vector<T>& BWoutH = bwh->second;
    
    std::unordered_map<std::string, std::vector<T> >::iterator bih = top.values.find("BONDS_INC_HYDROGEN");
        std::vector<T>& BIH = bih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bfc = top.values.find("BOND_FORCE_CONSTANT");
        std::vector<T>& BForceC = bfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bev = top.values.find("BOND_EQUIL_VALUE");
        std::vector<T>& BEQV = bev->second;


    //initialize angle force constants
    std::unordered_map<std::string, std::vector<T> >::iterator awh = top.values.find("ANGLES_WITHOUT_HYDROGEN");
        std::vector<T>& AWoutH = awh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aih = top.values.find("ANGLES_INC_HYDROGEN");
        std::vector<T>& AIH = aih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator afc = top.values.find("ANGLE_FORCE_CONSTANT");
        std::vector<T>& AForceC = afc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aev = top.values.find("ANGLE_EQUIL_VALUE");
        std::vector<T>& AEQV = aev->second;


    //initialize DIH force constants
    std::unordered_map<std::string, std::vector<T> >::iterator dih = top.values.find("DIHEDRALS_INC_HYDROGEN");
        std::vector<T>& DincH = dih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dwh = top.values.find("DIHEDRALS_WITHOUT_HYDROGEN");
        std::vector<T>& DWoutH = dwh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dfc = top.values.find("DIHEDRAL_FORCE_CONSTANT");
        std::vector<T>& DForceC = dfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dp = top.values.find("DIHEDRAL_PERIODICITY");
        std::vector<T>& DPeriod = dp->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dph = top.values.find("DIHEDRAL_PHASE");
        std::vector<T>& DPhase = dph->second;

    std::unordered_map<std::string, std::vector<T> >::iterator scee = top.values.find("SCEE_SCALE_FACTOR");
        std::vector<T>& SCEE_SF = scee->second;      

    std::unordered_map<std::string, std::vector<T> >::iterator scnb = top.values.find("SCNB_SCALE_FACTOR");
        std::vector<T>& SCNB_SF = scnb->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator lja = top.values.find("LENNARD_JONES_ACOEF");
        std::vector<T>& LJAC = lja->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator ljb = top.values.find("LENNARD_JONES_BCOEF");
        std::vector<T>& LJBC = ljb->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator ati = top.values.find("ATOM_TYPE_INDEX");
        std::vector<T>& ATI = ati->second;  

    std::unordered_map<std::string, std::vector<T> >::iterator nbpi = top.values.find("NONBONDED_PARM_INDEX");
        std::vector<T>& NBPIndex = nbpi->second;
  

    //IMPORT EXCLUDED ATOMS 
    std::unordered_map<std::string, std::vector<T> >::iterator nea = top.values.find("NUMBER_EXCLUDED_ATOMS");
        std::vector<T>& NEA = nea->second;
    std::unordered_map<std::string, std::vector<T> >::iterator eal = top.values.find("EXCLUDED_ATOMS_LIST");
        std::vector<T>& EAL = eal->second;
    //IMPORT EXCLUDED ATOMS

 atom_coords = new float[coord.Acoords.size()];
 BONDS_WITHOUT_HYDROGEN = new int[BWoutH.size()];
 BONDS_INC_HYDROGEN= new int[BIH.size()];
 BOND_EQUIL_VALUE= new float[BEQV.size()];
 BOND_FORCE_CONSTANT= new float[BForceC.size()];
 ANGLES_WITHOUT_HYDROGEN= new int[AWoutH.size()];
 ANGLES_INC_HYDROGEN= new int[AIH.size()];
 ANGLE_FORCE_CONSTANT= new float[AForceC.size()];
 ANGLE_EQUIL_VALUE= new float[AEQV.size()];
 DIHEDRALS_INC_HYDROGEN= new int[DincH.size()];
 DIHEDRALS_WITHOUT_HYDROGEN= new int[DWoutH.size()];
 DIHEDRAL_FORCE_CONSTANT= new float[DForceC.size()];
 DIHEDRAL_PERIODICITY= new float[DPeriod.size()];
 DIHEDRAL_PHASE= new float[DPhase.size()];
 SCEE_SCALE_FACTOR= new float[SCEE_SF.size()];
 SCNB_SCALE_FACTOR= new float[SCNB_SF.size()];
 LENNARD_JONES_ACOEF= new float[LJAC.size()];
 LENNARD_JONES_BCOEF= new float[LJBC.size()];
 ATOM_TYPE_INDEX= new int[ATI.size()];
 NONBONDED_PARM_INDEX= new int[NBPIndex.size()];
 NUMBER_EXCLUDED_ATOMS= new int[NEA.size()];
EXCLUDED_ATOMS_LIST= new int[EAL.size()];
MASS= new float[Mass.size()];




atom_coords_size = coord.Acoords.size();
for(int i = 0; i < coord.Acoords.size(); i++){
    atom_coords[i] = coord.Acoords[i];
}




BONDS_WITHOUT_HYDROGEN_size = BWoutH.size();
for(int i = 0; i < BWoutH.size(); i++){
    BONDS_WITHOUT_HYDROGEN[i] = std::get<int>(BWoutH[i]);

}



BONDS_INC_HYDROGEN_size = BIH.size();

for(int i = 0; i < BIH.size(); i++){
    BONDS_INC_HYDROGEN[i] = std::get<int>(BIH[i]);
}



BOND_FORCE_CONSTANT_size = BForceC.size();

for(int i = 0; i < BForceC.size(); i++){
    BOND_FORCE_CONSTANT[i] = std::get<double>(BForceC[i]);
}



BOND_EQUIL_VALUE_size = BEQV.size();
for(int i = 0; i < BEQV.size(); i++){
    BOND_EQUIL_VALUE[i] = std::get<double>(top.values.find("BOND_EQUIL_VALUE")->second[i]);
}



ANGLES_WITHOUT_HYDROGEN_size = AWoutH.size();
for(int i = 0; i < AWoutH.size(); i++){
    ANGLES_WITHOUT_HYDROGEN[i] = std::get<int>(top.values.find("ANGLES_WITHOUT_HYDROGEN")->second[i]);
}



ANGLES_INC_HYDROGEN_size = AIH.size();
for(int i = 0; i <AIH.size(); i++){
    ANGLES_INC_HYDROGEN[i] = std::get<int>(top.values.find("ANGLES_INC_HYDROGEN")->second[i]);
}



ANGLE_FORCE_CONSTANT_size = AForceC.size();
for(int i = 0; i < AForceC.size(); i++){
    ANGLE_FORCE_CONSTANT[i] = std::get<double>(top.values.find("ANGLE_FORCE_CONSTANT")->second[i]);
}




DIHEDRALS_INC_HYDROGEN_size = DincH.size();
for(int i = 0; i < DincH.size(); i++){
    DIHEDRALS_INC_HYDROGEN[i] = std::get<int>(top.values.find("DIHEDRALS_INC_HYDROGEN")->second[i]);
}



DIHEDRALS_WITHOUT_HYDROGEN_size = DWoutH.size();
for(int i = 0; i < DWoutH.size(); i++){
    DIHEDRALS_WITHOUT_HYDROGEN[i] = std::get<int>(top.values.find("DIHEDRALS_WITHOUT_HYDROGEN")->second[i]);
}


ANGLE_EQUIL_VALUE_size = AEQV.size();
for(int i = 0; i < AEQV.size(); i++){
    ANGLE_EQUIL_VALUE[i] = std::get<double>(top.values.find("ANGLE_EQUIL_VALUE")->second[i]);
}


DIHEDRAL_FORCE_CONSTANT_size = DForceC.size();
for(int i = 0; i < DForceC.size(); i++){
    DIHEDRAL_FORCE_CONSTANT[i] = std::get<double>(top.values.find("DIHEDRAL_FORCE_CONSTANT")->second[i]);
}

DIHEDRAL_PERIODICITY_size = DPeriod.size();
for(int i = 0; i < DPeriod.size(); i++){
    DIHEDRAL_PERIODICITY[i] = std::get<double>(top.values.find("DIHEDRAL_PERIODICITY")->second[i]);
}

DIHEDRAL_PHASE_size = DPhase.size();
for(int i = 0; i < top.values.find("DIHEDRAL_PHASE")->second.size(); i++){
    DIHEDRAL_PHASE[i] = std::get<double>(top.values.find("DIHEDRAL_PHASE")->second[i]);
}


SCEE_SCALE_FACTOR_size = SCEE_SF.size();
for(int i = 0; i < SCEE_SF.size(); i++){
    SCEE_SCALE_FACTOR[i] = std::get<double>(top.values.find("SCEE_SCALE_FACTOR")->second[i]);
}


SCNB_SCALE_FACTOR_size = SCNB_SF.size();
for(int i = 0; i < SCNB_SF.size(); i++){
    SCNB_SCALE_FACTOR[i] = std::get<double>(top.values.find("SCNB_SCALE_FACTOR")->second[i]);
}

LENNARD_JONES_ACOEF_size = LJAC.size();
for(int i = 0; i < LJAC.size(); i++){
    LENNARD_JONES_ACOEF[i] = std::get<double>(top.values.find("LENNARD_JONES_ACOEF")->second[i]);
}

LENNARD_JONES_BCOEF_size = LJBC.size();
for(int i = 0; i < LJBC.size(); i++){
    LENNARD_JONES_BCOEF[i] = std::get<double>(top.values.find("LENNARD_JONES_BCOEF")->second[i]);
}

ATOM_TYPE_INDEX_size = ATI.size();
for(int i = 0; i < ATI.size(); i++){
    ATOM_TYPE_INDEX[i] = std::get<int>(top.values.find("ATOM_TYPE_INDEX")->second[i]);
}

NONBONDED_PARM_INDEX_size = NBPIndex.size();
for(int i = 0; i < NBPIndex.size(); i++){
    NONBONDED_PARM_INDEX[i] = std::get<int>(top.values.find("NONBONDED_PARM_INDEX")->second[i]);
}

NUMBER_EXCLUDED_ATOMS_size = NEA.size();
for(int i = 0; i < NEA.size(); i++){
    NUMBER_EXCLUDED_ATOMS[i] = std::get<int>(top.values.find("NUMBER_EXCLUDED_ATOMS")->second[i]);
}

EXCLUDED_ATOMS_LIST_size = EAL.size();
for(int i = 0; i < EAL.size(); i++){
    EXCLUDED_ATOMS_LIST[i] = std::get<int>(top.values.find("EXCLUDED_ATOMS_LIST")->second[i]);
}


MASS_size = Mass.size();
for(int i = 0; i < Mass.size(); i++){
    MASS[i] = std::get<double>(Mass[i]);
}

};

Gdata::Gdata(){};

Gdata ready_data(Environment& coord, parm& top){

Gextract import(coord, top);
Gdata final;
final.atom_coords = import.atom_coords;
final.MASS = import.MASS;
final.BONDS_WITHOUT_HYDROGEN =import.BONDS_WITHOUT_HYDROGEN;
final.BONDS_INC_HYDROGEN =import.BONDS_INC_HYDROGEN;
final.BOND_EQUIL_VALUE =import.BOND_EQUIL_VALUE;
final.BOND_FORCE_CONSTANT = import.BOND_FORCE_CONSTANT;
final.ANGLES_WITHOUT_HYDROGEN = import.ANGLES_WITHOUT_HYDROGEN;
final.ANGLES_INC_HYDROGEN = import.ANGLES_INC_HYDROGEN;
final.ANGLE_FORCE_CONSTANT = import.ANGLE_FORCE_CONSTANT;
final.ANGLE_EQUIL_VALUE = import.ANGLE_EQUIL_VALUE;
final.DIHEDRALS_INC_HYDROGEN = import.DIHEDRALS_INC_HYDROGEN;
final.DIHEDRALS_WITHOUT_HYDROGEN = import.DIHEDRALS_WITHOUT_HYDROGEN;
final.DIHEDRAL_FORCE_CONSTANT = import.DIHEDRAL_FORCE_CONSTANT;
final.DIHEDRAL_PERIODICITY = import.DIHEDRAL_PERIODICITY;
final.DIHEDRAL_PHASE = import.DIHEDRAL_PHASE;
final.SCEE_SCALE_FACTOR = import.SCEE_SCALE_FACTOR;
final.SCNB_SCALE_FACTOR = import.SCNB_SCALE_FACTOR;
final.LENNARD_JONES_ACOEF = import.LENNARD_JONES_ACOEF;
final.LENNARD_JONES_BCOEF = import.LENNARD_JONES_BCOEF;
final.ATOM_TYPE_INDEX = import.ATOM_TYPE_INDEX;
final.NONBONDED_PARM_INDEX = import.NONBONDED_PARM_INDEX;
final.NUMBER_EXCLUDED_ATOMS = import.NUMBER_EXCLUDED_ATOMS;
final.MASS = import.MASS;
final.EXCLUDED_ATOMS_LIST = import.EXCLUDED_ATOMS_LIST;

//make sizes for each array
final.atom_coords_size = import.atom_coords_size;
final.MASS_size = import.MASS_size;
final.BONDS_WITHOUT_HYDROGEN_size =import.BONDS_WITHOUT_HYDROGEN_size;
final.BONDS_INC_HYDROGEN_size =import.BONDS_INC_HYDROGEN_size;
final.BOND_EQUIL_VALUE_size =import.BOND_EQUIL_VALUE_size;
final.BOND_FORCE_CONSTANT_size = import.BOND_FORCE_CONSTANT_size;
final.ANGLES_WITHOUT_HYDROGEN_size = import.ANGLES_WITHOUT_HYDROGEN_size;
final.ANGLES_INC_HYDROGEN_size = import.ANGLES_INC_HYDROGEN_size;
final.ANGLE_FORCE_CONSTANT_size = import.ANGLE_FORCE_CONSTANT_size;
final.ANGLE_EQUIL_VALUE_size = import.ANGLE_EQUIL_VALUE_size;
final.DIHEDRALS_INC_HYDROGEN_size = import.DIHEDRALS_INC_HYDROGEN_size;
final.DIHEDRALS_WITHOUT_HYDROGEN_size = import.DIHEDRALS_WITHOUT_HYDROGEN_size;
final.DIHEDRAL_FORCE_CONSTANT_size = import.DIHEDRAL_FORCE_CONSTANT_size;
final.DIHEDRAL_PERIODICITY_size = import.DIHEDRAL_PERIODICITY_size;
final.DIHEDRAL_PHASE_size = import.DIHEDRAL_PHASE_size;
final.SCEE_SCALE_FACTOR_size = import.SCEE_SCALE_FACTOR_size;
final.SCNB_SCALE_FACTOR_size = import.SCNB_SCALE_FACTOR_size;
final.LENNARD_JONES_ACOEF_size = import.LENNARD_JONES_ACOEF_size;
final.LENNARD_JONES_BCOEF_size = import.LENNARD_JONES_BCOEF_size;
final.ATOM_TYPE_INDEX_size = import.ATOM_TYPE_INDEX_size;
final.NONBONDED_PARM_INDEX_size = import.NONBONDED_PARM_INDEX_size;
final.NUMBER_EXCLUDED_ATOMS_size = import.NUMBER_EXCLUDED_ATOMS_size;
final.MASS_size = import.MASS_size;
final.EXCLUDED_ATOMS_LIST_size = import.EXCLUDED_ATOMS_LIST_size;
return final;



}