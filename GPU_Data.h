
#include "Env.h"
#include "parm.h"
#include <vector>
#ifndef GPUD
#define GPUD


class Gextract{
using T = std::variant<int, double, std::string>;
public:
Gextract(Environment& coord, parm& top);
float* atom_coords;
int* BONDS_WITHOUT_HYDROGEN;
int* BONDS_INC_HYDROGEN;
float* BOND_FORCE_CONSTANT;
float* BOND_EQUIL_VALUE;
int* ANGLES_WITHOUT_HYDROGEN;
int* ANGLES_INC_HYDROGEN;
float* ANGLE_FORCE_CONSTANT;
float* ANGLE_EQUIL_VALUE;
int* DIHEDRALS_INC_HYDROGEN;
int* DIHEDRALS_WITHOUT_HYDROGEN;
float* DIHEDRAL_FORCE_CONSTANT;
float* DIHEDRAL_PERIODICITY;
float* DIHEDRAL_PHASE;
float* SCEE_SCALE_FACTOR;
float* SCNB_SCALE_FACTOR;
float* LENNARD_JONES_ACOEF;
float* LENNARD_JONES_BCOEF;
int* ATOM_TYPE_INDEX;
int* NONBONDED_PARM_INDEX;
int* NUMBER_EXCLUDED_ATOMS;
int* EXCLUDED_ATOMS_LIST;
float* MASS;

int atom_coords_size;
int BONDS_WITHOUT_HYDROGEN_size;
int BONDS_INC_HYDROGEN_size;
int BOND_FORCE_CONSTANT_size;
int BOND_EQUIL_VALUE_size;
int ANGLES_WITHOUT_HYDROGEN_size;
int ANGLES_INC_HYDROGEN_size;
int ANGLE_FORCE_CONSTANT_size;
int ANGLE_EQUIL_VALUE_size;
int DIHEDRALS_INC_HYDROGEN_size;
int DIHEDRALS_WITHOUT_HYDROGEN_size;
int DIHEDRAL_FORCE_CONSTANT_size;
int DIHEDRAL_PERIODICITY_size;
int DIHEDRAL_PHASE_size;
int SCEE_SCALE_FACTOR_size;
int SCNB_SCALE_FACTOR_size;
int LENNARD_JONES_ACOEF_size;
int LENNARD_JONES_BCOEF_size;
int ATOM_TYPE_INDEX_size;
int NONBONDED_PARM_INDEX_size;
int NUMBER_EXCLUDED_ATOMS_size;
int EXCLUDED_ATOMS_LIST_size;
int MASS_size;


};

class Gdata{
    public:
    Gdata();
float* atom_coords;
int* BONDS_WITHOUT_HYDROGEN;
int* BONDS_INC_HYDROGEN;
float* BOND_FORCE_CONSTANT;
float* BOND_EQUIL_VALUE;
int* ANGLES_WITHOUT_HYDROGEN;
int* ANGLES_INC_HYDROGEN;
float* ANGLE_FORCE_CONSTANT;
float* ANGLE_EQUIL_VALUE;
int* DIHEDRALS_INC_HYDROGEN;
int* DIHEDRALS_WITHOUT_HYDROGEN;
float* DIHEDRAL_FORCE_CONSTANT;
float* DIHEDRAL_PERIODICITY;
float* DIHEDRAL_PHASE;
float* SCEE_SCALE_FACTOR;
float* SCNB_SCALE_FACTOR;
float* LENNARD_JONES_ACOEF;
float* LENNARD_JONES_BCOEF;
int* ATOM_TYPE_INDEX;
int* NONBONDED_PARM_INDEX;
int* NUMBER_EXCLUDED_ATOMS;
int* EXCLUDED_ATOMS_LIST;
float* MASS;


int atom_coords_size;
int BONDS_WITHOUT_HYDROGEN_size;
int BONDS_INC_HYDROGEN_size;
int BOND_FORCE_CONSTANT_size;
int BOND_EQUIL_VALUE_size;
int ANGLES_WITHOUT_HYDROGEN_size;
int ANGLES_INC_HYDROGEN_size;
int ANGLE_FORCE_CONSTANT_size;
int ANGLE_EQUIL_VALUE_size;
int DIHEDRALS_INC_HYDROGEN_size;
int DIHEDRALS_WITHOUT_HYDROGEN_size;
int DIHEDRAL_FORCE_CONSTANT_size;
int DIHEDRAL_PERIODICITY_size;
int DIHEDRAL_PHASE_size;
int SCEE_SCALE_FACTOR_size;
int SCNB_SCALE_FACTOR_size;
int LENNARD_JONES_ACOEF_size;
int LENNARD_JONES_BCOEF_size;
int ATOM_TYPE_INDEX_size;
int NONBONDED_PARM_INDEX_size;
int NUMBER_EXCLUDED_ATOMS_size;
int EXCLUDED_ATOMS_LIST_size;
int MASS_size;

};

Gdata ready_data(Environment& coord, parm& top);

#endif