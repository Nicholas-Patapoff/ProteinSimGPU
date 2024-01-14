#include "Vector_Math.h"
#include "Env.h"
#include <math.h>

void vect_add(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& product){

    for(int i = 0; i < v1.size(); i++){
        product[i] = v1[i] + v2[i]; 
    }

}

void resize(std::vector<double>& vect, double scale){
    for(int i = 0; i < vect.size(); i++){
        vect[i] *= scale;
    }
}

void magnitude(std::vector<double>& object, double& mag){
    double temp = 0;
    for(int i = 0; i < object.size(); i++){
        temp+= object[i] *object[i];
    }
    temp = std::sqrt(temp);
    mag = temp;
}

void displacement_vect(std::vector<double>& d, std::vector<double>& Acoords , int atom1, int atom2){
    for(int i = 0; i < 3; i++){
        d.push_back(Acoords[atom1 * 3 + i] - Acoords[atom2 * 3 + i]);
    }
}

void unit_vector(double& mag, std::vector<double> d, std::vector<double>& unitv){
    
    for(int i = 0; i < d.size(); i++){
        unitv.push_back(d[i]/mag);
    }
}

void dot(std::vector<double>& disp1, std::vector<double>& disp2, double& val){
    for(int i = 0; i < 3; i++){
        val += disp1[i] * disp2[i];
    }
}

void cross(std::vector<double>& vect1, std::vector<double>& vect2, std::vector<double>& cprod){

    cprod.push_back(vect1[1] * vect2[2] - vect1[2] * vect2[1]);
    cprod.push_back(-(vect1[0] * vect2[2] - vect1[2] * vect2[0]));
    cprod.push_back(vect1[0] * vect2[1] - vect1[1] * vect2[0]);

}
