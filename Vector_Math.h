#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#include <vector>
#include <memory>

void vect_add(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& product);
void resize(std::vector<double>& vect, double scale);
void magnitude(std::vector<double>& object, double& mag);
void displacement_vect(std::vector<double>& d, std::vector<double>& Acoords, int atom1, int atom2);
void unit_vector(double& mag, std::vector<double> d, std::vector<double>& unitv);
void dot(std::vector<double>& disp1, std::vector<double>& disp2, double& val);
void cross(std::vector<double>& vect1, std::vector<double>& vect2, std::vector<double>& cprod);



#endif