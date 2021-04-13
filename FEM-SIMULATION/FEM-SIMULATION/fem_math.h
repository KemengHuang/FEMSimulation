#pragma once
#ifndef FEM_MATH_H
#define FEM_MATH_H

#include "SVD.h"
#include<vector>

MatrixXd vec_double(MatrixXd F);
MatrixXf vec_float(MatrixXf F);
std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d);

#endif // !FEM_MATH_H

