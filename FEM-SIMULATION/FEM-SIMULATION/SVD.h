#pragma once
#ifndef FEM_SVD_H
#define FEM_SVD_H

#include "Eigen/Eigen"
using namespace Eigen;

struct SVDResult2D_double
{
    Matrix2d U;
    Matrix2d SIGMA;
    Matrix2d V;
};

struct SVDResult3D_float
{
    Matrix3f U;
    Matrix3f SIGMA;
    Matrix3f V;
};

struct SVDResult3D_double
{
    Eigen::Matrix3d U;
    Eigen::Matrix3d SIGMA;
    Eigen::Matrix3d V;
};


SVDResult2D_double SingularValueDecomposition2D_double(Matrix2d F);
SVDResult3D_float SingularValueDecomposition3D_float(Matrix3f F);
SVDResult3D_double SingularValueDecomposition3D_double(Matrix3d F);

#endif //FEM_SVD_H