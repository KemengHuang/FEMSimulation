#include "SVD.h"

SVDResult2D_double SingularValueDecomposition2D_double(Matrix2d F)
{
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SVDResult2D_double result;
    Matrix2d tempU = svd.matrixU();
    Matrix2d tempV = svd.matrixV();
    Vector2d singVals = svd.singularValues();
    Matrix2d tempSigma = Matrix2d::Zero();
    tempSigma(0, 0) = singVals(0);
    tempSigma(1, 1) = singVals(1);

    //Matrix2f R = tempU * tempV.transpose();
    //float sign = R.determinant();

    if (tempU.determinant() < 0 && tempV.determinant() > 0)
    {
        tempU(0, 1) *= -1;
        tempU(1, 1) *= -1;
        tempSigma(1, 1) *= -1;
    }
    else if (tempV.determinant() < 0 && tempU.determinant() > 0)
    {
        tempV(0, 1) *= -1;
        tempV(1, 1) *= -1;
        tempSigma(1, 1) *= -1;
    }

    result.U = tempU;
    result.V = tempV;
    result.SIGMA = tempSigma;

    return result;
}

SVDResult3D_float SingularValueDecomposition3D_float(Matrix3f F)
{
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SVDResult3D_float result;
    Matrix3f tempU = svd.matrixU();
    Matrix3f tempV = svd.matrixV();

    Vector3f singVals = svd.singularValues();
    Matrix3f tempSigma = Matrix3f::Zero();
    tempSigma(0, 0) = singVals(0);
    tempSigma(1, 1) = singVals(1);
    tempSigma(2, 2) = singVals(2);

    if (tempU.determinant() < 0 && tempV.determinant() > 0)
    {
        tempU(0, 2) *= -1;
        tempU(1, 2) *= -1;
        tempU(2, 2) *= -1;
        tempSigma(2, 2) *= -1;
    }
    else if (tempV.determinant() < 0 && tempU.determinant() > 0)
    {
        tempV(0, 2) *= -1;
        tempV(1, 2) *= -1;
        tempV(2, 2) *= -1;
        tempSigma(2, 2) *= -1;
    }

    //if (tempSigma(0, 0) < tempSigma(1, 1))
    //{
    //    float tempRecord = tempSigma(0, 0);
    //    tempSigma(0, 0) = tempSigma(1, 1);
    //    tempSigma(1, 1) = tempRecord;
    //}

    result.U = tempU;
    result.V = tempV;
    result.SIGMA = tempSigma;

    return result;
}

SVDResult3D_double SingularValueDecomposition3D_double(Matrix3d F)
{
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SVDResult3D_double result;
    Matrix3d tempU = svd.matrixU();
    Matrix3d tempV = svd.matrixV();

    Vector3d singVals = svd.singularValues();
    Matrix3d tempSigma = Matrix3d::Zero();
    tempSigma(0, 0) = singVals(0);
    tempSigma(1, 1) = singVals(1);
    tempSigma(2, 2) = singVals(2);

    if (tempU.determinant() < 0 && tempV.determinant() > 0)
    {
        tempU(0, 2) *= -1;
        tempU(1, 2) *= -1;
        tempU(2, 2) *= -1;
        tempSigma(2, 2) *= -1;
    }
    else if (tempV.determinant() < 0 && tempU.determinant() > 0)
    {
        tempV(0, 2) *= -1;
        tempV(1, 2) *= -1;
        tempV(2, 2) *= -1;
        tempSigma(2, 2) *= -1;
    }

    //if (tempSigma(0, 0) < tempSigma(1, 1))
    //{
    //    double tempRecord = tempSigma(0, 0);
    //    tempSigma(0, 0) = tempSigma(1, 1);
    //    tempSigma(1, 1) = tempRecord;
    //}

    result.U = tempU;
    result.V = tempV;
    result.SIGMA = tempSigma;
    return result;
}