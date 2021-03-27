#pragma once
#ifndef FEM2D_H
#define FEM2D_H

#include "mesh.h"

//Matrix4d project_ARAP_H_2D(Matrix2d A);
//Matrix2d calculateDms2D_double(double vertexes[][2], int index[3], int i);
//Matrix3d calculateDms3D_double(double vertexes[][3], int index[4], int i);
//MatrixXd computePFPX2D_double(const Matrix2d& InverseDm);
void initMesh2D(mesh2D& mesh, int type, double scale);
void fem_implicit2D(mesh2D& mesh);
void fem_explicit2D(mesh2D& mesh);
void Projected_Newton2D(mesh2D& mesh);

#endif // !FEM_H

