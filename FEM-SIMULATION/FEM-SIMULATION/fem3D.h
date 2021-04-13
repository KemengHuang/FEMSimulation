#pragma once
#ifndef FEM3D_H
#define FEM3D_H
#include "mesh.h"

void initMesh3D(mesh3D& mesh, int type, double scale);
void fem_implicit3D(mesh3D& mesh);
void fem_explicit3D(mesh3D& mesh);
void Projected_Newton3D(mesh3D& mesh, double& mfsum, int& total_cg_iterations, int& total_newton_iterations);

#endif // !FEM3D_H



