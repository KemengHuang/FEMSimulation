#pragma once
#ifndef FEM3D_H
#define FEM3D_H
#include "mesh.h"

void initMesh3D(mesh3D& mesh, int type, double scale);
void fem_implicit3D(mesh3D& mesh);
void fem_explicit3D(mesh3D& mesh);
void Projected_Newton3D(mesh3D& mesh);

#endif // !FEM3D_H



