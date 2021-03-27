#include "fem3D.h"
#include "fem_math.h"
#include<iostream>
#include <vector>
#include "fem_parameters.h"
using namespace std;
using namespace FEM;

void initMesh3D(mesh3D& mesh, int type, double scale) {
    mesh.InitMesh(type, scale);
}
