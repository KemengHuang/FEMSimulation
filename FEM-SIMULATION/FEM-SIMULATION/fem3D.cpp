#include "fem3D.h"
#include "fem_math.h"
#include<iostream>
#include <vector>
#include "fem_parameters.h"
using namespace std;
using namespace FEM;

void checkDegenerateElements() {

}

double calculateVolum(const vector<Vector3d>& vertexes, const vector<uint64_t>& index) {
    int id0 = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;
    double o1x = vertexes[index[id1]][0] - vertexes[index[id0]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[id0]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[id0]][2];
    Vector3d OA = Vector3d(o1x, o1y, o1z);

    double o2x = vertexes[index[id2]][0] - vertexes[index[id0]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[id0]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[id0]][2];
    Vector3d OB = Vector3d(o2x, o2y, o2z);

    double o3x = vertexes[index[id3]][0] - vertexes[index[id0]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[id0]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[id0]][2];
    Vector3d OC = Vector3d(o3x, o3y, o3z);

    Vector3d heightDir = OA.cross(OB);
    double bottomArea = heightDir.norm();
    heightDir.normalize();
    
    double volum = bottomArea * abs(heightDir.dot(OC));
    return volum;
}

Matrix3d calculateDms3D_double(const vector<Vector3d>& vertexes, const vector<uint64_t>& index, const int& i) {
    int id1 = (i + 1) % 4;
    int id2 = (i + 2) % 4;
    int id3 = (i + 3) % 4;
    double o1x = vertexes[index[id1]][0] - vertexes[index[i]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[i]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[i]][2];

    double o2x = vertexes[index[id2]][0] - vertexes[index[i]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[i]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[i]][2];

    double o3x = vertexes[index[id3]][0] - vertexes[index[i]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[i]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[i]][2];

    Matrix3d M;
    M(0, 0) = o1x; M(0, 1) = o2x; M(0, 2) = o3x;
    M(1, 0) = o1y; M(1, 1) = o2y; M(1, 2) = o3y;
    M(2, 0) = o1z; M(2, 1) = o2z; M(2, 2) = o3z;

    return M;
}

MatrixXd computePFPX3D_double(const Matrix3d& InverseDm) {
    MatrixXd PFPX = MatrixXd::Zero(9, 12);
    double m = InverseDm(0, 0), n = InverseDm(0, 1), o = InverseDm(0, 2);
    double p = InverseDm(1, 0), q = InverseDm(1, 1), r = InverseDm(1, 2);
    double s = InverseDm(2, 0), t = InverseDm(2, 1), u = InverseDm(2, 2);
    double t1 = -(m + p + s);
    double t2 = -(n + q + t);
    double t3 = -(o + r + u);
    PFPX(0, 0) = t1;  PFPX(0, 3) = m;  PFPX(0, 6) = p;  PFPX(0, 9) = s;
    PFPX(1, 1) = t1;  PFPX(1, 4) = m;  PFPX(1, 7) = p;  PFPX(1, 10) = s;
    PFPX(2, 2) = t1;  PFPX(2, 5) = m;  PFPX(2, 8) = p;  PFPX(2, 11) = s;
    PFPX(3, 0) = t2;  PFPX(3, 3) = n;  PFPX(3, 6) = q;  PFPX(3, 9) = t;
    PFPX(4, 1) = t2;  PFPX(4, 4) = n;  PFPX(4, 7) = q;  PFPX(4, 10) = t;
    PFPX(5, 2) = t2;  PFPX(5, 5) = n;  PFPX(5, 8) = q;  PFPX(5, 11) = t;
    PFPX(6, 0) = t3;  PFPX(6, 3) = o;  PFPX(6, 6) = r;  PFPX(6, 9) = u;
    PFPX(7, 1) = t3;  PFPX(7, 4) = o;  PFPX(7, 7) = r;  PFPX(7, 10) = u;
    PFPX(8, 2) = t3;  PFPX(8, 5) = o;  PFPX(8, 8) = r;  PFPX(8, 11) = u;
    return PFPX;
}

void initMesh3D(mesh3D& mesh, int type, double scale) {
    //mesh.InitMesh(type, scale);
    string path = "tetrahedraMesh/cube5.msh";
    mesh.load_tetrahedraMesh(path, scale);
    
    //double angle = (PI / 2);// *explicit_time_step;

    //rotationTop << cos(angleYVelocity), 0, -sin(angleYVelocity), 0, 1, 0, sin(angleYVelocity), 0, cos(angleYVelocity);
    //Matrix3d rotation;
    //rotation << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
    if (type == 1) {
        for (int i = 0; i < mesh.vertexNum; i++) {
            mesh.vertexes[i] = mesh.vertexes[i] * 2 -Vector3d(scale, scale, scale);
            //mesh.vertexes[i][1] -= 0.8;
        }
    }
    double min = 10000000000000;
    for (int i = 0; i < mesh.tetrahedraNum; i++) {
        mesh.DM_triangle_inverse.push_back(calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[i], 0).inverse());
        double vlm = calculateVolum(mesh.vertexes, mesh.tetrahedras[i]) / 6;

        if (vlm < min) {
            min = vlm;
        }

        mesh.volum.push_back(vlm);
        mesh.masses[mesh.tetrahedras[i][0]] += vlm * density / 4;
        mesh.masses[mesh.tetrahedras[i][1]] += vlm * density / 4;
        mesh.masses[mesh.tetrahedras[i][2]] += vlm * density / 4;
        mesh.masses[mesh.tetrahedras[i][3]] += vlm * density / 4;
    }


    Matrix3d rotationTop, rotationDown;
    double angleYVelocity = (PI/2);// *explicit_time_step;

    rotationTop << cos(angleYVelocity), 0, -sin(angleYVelocity), 0, 1, 0, sin(angleYVelocity), 0, cos(angleYVelocity);
    rotationDown << cos(-angleYVelocity), 0, -sin(-angleYVelocity), 0, 1, 0, sin(-angleYVelocity), 0, cos(-angleYVelocity);


    //for (int ii = 0; ii < mesh.vertexNum; ii++) {
    //    if (mesh.vertexes[ii][0] > 0.39999) {
    //        //mesh.velocities[ii] *= 0;
    //        //mesh.velocities[ii][0] = 1;
    //        mesh.Constraints[ii].setZero();
    //        //mesh.d_velocities[ii][1] = 1;
    //        //continue;
    //        mesh.vertexes[ii] = rotationTop * mesh.vertexes[ii];
    //    }

    //    if (mesh.vertexes[ii][0] < -0.39999) {
    //        //mesh.velocities[ii] *= 0;
    //        //mesh.velocities[ii][0] = -1;
    //        mesh.Constraints[ii].setZero();
    //        //mesh.d_velocities[ii][1] = 1;
    //        //continue;
    //        mesh.vertexes[ii] = rotationDown * mesh.vertexes[ii];
    //    }
    //}

    cout << "   " << min;
    cout << endl;
}

Matrix3d computePEPF_StableNHK3D_double(const Matrix3d& F, const double& lengthRate, const double& volumRate) {
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();

    double u = lengthRate, r = volumRate;
    Matrix3d pI3pF;
    pI3pF(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    pI3pF(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    pI3pF(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
    
    pI3pF(1, 0) = F(2, 1) * F(0, 2) - F(2, 2) * F(0, 1);
    pI3pF(1, 1) = F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
    pI3pF(1, 2) = F(2, 0) * F(0, 1) - F(2, 1) * F(0, 0);

    pI3pF(2, 0) = F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2);
    pI3pF(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    pI3pF(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
    Matrix3d PEPF = u * F + (r * (S.determinant() - 1) - u) * pI3pF;
    return PEPF;
}

Matrix3d computePEPF_StableNHK3D_2_double(const Matrix3d& F, const double& lengthRate, const double& volumRate) {
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();

    double I2 = (S * S).trace();
    double I3 = S.determinant();

    double u = lengthRate, r = volumRate;
    Matrix3d pI3pF;
    pI3pF(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    pI3pF(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    pI3pF(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);

    pI3pF(1, 0) = F(2, 1) * F(0, 2) - F(2, 2) * F(0, 1);
    pI3pF(1, 1) = F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
    pI3pF(1, 2) = F(2, 0) * F(0, 1) - F(2, 1) * F(0, 0);

    pI3pF(2, 0) = F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2);
    pI3pF(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    pI3pF(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
    Matrix3d PEPF = u * (1 - 1 / (I2 + 1)) * F + (r * (I3 - 1 - u * 3 / (r * 4))) * pI3pF;
    return PEPF;
}

Matrix3d computePEPF_Aniostropic3D_double(const Matrix3d& F, Vector3d direction, const double& scale) {

    direction.normalize();
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, V, R, S, sigma;
    
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();

    double I4 = direction.transpose() * S * direction;
    double I5 = direction.transpose() * S.transpose() * S * direction;

    if (I4 == 0) {
        system("pause");
    }

    double s = 0;
    if (I4 < 0) {
        s = -1;
    }
    else if (I4 > 0) {
        s = 1;
    }

    Matrix3d PEPF = scale * (1 - s / sqrt(I5)) * F * direction * direction.transpose();
    return PEPF;
}

void boundary_process(mesh3D& mesh, const int& index, bool resetConstraint) {
    int constraintX = 0, constraintY = 0, constraintZ = 0;
    double d_x = 0, d_y = 0, d_z = 0;
    mesh.d_velocities[index] = Vector3d(0, 0, 0);
    if (resetConstraint) {
       mesh.Constraints[index].setIdentity();
    }
    if (mesh.vertexes[index][1] < -1) {
        d_y = -1 - mesh.vertexes[index][1];
        mesh.vertexes[index][1] = -0.9999;
        mesh.d_velocities[index][1] = mesh.velocities[index][1] * -0.1;
        mesh.velocities[index][1] = 0;
        constraintY = 1;
    }

    if (mesh.vertexes[index][0] < -1) {
        d_x = -1 - mesh.vertexes[index][0];
        mesh.vertexes[index][0] = -0.9999;
        mesh.d_velocities[index][0] = mesh.velocities[index][0] * -0.1;
        mesh.velocities[index][0] = 0;
        constraintX = 1;
    }

    if (mesh.vertexes[index][2] < -1) {
        d_z = -1 - mesh.vertexes[index][2];
        mesh.vertexes[index][2] = -0.9999;
        mesh.d_velocities[index][2] = mesh.velocities[index][2] * -0.1;
        mesh.velocities[index][2] = 0;
        constraintZ = 1;
    }

    if (mesh.vertexes[index][1] > 1) {
        d_y = 1 - mesh.vertexes[index][1];
        mesh.vertexes[index][1] = 0.9999;
        mesh.d_velocities[index][1] = mesh.velocities[index][1] * -0.1;
        mesh.velocities[index][1] = 0;
        constraintY = 1;
    }

    if (mesh.vertexes[index][0] > 1) {
        d_x = 1 - mesh.vertexes[index][0];
        mesh.vertexes[index][0] = 0.9999;
        mesh.d_velocities[index][0] = mesh.velocities[index][0] * -0.1;
        mesh.velocities[index][0] = 0;
        constraintX = 1;
    }

    if (mesh.vertexes[index][2] > 1) {
        d_z = 1 - mesh.vertexes[index][2];
        mesh.vertexes[index][2] = 0.9999;
        mesh.d_velocities[index][2] = mesh.velocities[index][2] * -0.1;
        mesh.velocities[index][2] = 0;
        constraintZ = 1;
    }
    Vector3d Dx(d_x, d_y, d_z);
    Vector3d Px(constraintX, 0, 0), Py(0, constraintY, 0), Pz(0, 0, constraintZ);
    mesh.d_positions[index] = Dx;
    mesh.Constraints[index] = mesh.Constraints[index] - Px * Px.transpose() - Py * Py.transpose() - Pz * Pz.transpose();
}

MatrixXd project_StabbleNHK_H_3D(const Matrix3d& F, const double& lengthRate, const double& volumRate) {
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, sigma, V, A;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    double I3 = sigma(0, 0) * sigma(1, 1) * sigma(2, 2);

    double u = lengthRate, r = volumRate;
    A.setZero();
    A(0, 0) = u + r * sigma(1, 1) * sigma(2, 2) * sigma(1, 1) * sigma(2, 2);
    A(0, 1) = sigma(2, 2) * (r * (2 * I3 - 1) - u); 
    A(0, 2) = sigma(1, 1) * (r * (2 * I3 - 1) - u);
    A(1, 0) = sigma(2, 2) * (r * (2 * I3 - 1) - u);
    A(1, 1) = u + r * sigma(0, 0) * sigma(2, 2) * sigma(0, 0) * sigma(2, 2);
    A(1, 2) = sigma(0, 0) * (r * (2 * I3 - 1) - u);
    A(2, 0) = sigma(1, 1) * (r * (2 * I3 - 1) - u);
    A(2, 1) = sigma(0, 0) * (r * (2 * I3 - 1) - u); 
    A(2, 2) = u + r * sigma(1, 1) * sigma(0, 0) * sigma(1, 1) * sigma(0, 0);

    double lamda[9];
    lamda[0] = A.eigenvalues()(0).real();
    lamda[1] = A.eigenvalues()(1).real();
    lamda[2] = A.eigenvalues()(2).real();
    
    lamda[3] = u + sigma(2, 2) * (r * (I3 - 1) - u);
    lamda[4] = u + sigma(0, 0) * (r * (I3 - 1) - u);
    lamda[5] = u + sigma(1, 1) * (r * (I3 - 1) - u);
    
    lamda[6] = u - sigma(2, 2) * (r * (I3 - 1) - u);
    lamda[7] = u - sigma(0, 0) * (r * (I3 - 1) - u);
    lamda[8] = u - sigma(1, 1) * (r * (I3 - 1) - u);
    Matrix3d D0, D1, D2;
    D0 << 1, 0, 0, 0, 0, 0, 0, 0, 0;
    D1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
    D2 << 0, 0, 0, 0, 0, 0, 0, 0, 1;

    D0 = U * D0 * V.transpose();
    D1 = U * D1 * V.transpose();
    D2 = U * D2 * V.transpose();

    Matrix3d Q[9];

    for (int i = 0; i < 3; i++) {
        double z0 = sigma(1, 1) * lamda[i] + sigma(0, 0) * sigma(2, 2);
        double z1 = sigma(0, 0) * lamda[i] + sigma(1, 1) * sigma(2, 2);
        double z2 = lamda[i] * lamda[i] - sigma(2, 2) * sigma(2, 2);
        Q[i] = z0 * D0 + z1 * D1 + z2 * D2;
    }
    
    Q[3] << 0, -1, 0, 1, 0, 0, 0, 0, 0;
    Q[4] << 0, 0, 0, 0, 0, 1, 0, -1, 0;
    Q[5] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    Q[6] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
    Q[7] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
    Q[8] << 0, 0, 1, 0, 0, 0, 1, 0, 0;

    double ml = 1 / sqrt(2);
    Q[3] = ml * U * Q[3] * V.transpose();
    Q[4] = ml * U * Q[4] * V.transpose();
    Q[5] = ml * U * Q[5] * V.transpose();
    Q[6] = ml * U * Q[6] * V.transpose();
    Q[7] = ml * U * Q[7] * V.transpose();
    Q[8] = ml * U * Q[8] * V.transpose();
    MatrixXd H = MatrixXd::Zero(9, 9);
    //H.setZero(9, 9);
    
    for (int i = 0; i < 9; i++) {
        if (lamda[i] > 0) {
            H += lamda[i] * vec_double(Q[i]) * vec_double(Q[i]).transpose();
        }
    }
    return H;
}

MatrixXd project_StabbleNHK_2_H_3D(const Matrix3d& F, const double& lengthRate, const double& volumRate) {
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, sigma, V, A;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    double I3 = sigma(0, 0) * sigma(1, 1) * sigma(2, 2);
    double I2 = sigma(0, 0) * sigma(0, 0) + sigma(1, 1) * sigma(1, 1) + sigma(2, 2) * sigma(2, 2);
    double g2 = sigma(0, 0) * sigma(1, 1) * sigma(0, 0) * sigma(1, 1) +
        sigma(0, 0) * sigma(2, 2) * sigma(0, 0) * sigma(2, 2) +
        sigma(2, 2) * sigma(1, 1) * sigma(2, 2) * sigma(1, 1);

    double u = lengthRate, r = volumRate;

    double n = 2 * u / ((I2 + 1) * (I2 + 1) * (r * (I3 - 1) - 3 * u / 4));
    double p = r / (r * (I3 - 1) - 3 * u / 4);
    double c2 = -g2 * p - I2 * n;
    double c1 = -(1 + 2 * I3 * p) * I2 - 6 * I3 * n + (g2 * I2 - 9 * I3 * I3) * p * n;
    double c0 = -(2 + 3 * I3 * p) * I3 + (I2 * I2 - 4 * g2) * n + 2 * I3 * p * n * (I2 * I2 - 3 * g2);

    vector<double> roots = NewtonSolverForCubicEquation(1, c2, c1, c0);
    if (roots.size() < 3) {
        system("pause");
    }

    Matrix3d D[3];
    double q[3];
    Matrix3d Q[9];
    double lamda[9];
    double Ut = u * (1 - 1 / (I2 + 1));
    double alpha = 1 + 3 * u / r / 4;

    for (int i = 0; i < 3; i++) {
        double alpha0 = roots[i] * (sigma(1, 1) + sigma(0, 0) * sigma(2, 2) * n + I3 * sigma(1, 1) * p) +
            sigma(0, 0) * sigma(2, 2) + sigma(1, 1) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1) + sigma(2, 2) * sigma(2, 2)) * n +
            I3 * sigma(0, 0) * sigma(2, 2) * p +
            sigma(0, 0) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2) *
            (sigma(1, 1) * sigma(1, 1) - sigma(2, 2) * sigma(2, 2)) * p * n;

        double alpha1 = roots[i] * (sigma(0, 0) + sigma(1, 1) * sigma(2, 2) * n + I3 * sigma(0, 0) * p) +
            sigma(1, 1) * sigma(2, 2) - sigma(0, 0) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1) - sigma(2, 2) * sigma(2, 2)) * n +
            I3 * sigma(1, 1) * sigma(2, 2) * p -
            sigma(1, 1) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2) *
            (sigma(0, 0) * sigma(0, 0) - sigma(2, 2) * sigma(2, 2)) * p * n;

        double alpha2 = roots[i] * roots[i] - roots[i] * (sigma(0, 0) * sigma(0, 0) + sigma(1, 1) * sigma(1, 1)) * (n + sigma(2, 2) * sigma(2, 2) * p) -
            sigma(2, 2) * sigma(2, 2) - 2 * I3 * n - 2 * I3 * sigma(2, 2) * sigma(2, 2) * p +
            ((sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2)) * ((sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2)) * p * n;

        q[i] = 1 / sqrt(alpha0 * alpha0 + alpha1 * alpha1 + alpha2 * alpha2);
        D[i] << alpha0, 0, 0, 0, alpha1, 0, 0, 0, alpha2;
        Q[i] = q[i] * U * D[i] * V.transpose();
        lamda[i] = r * (I3 - alpha) * roots[i] + Ut;
    }

    lamda[3] = Ut + sigma(2, 2) * r * (I3 - alpha);
    lamda[4] = Ut + sigma(0, 0) * r * (I3 - alpha);
    lamda[5] = Ut + sigma(1, 1) * r * (I3 - alpha);

    lamda[6] = Ut - sigma(2, 2) * r * (I3 - alpha);
    lamda[7] = Ut - sigma(0, 0) * r * (I3 - alpha);
    lamda[8] = Ut - sigma(1, 1) * r * (I3 - alpha);


    Q[3] << 0, -1, 0, 1, 0, 0, 0, 0, 0;
    Q[4] << 0, 0, 0, 0, 0, 1, 0, -1, 0;
    Q[5] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    Q[6] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
    Q[7] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
    Q[8] << 0, 0, 1, 0, 0, 0, 1, 0, 0;

    //for (int i = 0; i < 9; i++) {
    //    cout << "Q" << i << endl;
    //    cout << Q[i] << endl;
    //}

    double ml = 1 / sqrt(2);
    Q[3] = ml * U * Q[3] * V.transpose();
    Q[4] = ml * U * Q[4] * V.transpose();
    Q[5] = ml * U * Q[5] * V.transpose();
    Q[6] = ml * U * Q[6] * V.transpose();
    Q[7] = ml * U * Q[7] * V.transpose();
    Q[8] = ml * U * Q[8] * V.transpose();
    MatrixXd H = MatrixXd::Zero(9, 9);
    //H.setZero(9, 9);

    for (int i = 0; i < 9; i++) {
        if (lamda[i] > 0) {
            H += lamda[i] * vec_double(Q[i]) * vec_double(Q[i]).transpose();
        }
    }
    //cout << "H:" << endl;
    //cout << H << endl;
    return H;
}

MatrixXd project_ANIOSI5_H_3D(const Matrix3d& F, Vector3d direction, const double& scale) {
    direction.normalize();
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, sigma, V, S;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();
    double I4 = direction.transpose() * S * direction;
    double I5 = direction.transpose() * S.transpose() * S * direction;

    if (abs(I5) < 1e-15) return MatrixXd::Zero(9, 9);

    double s = 0;
    if (I4 < 0) {
        s = -1;
    }
    else if (I4 > 0) {
        s = 1;
    }

    double lamda0 = scale;
    double lamda1 = scale * (1 - s / sqrt(I5));
    double lamda2 = lamda1;
    //double lamda2 = lamda1;
    Matrix3d Q0, Q1, Q2, A;
    A = direction * direction.transpose();
    Q0 = (1 / sqrt(I5)) * F * A;

    Matrix3d Tx, Ty, Tz;
    
    Tx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
    Ty << 0, 0, -1, 0, 0, 0, 1, 0, 0;
    Tz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
    Vector3d directionM = V.transpose() * direction;

    Tx *= 1.f / sqrt(2.f);
    Ty *= 1.f / sqrt(2.f);
    Tz *= 1.f / sqrt(2.f);

    Q1 = U * Tx * sigma * V.transpose() * A;
    Q2 = (sigma(1, 1) * directionM[1]) * U * Tz * sigma * V.transpose() * A - (sigma(2, 2) * directionM[2]) * U * Ty * sigma * V.transpose() * A;

    MatrixXd H = lamda0 * vec_double(Q0) * vec_double(Q0).transpose();
    if (lamda1 > 0) {
        H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();
        H += lamda2 * vec_double(Q2) * vec_double(Q2).transpose();
    }
    
    return H;
}

double angle = 0;

void fem_explicit3D(mesh3D& mesh) {
    Vector3d direction = Vector3d(0, 1, 0);
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
        //double volum = calculateVolum(mesh.vertexes, mesh.tetrahedras[ii]);
        //cout << "F:\n" << F << endl;
        Matrix3d PEPF; PEPF.setZero();
        //PEPF += computePEPX_ARAP2D_double(F);

        PEPF += computePEPF_StableNHK3D_2_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic3D_double(F, direction, aniosScale);
        //PEPF += computePEPX_Aniostropic_double(F, Vector2d(1, 0), 0.3);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);

        //cout << "pepf:\n" << pepf << endl;
        //cout << "PFPX:\n" << PFPX << endl;
        MatrixXd f = -PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-9) {
                f(i, 0) = 0;
            }
            else {
                f(i, 0) = mesh.volum[ii] * f(i, 0);
            }
        }

        mesh.forces[mesh.tetrahedras[ii][0]][0] += f(0, 0);
        mesh.forces[mesh.tetrahedras[ii][0]][1] += f(1, 0);
        mesh.forces[mesh.tetrahedras[ii][0]][2] += f(2, 0);

        mesh.forces[mesh.tetrahedras[ii][1]][0] += f(3, 0);
        mesh.forces[mesh.tetrahedras[ii][1]][1] += f(4, 0);
        mesh.forces[mesh.tetrahedras[ii][1]][2] += f(5, 0);

        mesh.forces[mesh.tetrahedras[ii][2]][0] += f(6, 0);
        mesh.forces[mesh.tetrahedras[ii][2]][1] += f(7, 0);
        mesh.forces[mesh.tetrahedras[ii][2]][2] += f(8, 0);

        mesh.forces[mesh.tetrahedras[ii][3]][0] += f(9, 0);
        mesh.forces[mesh.tetrahedras[ii][3]][1] += f(10, 0);
        mesh.forces[mesh.tetrahedras[ii][3]][2] += f(11, 0);
    }

    //Matrix3d rotationTop, rotationDown;
    //double angleYVelocity = (PI / 12)* explicit_time_step;
    //angle += angleYVelocity;
    //if (angle < PI/2) {
    //    rotationTop << cos(angleYVelocity), 0, -sin(angleYVelocity), 0, 1, 0, sin(angleYVelocity), 0, cos(angleYVelocity);
    //    rotationDown << cos(-angleYVelocity), 0, -sin(-angleYVelocity), 0, 1, 0, sin(-angleYVelocity), 0, cos(-angleYVelocity);
    //}
    //else {
    //    rotationTop << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    //    rotationDown = rotationTop;
    //}


    //float sum = 0.f;
    //n++;
    //if (sqrt(sum) < 1e-4) return;
    //double gravity = -9.8;
    double gravity = -9.8;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {

        double acx = mesh.forces[ii][0] / mesh.masses[ii];
        double acy = mesh.forces[ii][1] / mesh.masses[ii] + gravity;
        double acz = mesh.forces[ii][2] / mesh.masses[ii];

        mesh.forces[ii][0] = 0;
        mesh.forces[ii][1] = 0;
        mesh.forces[ii][2] = 0;

        mesh.velocities[ii][0] += acx * explicit_time_step;
        mesh.velocities[ii][1] += acy * explicit_time_step;
        mesh.velocities[ii][2] += acz * explicit_time_step;

        //if (mesh.vertexes[ii][0] > 0.19999) {
        //    mesh.velocities[ii] *= 0;
        //    mesh.velocities[ii][0] = 1;
        //}

        //if (mesh.vertexes[ii][0] < -0.19999) {
        //    mesh.velocities[ii] *= 0;
        //    mesh.velocities[ii][0] = -1;
        //}




        //if (mesh.vertexes[ii][1] > 0.19999) {
        //    mesh.vertexes[ii] = rotationTop * mesh.vertexes[ii];
        //    mesh.velocities[ii] *= 0;
        //    //mesh.velocities[ii][1] = 1;
        //}

        //if (mesh.vertexes[ii][1] < -0.19999) {
        //    mesh.vertexes[ii] = rotationDown * mesh.vertexes[ii];
        //    mesh.velocities[ii] *= 0;
        //    //mesh.velocities[ii][1] = -1;
        //}



        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * explicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * explicit_time_step;
        mesh.vertexes[ii][2] += mesh.velocities[ii][2] * explicit_time_step;

        boundary_process(mesh, ii, false);
    }
}



void fem_implicit3D(mesh3D& mesh) {

    Vector3d direction = Vector3d(0, 1, 0);
    float tolerance = 1e-9;
    VectorXd P(mesh.vertexNum * 3);
    vector<Vector3d> b0(mesh.vertexNum, Vector3d(0, 0, 0));
    P.setZero();
    double fsum = 0;


    //Matrix3d rotationTop, rotationDown;
    //double angleYVelocity = (PI) * implicit_time_step;
    //angle += angleYVelocity;
    //if (angle < PI / 2) {
    //    rotationTop << cos(angleYVelocity), 0, -sin(angleYVelocity), 0, 1, 0, sin(angleYVelocity), 0, cos(angleYVelocity);
    //    rotationDown << cos(-angleYVelocity), 0, -sin(-angleYVelocity), 0, 1, 0, sin(-angleYVelocity), 0, cos(-angleYVelocity);
    //}
    //else {
    //    rotationTop << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    //    rotationDown = rotationTop;
    //}

    //for (int ii = 0; ii < mesh.vertexNum; ii++) {
    //    if (mesh.vertexes[ii][1] > 0.19999) {
    //        mesh.velocities[ii] *= 0;
    //        //mesh.velocities[ii][1] = 1;
    //        mesh.Constraints[ii].setZero();
    //        //mesh.d_velocities[ii][1] = 1;
    //        //continue;
    //        mesh.vertexes[ii] = rotationTop * mesh.vertexes[ii];
    //    }

    //    if (mesh.vertexes[ii][1] < -0.19999) {
    //        mesh.velocities[ii] *= 0;
    //        //mesh.velocities[ii][1] = -1;
    //        mesh.Constraints[ii].setZero();
    //        mesh.vertexes[ii] = rotationDown * mesh.vertexes[ii];
    //        //mesh.d_velocities[ii][1] = 1;
    //        //continue;
    //    }
    //}



    for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic3D_double(F, direction, aniosScale);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);
        //cout << "pepf:\n" << pepf << endl;
        //cout << "PFPX:\n" << PFPX << endl;
        MatrixXd f = -PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < tolerance) {
                f(i, 0) = 0;
            }
            else {
                f(i, 0) = mesh.volum[ii] * f(i, 0);
            }
            fsum += abs(f(i, 0));
        }

        for (int i = 0; i < 12; i++) {
            b0[mesh.tetrahedras[ii][i / 3]][i % 3] += f(i, 0) * implicit_time_step;
        }

        //Init Precondition Matrix
        MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumRate);

        //cout << "Hq:" << endl;
        //cout << mesh.volum[ii] * PFPX.transpose() * Hq * PFPX <<endl;

        Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
        MatrixXd H = -mesh.volum[ii] * PFPX.transpose() * Hq * PFPX;


        VectorXd tempC(12);
        for (int i = 0; i < 12; i++) {
            tempC(i) = mesh.velocities[mesh.tetrahedras[ii][i / 3]][i % 3];
        }

        VectorXd tempQ = implicit_time_step * implicit_time_step * H * tempC;

        for (int i = 0; i < 12; i++) {
            b0[mesh.tetrahedras[ii][i / 3]][i % 3] += tempQ(i);
        }


        for (int i = 0; i < 12; i++) {
            tempC(i) = mesh.d_positions[mesh.tetrahedras[ii][i / 3]][i % 3];
        }

        tempQ = implicit_time_step * H * tempC;

        for (int i = 0; i < 12; i++) {
            b0[mesh.tetrahedras[ii][i / 3]][i % 3] += tempQ(i);
        }

        for (int i = 0; i < 12; i++) {
            P[mesh.tetrahedras[ii][i / 3] * 3 + i % 3] += mesh.volum[ii] * density / 4 - implicit_time_step * implicit_time_step * H(i, i);
        }
    }

    //vector<Vector3d> dV(mesh.vertexNum, Vector3d(0, 0, 0));
    std::cout << "fsum: " << fsum << endl;

    vector<Vector3d> tempDeltaVelocites(mesh.vertexNum, Vector3d(0, 0, 0));
    /*for (int i = 0; i < mesh.vertexNum; i++) {
        tempDeltaVelocites.push_back(Vector3d(mesh.d_velocities[i]));
    }*/
    double deltaN = 0;
    double localOptimal = DBL_MAX;
    bool getLocalOpt = false;

    if (fsum>tolerance) {
        
        //cout << "P\n" << P << endl;
        double delta0 = 0;
        double deltaO = 0;
        
        vector<Vector3d> r(mesh.vertexNum, Vector3d(0, 0, 0));
        vector<Vector3d> c(mesh.vertexNum, Vector3d(0, 0, 0));

        for (int i = 0; i < mesh.vertexNum; i++) {
            //P(i * 3) = 1;
            //compute delta0
            //P(i * 3) = 1; P(i * 3 + 1) = 1; P(i * 3 + 2) = 1;
            double vx = 1 / (P(i * 3));// abs(P(i * 3)) > 1e-5 ? 1 / (P(i * 3)) : 1;
            double vy = 1 / (P(i * 3 + 1));// abs(P(i * 3) + 1) > 1e-5 ? 1 / (P(i * 3) + 1) : 1;
            double vz = 1 / (P(i * 3 + 2));// abs(P(i * 3) + 2) > 1e-5 ? 1 / (P(i * 3) + 2) : 1;
            Vector3d filter_b = mesh.Constraints[i] * b0[i];
            delta0 += filter_b[0] * filter_b[0] * vx;
            delta0 += filter_b[1] * filter_b[1] * vy;
            delta0 += filter_b[2] * filter_b[2] * vz;
        }

        for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
            MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
            MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
            MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumRate);
            Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
            MatrixXd H = -mesh.volum[ii] * PFPX.transpose() * Hq * PFPX;
            MatrixXd M(12, 12);// = mass * MatrixXd
            M.setIdentity();
            M = mesh.volum[ii] * density / 4 * M;

            MatrixXd A = M - implicit_time_step * implicit_time_step * H;

            VectorXd tempC(12);
            for (int i = 0; i < 12; i++) {
                tempC(i) = mesh.d_velocities[mesh.tetrahedras[ii][i / 3]][i % 3];
            }

            VectorXd tempQ = A * tempC;

            for (int i = 0; i < 12; i++) {
                r[mesh.tetrahedras[ii][i / 3]][i % 3] += -tempQ(i);
            }
        }

        for (int i = 0; i < mesh.vertexNum; i++) {
            r[i] += b0[i];
            r[i] = mesh.Constraints[i] * r[i];
            c[i][0] = r[i][0] * P(i * 3);
            c[i][1] = r[i][1] * P(i * 3 + 1);
            c[i][2] = r[i][2] * P(i * 3 + 2);
            c[i] = mesh.Constraints[i] * c[i];

            deltaN += c[i][0] * r[i][0];
            deltaN += c[i][1] * r[i][1];
            deltaN += c[i][2] * r[i][2];
        }
        //localOptimal = deltaN;
        
        double errorRate = 1e-6;
        std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN << endl;
        //PCG main loop
        int cgCounts = 0;
        while (deltaN > errorRate* errorRate* delta0) {
            cgCounts++;
            std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN <<"      iteration_counts:      "<<cgCounts<<endl;
            vector<Vector3d> q(mesh.vertexNum, Vector3d(0, 0, 0));
            for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
                MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
                MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
                MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumRate);
                Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
                MatrixXd H = -mesh.volum[ii] * PFPX.transpose() * Hq * PFPX;
                MatrixXd M(12, 12);// = mass * MatrixXd
                M.setIdentity();
                M = mesh.volum[ii] * density / 4 * M;

                MatrixXd A = M - implicit_time_step * implicit_time_step * H;
                VectorXd tempC(12);

                for (int i = 0; i < 12; i++) {
                    tempC(i) = c[mesh.tetrahedras[ii][i / 3]][i % 3];
                }
                VectorXd tempQ = A * tempC;

                for (int i = 0; i < 12; i++) {
                    q[mesh.tetrahedras[ii][i / 3]][i % 3] += tempQ(i);
                }
            }

            double tempSum = 0;
            for (int i = 0; i < mesh.vertexNum; i++) {
                q[i] = mesh.Constraints[i] * q[i];
                tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1] + c[i][2] * q[i][2]);
            }
            double alpha = deltaN / tempSum;
            //cout << "tempSum:------------------"<<tempSum << endl;
            //if(tempSum)
            deltaO = deltaN;
            deltaN = 0;
            vector<Vector3d> s(mesh.vertexNum, Vector3d(0, 0, 0));
            for (int i = 0; i < mesh.vertexNum; i++) {
                mesh.d_velocities[i][0] = mesh.d_velocities[i][0] + alpha * c[i][0];
                mesh.d_velocities[i][1] = mesh.d_velocities[i][1] + alpha * c[i][1];
                mesh.d_velocities[i][2] = mesh.d_velocities[i][2] + alpha * c[i][2];

                r[i][0] = r[i][0] - alpha * q[i][0];
                r[i][1] = r[i][1] - alpha * q[i][1];
                r[i][2] = r[i][2] - alpha * q[i][2];

                s[i][0] = r[i][0] * P(i * 3);
                s[i][1] = r[i][1] * P(i * 3 + 1);
                s[i][2] = r[i][2] * P(i * 3 + 2);

                deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2]);
            }

            if (deltaN < localOptimal) {
                localOptimal = deltaN;
                getLocalOpt = true;
                for (int j = 0; j < mesh.vertexNum; j++) {
                    tempDeltaVelocites[j] = Vector3d(mesh.d_velocities[j]);
                }
            }

            for (int i = 0; i < mesh.vertexNum; i++) {
                c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
                c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
                c[i][2] = s[i][2] + (deltaN / deltaO) * c[i][2];
                c[i] = mesh.Constraints[i] * c[i];
            }
        }
    }
    //cout << "deltaN:  " << deltaN << endl;
    double gravity = -9.8;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {
        //clear float numerical error
        //if (getLocalOpt) {
        //    mesh.d_velocities[ii] = tempDeltaVelocites[ii];
        //}

        //if (abs(mesh.d_velocities[ii][0]) < 1e-10) {
        //    mesh.d_velocities[ii][0] = 0;
        //}
        //if (abs(mesh.d_velocities[ii][1]) < 1e-10) {
        //    mesh.d_velocities[ii][1] = 0;
        //}
        //if (abs(mesh.d_velocities[ii][2]) < 1e-10) {
        //    mesh.d_velocities[ii][2] = 0;
        //}

        mesh.velocities[ii][0] += mesh.d_velocities[ii][0];
        mesh.velocities[ii][1] += mesh.d_velocities[ii][1] + gravity * implicit_time_step;
        mesh.velocities[ii][2] += mesh.d_velocities[ii][2];

        //if (mesh.vertexes[ii][1] > 0.19999) {
        //    mesh.velocities[ii] *= 0;
        //    mesh.velocities[ii][1] = 1;
        //}

        //if (mesh.vertexes[ii][1] < -0.19999) {
        //    mesh.velocities[ii] *= 0;
        //    mesh.velocities[ii][1] = -1;
        //}

        //if (mesh.vertexes[ii][1] > 0.19999) {
        //    mesh.vertexes[ii] = rotationTop * mesh.vertexes[ii];
        //    mesh.velocities[ii] *= 0;
        //    //mesh.velocities[ii][1] = 1;
        //}

        //if (mesh.vertexes[ii][1] < -0.19999) {
        //    mesh.vertexes[ii] = rotationDown * mesh.vertexes[ii];
        //    mesh.velocities[ii] *= 0;
        //    //mesh.velocities[ii][1] = -1;
        //}
        bool resetContraint = true;
        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * implicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * implicit_time_step;
        mesh.vertexes[ii][2] += mesh.velocities[ii][2] * implicit_time_step;

        //if (abs(mesh.vertexes[ii][0]) > 0.5) {
        //    resetContraint = true;
        //}

        boundary_process(mesh, ii, resetContraint);
    }
}

void Projected_Newton3D(mesh3D& mesh, double &mfsum, int& total_cg_iterations, int& total_newton_iterations) {
    Vector3d direction = Vector3d(0, 1, 0);
    float tolerance = 1e-8;
    VectorXd P(mesh.vertexNum * 3);
    vector<Vector3d> b0(mesh.vertexNum, Vector3d(0, 0, 0));
    P.setZero();
    double fsum = 0;
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic3D_double(F, direction, aniosScale);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);
        //cout << "pepf:\n" << pepf << endl;
        //cout << "PFPX:\n" << PFPX << endl;
        MatrixXd f = PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-8) {
                f(i, 0) = 0;
            }
            else {
                f(i, 0) = mesh.volum[ii] * f(i, 0);
            }
            fsum += abs(f(i, 0));
        }

        for (int i = 0; i < 12; i++) {
            b0[mesh.tetrahedras[ii][i / 3]][i % 3] += f(i, 0);
        }

        MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumRate);
        Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
        MatrixXd H = mesh.volum[ii] * PFPX.transpose() * Hq * PFPX;

        for (int i = 0; i < 12; i++) {
            P[mesh.tetrahedras[ii][i / 3] * 3 + i % 3] += H(i, i);
        }

    }
    mfsum = fsum;
    //std::cout << "fsum: " << fsum << endl;

    vector<Vector3d> tempDeltaX(mesh.vertexNum, Vector3d(0, 0, 0));

    double deltaN = 0;
    double localOptimal = DBL_MAX;
    bool getLocalOpt = false;
    vector<Vector3d> dX(mesh.vertexNum, Vector3d(0, 0, 0));
    if (fsum > tolerance) {
        double delta0 = 0;
        double deltaO = 0;

        vector<Vector3d> r(mesh.vertexNum, Vector3d(0, 0, 0));
        vector<Vector3d> c(mesh.vertexNum, Vector3d(0, 0, 0));
        
        for (int i = 0; i < mesh.vertexNum; i++) {
            double vx = 1 / (P(i * 3));
            double vy = 1 / (P(i * 3 + 1));
            double vz = 1 / (P(i * 3 + 2));
            Vector3d filter_b = mesh.Constraints[i] * b0[i];
            delta0 += filter_b[0] * filter_b[0] * vx;
            delta0 += filter_b[1] * filter_b[1] * vy;
            delta0 += filter_b[2] * filter_b[2] * vz;
        }

        for (int i = 0; i < mesh.vertexNum; i++) {
            r[i] = b0[i];
            r[i] = mesh.Constraints[i] * r[i];
            c[i][0] = r[i][0] * P(i * 3);
            c[i][1] = r[i][1] * P(i * 3 + 1);
            c[i][2] = r[i][2] * P(i * 3 + 2);
            c[i] = mesh.Constraints[i] * c[i];

            deltaN += c[i][0] * r[i][0];
            deltaN += c[i][1] * r[i][1];
            deltaN += c[i][2] * r[i][2];
        }
        //localOptimal = deltaN;

        double errorRate = 1e-1;
        //std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN << endl;
        //PCG main loop
        int cgCounts = 0;
        while (cgCounts<100000 && deltaN > errorRate* errorRate* delta0) {
            cgCounts++;
            //std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN << "      iteration_counts:      " << cgCounts << endl;
            vector<Vector3d> q(mesh.vertexNum, Vector3d(0, 0, 0));
            for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
                MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
                MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
                MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumRate);
                Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
                MatrixXd H = mesh.volum[ii] * PFPX.transpose() * Hq * PFPX;
                //MatrixXd M(12, 12);// = mass * MatrixXd
                //M.setIdentity();
                //M = mesh.volum[ii] * density / 4 * M;

                //MatrixXd A = M - implicit_time_step * implicit_time_step * H;
                VectorXd tempC(12);

                for (int i = 0; i < 12; i++) {
                    tempC(i) = c[mesh.tetrahedras[ii][i / 3]][i % 3];
                }
                VectorXd tempQ = H * tempC;

                for (int i = 0; i < 12; i++) {
                    q[mesh.tetrahedras[ii][i / 3]][i % 3] += tempQ(i);
                }
            }

            double tempSum = 0;
            for (int i = 0; i < mesh.vertexNum; i++) {
                q[i] = mesh.Constraints[i] * q[i];
                tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1] + c[i][2] * q[i][2]);
            }
            double alpha = deltaN / tempSum;
            //cout << "tempSum:------------------"<<tempSum << endl;
            //if(tempSum)
            deltaO = deltaN;
            deltaN = 0;
            vector<Vector3d> s(mesh.vertexNum, Vector3d(0, 0, 0));
            for (int i = 0; i < mesh.vertexNum; i++) {
                dX[i][0] = dX[i][0] + alpha * c[i][0];
                dX[i][1] = dX[i][1] + alpha * c[i][1];
                dX[i][2] = dX[i][2] + alpha * c[i][2];

                r[i][0] = r[i][0] - alpha * q[i][0];
                r[i][1] = r[i][1] - alpha * q[i][1];
                r[i][2] = r[i][2] - alpha * q[i][2];

                s[i][0] = r[i][0] * P(i * 3);
                s[i][1] = r[i][1] * P(i * 3 + 1);
                s[i][2] = r[i][2] * P(i * 3 + 2);

                deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2]);
            }

            if (deltaN < localOptimal) {
                localOptimal = deltaN;
                getLocalOpt = true;
                for (int j = 0; j < mesh.vertexNum; j++) {
                    tempDeltaX[j] = Vector3d(dX[j]);
                }
            }

            for (int i = 0; i < mesh.vertexNum; i++) {
                c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
                c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
                c[i][2] = s[i][2] + (deltaN / deltaO) * c[i][2];
                c[i] = mesh.Constraints[i] * c[i];
            }
        }
        cout << "cg_counts:   " << cgCounts << endl;
        total_cg_iterations += cgCounts;
    }
    total_newton_iterations++;
    for (int i = 0; i < mesh.vertexNum; i++) {
        if (getLocalOpt) {
            dX[i] = tempDeltaX[i];
        }
        double dt = 0.1;
        mesh.vertexes[i][0] -= dX[i][0] * dt;
        mesh.vertexes[i][1] -= dX[i][1] * dt;
        mesh.vertexes[i][2] -= dX[i][2] * dt;
    }
}
