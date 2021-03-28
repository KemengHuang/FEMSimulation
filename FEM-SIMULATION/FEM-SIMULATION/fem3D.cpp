#include "fem3D.h"
#include "fem_math.h"
#include<iostream>
#include <vector>
#include "fem_parameters.h"
using namespace std;
using namespace FEM;

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
    mesh.InitMesh(type, scale);

    Matrix3d rotation;
    rotation << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if (type == 1) {
        for (int i = 0; i < mesh.vertexNum; i++) {
            mesh.vertexes[i] = rotation * mesh.vertexes[i];
        }
    }

    for (int i = 0; i < mesh.tetrahedraNum; i++) {
        mesh.DM_triangle_inverse.push_back(calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[i], 0).inverse());
    }
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

Matrix3d computePEPF_Aniostropic3D_double(const Matrix3d& F, Vector3d direction, const double& scale) {
    //double x = 0, y = 1;
    //double rate = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    //direction /= rate;
    direction.normalize();
    SVDResult3D_double svdResult = SingularValueDecomposition3D_double(F);
    Matrix3d U, V, R, S, sigma;
    //U = svdResult.U;
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

    //double CC4 = direction.transpose() * F.transpose() * F * direction;
    Matrix3d PEPF = scale * (1 - s / sqrt(I5)) * F * direction * direction.transpose();
    return PEPF;
}

void boundary_process(mesh3D& mesh, const int& index) {
    if (mesh.vertexes[index][1] < -1) {
        mesh.vertexes[index][1] = -1; mesh.velocities[index][1] *= -0.1;
    }

    if (mesh.vertexes[index][0] < -1) {
        mesh.vertexes[index][0] = -1; mesh.velocities[index][0] *= -0.1;
    }

    if (mesh.vertexes[index][2] < -1) {
        mesh.vertexes[index][2] = -1; mesh.velocities[index][2] *= -0.1;
    }

    if (mesh.vertexes[index][1] > 1) {
        mesh.vertexes[index][1] = 1; mesh.velocities[index][1] *= -0.1;
    }

    if (mesh.vertexes[index][0] > 1) {
        mesh.vertexes[index][0] = 1; mesh.velocities[index][0] *= -0.1;
    }

    if (mesh.vertexes[index][2] > 1) {
        mesh.vertexes[index][2] = 1; mesh.velocities[index][2] *= -0.1;
    }
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
    if (lamda1 >= 0) {
        H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();
        H += lamda2 * vec_double(Q2) * vec_double(Q2).transpose();
    }
    
    return H;
}

void fem_explicit3D(mesh3D& mesh) {
    Vector3d direction = Vector3d(1, 1, 1);
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix3d PEPF; PEPF.setZero();
        //PEPF += computePEPX_ARAP2D_double(F);

        PEPF += computePEPF_StableNHK3D_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic3D_double(F, direction, aniosScale);
        //PEPF += computePEPX_Aniostropic_double(F, Vector2d(1, 0), 0.3);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);

        //cout << "pepf:\n" << pepf << endl;
        //cout << "PFPX:\n" << PFPX << endl;
        MatrixXd f = -PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-15) f(i, 0) = 0;
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
    //float sum = 0.f;
    //n++;
    //if (sqrt(sum) < 1e-4) return;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {

        double acx = mesh.forces[ii][0] * mass_inverse;
        double acy = mesh.forces[ii][1] * mass_inverse - 9.8;
        double acz = mesh.forces[ii][2] * mass_inverse;

        mesh.forces[ii][0] = 0;
        mesh.forces[ii][1] = 0;
        mesh.forces[ii][2] = 0;

        mesh.velocities[ii][0] += acx * explicit_time_step;
        mesh.velocities[ii][1] += acy * explicit_time_step;
        mesh.velocities[ii][2] += acz * explicit_time_step;

        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * explicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * explicit_time_step;
        mesh.vertexes[ii][2] += mesh.velocities[ii][2] * explicit_time_step;

        boundary_process(mesh, ii);
    }
}

void fem_implicit3D(mesh3D& mesh) {

    Vector3d direction = Vector3d(1, 1, 1);

    VectorXd P(mesh.vertexNum * 3), b0(mesh.vertexNum * 3);
    P.setZero(); b0.setZero();
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix3d PEPF = computePEPF_StableNHK3D_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic3D_double(F, direction, aniosScale);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);
        //cout << "pepf:\n" << pepf << endl;
        //cout << "PFPX:\n" << PFPX << endl;
        MatrixXd f = -PFPX.transpose() * pepf;
       //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-15) f(i, 0) = 0;
        }
        
        //compute the b for formula: Ax = b
        b0[mesh.tetrahedras[ii][0] * 3] += f(0, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][0] * 3 + 1] += f(1, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][0] * 3 + 2] += f(2, 0) * implicit_time_step;

        b0[mesh.tetrahedras[ii][1] * 3] += f(3, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][1] * 3 + 1] += f(4, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][1] * 3 + 2] += f(5, 0) * implicit_time_step;

        b0[mesh.tetrahedras[ii][2] * 3] += f(6, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][2] * 3 + 1] += f(7, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][2] * 3 + 2] += f(8, 0) * implicit_time_step;

        b0[mesh.tetrahedras[ii][3] * 3] += f(9, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][3] * 3 + 1] += f(10, 0) * implicit_time_step;
        b0[mesh.tetrahedras[ii][3] * 3 + 2] += f(11, 0) * implicit_time_step;

        //Init Precondition Matrix
        MatrixXd Hq = project_StabbleNHK_H_3D(F, lengthRate, volumRate);
        Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
        MatrixXd H = -PFPX.transpose() * Hq * PFPX;


        VectorXd tempC(12);
        tempC(0) = mesh.velocities[mesh.tetrahedras[ii][0]][0];
        tempC(1) = mesh.velocities[mesh.tetrahedras[ii][0]][1];
        tempC(2) = mesh.velocities[mesh.tetrahedras[ii][0]][2];
        tempC(3) = mesh.velocities[mesh.tetrahedras[ii][1]][0];
        tempC(4) = mesh.velocities[mesh.tetrahedras[ii][1]][1];
        tempC(5) = mesh.velocities[mesh.tetrahedras[ii][1]][2];
        tempC(6) = mesh.velocities[mesh.tetrahedras[ii][2]][0];
        tempC(7) = mesh.velocities[mesh.tetrahedras[ii][2]][1];
        tempC(8) = mesh.velocities[mesh.tetrahedras[ii][2]][2];
        tempC(9) = mesh.velocities[mesh.tetrahedras[ii][3]][0];
        tempC(10) = mesh.velocities[mesh.tetrahedras[ii][3]][1];
        tempC(11) = mesh.velocities[mesh.tetrahedras[ii][3]][2];

        VectorXd tempQ = implicit_time_step * implicit_time_step * H * tempC;

        b0[mesh.tetrahedras[ii][0] * 3] += tempQ(0);
        b0[mesh.tetrahedras[ii][0] * 3 + 1] += tempQ(1);
        b0[mesh.tetrahedras[ii][0] * 3 + 2] += tempQ(2);

        b0[mesh.tetrahedras[ii][1] * 3] += tempQ(3);
        b0[mesh.tetrahedras[ii][1] * 3 + 1] += tempQ(4);
        b0[mesh.tetrahedras[ii][1] * 3 + 2] += tempQ(5);

        b0[mesh.tetrahedras[ii][2] * 3] += tempQ(6);
        b0[mesh.tetrahedras[ii][2] * 3 + 1] += tempQ(7);
        b0[mesh.tetrahedras[ii][2] * 3 + 2] += tempQ(8);

        b0[mesh.tetrahedras[ii][3] * 3] += tempQ(9);
        b0[mesh.tetrahedras[ii][3] * 3 + 1] += tempQ(10);
        b0[mesh.tetrahedras[ii][3] * 3 + 2] += tempQ(11);
        //cout << "H:\n" << H << endl;

        //double threshold = 1e-10;
        P[mesh.tetrahedras[ii][0] * 3] += mass - implicit_time_step * implicit_time_step * H(0, 0);
        P[mesh.tetrahedras[ii][0] * 3 + 1] += mass - implicit_time_step * implicit_time_step * H(1, 1);
        P[mesh.tetrahedras[ii][0] * 3 + 2] += mass - implicit_time_step * implicit_time_step * H(2, 2);

        P[mesh.tetrahedras[ii][1] * 3] += mass - implicit_time_step * implicit_time_step * H(3, 3);
        P[mesh.tetrahedras[ii][1] * 3 + 1] += mass - implicit_time_step * implicit_time_step * H(4, 5);
        P[mesh.tetrahedras[ii][1] * 3 + 2] += mass - implicit_time_step * implicit_time_step * H(5, 5);

        P[mesh.tetrahedras[ii][2] * 3] += mass - implicit_time_step * implicit_time_step * H(6, 6);
        P[mesh.tetrahedras[ii][2] * 3 + 1] += mass - implicit_time_step * implicit_time_step * H(7, 7);
        P[mesh.tetrahedras[ii][2] * 3 + 2] += mass - implicit_time_step * implicit_time_step * H(8, 8);

        P[mesh.tetrahedras[ii][3] * 3] += mass - implicit_time_step * implicit_time_step * H(9, 9);
        P[mesh.tetrahedras[ii][3] * 3 + 1] += mass - implicit_time_step * implicit_time_step * H(10, 10);
        P[mesh.tetrahedras[ii][3] * 3 + 2] += mass - implicit_time_step * implicit_time_step * H(11, 11);
    }
    //cout << "P\n" << P << endl;
    double delta0 = 0;
    double deltaO = 0;
    double deltaN = 0;
    vector<Vector3d> r(mesh.vertexNum, Vector3d(0, 0, 0));
    vector<Vector3d> c(mesh.vertexNum, Vector3d(0, 0, 0));
    vector<Vector3d> dV(mesh.vertexNum, Vector3d(0, 0, 0));
    //double r[9][2];
    //double c[9][2];
    //double dX[9][2] = { 0 };
    for (int i = 0; i < mesh.vertexNum; i++) {
        //compute delta0
        double vx = abs(P(i * 3)) > 1e-5 ? 1 / P(i * 3) : 1;
        double vy = abs(P(i * 3) + 1) > 1e-5 ? 1 / P(i * 3 + 1) : 1;
        double vz = abs(P(i * 3) + 2) > 1e-5 ? 1 / P(i * 3 + 2) : 1;
        delta0 += b0[i * 3] * b0[i * 3] * vx;
        delta0 += b0[i * 3 + 1] * b0[i * 3 + 1] * vy;
        delta0 += b0[i * 3 + 2] * b0[i * 3 + 2] * vz;

        //init r
        r[i][0] = b0[i * 3];
        r[i][1] = b0[i * 3 + 1];
        r[i][2] = b0[i * 3 + 2];

        //mesh.forces[i][0] = 0;
        //mesh.forces[i][1] = 0;
        //init c
        c[i][0] = r[i][0] * P(i * 3);
        c[i][1] = r[i][1] * P(i * 3 + 1);
        c[i][2] = r[i][2] * P(i * 3 + 2);

        deltaN += c[i][0] * r[i][0];
        deltaN += c[i][1] * r[i][1];
        deltaN += c[i][2] * r[i][2];
    }

    double errorRate = 1e-15;
    //PCG main loop
    while (deltaN > errorRate* errorRate* delta0) {
        vector<Vector3d> q(mesh.vertexNum, Vector3d(0, 0, 0));
        for (int ii = 0; ii < mesh.tetrahedraNum; ii++) {
            MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
            MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];
            MatrixXd Hq = project_StabbleNHK_H_3D(F, lengthRate, volumRate);
            Hq += project_ANIOSI5_H_3D(F, direction, aniosScale);
            MatrixXd H = -PFPX.transpose() * Hq * PFPX;
            MatrixXd M(12, 12);// = mass * MatrixXd
            M.setIdentity();
            M = mass * M;
            MatrixXd A = M - implicit_time_step * implicit_time_step * H;
            VectorXd tempC(12);
            tempC(0) = c[mesh.tetrahedras[ii][0]][0];
            tempC(1) = c[mesh.tetrahedras[ii][0]][1];
            tempC(2) = c[mesh.tetrahedras[ii][0]][2];
            tempC(3) = c[mesh.tetrahedras[ii][1]][0];
            tempC(4) = c[mesh.tetrahedras[ii][1]][1];
            tempC(5) = c[mesh.tetrahedras[ii][1]][2];
            tempC(6) = c[mesh.tetrahedras[ii][2]][0];
            tempC(7) = c[mesh.tetrahedras[ii][2]][1];
            tempC(8) = c[mesh.tetrahedras[ii][2]][2];
            tempC(9) = c[mesh.tetrahedras[ii][3]][0];
            tempC(10) = c[mesh.tetrahedras[ii][3]][1];
            tempC(11) = c[mesh.tetrahedras[ii][3]][2];
            VectorXd tempQ = A * tempC;
            q[mesh.tetrahedras[ii][0]][0] += tempQ(0);
            q[mesh.tetrahedras[ii][0]][1] += tempQ(1);
            q[mesh.tetrahedras[ii][0]][2] += tempQ(2);
            q[mesh.tetrahedras[ii][1]][0] += tempQ(3);
            q[mesh.tetrahedras[ii][1]][1] += tempQ(4);
            q[mesh.tetrahedras[ii][1]][2] += tempQ(5);
            q[mesh.tetrahedras[ii][2]][0] += tempQ(6);
            q[mesh.tetrahedras[ii][2]][1] += tempQ(7);
            q[mesh.tetrahedras[ii][2]][2] += tempQ(8);
            q[mesh.tetrahedras[ii][3]][0] += tempQ(9);
            q[mesh.tetrahedras[ii][3]][1] += tempQ(10);
            q[mesh.tetrahedras[ii][3]][2] += tempQ(11);
        }

        double tempSum = 0;
        for (int i = 0; i < mesh.vertexNum; i++) {
            tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1] + c[i][2] * q[i][2]);
        }
        double alpha = deltaN / tempSum;

        deltaO = deltaN;
        deltaN = 0;
        vector<Vector3d> s(mesh.vertexNum, Vector3d(0, 0, 0));
        for (int i = 0; i < mesh.vertexNum; i++) {
            dV[i][0] = dV[i][0] + alpha * c[i][0];
            dV[i][1] = dV[i][1] + alpha * c[i][1];
            dV[i][2] = dV[i][2] + alpha * c[i][2];

            r[i][0] = r[i][0] - alpha * q[i][0];
            r[i][1] = r[i][1] - alpha * q[i][1];
            r[i][2] = r[i][2] - alpha * q[i][2];

            s[i][0] = r[i][0] * P(i * 3);
            s[i][1] = r[i][1] * P(i * 3 + 1);
            s[i][2] = r[i][2] * P(i * 3 + 2);

            deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2]);
        }


        for (int i = 0; i < mesh.vertexNum; i++) {
            c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
            c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
            c[i][2] = s[i][2] + (deltaN / deltaO) * c[i][2];
        }
    }
    //cout << "deltaN:  " << deltaN << endl;
    double gravity = -9.8;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {
        //clear float numerical error
        if (abs(dV[ii][0]) < 1e-15) {
            dV[ii][0] = 0;
        }
        if (abs(dV[ii][1]) < 1e-15) {
            dV[ii][1] = 0;
        }
        if (abs(dV[ii][2]) < 1e-15) {
            dV[ii][2] = 0;
        }

        mesh.velocities[ii][0] += dV[ii][0];
        mesh.velocities[ii][1] += dV[ii][1] + gravity * implicit_time_step;
        mesh.velocities[ii][2] += dV[ii][2];

        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * implicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * implicit_time_step;
        mesh.vertexes[ii][2] += mesh.velocities[ii][2] * implicit_time_step;
        boundary_process(mesh, ii);
    }
}
