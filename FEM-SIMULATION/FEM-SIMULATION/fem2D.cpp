#include "fem2D.h"
#include "fem_math.h"
#include<iostream>
#include <vector>
#include "fem_parameters.h"
using namespace std;
using namespace FEM;

Matrix4d project_ARAP_H_2D(const Matrix2d& F) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d Ue, sigma, Ve;
    Ue = svdResult.U;
    sigma = svdResult.SIGMA;
    Ve = svdResult.V;

    Matrix2d twist;

    twist << 0, -1, 1, 0;
    twist *= 1.f / sqrt(2.f);
    const MatrixXd e = vec_double(Ue * twist * Ve.transpose());
    const float I_1 = sigma.trace();
    const float filtered = (I_1 >= 2.0) ? 2.0 / I_1 : 1.0;

    Matrix4d H;
    H.setIdentity();
    H -= filtered * (e * e.transpose());
    H *= 2.0;
    return H;
}

Matrix4d project_ANIOSI5_H_2D(const Matrix2d& F, Vector2d direction, const double& scale) {
    direction.normalize();
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, sigma, V, S;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();
    double I4 = direction.transpose() * S * direction;
    double I5 = direction.transpose() * S.transpose() * S * direction;

    if (abs(I5) < 1e-10) return MatrixXd::Zero(4, 4);

    double s = 0;
    if (I4 < 0) {
        s = -1;
    }
    else if (I4 > 0) {
        s = 1;
    }

    double lamda0 = scale;
    double lamda1 = scale * (1 - s / sqrt(I5));
    //double lamda2 = lamda1;
    Matrix2d Q0, Q1, A;
    A = direction * direction.transpose();
    Q0 = (1 / sqrt(I5)) * F * A;
    
    Matrix2d twist;
    twist << 0, -1, 1, 0;
    twist *= 1.f / sqrt(2.f);
    Q1 = U * twist * sigma * V.transpose() * A;

    Matrix4d H = lamda0 * vec_double(Q0) * vec_double(Q0).transpose();
    if (lamda1 >= 0) {
        H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();
    }

    return H;
}

Matrix4d project_StabbleNHK_H_2D(const Matrix2d& F, const double& lengthRate, const double& volumRate) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, sigma, V, A;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    double u = lengthRate, r = volumRate;
    A.setZero();
    A(0, 0) = u;
    A(1, 1) = u;

    double lamda0 = A.eigenvalues()(0).real();
    double lamda1 = A.eigenvalues()(1).real();

    Matrix2d D0, D1;
    D0 << 1, 0, 0, 0;
    D1 << 0, 0, 0, 1;

    D0 = U * D0 * V.transpose();
    D1 = U * D1 * V.transpose();

    double z0 = sigma(1, 1) * lamda0;
    double z1 = sigma(0, 0) * lamda0;

    Matrix2d Q0 = z0 * D0 + z1 * D1;

    z0 = sigma(1, 1) * lamda1;
    z1 = sigma(0, 0) * lamda1;

    Matrix2d Q1 = z0 * D0 + z1 * D1;

    Matrix4d H;
    H.setZero();
    H += lamda0 * vec_double(Q0) * vec_double(Q0).transpose();
    H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();

    Matrix2d twist0, twist1;
    twist0 << 0, -1, 1, 0;
    twist1 << 0, 1, 1, 0;
    twist0 *= 1.f / sqrt(2.f);
    twist1 *= 1.f / sqrt(2.f);
    twist0 = U * twist0 * V.transpose();
    twist1 = U * twist1 * V.transpose();

    double lamda2 = u;
    double lamda3 = u;

    H += lamda2 * vec_double(twist0) * vec_double(twist0).transpose();
    H += lamda3 * vec_double(twist1) * vec_double(twist1).transpose();
    //cout << H << endl;
    return H;
}

Matrix2d calculateDms2D_double(const vector<Vector2d>& vertexes, const vector<uint64_t>& index, const int& i) {
    int id1 = (i + 1) % 3;
    int id2 = (i + 2) % 3;
    double o1x = vertexes[index[id1]][0] - vertexes[index[i]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[i]][1];

    double o2x = vertexes[index[id2]][0] - vertexes[index[i]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[i]][1];

    Matrix2d M; 
    M(0, 0) = o1x; M(0, 1) = o2x;
    M(1, 0) = o1y; M(1, 1) = o2y;
    
    return M;
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

MatrixXd computePFPX2D_double(const Matrix2d& InverseDm) {
    MatrixXd PFPX = MatrixXd::Zero(4, 6);
    float r00 = InverseDm(0, 0), r01 = InverseDm(0, 1);
    float r10 = InverseDm(1, 0), r11 = InverseDm(1, 1);
    float s0 = r00 + r10;
    float s1 = r01 + r11;
    PFPX(0, 0) = -s0;  PFPX(2, 0) = -s1;  PFPX(1, 1) = -s0;  PFPX(3, 1) = -s1;
    PFPX(0, 2) = r00;  PFPX(2, 2) = r01;  PFPX(1, 3) = r00;  PFPX(3, 3) = r01;
    PFPX(0, 4) = r10;  PFPX(2, 4) = r11;  PFPX(1, 5) = r10;  PFPX(3, 5) = r11;
    return PFPX;
}

Matrix2d computePEPF_ARAP2D_double(const Matrix2d& F) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R;
    U = svdResult.U;
    V = svdResult.V;
    R = U * V.transpose();
    Matrix2d PEPF = F - R;
    return PEPF;
}

Matrix2d computePEPF_StableNHK_double(const Matrix2d& F, const double& lengthRate, const double& volumRate) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();
    
    double u = lengthRate, r = volumRate;
    Matrix2d pI3pF;
    pI3pF(0, 0) = F(1, 1); pI3pF(0, 1) = -F(1, 0);
    pI3pF(1, 0) = -F(0, 1); pI3pF(1, 1) = F(0, 0);
    Matrix2d PEPF = u * F + (r * (S.determinant() - 1) - u) * pI3pF;
    return PEPF;
}

Matrix2d computePEPF_Aniostropic_double(const Matrix2d& F, Vector2d direction, const double& scale) {
    //double x = 0, y = 1;
    //double rate = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    //direction /= rate;
    direction.normalize();
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();

    double I4 = direction.transpose() * S * direction;
    double I5 = direction.transpose() * S.transpose() * S * direction;

    if(I4==0){
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
    Matrix2d PEPF = scale * (1 - s / sqrt(I5)) * F * direction * direction.transpose();
    return PEPF;
}

void boundary_process(mesh2D& mesh, const int& index) {
    if (mesh.vertexes[index][1] < -1) {
        mesh.vertexes[index][1] = -1; mesh.velocities[index][1] *= -0.1;
    }

    if (mesh.vertexes[index][0] < -1) {
        mesh.vertexes[index][0] = -1; mesh.velocities[index][0] *= -0.1;
    }

    if (mesh.vertexes[index][1] > 1) {
        mesh.vertexes[index][1] = 1; mesh.velocities[index][1] *= -0.1;
    }

    if (mesh.vertexes[index][0] > 1) {
        mesh.vertexes[index][0] = 1; mesh.velocities[index][0] *= -0.1;
    }
}

///////////////////////FEM SIMULATIONS//////////////////////////////////

void fem_explicit2D(mesh2D& mesh) {
    Vector2d direction = Vector2d(1, 0);
    for (int ii = 0; ii < mesh.triangleNum; ii++) {
        MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix2d PEPF; PEPF.setZero();
        //PEPF += computePEPX_ARAP2D_double(F);

        PEPF += computePEPF_StableNHK_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic_double(F, direction, aniosScale);
        //PEPF += computePEPX_Aniostropic_double(F, Vector2d(1, 0), 0.3);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);

        /*cout << "pepf:\n" << pepf << endl;
        cout << "PFPX:\n" << PFPX << endl;*/
        MatrixXd f = -PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-10) f(i, 0) = 0;
        }

        mesh.forces[mesh.triangles[ii][0]][0] += f(0, 0);
        mesh.forces[mesh.triangles[ii][0]][1] += f(1, 0);
        mesh.forces[mesh.triangles[ii][1]][0] += f(2, 0);
        mesh.forces[mesh.triangles[ii][1]][1] += f(3, 0);
        mesh.forces[mesh.triangles[ii][2]][0] += f(4, 0);
        mesh.forces[mesh.triangles[ii][2]][1] += f(5, 0);
    }
    //float sum = 0.f;
    //n++;
    //if (sqrt(sum) < 1e-4) return;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {

        double acx = mesh.forces[ii][0] * mass_inverse;
        double acy = mesh.forces[ii][1] * mass_inverse - 9.8;

        mesh.forces[ii][0] = 0;
        mesh.forces[ii][1] = 0;

        mesh.velocities[ii][0] += acx * explicit_time_step;
        mesh.velocities[ii][1] += acy * explicit_time_step;

            //mesh.velocities[2][1] = 0;
            //mesh.velocities[3][1] = 0;
            //mesh.velocities[4][1] = 0;
            //mesh.velocities[2][0] = 0;
            //mesh.velocities[3][0] = 0;
            //mesh.velocities[4][0] = 0;


            //mesh.velocities[12][1] = 0;
            //mesh.velocities[13][1] = 0;
            //mesh.velocities[14][1] = 0;
            //mesh.velocities[12][0] = 0;
            //mesh.velocities[13][0] = 0;
            //mesh.velocities[14][0] = 0;

        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * explicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * explicit_time_step;

        boundary_process(mesh, ii);
    }
}

void initMesh2D(mesh2D &mesh, int type, double scale) {

    mesh.InitMesh(type, scale);
    Matrix2d rotation;
    rotation << 0, -1, 1, 0;
    if (type == 1) {
        for (int i = 0; i < mesh.vertexNum; i++) {
            mesh.vertexes[i] = rotation * mesh.vertexes[i];
        }
    }

    for (int i = 0; i < mesh.triangleNum; i++) {
        mesh.DM_triangle_inverse.push_back(calculateDms2D_double(mesh.vertexes, mesh.triangles[i], 0).inverse());
    }
    
    //if (true) {
    //    mesh.vertexes[3][0] -= 0.3;
    //    mesh.vertexes[3][1] -= 0.3;
    //}
}
 
void Projected_Newton2D(mesh2D& mesh) {
    
    VectorXd P(mesh.vertexNum * 2);
    P.setZero();
    for (int ii = 0; ii < mesh.triangleNum; ii++) {
        MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix2d PEPF = computePEPF_StableNHK_double(F, lengthRate, volumRate);

        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);
        /*cout << "pepf:\n" << pepf << endl;
        cout << "PFPX:\n" << PFPX << endl;*/
        MatrixXd f = PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        //for (int i = 0; i < f.rows(); i++) {
        //    if (abs(f(i, 0)) < 1e-10) f(i, 0) = 0;
        //}

        //compute the b for formula: Ax = b
        mesh.forces[mesh.triangles[ii][0]][0] += f(0, 0);
        mesh.forces[mesh.triangles[ii][0]][1] += f(1, 0);
        mesh.forces[mesh.triangles[ii][1]][0] += f(2, 0);
        mesh.forces[mesh.triangles[ii][1]][1] += f(3, 0);
        mesh.forces[mesh.triangles[ii][2]][0] += f(4, 0);
        mesh.forces[mesh.triangles[ii][2]][1] += f(5, 0);

        //Init Precondition Matrix
        Matrix4d Hq = project_StabbleNHK_H_2D(F, lengthRate, volumRate);
        MatrixXd H = PFPX.transpose() * Hq * PFPX;

        //cout << "H:\n" << H << endl;

        //double threshold = 1e-10;
        P(mesh.triangles[ii][0] * 2) += H(0, 0);
        P(mesh.triangles[ii][0] * 2 + 1) += H(1, 1);
        P(mesh.triangles[ii][1] * 2) += H(2, 2);
        P(mesh.triangles[ii][1] * 2 + 1) += H(3, 3);
        P(mesh.triangles[ii][2] * 2) += H(4, 4);
        P(mesh.triangles[ii][2] * 2 + 1) += H(5, 5);
    }
    //cout << "P\n" << P << endl;
    double delta0 = 0;
    double deltaO = 0;
    double deltaN = 0;
    vector<Vector2d> r(mesh.vertexNum, Vector2d(0, 0));
    vector<Vector2d> c(mesh.vertexNum, Vector2d(0, 0));
    vector<Vector2d> dX(mesh.vertexNum, Vector2d(0, 0));
    //double r[9][2];
    //double c[9][2];
    //double dX[9][2] = { 0 };
    for (int i = 0; i < mesh.vertexNum; i++) {
        //compute delta0
        double vx = abs(P(i * 2)) > 1e-5 ? 1 / P(i * 2) : 1;
        double vy = abs(P(i * 2) + 1) > 1e-5 ? 1 / P(i * 2 + 1) : 1;
        delta0 += mesh.forces[i][0] * mesh.forces[i][0] * vx;
        delta0 += mesh.forces[i][1] * mesh.forces[i][1] * vy;
        
        //init r
        r[i][0] = mesh.forces[i][0];
        r[i][1] = mesh.forces[i][1];

        mesh.forces[i][0] = 0;
        mesh.forces[i][1] = 0;
        //init c
        c[i][0] = r[i][0] * P(i * 2);
        c[i][1] = r[i][1] * P(i * 2 + 1);

        deltaN += c[i][0] * r[i][0];
        deltaN += c[i][1] * r[i][1];
    }

    double errorRate = 1e-15;

    //PCG main loop
    while (deltaN > errorRate * errorRate * delta0) {
        vector<Vector2d> q(mesh.vertexNum, Vector2d(0, 0));
        for (int ii = 0; ii < mesh.triangleNum; ii++) {
            MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
            MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
            Matrix4d Hq = project_StabbleNHK_H_2D(F, lengthRate, volumRate);
            MatrixXd H = PFPX.transpose() * Hq * PFPX;
            VectorXd tempC(6);
            tempC(0) = c[mesh.triangles[ii][0]][0];
            tempC(1) = c[mesh.triangles[ii][0]][1];
            tempC(2) = c[mesh.triangles[ii][1]][0];
            tempC(3) = c[mesh.triangles[ii][1]][1];
            tempC(4) = c[mesh.triangles[ii][2]][0];
            tempC(5) = c[mesh.triangles[ii][2]][1];
            VectorXd tempQ = H * tempC;
            q[mesh.triangles[ii][0]][0] += tempQ(0);
            q[mesh.triangles[ii][0]][1] += tempQ(1);
            q[mesh.triangles[ii][1]][0] += tempQ(2);
            q[mesh.triangles[ii][1]][1] += tempQ(3);
            q[mesh.triangles[ii][2]][0] += tempQ(4);
            q[mesh.triangles[ii][2]][1] += tempQ(5);
        }

        double tempSum = 0;
        for (int i = 0; i < mesh.vertexNum; i++) {
            tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1]);
        }
        double alpha = deltaN / tempSum;

        deltaO = deltaN;
        deltaN = 0;
        vector<Vector2d> s(mesh.vertexNum, Vector2d(0, 0));
        for (int i = 0; i < mesh.vertexNum; i++) {
            dX[i][0] = dX[i][0] + alpha * c[i][0];
            dX[i][1] = dX[i][1] + alpha * c[i][1];

            r[i][0] = r[i][0] - alpha * q[i][0];
            r[i][1] = r[i][1] - alpha * q[i][1];

            s[i][0] = r[i][0] * P(i * 2);
            s[i][1] = r[i][1] * P(i * 2 + 1);

            deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1]);
        }

        
        for (int i = 0; i < mesh.vertexNum; i++) {
            c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
            c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
        }
    }

    for (int i = 0; i < mesh.vertexNum; i++) {
        mesh.vertexes[i][0] -= dX[i][0] * implicit_time_step;
        mesh.vertexes[i][1] -= dX[i][1] * implicit_time_step;
    }
}


void fem_implicit2D(mesh2D& mesh) {

    Vector2d direction = Vector2d(1, 2);
    
    VectorXd P(mesh.vertexNum * 2), b0(mesh.vertexNum * 2);
    P.setZero(); b0.setZero();
    for (int ii = 0; ii < mesh.triangleNum; ii++) {
        MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix2d PEPF = computePEPF_StableNHK_double(F, lengthRate, volumRate);
        PEPF += computePEPF_Aniostropic_double(F, direction, aniosScale);
        //cout << "PEPF:\n" << PEPF << endl;
        MatrixXd pepf = vec_double(PEPF);
        /*cout << "pepf:\n" << pepf << endl;
        cout << "PFPX:\n" << PFPX << endl;*/
        MatrixXd f = -PFPX.transpose() * pepf;
        //cout << "f:\n" << f << endl;
        for (int i = 0; i < f.rows(); i++) {
            if (abs(f(i, 0)) < 1e-16) f(i, 0) = 0;
        }

        //compute the b for formula: Ax = b
        b0[mesh.triangles[ii][0] * 2] += f(0, 0) * implicit_time_step;
        b0[mesh.triangles[ii][0] * 2 + 1] += f(1, 0) * implicit_time_step;
        b0[mesh.triangles[ii][1] * 2] += f(2, 0) * implicit_time_step;
        b0[mesh.triangles[ii][1] * 2 + 1] += f(3, 0) * implicit_time_step;
        b0[mesh.triangles[ii][2] * 2] += f(4, 0) * implicit_time_step;
        b0[mesh.triangles[ii][2] * 2 + 1] += f(5, 0) * implicit_time_step;

        //Init Precondition Matrix
        Matrix4d Hq = project_StabbleNHK_H_2D(F, lengthRate, volumRate);
        Hq += project_ANIOSI5_H_2D(F, direction, aniosScale);
        MatrixXd H = -PFPX.transpose() * Hq * PFPX;


        VectorXd tempC(6);
        tempC(0) = mesh.velocities[mesh.triangles[ii][0]][0];
        tempC(1) = mesh.velocities[mesh.triangles[ii][0]][1];
        tempC(2) = mesh.velocities[mesh.triangles[ii][1]][0];
        tempC(3) = mesh.velocities[mesh.triangles[ii][1]][1];
        tempC(4) = mesh.velocities[mesh.triangles[ii][2]][0];
        tempC(5) = mesh.velocities[mesh.triangles[ii][2]][1];

        VectorXd tempQ = implicit_time_step * implicit_time_step * H * tempC;
        b0[mesh.triangles[ii][0] * 2] += tempQ(0);
        b0[mesh.triangles[ii][0] * 2 + 1] += tempQ(1);
        b0[mesh.triangles[ii][1] * 2] += tempQ(2);
        b0[mesh.triangles[ii][1] * 2 + 1] += tempQ(3);
        b0[mesh.triangles[ii][2] * 2] += tempQ(4);
        b0[mesh.triangles[ii][2] * 2 + 1] += tempQ(5);
        //cout << "H:\n" << H << endl;

        //double threshold = 1e-10;
        P(mesh.triangles[ii][0] * 2) += mass - implicit_time_step * implicit_time_step * H(0, 0);
        P(mesh.triangles[ii][0] * 2 + 1) += mass - implicit_time_step * implicit_time_step * H(1, 1);
        P(mesh.triangles[ii][1] * 2) += mass - implicit_time_step * implicit_time_step * H(2, 2);
        P(mesh.triangles[ii][1] * 2 + 1) += mass - implicit_time_step * implicit_time_step * H(3, 3);
        P(mesh.triangles[ii][2] * 2) += mass - implicit_time_step * implicit_time_step * H(4, 4);
        P(mesh.triangles[ii][2] * 2 + 1) += mass - implicit_time_step * implicit_time_step * H(5, 5);
    }
    //cout << "P\n" << P << endl;
    double delta0 = 0;
    double deltaO = 0;
    double deltaN = 0;
    vector<Vector2d> r(mesh.vertexNum, Vector2d(0, 0));
    vector<Vector2d> c(mesh.vertexNum, Vector2d(0, 0));
    vector<Vector2d> dV(mesh.vertexNum, Vector2d(0, 0));
    //double r[9][2];
    //double c[9][2];
    //double dX[9][2] = { 0 };
    for (int i = 0; i < mesh.vertexNum; i++) {
        //compute delta0
        double vx = abs(P(i * 2)) > 1e-5 ? 1 / P(i * 2) : 1;
        double vy = abs(P(i * 2) + 1) > 1e-5 ? 1 / P(i * 2 + 1) : 1;
        delta0 += b0[i * 2] * b0[i * 2] * vx;
        delta0 += b0[i * 2 + 1] * b0[i * 2 + 1] * vy;

        //init r
        r[i][0] = b0[i * 2];
        r[i][1] = b0[i * 2 + 1];

        //mesh.forces[i][0] = 0;
        //mesh.forces[i][1] = 0;
        //init c
        c[i][0] = r[i][0] * P(i * 2);
        c[i][1] = r[i][1] * P(i * 2 + 1);

        deltaN += c[i][0] * r[i][0];
        deltaN += c[i][1] * r[i][1];
    }

    double errorRate = 1e-15;
    //PCG main loop
    while (deltaN > errorRate* errorRate* delta0) {
        vector<Vector2d> q(mesh.vertexNum, Vector2d(0, 0));
        for (int ii = 0; ii < mesh.triangleNum; ii++) {
            MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
            MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
            Matrix4d Hq = project_StabbleNHK_H_2D(F, lengthRate, volumRate);
            Hq += project_ANIOSI5_H_2D(F, direction, aniosScale);
            MatrixXd H = -PFPX.transpose() * Hq * PFPX;
            MatrixXd M(6, 6);// = mass * MatrixXd
            M.setIdentity();
            M = mass * M;
            MatrixXd A = M - implicit_time_step * implicit_time_step * H;
            VectorXd tempC(6);
            tempC(0) = c[mesh.triangles[ii][0]][0];
            tempC(1) = c[mesh.triangles[ii][0]][1];
            tempC(2) = c[mesh.triangles[ii][1]][0];
            tempC(3) = c[mesh.triangles[ii][1]][1];
            tempC(4) = c[mesh.triangles[ii][2]][0];
            tempC(5) = c[mesh.triangles[ii][2]][1];
            VectorXd tempQ = A * tempC;
            q[mesh.triangles[ii][0]][0] += tempQ(0);
            q[mesh.triangles[ii][0]][1] += tempQ(1);
            q[mesh.triangles[ii][1]][0] += tempQ(2);
            q[mesh.triangles[ii][1]][1] += tempQ(3);
            q[mesh.triangles[ii][2]][0] += tempQ(4);
            q[mesh.triangles[ii][2]][1] += tempQ(5);
        }
        
        double tempSum = 0;
        for (int i = 0; i < mesh.vertexNum; i++) {
            tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1]);
        }
        double alpha = deltaN / tempSum;

        deltaO = deltaN;
        deltaN = 0;
        vector<Vector2d> s(mesh.vertexNum, Vector2d(0, 0));
        for (int i = 0; i < mesh.vertexNum; i++) {
            dV[i][0] = dV[i][0] + alpha * c[i][0];
            dV[i][1] = dV[i][1] + alpha * c[i][1];

            r[i][0] = r[i][0] - alpha * q[i][0];
            r[i][1] = r[i][1] - alpha * q[i][1];

            s[i][0] = r[i][0] * P(i * 2);
            s[i][1] = r[i][1] * P(i * 2 + 1);

            deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1]);
        }


        for (int i = 0; i < mesh.vertexNum; i++) {
            c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
            c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
        }
    }
    //cout << "deltaN:  " << deltaN << endl;
    double gravity = -9.8;
    for (int ii = 0; ii < mesh.vertexNum; ii++) {
        //clear float numerical error
        if (abs(dV[ii][0]) < 1e-16) {
            dV[ii][0] = 0;
        }
        if (abs(dV[ii][1]) < 1e-16) {
            dV[ii][1] = 0;
        }

        mesh.velocities[ii][0] += dV[ii][0];
        mesh.velocities[ii][1] += dV[ii][1] + gravity * implicit_time_step;

       /* mesh.velocities[2][1] = 0;
        mesh.velocities[3][1] = 0;
        mesh.velocities[4][1] = 0;
        mesh.velocities[2][0] = 0;
        mesh.velocities[3][0] = 0;
        mesh.velocities[4][0] = 0;


        mesh.velocities[12][1] = 0;
        mesh.velocities[13][1] = 0;
        mesh.velocities[14][1] = 0;
        mesh.velocities[12][0] = 0;
        mesh.velocities[13][0] = 0;
        mesh.velocities[14][0] = 0;*/

        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * implicit_time_step;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * implicit_time_step;

        boundary_process(mesh, ii);
    }
}

void PCG_Test() {
    Matrix3d A, P;
    P.setZero();
    A << 1, 1, -2, 1, 2, -3, -2, -3, 5;
    cout << A << endl;
    Vector3d b;
    b << 2, 3, -5;
    for (int i = 0; i < 3; i++) {
        P(i, i) = 1;// / A(i, i);
    }
    float delta = b.transpose() * P * b;
    cout << delta << endl;
    double alpha = 0, beta = 0;
    Vector3d x, s, r, c, q;
    x.setZero(); s.setZero(); r.setZero(); c.setZero(); q.setZero();
    
    cout << "P\n" << P << endl;

    r = b - A * x;
    //cout << r << endl;
    cout << "r\n" << r << endl;
    c = P.inverse() * r;
    cout << "c\n" << c << endl;
    //d = z;
    float error_n = r.transpose() * c;
    cout << "error\n"<<error_n<<endl;
    float error_o = error_n;
    int iteration_times = 1000;
    while (error_n > 0.000001*delta) {
        q = A * c;
        cout << "q\n" << q << endl;
        alpha = error_n / (c.transpose() * q);
        cout << "alpha\n" << alpha << endl;
        x = x + alpha * c;
        cout << "x\n" << x << endl;
        r = r - alpha * q;
        s = P.inverse() * r;
        error_o = error_n;
        error_n = r.transpose() * s; 
        c = s + (error_n / error_o) * c;
    }
    cout << "x\n" << x << endl;
}
