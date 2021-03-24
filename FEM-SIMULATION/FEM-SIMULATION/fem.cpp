#include "fem.h"
#include<iostream>
#include <vector>
using namespace std;
Matrix4d project_ARAP_H_2D(Matrix2d A) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(A);
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

Matrix2d calculateDms2D_double(vector<Vector2d> vertexes, vector<uint64_t> index, int i) {
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

Matrix3d calculateDms3D_double(vector<Vector3d> vertexes, vector<uint64_t> index, int i) {
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

Matrix2d computePEPX_ARAP2D_double(const Matrix2d& F) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R;
    U = svdResult.U;
    V = svdResult.V;
    R = U * V.transpose();
    Matrix2d PEPF = F - R;
    return PEPF;
}

Matrix2d computePEPX_StableNHK_double(const Matrix2d& F) {
    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();
    
    double u = 0.3, r = 0.7;
    Matrix2d pI3pF;
    pI3pF(0, 0) = F(1, 1); pI3pF(0, 1) = -F(1, 0);
    pI3pF(1, 0) = -F(0, 1); pI3pF(1, 1) = F(0, 0);
    Matrix2d PEPF = u * F + (r * (S.determinant() - 1) - u) * pI3pF;
    return PEPF;
}

Matrix2d computePEPX_Aniostropic_double(const Matrix2d& F, Vector2d direction, const double& scale) {
    //double x = 0, y = 1;
    double rate = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    direction /= rate;

    SVDResult2D_double svdResult = SingularValueDecomposition2D_double(F);
    Matrix2d U, V, R, S, sigma;
    //U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    S = V * sigma * V.transpose();

    double I4 = direction.transpose() * S * direction;
    double I5 = direction.transpose() * S.transpose() * S * direction;

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
//int n = 0;
void fem_explicit2D(mesh2D& mesh) {
    for (int ii = 0; ii < mesh.triangleNum; ii++) {
        MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
        //cout << "F:\n" << F << endl;
        Matrix2d PEPF; PEPF.setZero();
        //Matrix2d PEPF = computePEPX_ARAP2D_double(F);
        PEPF += computePEPX_StableNHK_double(F);
        PEPF += computePEPX_Aniostropic_double(F, Vector2d(1, 0), 1);
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
        double mass = 0.01;
        mass = 1.0 / mass;

        double acx = mesh.forces[ii][0] * mass;
        double acy = mesh.forces[ii][1] * mass - 9.8;

        mesh.forces[ii][0] = 0;
        mesh.forces[ii][1] = 0;

        double time = 0.0001;

        mesh.velocities[ii][0] += acx * time;
        mesh.velocities[ii][1] += acy * time;
        //if (n < 50000) {
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
        //}
        //else {
        //    return;
        //}
        //n++;
        mesh.vertexes[ii][0] += mesh.velocities[ii][0] * time;
        mesh.vertexes[ii][1] += mesh.velocities[ii][1] * time;

        if (mesh.vertexes[ii][1] < -1) {
            mesh.vertexes[ii][1] = -1; mesh.velocities[ii][1] *= -0.1;
        }

        if (mesh.vertexes[ii][0] < -1) {
            mesh.vertexes[ii][0] = -1; mesh.velocities[ii][0] *= -0.1;
        }

        if (mesh.vertexes[ii][1] > 1) {
            mesh.vertexes[ii][1] = 1; mesh.velocities[ii][1] *= -0.1;
        }

        if (mesh.vertexes[ii][0] > 1) {
            mesh.vertexes[ii][0] = 1; mesh.velocities[ii][0] *= -0.1;
        }
    }
}

void initMesh(mesh2D &mesh, int type, double scale) {

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
        Matrix2d PEPF = computePEPX_ARAP2D_double(F);

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
        Matrix4d Hq = project_ARAP_H_2D(F);
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

    double errorRate = 0.001;

    //PCG main loop
    while (deltaN > errorRate * errorRate * delta0) {
        vector<Vector2d> q(mesh.vertexNum, Vector2d(0, 0));
        for (int ii = 0; ii < mesh.triangleNum; ii++) {
            MatrixXd PFPX = computePFPX2D_double(mesh.DM_triangle_inverse[ii]);
            MatrixXd F = calculateDms2D_double(mesh.vertexes, mesh.triangles[ii], 0) * mesh.DM_triangle_inverse[ii];
            Matrix4d Hq = project_ARAP_H_2D(F);
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

    float stepSize = 0.001;
    for (int i = 0; i < mesh.vertexNum; i++) {
        mesh.vertexes[i][0] -= dX[i][0] * stepSize;
        mesh.vertexes[i][1] -= dX[i][1] * stepSize;
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
