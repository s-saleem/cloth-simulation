#include <dV_membrane_corotational_dq.h>
#include <iostream>

static void getF(Eigen::Matrix3d &F, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element) {
    // Get reference coords X
    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    // Get reference normal N
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0);
    N /= N.norm();

    // Get world coords x and normal n
    Eigen::Matrix3d x;
    for(int i = 0; i < element.size(); i++) {
        x.col(i) = q.segment<3>(3*element(i));
    }

    Eigen::Vector3d n = (x.col(1) - x.col(0)).cross(x.col(2) - x.col(0));
    n /= n.norm();

    Eigen::Matrix34d F1;
    F1 << x, n;
    Eigen::Matrix43d F2;
    F2 << dX, 
        N.transpose();

    F = F1*F2;
}

static void get_dF_dq(Eigen::Matrix99d &dF, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element) {
    // Get reference coords X
    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    // Get reference normal N
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0);
    N /= N.norm();

    // Get world coords x and normal n
    Eigen::Matrix3d x;
    for(int i = 0; i < element.size(); i++) {
        x.col(i) = q.segment<3>(3*element(i));
    }

    Eigen::Vector3d n = (x.col(1) - x.col(0)).cross(x.col(2) - x.col(0));
    n /= n.norm();
    
    // Calcualte dF/dq = B + N * dn/dq
    // B 
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 9);
    for(int i = 0; i < dX.rows(); i++) {
        for(int j = 0; j < dX.cols(); j++) {
            B(j, 3*i) = dX(i, j);
            B(j + 3, 3*i + 1) = dX(i, j);
            B(j + 6, 3*i + 2) = dX(i, j);
        }
    }

    // N
    Eigen::MatrixXd N_extended = Eigen::MatrixXd::Zero(9, 3);
    for(int i = 0; i < N.size(); i++) {
        N_extended(i, 0) = N(i);
        N_extended(i + 3, 1) = N(i);
        N_extended(i + 6, 2) = N(i);
    }

    // dn/dq
    Eigen::Vector3d deltaX1, deltaX2;
    deltaX1 = x.col(1) - x.col(0);
    deltaX2 = x.col(2) - x.col(0);
    Eigen::Matrix3d dx1_cross, dx2_cross;
    dx1_cross << 0, -deltaX1.z(), deltaX1.y(),
        deltaX1.z(), 0, -deltaX1.x(),
        -deltaX1.y(), deltaX1.x(), 0;

    dx2_cross << 0, -deltaX2.z(), deltaX2.y(),
        deltaX2.z(), 0, -deltaX2.x(),
        -deltaX2.y(), deltaX2.x(), 0;

    Eigen::MatrixXd d_delta_x1(3, 9), d_delta_x2(3,9);
    d_delta_x1 << -Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero();
    d_delta_x2 << -Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

    double norm = deltaX1.cross(deltaX2).norm();

    Eigen::MatrixXd dn_dq(3,9);
    dn_dq = 1 / norm * (Eigen::Matrix3d::Identity() - n * n.transpose()) *
                            (dx1_cross * d_delta_x2 - dx2_cross * d_delta_x1);

    dF = B + N_extended * dn_dq;
}

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
    Eigen::Matrix3d F;
    getF(F, q, dX, V, element);

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient 
    S(0) = sqrt(S(0));
    S(1) = sqrt(S(1));
    S(2) = sqrt(S(2));

    // if(S(0) < 0.9 || S(1) < 0.9 || S(2) < 0.9 || S(0) > 1.01 || S(1) > 1.01 || S(2) > 1.01) {
    //     std::cout << "V S: " << S << std::endl;
    //     exit(0);
    // }

    // Calcualte dpsi/dF
    Eigen::Matrix3d dpsi_dS;
    dpsi_dS = 2 * mu * (Eigen::Matrix3d(S.asDiagonal()) - Eigen::Matrix3d::Identity()); 
    dpsi_dS += (lambda * (S(0) + S(1) + S(2) - 3) * Eigen::Matrix3d::Identity());

    Eigen::Matrix3d dpsi_dF;
    dpsi_dF = U * dpsi_dS * W.transpose();

    // Calcualte dF/dq = B + N * dn/dq
    Eigen::Matrix99d dF_dq;
    get_dF_dq(dF_dq, q, dX, V, element);


    // dV = dF^T/dq * dpsi/dF
    Eigen::Vector9d dpsi_flat;
    dpsi_flat << dpsi_dF.row(0).transpose(), dpsi_dF.row(1).transpose(), dpsi_dF.row(2).transpose();
    dV = area * dF_dq.transpose() * dpsi_flat;

    // std::cout << "dF/dq: \n" << dF_dq << std::endl;
    // std::cout << "dn/dq: \n" << dn_dq << std::endl;
    // std::cout << "dV: \n" << dV << std::endl;
    // exit(0);
    // if (n.transpose() * N < 0.99) {
    //     std::cout << "dV : n dot N: " << n.transpose() * N << std::endl;
    //     exit(0);
    // }
}
