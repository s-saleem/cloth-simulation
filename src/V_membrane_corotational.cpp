#include <V_membrane_corotational.h>
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

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    
   

    Eigen::Matrix3d F;
    getF(F, q, dX, V, element);

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F);
    Eigen::Vector3d s_values = svd.singularValues();
    s_values(0) = sqrt(s_values(0));
    s_values(1) = sqrt(s_values(1));
    s_values(2) = sqrt(s_values(2));

    double deformation = 0;
    double volume_preservation = 0;
    for(int i = 0; i < s_values.size(); i++) {
        deformation += (s_values(i) - 1)*(s_values(i) - 1);
        volume_preservation += s_values(i);
    }
    deformation *= mu;
    volume_preservation = lambda/2 * (volume_preservation - 3)*(volume_preservation - 3);

    
    double psi = deformation + volume_preservation;
    energy = area * psi;

    // std::cout << "s_values: " << s_values << std::endl;
    // if (n.transpose() * N < 0.99) {
    //     std::cout << "V : n dot N: " << n.transpose() * N << std::endl;
    //     exit(0);
    // }
}
