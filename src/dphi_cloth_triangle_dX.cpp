#include <dphi_cloth_triangle_dX.h>

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    // Get vertex coords X
    Eigen::Vector3d X0, X1, X2, X3;

    X0 = V.row(element(0)).transpose();
    X1 = V.row(element(1)).transpose();
    X2 = V.row(element(2)).transpose();   
    
    // get point differences from X0
    Eigen::Matrix32d T;
    T << X1 - X0, X2 - X0;

    Eigen::RowVector2d neg_one;
    neg_one << -1,-1;

    dphi << neg_one * (T.transpose()*T).inverse() * T.transpose(), (T.transpose()*T).inverse() * T.transpose();
}