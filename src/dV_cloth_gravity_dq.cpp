#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {

    Eigen::VectorXd G(M.cols()); 
    for(int i = 0; i < G.rows()/3; i++) {
        G.segment<3>(3*i) = g;
    }

    fg = -M*G;
}
