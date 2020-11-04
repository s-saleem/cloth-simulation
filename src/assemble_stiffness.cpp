#include <assemble_stiffness.h>
#include <iostream>

typedef Eigen::Triplet<double> T;
static void copyMatrixToTripletList(std::vector<T> &tripletList, int startRow, int startCol, Eigen::MatrixXd K) {
    for (int i = 0; i < K.rows(); i++) {
        for (int j = 0; j < K.cols(); j++) {
            tripletList.push_back(T(startRow + i, startCol + j, K(i, j)));
        }
    }
}

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 
        
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows());

    for(int i = 0; i < F.rows(); i++) {
        // tetrahedron vertices
        Eigen::RowVector3i element;
        element = F.row(i);

        // dphi/dx
        Eigen::Matrix<double, 1,9> tmp_row;
        tmp_row = dX.row(i);

        Eigen::Matrix3d dphi_dX;
        dphi_dX = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());

        Eigen::Matrix99d K_e;
        d2V_membrane_corotational_dq2(K_e, q, dphi_dX, V, element, a0[i], mu, lambda);
        K_e *= -1;

        for(int j = 0; j < element.size(); j++) {
            for(int k = 0; k < element.size(); k++) {
                copyMatrixToTripletList(tripletList, element(j)*3, element(k)*3, K_e.block<3, 3>(3*j, 3*k));
            }
        }
    }

    K.resize(q.size(), q.size());
    K.setFromTriplets(tripletList.begin(), tripletList.end());
};
