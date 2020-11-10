#include <mass_matrix_mesh.h>
#include <iostream>

typedef Eigen::Triplet<double> T;
static void copyMatrixToTripletList(std::vector<T> &tripletList, int startRow, int startCol, Eigen::MatrixXd K) {
    for (int i = 0; i < K.rows(); i++) {
        for (int j = 0; j < K.cols(); j++) {
            tripletList.push_back(T(startRow + i, startCol + j, K(i, j)));
        }
    }
}

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {

    /** imported from MATLAB using symbolic integrator
    * N = [(1-phi1-phi2)*eye(3) phi1*eye(3) phi2*eye(3)]
    * int(int(int(N'*N, phi2, 0, 1-phi1), phi1, 0, 1)
    */
    Eigen::MatrixXd M_triangle(9,9);
    M_triangle << 1/12.0,    0,    0, 1/24.0,    0,    0, 1/24.0,    0,    0,
                    0, 1/12.0,    0,    0, 1/24.0,    0,    0, 1/24.0,    0,
                    0,    0, 1/12.0,    0,    0, 1/24.0,    0,    0, 1/24.0,
                    1/24.0,    0,    0, 1/12.0,    0,    0, 1/24.0,    0,    0,
                    0, 1/24.0,    0,    0, 1/12.0,    0,    0, 1/24.0,    0,
                    0,    0, 1/24.0,    0,    0, 1/12.0,    0,    0, 1/24.0,
                    1/24.0,    0,    0, 1/24.0,    0,    0, 1/12.0,    0,    0,
                    0, 1/24.0,    0,    0, 1/24.0,    0,    0, 1/12.0,    0,
                    0,    0, 1/24.0,    0,    0, 1/24.0,    0,    0, 1/12.0;

    M_triangle *= 2;
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows());

    for(int i = 0; i < F.rows(); i++) {
        // triangle vertices
        Eigen::RowVector3i element;
        element = F.row(i);

        Eigen::MatrixXd M_i = density * areas[i] * M_triangle;

        for(int j = 0; j < element.size(); j++) {
            for(int k = 0; k < element.size(); k++) {
                copyMatrixToTripletList(tripletList, element(j)*3, element(k)*3, M_i.block<3, 3>(3*j, 3*k));
            }
        }
    }

    M.resize(q.size(), q.size());
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}
 