#include <mass_vector.h>

void mass_vector(Eigen::VectorXd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                double density) {

    M.resize(q.size() / 3.0);
    for(int i = 0; i < F.rows(); i++) {
         // triangle vertices
        Eigen::RowVector3i element;
        element = F.row(i);

        for(int j = 0; j < element.size(); j++) {
            int vertex = element[j];
            M[vertex] += density / 3.0;
        }
    }
}