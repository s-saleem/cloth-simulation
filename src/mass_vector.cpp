#include <mass_vector.h>

void mass_vector(Eigen::VectorXd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                double density) {

    M.resize(q.size() / 3.0);
    for(int i = 0; i < F.rows(); i++) {
         // triangle vertices
        Eigen::RowVector3i element;
        element = F.row(i);

        Eigen::Vector3d p0 = q.segment<3>(3 * element(0));
        Eigen::Vector3d p1 = q.segment<3>(3 * element(1));
        Eigen::Vector3d p2 = q.segment<3>(3 * element(2));
        double area = (p1 - p0).cross(p2 - p0).norm() / 2.0;
        for(int j = 0; j < element.size(); j++) {
            int vertex = element[j];
            M[vertex] += density * area / 3.0;
        }
    }
}