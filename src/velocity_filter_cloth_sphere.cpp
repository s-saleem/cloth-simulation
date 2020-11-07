#include <velocity_filter_cloth_sphere.h>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {

    for(int i = 0; i < indices.size(); ++i) {
        Eigen::Vector3d velocity = qdot.segment<3>(3*i);
        Eigen::Vector3d normal = normals[i];
        
        double ntv = normal.transpose() * velocity;
        if ( ntv < 0) {
            velocity -= ntv * normal;
        }
    }

}