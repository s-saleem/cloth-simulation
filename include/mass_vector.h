#include <EigenTypes.h>

//Input: 
//  V - the nx3 matrix of undeformed vertex positions
//  F - the mx3 matrix of triangle-vertex indices
//  density - the density of the cloth material
//Output:
//  M - mass vector for the entire mesh
void mass_vector(Eigen::VectorXd &M, Eigen::Ref<const Eigen::VectorXd> q,
                Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                double density);