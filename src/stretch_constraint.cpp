#include <stretch_constraint.h>

double stretch_constraint(Eigen::Vector3d p1, Eigen::Vector3d p2, double l0) {
    return (p1 - p2).norm() - l0;
}