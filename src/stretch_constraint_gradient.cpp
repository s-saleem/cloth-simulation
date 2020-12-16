#include <stretch_constraint_gradient.h>

Eigen::Vector6d stretch_constraint_gradient(Eigen::Vector3d p1, Eigen::Vector3d p2) {
    Eigen::Vector6d gradient;

    double norm = (p1 - p2).norm();

    gradient.segment<3>(0) = (p1 - p2)/norm;
    gradient.segment<3>(3) = -(p1 - p2)/norm;

    return gradient;
}