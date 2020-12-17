#include <EigenTypes.h>
#include <iostream>

double bend_constraint(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d p4) {
    Eigen::Vector3d n1 = (p2 - p1).cross(p3 - p1);
    n1 /= n1.norm();

    Eigen::Vector3d n2 = (p2 - p1).cross(p4 - p1);
    n2 /= n2.norm();

    return acos(n1.dot(n2)) - M_PI;
}