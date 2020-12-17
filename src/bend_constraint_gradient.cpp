#include <EigenTypes.h>
#include <iostream>

Eigen::Vector12d bend_constraint_gradient(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d p4) {
    p2 = p2 - p1;
    p3 = p3 - p1;
    p4 = p4 - p1;
    Eigen::Vector3d n1 = (p2).cross(p3);
    n1 /= n1.norm();

    Eigen::Vector3d n2 = (p2).cross(p4);
    n2 /= n2.norm();

    double d = n1.dot(n2);
    double s = -1 / sqrt(1 - d*d);

    // std::cout << "d: " << d << std::endl;
    // std::cout << "s: " << s << std::endl;

    // Test
    Eigen::Vector3d grad_p3 = (p2.cross(n2) + (n1.cross(p2) * d)) / p2.cross(p3).norm();
    Eigen::Vector3d grad_p4 = (p2.cross(n1) + (n2.cross(p2) * d)) / p2.cross(p4).norm();

    Eigen::Vector3d grad_p2 = -(p3.cross(n2) + (n1.cross(p3) * d)) / p2.cross(p3).norm() - (p4.cross(n1) + (n2.cross(p4) * d)) / p2.cross(p4).norm();


    Eigen::Vector3d grad_p1 = -grad_p2 - grad_p3 - grad_p4;

    Eigen::Vector12d grad;
    grad << grad_p1, grad_p2, grad_p3, grad_p4;

    // std::cout << "grad: \n" << grad << std::endl;

    return grad;
}