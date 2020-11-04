#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    Eigen::Matrix36d B;
    B << -Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity();

    Eigen::Vector6d q;
    q << q0, q1;

    Eigen::Vector3d Bq = B * q;

    f << stiffness * (1 - l0 / Bq.norm()) * B.transpose() * Bq;
    
}