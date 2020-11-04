#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 
        
        f = Eigen::VectorXd::Zero(q.size());
        
        for (int i = 0; i < F.rows(); i++) {
            // tetrahedron vertices
            Eigen::RowVector3i element;
            element = F.row(i);

            // dphi/dx
            Eigen::Matrix<double, 1,9> tmp_row;
            tmp_row = dX.row(i);

            Eigen::Matrix3d dphi_dX;
            dphi_dX = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());

            // compute internal tertrahedron force
            Eigen::Vector9d f_e;
            dV_membrane_corotational_dq(f_e, q, dphi_dX, V, element, a0[i], mu, lambda);
            f_e *= -1;

            for(int i = 0; i < element.size() ; i++) {
                f.segment<3>(element(i)*3) += f_e.segment<3>(i*3);
            }
        }
       
    };
