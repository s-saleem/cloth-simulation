#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

//Assignment 4 code 
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <igl/triangle_triangle_adjacency.h>
#include <init_state.h>
#include <find_max_vertices.h>
#include <fixed_point_constraints.h>

#include <mass_matrix_mesh.h>
#include <linearly_implicit_euler.h>
#include <dsvd.h>

#include <dphi_cloth_triangle_dX.h>
#include <T_cloth.h>
#include <V_membrane_corotational.h>
#include <dV_membrane_corotational_dq.h>
#include <d2V_membrane_corotational_dq2.h>
#include <dV_cloth_gravity_dq.h>
#include <V_spring_particle_particle.h>
#include <dV_spring_particle_particle_dq.h>
#include <assemble_forces.h>
#include <assemble_stiffness.h>

//collision detection stuff
#include <collision_detection_cloth_sphere.h>
#include <velocity_filter_cloth_sphere.h>

#include <mass_vector.h>
#include <stretch_constraint.h>
#include <stretch_constraint_gradient.h>
#include <bend_constraint.h>
#include <bend_constraint_gradient.h>

//Variable for geometry
Eigen::MatrixXd V, V_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F, F_skin; //faces of simulation mesh
Eigen::MatrixXi E; //edges of simulation mesh (which will become springs)
Eigen::VectorXd l0; //original length of all edges in the mesh

//   TT   #F by #3 adjacent matrix, the element i,j is the id of the triangle
//        adjacent to the j edge of triangle i
//   TTi  #F by #3 adjacent matrix, the element i,j is the id of edge of the
//        triangle TT(i,j) that is adjacent with triangle i
Eigen::MatrixXi TT;
Eigen::MatrixXi TTi;

Eigen::MatrixXd V_sphere, V_sphere_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F_sphere, F_sphere_skin; //faces of simulation mesh

Eigen::SparseMatrixd N;

//material parameters
double density = 1;
double YM = 1e6; //young's modulus
double mu = 0.4; //poissons ratio
double C = (YM*mu)/((1.0+mu)*(1.0-2.0*mu));
double D = YM/(2.0*(1.0+mu));

//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::SparseMatrixd P;
Eigen::VectorXd x0; 

//mass matrix
Eigen::VectorXd M; //mass vector
Eigen::VectorXd a0; //areas
Eigen::MatrixXd dX; //dX matrices for computing deformation gradients

//scratch memory for assembly
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::VectorXd gravity;
Eigen::VectorXd wind;
Eigen::VectorXd qtmp;
Eigen::SparseMatrixd tmp_stiffness;

std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points; //need this for interaction 

//collision detection stuff
bool collision_detection_on = false; 
std::vector<unsigned int> collision_indices;
std::vector<Eigen::Vector3d> collision_normals;

bool fully_implicit = false;

//selection spring
double k_selected = 1e5;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {  

    double V_ele, T_ele, KE,PE;

    spring_points.clear();

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
    
    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        // spring_points.push_back(std::make_pair((P.transpose()*q+x0).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
        spring_points.push_back(std::make_pair(q.segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
    }


    /**
     * Step 1 - update velocity for external force not described in constraints
     * gravity - v_i += dt * g
     */
    Eigen::VectorXd G(qdot.size());
    
    for(int i = 0; i < qdot.size(); i++) {
        G[i] = gravity[i % 3];
    }

    Eigen::VectorXd W(qdot.size());
    
    for(int i = 0; i < qdot.size(); i++) {
        W[i] = wind[i % 3];// * sin(dt * M_PI) * sin(dt * M_PI);
    }

    Eigen::VectorXd f = G;// + W;
    for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
        dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, q.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
        f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
    }

    qdot += dt * f;// * 0;

    /**
     * Damp velocities
     */
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    Eigen::Vector3d cov = Eigen::Vector3d::Zero();
    double mass = 0;
    for(int i = 0; i < V.rows(); i++) {
        mass += M(i);
        com += q.segment<3>(3*i) * M(i);
        cov += qdot.segment<3>(3*i) * M(i);
    }
    com /= mass;
    cov /= mass;

    Eigen::Vector3d L = Eigen::Vector3d::Zero();
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    for(int i = 0; i < V.rows(); i++) {
        Eigen::Vector3d r = q.segment<3>(3*i) - com;
        Eigen::Vector3d v = qdot.segment<3>(3*i);

        L += r.cross(M(i) * v);

        Eigen::Matrix3d r_skew;
        r_skew << 0, -r.z(), r.y(),
                r.z(), 0, -r.x(),
                -r.y(), r.x(), 0;

        I += r_skew * r_skew.transpose() * M(i);
    }

    Eigen::Vector3d omega = I.inverse() * L;

    for(int i = 0; i < V.rows(); i++) {
        Eigen::Vector3d r = q.segment<3>(3*i) - com;
        Eigen::Vector3d v = qdot.segment<3>(3*i);
        Eigen::Vector3d delta_v = cov + omega.cross(r) - v;

        qdot.segment<3>(3*i) += 0.005 * delta_v;
    }

    qtmp = q + dt * qdot;

    /**
     * Step 2 - generate collison constraints
     */
    // std::vector<int> sphere_collisions;
    // Eigen::Vector3d sphere_center = Eigen::Vector3d (0.0, 0.0, 0.4);
    // double radius = 0.22;
    // for(int i = 0; i < V.rows(); i++) {
    //     Eigen::Vector3d p = qtmp.segment<3>(3*i);

    //     if((p - sphere_center).norm() <= radius) {
    //         sphere_collisions.push_back(i);
    //     }
    // }


    /**
     * Step 3 - for X iterations project constraints
     */
    int solver_iterations = 1;
    for(int i = 0; i < solver_iterations; i++) {
        // //sphere collisions
        // for(int j = 0; j < sphere_collisions.size(); j++) {
        //     int idx = sphere_collisions[j];
        //     Eigen::Vector3d p1 = qtmp.segment<3>(3 * idx);

        //     Eigen::Vector3d n = (p1 - sphere_center);
        //     n /= n.norm();
        //     Eigen::Vector3d p2 = n * radius;

        //     double l = 0;
    
        //     double constraint = stretch_constraint(p1, p2, l);

        //     if(constraint < 0) {
        //         Eigen::Vector6d gradient = stretch_constraint_gradient(p1, p2);

        //         // -2 because only p1 will move
        //         Eigen::Vector6d delta_p = -2 * constraint / (gradient.norm() * gradient.norm()) * gradient;
                
        //         qtmp.segment<3>(3 * idx) += delta_p.segment<3>(0);
        //     }
        // }

        // fixed point constraints
        for(int j = 0; j < fixed_point_indices.size(); j++) {

            int idx = fixed_point_indices[j];
            Eigen::Vector3d p1 = qtmp.segment<3>(3 * idx);
            Eigen::Vector3d p2 = V.row(idx);

            double l = 0.0;
    
            double constraint = stretch_constraint(p1, p2, l);
            if(constraint > 0) {
                Eigen::Vector6d gradient = stretch_constraint_gradient(p1, p2);

                // -2 because only p1 will move
                Eigen::Vector6d delta_p = -2 * constraint / (gradient.norm() * gradient.norm()) * gradient;
                
                qtmp.segment<3>(3 * idx) += delta_p.segment<3>(0);
            }
        }

        // loop through stretch constraints
        double stretch_stiffness = 0.7;
        for(int j = 0; j < E.rows(); j++) {
            int idx1 = E(j, 0);
            int idx2 = E(j, 1);

            Eigen::Vector3d p1 = qtmp.segment<3>(3 * idx1);
            Eigen::Vector3d p2 = qtmp.segment<3>(3 * idx2);

            double l = l0(j);
    
            double constraint = stretch_constraint(p1, p2, l);
            if(abs(constraint) > 1e-8) {
                Eigen::Vector6d gradient = stretch_constraint_gradient(p1, p2);

                Eigen::Vector6d delta_p = -constraint / (gradient.norm() * gradient.norm()) * gradient;

                // keep k linear with iterations (only 1 iteration works right now)
                double k = 1 - pow(1 - stretch_stiffness, 1 / solver_iterations);

                double w1 = 1 / M(idx1);
                double w2 = 1 / M(idx2);

                qtmp.segment<3>(3 * idx1) += w1 / (w1 + w2) * k * delta_p.segment<3>(0);
                qtmp.segment<3>(3 * idx2) += w2 / (w1 + w2) * k * delta_p.segment<3>(3);
            }
        }

        double bend_stiffness = 0.8;
        for(int j = 0; j < TT.rows(); j++) {
            for(int k = 0; k < TT.row(j).size(); k++) {
                if(TT(j, k) > j) {
                    Eigen::Vector3i verts_T1 = F.row(j);
                    Eigen::Vector3i verts_T2 = F.row(TT(j, k));

                    int idx1 = verts_T1(k);
                    int idx2 = verts_T1((k + 1) % 3);
                    int idx3 = verts_T1((k + 2) % 3);
                    int idx4 = verts_T2((TTi(j, k) + 2) % 3);

                    Eigen::Vector3d p1 = qtmp.segment<3>(3 * idx1);
                    Eigen::Vector3d p2 = qtmp.segment<3>(3 * idx2);
                    Eigen::Vector3d p3 = qtmp.segment<3>(3 * idx3);
                    Eigen::Vector3d p4 = qtmp.segment<3>(3 * idx4);

                    Eigen::Vector3d n1 = (p2 - p1).cross(p3 - p1)/ (p2 - p1).cross(p3 - p1).norm();
                    Eigen::Vector3d n2 = (p2 - p1).cross(p4 - p1)/ (p2 - p1).cross(p4 - p1).norm();
                    double d = n1.dot(n2);

                    double constraint = bend_constraint(p1, p2, p3, p4);

                    if(d*d < 1) {
                        Eigen::Vector12d gradient = bend_constraint_gradient(p1, p2, p3, p4);
                        Eigen::Vector12d delta_p = -sqrt(1 - d*d) *constraint / (gradient.norm() * gradient.norm()) * gradient;

                        // keep k linear with iterations (only 1 iteration works right now)
                        double k = 1 - pow(1 - bend_stiffness, 1 / solver_iterations);

                        double w1 = 1 / M(idx1);
                        double w2 = 1 / M(idx2);
                        double w3 = 1 / M(idx3);
                        double w4 = 1 / M(idx4);

                        qtmp.segment<3>(3 * idx1) += 4 * w1 / (w1 + w2 + w3 + w4) * k * delta_p.segment<3>(0);
                        qtmp.segment<3>(3 * idx2) += 4 * w2 / (w1 + w2 + w3 + w4) * k * delta_p.segment<3>(3);
                        qtmp.segment<3>(3 * idx3) += 4 * w3 / (w1 + w2 + w3 + w4) * k * delta_p.segment<3>(6);
                        qtmp.segment<3>(3 * idx4) += 4 * w4 / (w1 + w2 + w3 + w4) * k * delta_p.segment<3>(9);
                    }
                }
            }
        }
    }


    /**
     * Step 4 - update velocity with new postions (after constraints)
     */

    qdot = (qtmp - q) / dt;
    q = qtmp;

}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    // Visualize::update_vertex_positions(0, P.transpose()*q + x0);
    Visualize::update_vertex_positions(0, q);

}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    if(key =='N') {
        std::cout<<"toggle integrators \n";
        fully_implicit = !fully_implicit;
    } 
    if(key == 'C') {
        collision_detection_on = !collision_detection_on;
        Visualize::set_visible(1, collision_detection_on);
    }

    return false;
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 
    igl::readOBJ("../data/square_cloth.obj", V, F);
    igl::readOBJ("../data/sphere.obj", V_sphere, F_sphere);

    //setup simulation 
    init_state(q,qdot,V);

    // Set up edges and adjacent faces
    igl::edges(F, E);
    igl::edge_lengths(V, E, l0);
    igl::triangle_triangle_adjacency(F, TT, TTi);

    //add geometry to scene
    V_skin = V;
    F_skin = F;
    N.resize(V.rows(), V.rows());
    N.setIdentity();

    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);

    //add collision sphere to scene 
    V_sphere_skin = V_sphere;
    F_sphere_skin = F_sphere;
    N.resize(V.rows(), V.rows());
    N.setIdentity();


    Visualize::add_object_to_scene(V_sphere,F_sphere, V_sphere_skin, F_sphere_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::set_visible(1, collision_detection_on);
    
    //Mass Matrix
    // mass_matrix_mesh(M, q,  V, F, density, a0);
    mass_vector(M, q, V, F, density);

    if(M.size() == 0) {
        std::cout<<"Mass vector not implemented, exiting.\n";
        exit(1);
    }
    
    //should be max verts for cloth simulation
    find_max_vertices(fixed_point_indices, V, 0.001);
    // int edges = E.rows();
    // E.conservativeResize(E.rows() + fixed_point_indices.size(), E.cols());
    // l0.conservativeResize(l0.rows() + fixed_point_indices.size());
    // for(int i = 0; i < fixed_point_indices.size(); i++) {
    //     E(edges + i, 0) = fixed_point_indices[i];
    //     E(edges + i, 1) = fixed_point_indices[i];
    //     l0(edges + i) = 0.0;
        
    // }
    
    P.resize(q.rows(),q.rows());
    P.setIdentity();
    fixed_point_constraints(P, q.rows(), fixed_point_indices);

    for(int i = 0; i < fixed_point_indices.size(); i++) {
        M(fixed_point_indices[i]) = 999.0;
    }
    
    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else    
    
    //constant gravity vector
    gravity.resize(3, 1);
    gravity << 0, -9.5, 0;

    //constant wind vector
    wind.resize(3, 1);
    wind << 0, 0, 0.7;


    //std::cout<<"Gravity "<<gravity.transpose()<<"\n";
    
    //correct M, q and qdot so they are the right size
    // q = P*q;
    // qdot = P*qdot;
    // M = P*M;
    
    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

