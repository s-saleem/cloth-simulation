#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(q_size - 3 * indices.size());

    P.resize(q_size - 3 * indices.size(), q_size);
    int j = 0;
    for(int i = 0; i < q_size / 3; i++) {
        if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
            tripletList.push_back(T(3*j, 3*i, 1));
            tripletList.push_back(T(3*j + 1, 3*i + 1, 1));
            tripletList.push_back(T(3*j + 2, 3*i + 2, 1));
            j++;
        }
    }

    P.setFromTriplets(tripletList.begin(), tripletList.end());

}