#ifndef DOUBLEBEAM_CPP_EIGEN_HELPERS_H
#define DOUBLEBEAM_CPP_EIGEN_HELPERS_H


#include <complex>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


// typedefs for Eigen::Tensor inkeeping with the Eigen::Matrix conventions
namespace Eigen {
    using Tensor1d = Eigen::Tensor<double, 1>;
    using Tensor2cd = Eigen::Tensor<std::complex<double>, 2>;
    using Tensor3cd = Eigen::Tensor<std::complex<double>, 3>;
} // namespace Eigen


template <typename Scalar>
auto as_matrix(Eigen::Tensor<Scalar, 2>& tensor) {
    return Eigen::Map<Eigen::Matrix<Scalar, 2, 2>>(tensor.data(), 2, 2);
}


/**
 * Helper function to extract nth 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return nth element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
Eigen::Matrix<Scalar, 2, 2> nth_matrix(Eigen::Tensor<Scalar, 3>& in,
                                       typename Eigen::Tensor<Scalar, 3>::Index n) {
    // assignment evals the tensor which is required because only evaluated tensors have data()
    // member.
    Eigen::Tensor<Scalar, 2> a = in.chip(n, N_index).reshape(std::array<long, 2>{2, 2});
    return Eigen::Map<Eigen::Matrix<Scalar, 2, 2>>(a.data(), 2, 2);
}


/**
 * Helper function to extract last 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return Last element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
Eigen::Matrix<Scalar, 2, 2> last_element(Eigen::Tensor<Scalar, 3>& in) {
    return nth_matrix<Scalar, N_index>(in, in.dimension(N_index) - 1);
}

/**
 * Helper function to extract first 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return first element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
auto first_element(Eigen::Tensor<Scalar, 3>& in) {
   return nth_matrix<Scalar, N_index>(in, 0);
}

/**
 * Helper function to extract nth 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return nth element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
Eigen::Matrix<Scalar, 2, 2> nth_matrix(const Eigen::Tensor<Scalar, 3>& in,
                                       typename Eigen::Tensor<Scalar, 3>::Index n) {
    // assignment evals the tensor which is required because only evaluated tensors have data()
    // member.
    Eigen::Tensor<Scalar, 2> a = in.chip(n, N_index).reshape(std::array<long, 2>{2, 2});
    return Eigen::Map<Eigen::Matrix<Scalar, 2, 2>>(a.data(), 2, 2);
}


/**
 * Helper function to extract last 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return Last element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
Eigen::Matrix<Scalar, 2, 2> last_element(const Eigen::Tensor<Scalar, 3>& in) {
    return nth_matrix<Scalar, N_index>(in, in.dimension(N_index) - 1);
}

/**
 * Helper function to extract first 2x2 matrix from a Nx2x2 tensor.
 * @tparam Scalar Scalar data type of tensor elements.
 * @tparam N_index Index of the "long" axis of the tensor which contains N elements.
 * @param in Nx2x2 tensor
 * @return first element (matrix) of the tensor.
 */
template <typename Scalar, int N_index = 0>
auto first_element(const Eigen::Tensor<Scalar, 3>& in) {
    return nth_matrix<Scalar, N_index>(in, 0);
}




#endif // DOUBLEBEAM_CPP_EIGEN_HELPERS_H