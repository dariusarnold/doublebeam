/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DOUBLEBEAM_CPP_EIGEN_HELPERS_H
#define DOUBLEBEAM_CPP_EIGEN_HELPERS_H


#include <complex>
#include <filesystem>
#include <fstream>
#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <cnpy.h>

#include "utils.hpp"


// typedefs for Eigen::Tensor inkeeping with the Eigen::Matrix conventions
namespace Eigen {
    using Tensor1d = Eigen::Tensor<double, 1>;
    using Tensor2cd = Eigen::Tensor<std::complex<double>, 2>;
    using Tensor3cd = Eigen::Tensor<std::complex<double>, 3>;
    using Tensor4cd = Eigen::Tensor<std::complex<double>, 4>;
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
    return Eigen::Map<Eigen::Matrix<Scalar, 2, 2>>(a.data(), 2, 2).transpose();
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
    return Eigen::Map<Eigen::Matrix<Scalar, 2, 2>>(a.data(), 2, 2).transpose();
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

namespace impl {
    /**
     * Create Eigen::Tensor with sizes of dimensions specified in a vector.
     * @tparam Scalar Scalar type of tensor elements.
     * @tparam Dimensions Rank of tensor (number of dimensions).
     * @tparam Layout Layout of returned tensor, either column major or row major.
     * @tparam I Index sequence for vector.
     * @param data Pointer to raw array data of tensor.
     * @param shape Array containing shape of tensor (size of dimensions).
     * @return Tensor with data given from pointer in given shape.
     */
    template <typename Scalar, int Dimensions, Eigen::StorageOptions Layout = Eigen::ColMajor,
              size_t... I>
    Eigen::Tensor<Scalar, Dimensions, Layout> make_tensor(Scalar* data, std::vector<size_t> shape,
                                                          std::index_sequence<I...>) {
        return Eigen::TensorMap<Eigen::Tensor<Scalar, Dimensions, Layout>>(data, shape[I]...);
    }
} // namespace impl


/**
 * Load npy file and return it as tensor.
 * While the datatype and the dimensions are specified in the npy file, they are required to declare
 * the return type of this function, so they must be given by the user.
 * @tparam Scalar Scalar element type of npy file.
 * @tparam Dimensions Rank of tensor/number of dimensions of the numpy array in the file.
 * @param path Path to .npy file.
 * @return
 */
template <typename Scalar, int Dimensions>
Eigen::Tensor<Scalar, Dimensions> load_npy(const std::filesystem::path& path) {

    std::ifstream file{path, std::ios::binary};
    if (not file.is_open()) {
        throw std::runtime_error(("Failed to open file " + path.string()));
    }
    const int initial_buffer_size = 250;
    std::vector<unsigned char> buffer(initial_buffer_size);
    // size of file content in bytes
    const int SIZE_MAGIC_STRING = 6;
    const int SIZE_VERSION_NUMBER = 2;
    const int SIZE_HEADER_LEN = 2;
    int offset = 0;
    // ifstream::read expects pointer to signed char, while the file content is unsigned
    auto buffer_as_signed = [&]() { return reinterpret_cast<char*>(buffer.data()); };
    file.read(buffer_as_signed(), SIZE_MAGIC_STRING);
    // TODO magic string should be \x93NUMPy, why does this read \x23????
    if (not std::strcmp(buffer_as_signed(), "\x23NUMPY")) {
        throw std::runtime_error(path.string() +
                                 " is not a valid npy file. "
                                 "Read invalid magic string " +
                                 std::string(buffer.begin(), buffer.end() + SIZE_MAGIC_STRING) +
                                 ". Expected \x93NUMPY");
    }
    offset += SIZE_MAGIC_STRING;
    file.get(buffer_as_signed() + offset, SIZE_VERSION_NUMBER);
    if (not std::strcmp(buffer_as_signed() + offset, "\x00\x01")) {
        throw std::runtime_error(
            "Invalid version number " +
            std::string(buffer.begin() + offset, buffer.end() + offset + SIZE_VERSION_NUMBER));
    }
    offset += SIZE_VERSION_NUMBER;
    file.read(buffer_as_signed() + offset, SIZE_HEADER_LEN);
    // extract header length as little endian 2 byte unsigned int
    std::uint16_t header_length = *reinterpret_cast<uint16_t*>(buffer.data());
    offset += SIZE_HEADER_LEN;
    if (header_length > initial_buffer_size - offset) {
        buffer.resize(offset + header_length);
    }
    file.read(buffer_as_signed() + offset, header_length);

    size_t word_size;
    std::vector<size_t> shape;
    bool col_wise_order;
    // delegate header parsing
    cnpy::parse_npy_header(reinterpret_cast<unsigned char*>(buffer.data()), word_size, shape,
                           col_wise_order);
    cnpy::NpyArray P = cnpy::npy_load(path.string());
    if (Dimensions != P.shape.size()) {
        throw std::runtime_error(impl::Formatter()
                                 << "Got " << Dimensions
                                 << " dimensions as template argument, but file contains "
                                 << P.shape.size() << ".");
    }
    auto is = std::make_index_sequence<Dimensions>();
    if (not col_wise_order) {
        // load the tensor in row wise order and shuffle it
        auto tensor =
            impl::make_tensor<Scalar, Dimensions, Eigen::RowMajor>(P.data<Scalar>(), P.shape, is);
        auto shuffle = std::vector<int>(Dimensions);
        std::iota(shuffle.begin(), shuffle.end(), 0);
        std::reverse(shuffle.begin(), shuffle.end());
        return tensor.swap_layout().shuffle(shuffle).eval();
    }
    return impl::make_tensor<Scalar, Dimensions>(P.data<Scalar>(), P.shape, is);
}

#endif // DOUBLEBEAM_CPP_EIGEN_HELPERS_H