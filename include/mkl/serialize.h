#pragma once

// Serialization utilities for sparse matrices and solver components.

#include <fstream>
#include <string>
#include <vector>
#include <cstring>

#include "mkl.h"
#include "mkl_spblas.h"
#include "csr.h"

namespace serialize {

// Write a CSR matrix to a binary file.
template<typename T>
bool write_csr(const std::string &path, const csr<T> &matrix) {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    // Write dimensions
    file.write(reinterpret_cast<const char*>(&matrix.n), sizeof(matrix.n));
    file.write(reinterpret_cast<const char*>(&matrix.m), sizeof(matrix.m));

    // Write values
    std::size_t v_size = matrix.v.size();
    file.write(reinterpret_cast<const char*>(&v_size), sizeof(v_size));
    file.write(reinterpret_cast<const char*>(matrix.v.data()), v_size * sizeof(T));

    // Write row pointers
    std::size_t r_size = matrix.r.size();
    file.write(reinterpret_cast<const char*>(&r_size), sizeof(r_size));
    file.write(reinterpret_cast<const char*>(matrix.r.data()), r_size * sizeof(MKL_INT));

    // Write column indices
    std::size_t c_size = matrix.c.size();
    file.write(reinterpret_cast<const char*>(&c_size), sizeof(c_size));
    file.write(reinterpret_cast<const char*>(matrix.c.data()), c_size * sizeof(MKL_INT));

    return file.good();
}

// Read a CSR matrix from a binary file.
template<typename T>
bool read_csr(const std::string &path, csr<T> &matrix) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;

    // Read dimensions
    std::size_t n, m;
    file.read(reinterpret_cast<char*>(&n), sizeof(n));
    file.read(reinterpret_cast<char*>(&m), sizeof(m));

    matrix = csr<T>(n, m);

    // Read values
    std::size_t v_size;
    file.read(reinterpret_cast<char*>(&v_size), sizeof(v_size));
    matrix.v.resize(v_size);
    file.read(reinterpret_cast<char*>(matrix.v.data()), v_size * sizeof(T));

    // Read row pointers
    std::size_t r_size;
    file.read(reinterpret_cast<char*>(&r_size), sizeof(r_size));
    matrix.r.resize(r_size);
    file.read(reinterpret_cast<char*>(matrix.r.data()), r_size * sizeof(MKL_INT));

    // Read column indices
    std::size_t c_size;
    file.read(reinterpret_cast<char*>(&c_size), sizeof(c_size));
    matrix.c.resize(c_size);
    file.read(reinterpret_cast<char*>(matrix.c.data()), c_size * sizeof(MKL_INT));

    return file.good();
}

// Write a dense matrix to a binary file.
template<typename T>
bool write_dense(const std::string &path, const T *data,
                 std::size_t rows, std::size_t cols) {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    file.write(reinterpret_cast<const char*>(data), rows * cols * sizeof(T));

    return file.good();
}

// Read a dense matrix from a binary file.
template<typename T>
bool read_dense(const std::string &path, T **data,
                std::size_t &rows, std::size_t &cols) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;

    file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    file.read(reinterpret_cast<char*>(&cols), sizeof(cols));

    *data = (T*)mkl_malloc(sizeof(T) * rows * cols, 64);
    file.read(reinterpret_cast<char*>(*data), rows * cols * sizeof(T));

    return file.good();
}

// Write problem parameters to a text file for reference.
inline bool write_params(const std::string &path,
                         std::size_t n_points,
                         std::size_t n_blocks_dim,
                         double tau,
                         double beta) {
    std::ofstream file(path);
    if (!file) return false;

    file << "n_points=" << n_points << "\n";
    file << "n_blocks_dim=" << n_blocks_dim << "\n";
    file << "n_blocks=" << n_blocks_dim * n_blocks_dim << "\n";
    file << "tau=" << tau << "\n";
    file << "beta=" << beta << "\n";

    return file.good();
}

} // namespace serialize
