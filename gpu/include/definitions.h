#pragma once 

#include <vector>
#include <tuple>
#include <numbers>

#include "ranges.h"
#include "csr.h"

namespace type {
    using real_t = double;
};

template <typename T>
using type_p = std::tuple<std::size_t, T*>;
using real_p = type_p<type::real_t>;

template <typename T>
using vv = std::vector<std::vector<T>>;

enum class solver {
    sparse_qr, 
    dense_cholesky,
    rfp_cholesky,
    pardiso
};

constexpr double π = std::numbers::pi_v<double>;


using range_t = numerics::linrange<double>;
using real_t = type::real_t; 
using vector_t = std::vector<real_t>;
auto static to_real_t = [](std::size_t n){return static_cast<real_t>(n);};

using csr_t = csr<type::real_t>;
using csr_v = std::vector<csr<type::real_t>>;