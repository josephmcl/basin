#pragma once 

#include <vector>
#include <tuple>
#include <numbers>
#include <type_traits>

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

enum class storage_type {
    dense, 
    csr
};

using storage = storage_type;
struct solver_info {
    bool         overwrite_b;
    storage_type storage;
};

constexpr bool overwrite = true;
constexpr bool pass_both = false;

constexpr solver_info solver_info_key[] = {
    {overwrite, storage::csr  },
    {overwrite, storage::dense}}; // Used by dense cholesky.
    
enum class solver_t : std::size_t {
    sparse_qr      = 0,
    dense_cholesky = 1
};

/* */
namespace solver {
    storage_type storage(solver_t s);
    bool overwrite(solver_t s);
};

constexpr double π = std::numbers::pi_v<double>;


using range_t = numerics::linrange<double>;
using real_t = type::real_t; 
using vector_t = std::vector<real_t>;
auto static to_real_t = [](std::size_t n){return static_cast<real_t>(n);};

using csr_t = csr<type::real_t>;
using csr_v = std::vector<csr<type::real_t>>;

