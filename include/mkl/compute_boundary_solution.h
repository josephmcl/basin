// std
#include <functional>
// you
#include "ranges.h"
#include "definitions.h"
// lib
#include "mkl.h"


using real_t = type::real_t; 
using range_t = numerics::linrange<real_t>;
using boundary_functions = std::vector<
    std::function<real_t(real_t, real_t)>>;
using boundary_vectors = std::array<real_t, 4>;

void compute_boundary_solution(
    double                     **g, 
    std::vector<range_t>  const &ranges,
    boundary_functions    const  bf, 
    boundary_vectors      const  b);