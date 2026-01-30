// std
#include <functional>
// you
#include "ranges.h"
#include "definitions.h"
// lib
#include "mkl.h"

using real_t = type::real_t;
using range_t = numerics::linrange<real_t>;

void compute_sources(
    real_t                                      **F, 
    std::vector<range_t>                  const  &grids,
    std::function<real_t(real_t, real_t)> const   f);