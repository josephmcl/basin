#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"

#include "rocsparse.h"

using real_t = type::real_t;

void compute_g(
    real_t                            **g, 
    std::vector<csr_t>                 &boundaries,
    real_t                             *solutions,
    real_t                             *sources, 
    vv<std::size_t>                    &boundary_type_map,
    vv<std::size_t>                    &boundary_data_map,
    components                         &sbp);