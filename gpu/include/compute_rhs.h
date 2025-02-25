#pragma once

#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

void compute_rhs(
    real_t             *rhs,
    real_t             *g, 
    std::vector<csr_t> &F,
    real_t             *λ,
    vv<std::size_t>    &F_symbols,
    components         &sbp);