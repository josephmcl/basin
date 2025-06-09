#pragma once

#include <vector>
#include <tuple>

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "rocblas.h"

void compute_lambda_b(
    real_t             *LAMBDAb, 
    std::vector<csr_t> &Fsparse, 
    real_t             *Mg, 
    vv<std::size_t>    &FT_symbols, 
    components         &sbp);

    void compute_lambda_b_csr(
        real_t     *λb, 
        csr_t      *F, 
        real_t     *Mg, 
        components &sbp);