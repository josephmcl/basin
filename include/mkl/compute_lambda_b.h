#pragma once

#include <vector>
#include <tuple>

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"

void compute_lambda_b(
    real_t *λb, 
    std::vector<sparse_matrix_t> &Fsparse, 
    real_t *Mg, 
    vv<std::size_t> &FT_symbols, 
    components &sbp);