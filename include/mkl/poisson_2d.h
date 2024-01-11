#pragma once

#include "measure.h"
#include "components.h"
#include "compute_lambda_matrix.h"
#include "definitions.h"
#include "connect.h"
#include "compute_b.h"
#include "compute_boundary_solution.h"
#include "compute_sources.h"
#include "compute_g.h"

#include "omp.h"
#include "mkl.h"
#include "mkl_spblas.h"

#include <iomanip>
#include <string>

namespace poisson_2d {

    void problem(std::size_t vln, std::size_t eln);

} /* namespace poisson_2d */