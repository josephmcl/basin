#pragma once

#include "measure.h"
#include "components.h"
#include "compute_lambda_a.h"
#include "definitions.h"
#include "connect.h"
#include "compute_b.h"
#include "compute_boundary_solution.h"
#include "compute_sources.h"
#include "compute_g.h"
#include "compute_f_symbols.h"
#include "compute_f.h"
#include "compute_m.h"
#include "compute_mf.h"
#include "compute_d.h"
#include "compute_mg.h"
#include "compute_lambda_b.h"
#include "compute_lambda.h"
#include "compute_rhs.h"
#include "compute_u.h"


#include "omp.h"
#include "mkl.h"
#include "mkl_spblas.h"

#include <iomanip>
#include <string>
#include <chrono>
#include <thread>

namespace poisson_2d {

    void problem(std::size_t vln, std::size_t eln);

} /* namespace poisson_2d */