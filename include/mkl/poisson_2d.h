
#include "measure.h"
#include "components.h"
#include "compute_lambda_matrix.h"
#include "definitions.h"
#include "connect.h"
#include "compute_b.h"

#include "omp.h"

#include <iomanip>
#include <string>

namespace poisson_2d {

    void problem(std::size_t vln, std::size_t eln);

} /* namespace poisson_2d */