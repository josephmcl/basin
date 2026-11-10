
#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include <cuda_runtime.h>



static __global__ void vector_csr_kernel(const int m,
                                const int *__restrict__ row_offsets,
                                const int *__restrict__ cols,
                                const double *__restrict__ vals,
                                const double *__restrict__ x,
                                double *__restrict__ y,
                                const double alpha,
                                const double beta)