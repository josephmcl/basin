#include "compute_boundary_solution.h"
#include <hip/hip_runtime.h>
#include <vector>
#include <iostream>

// HIP Kernel for computing boundary solutions
__global__ void compute_boundary_solution_kernel(
    double *g, 
    range_t *ranges, 
    boundary_functions bf, 
    boundary_vectors b,
    std::size_t I, std::size_t J, std::size_t num_ranges) 
{
    std::size_t i = blockIdx.x;  // Boundary index (0 to 3)
    std::size_t j = blockIdx.y;  // Range index

    if (j >= num_ranges) return;

    range_t range = ranges[j];

    for (auto e = range.begin(); e != range.end(); ++e) {
        std::size_t index = (I * i) + (J * j) + e.index;
        g[index] = (i < 2) 
            ? bf[i](b[i], *e)  // E/W boundary
            : bf[i](*e, b[i]);  // N/S boundary
    }
}

// HIP Version of `compute_boundary_solution`
void compute_boundary_solution(
    double **g, 
    std::vector<range_t> const &ranges,
    boundary_functions const bf, 
    boundary_vectors const b) 
{
    std::size_t size = 0;
    for (auto &e: ranges) 
        size += e.size();
    
    std::size_t I = size;
    std::size_t J = ranges[0].size();
    size *= 4;

    // Allocate device memory for g
    hipMalloc((void**)g, sizeof(double) * size);

    // Allocate device memory for ranges (copy from host)
    range_t *d_ranges;
    hipMalloc((void**)&d_ranges, sizeof(range_t) * ranges.size());
    hipMemcpy(d_ranges, ranges.data(), sizeof(range_t) * ranges.size(), hipMemcpyHostToDevice);

    // Launch HIP kernel
    dim3 gridDim(4, ranges.size());  // (4 boundaries, num_ranges)
    hipLaunchKernelGGL(compute_boundary_solution_kernel, gridDim, 1, 0, 0, *g, d_ranges, bf, b, I, J, ranges.size());

    // Free device memory for ranges
    hipFree(d_ranges);
}
