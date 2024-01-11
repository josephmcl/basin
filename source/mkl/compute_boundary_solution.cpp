#include "compute_boundary_solution.h"

// NOTE: Assumes number of grid points in x and y dims are equal. 
void compute_boundary_solution(
    double                     **g, 
    std::vector<range_t>  const &ranges,
    boundary_functions    const  bf, 
    boundary_vectors      const  b) {

    std::size_t size = 0;
    for (auto &e: ranges) 
        size += e.size();
    std::size_t I = size;
    std::size_t J = ranges[0].size();
    size *= 4;

    (*g) = (double *) mkl_malloc(sizeof(double) * size, 64);

    for (std::size_t i = 0; i != 4; ++i) {
        for (std::size_t j = 0; j != ranges.size(); ++j) {
            auto range = ranges[j];            
            for (auto e = range.begin(); e != range.end(); ++e) {
                (*g)[(I * i) + (J * j) + e.index] = (i < 2) 
                    ? bf[i](b[i], *e) 
                    : bf[i](*e, b[i]);
            }
        }
    }
}