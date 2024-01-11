#include "compute_sources.h"

// Given a list of n ranges, generate n^2 matrices of source data: 
// (0, 0), (0, 1), ... (1, 0), (1, 1), ... (n-1, 0) ... (n-1, n-1)  
// assumed that f(x, y) = source_x(x, y) + source_y(x, y)
void compute_sources(
    real_t                                      **F, 
    std::vector<range_t>                  const  &grids,
    std::function<real_t(real_t, real_t)> const   f) {

    (*F) = (real_t *) mkl_malloc(sizeof(real_t) * grids.size() * grids.size() * grids[0].size() * grids[0].size(), 64);
    std::size_t i = 0;
    std::size_t v;
    for (auto &x : grids) { for (auto &y : grids) {
        //MatCreateSeqDense(PETSC_COMM_SELF, x.size(), y.size(), NULL, 
        //    &F[i]);
        
        for (auto xi = x.begin(); xi != x.end(); ++xi) { 
            for (auto yi = y.begin(); yi != y.end(); ++yi) {
                v = ((x.size() + y.size()) * i) + (y.size() * xi.index) + yi.index;
                (*F)[v] = f(*xi, *yi); 
                //MatSetValue(F[i], xi.index, yi.index, f(*xi, *yi), 
                //    INSERT_VALUES);
            }
        }
        // finalize<fw>(F[i]);
        ++i;
    }}
}