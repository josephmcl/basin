#include "compute_g.h"

void compute_g(
    real_t                            **g, 
    std::vector<sparse_matrix_t>       &boundaries,
    real_t                             *solutions,
    real_t                             *sources, 
    vv<std::size_t>                    &boundary_type_map,
    vv<std::size_t>                    &boundary_data_map,
    components                         &sbp) {

    // compute ith: (boundary * solution ...) - H_tilde sources[i]
    //                       bm       bm         stat     m2 block 
    //    v         v        v        v         v         v      
    // (b_LB_W * g_LB_W + b_LB_S * g_LB_S) - H_tilde * F_LB[:]

    (*g) = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u, 64);
    memset(*g, 0, sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

    std::size_t n2 = sbp.n * sbp.n;
    std::size_t face_size = sbp.n * sbp.n_blocks_dim;
    real_t *gi, *solution, *source;
    sparse_status_t status;

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;

    sparse_matrix_t h;
    sbp.hl.mkl(&h);

    std::size_t relblock; 
    for (std::size_t block = 0; block != sbp.rank_limit_u; ++block) {
        
        relblock = sbp.rank_index_u[block];
        // block is the local index local -> 0, 1, 2, ... 
        // but relblock is element index --> k, k+1, k+2, ... 
        gi = &(*g)[0] + (n2 * block);

        // 4 is known quantity as our blocks are rect. and orth.
        constexpr std::size_t faces = 4;
        for (std::size_t face = 0; face != faces; ++face) {
            
            // b_type_index gives the types (direction) of a boundary
            auto b_type_index = boundary_type_map[relblock][face];
            if (b_type_index != 0) {
                
                // ti gives face + b_type_index 
                auto ti = face + ((b_type_index - 1) * faces);
                auto boundary = boundaries[ti];

                // di gives which index of boundary data on a given face 
                // [ 0][ 1][ 2] 
                // each index is length n
                // [..face 0..][..face 1..][..face 2..][..face 3..] 
                // each face is length n * n_block_dim
                // [...............boundary solution..............]
                auto di = boundary_data_map[relblock][face] - 1;
                solution = &solutions[0] + (face * face_size) + (di * sbp.n);

                
                status = mkl_sparse_d_mv(
                    SPARSE_OPERATION_NON_TRANSPOSE, 1., boundary, da,
                    solution, 
                    1., gi);
                mkl_sparse_status(status);
            }
        }

        source = &sources[0] + (n2 + block);
        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE, -1., h, da,
            source, 
            1., gi);
        mkl_sparse_status(status);
    }

}