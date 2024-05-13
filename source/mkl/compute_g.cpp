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

    real_t *gtemp = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u, 64);
    memset(gtemp, 0, sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

    std::size_t n2 = sbp.n * sbp.n;
    std::size_t face_size = sbp.n * sbp.n_blocks_dim;
    real_t *gi, *gti, *solution, *source;
    sparse_status_t status;

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;

    sparse_matrix_t h;
    sbp.hl.mkl(&h);

    sparse_matrix_t tmat;
    csr<real_t> eye(sbp.n*sbp.n, sbp.n*sbp.n);
    for (std::size_t i = 0; i < sbp.n * sbp.n; ++i) {
      eye(1., i, i);
    }
    eye.mkl(&tmat);
    
    double *mtemp = (double *) mkl_malloc(sizeof(double) * sbp.n*sbp.n * sbp.n*sbp.n, 64);
    for (std::size_t i = 0; i < sbp.n * sbp.n * sbp.n * sbp.n; ++i) {
      mtemp[i] = 0.;
    }
    mkl_sparse_d_spmmd(SPARSE_OPERATION_NON_TRANSPOSE, h, tmat, SPARSE_LAYOUT_ROW_MAJOR, mtemp, sbp.n*sbp.n);
    for (std::size_t i = 0; i < sbp.n * sbp.n; ++i) {
      for (std::size_t j = 0; j < sbp.n * sbp.n; ++j) {
        std::cout << mtemp[i * sbp.n * sbp.n + j] << " ";
      }
      std::cout << std::endl;
    }

    std::size_t relblock; 
    for (std::size_t block = 0; block != sbp.rank_limit_u; ++block) {
        
        //std::cout << "block " << block << std::endl;
        relblock = sbp.rank_index_u[block];
        std::cout << "block " << relblock << std::endl;
        // block is the local index local -> 0, 1, 2, ... 
        // but relblock is element index --> k, k+1, k+2, ... 
        gi = &(*g)[0] + (n2 * relblock);
        gti = &(gtemp)[0] + (n2 * relblock);

        
        // 4 is known quantity as our blocks are rect. and orth.
        constexpr std::size_t faces = 4;
        for (std::size_t face = 0; face != faces; ++face) {
            
            // b_type_index gives the types (direction) of a boundary
            auto b_type_index = boundary_type_map[relblock][face];
            if (b_type_index != 0) {
                
                // std::cout << "face - " << face + ((b_type_index - 1) * faces) << std::endl;
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
                // gi += boundary:matrix * solution:vector
                status = mkl_sparse_d_mv(
                    SPARSE_OPERATION_NON_TRANSPOSE, 1., 
                    boundary, da,
                    solution, 
                    1., gi);
                mkl_sparse_status(status);
            }
        }
        

        source = &sources[0] + (sbp.n * sbp.n * relblock);

        // std::size_t source_offset = relblock * sbp.n * sbp.n;

        for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
            gti[i] = -source[i];
        }

        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE, 1., 
            h, da,
            gti, 
            1., gi);
        mkl_sparse_status(status);

    }

}