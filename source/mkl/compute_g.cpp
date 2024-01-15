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

    /// petsc_vector temp2, temp3, temp5, temp6;
    //auto v = std::vector<double>(sbp.n * sbp.n);
    //auto mi = std::vector<int>(sbp.n * sbp.n);
    //for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
    //    mi[i] = i;
    //} 

    /*
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp2);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp3);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp5);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp6);
    g.resize(sbp.n_blocks);
    */

    (*g) = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks, 64);
    memset(*g, 0, sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks);

    std::size_t n2 = sbp.n * sbp.n;
    std::size_t face_size = sbp.n * sbp.n_blocks_dim;
    real_t *gi, *solution, *source;
    sparse_status_t status;

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    da.mode = SPARSE_FILL_MODE_UPPER;
	da.diag = SPARSE_DIAG_NON_UNIT;

    sparse_matrix_t h;
    sbp.hl.mkl(&h);

    for (std::size_t block = 0; block != sbp.n_blocks; ++block) {

        // VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &g[block]);

        gi = &(*g)[0] + (n2 * block);

        // 4 is known quantity as our blocks are rect. and orth.
        constexpr std::size_t faces = 4;
        for (std::size_t face = 0; face != faces; ++face) {
            
            // b_type_index gives the types (direction) of a boundary
            auto b_type_index = boundary_type_map[block][face];
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
                auto di = boundary_data_map[block][face] - 1;
                solution = &solutions[0] + (face * face_size) + (di * sbp.n);

                //sparse_status_t mkl_sparse_d_mv (const sparse_operation_t operation, const double alpha, const sparse_matrix_t A, const struct matrix_descr descr, const double *x, const double beta, double *y);
                
                status = mkl_sparse_d_mv(
                    SPARSE_OPERATION_NON_TRANSPOSE, 1., boundary, da,
                    solution, 
                    1., gi);
                mkl_sparse_status(status);

                //MatMult(boundary, solution, temp2);
                //VecAXPY(g[block], 1., temp2); // g[block] += temp2 * 1;
            }
        }

        source = &sources[0] + (n2 + block);
        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE, -1., h, da,
            source, 
            1., gi);
        mkl_sparse_status(status);

        // Get the raw dense matrix data, effectively reshape to a vector
        /*
        MatGetValues(
            sources[block], sbp.n, &mi[0], sbp.n, &mi[0], &v[0]);
        VecSetValues(temp5, sbp.n * sbp.n, &mi[0], &v[0], INSERT_VALUES);
        finalize<fw>(temp5);
        MatMult(sbp.hl, temp5, temp6);
        VecAXPY(g[block], -1., temp6); // g[block] -= temp5 
        finalize<fw>(g[block]);
        */
        // VecView(g[block], PETSC_VIEWER_STDOUT_SELF);
    }

}