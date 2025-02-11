#include <rocsparse.h>
#include <hip/hip_runtime.h>
#include <vector>
#include <iostream>

typedef double real_t;

void compute_g(
    real_t** g, 
    std::vector<rocsparse_mat_descr> &boundaries,  // rocsparse descriptors for boundary matrices
    real_t* solutions,  // Device pointer
    real_t* sources,    // Device pointer
    std::vector<std::vector<std::size_t>> &boundary_type_map,
    std::vector<std::vector<std::size_t>> &boundary_data_map,
    components &sbp)  
{
    // 1) Create a ROCm rocsparse handle
    rocsparse_handle handle;
    rocsparse_create_handle(&handle);

    // 2) Allocate memory for g and gtemp on the device
    std::size_t total_size = sbp.n * sbp.n * sbp.rank_limit_u;
    hipMalloc((void**)g, total_size * sizeof(real_t));
    hipMemset(*g, 0, total_size * sizeof(real_t));

    real_t* gtemp;
    hipMalloc((void**)&gtemp, total_size * sizeof(real_t));
    hipMemset(gtemp, 0, total_size * sizeof(real_t));

    // 3) Get the H matrix (assumed already in CSR format)
    rocsparse_mat_descr h_descr;
    rocsparse_create_mat_descr(&h_descr);
    rocsparse_set_mat_type(h_descr, rocsparse_matrix_type_general);
    rocsparse_set_mat_index_base(h_descr, rocsparse_index_base_zero);

    // Assume sbp.hl contains the CSR matrix for H.
    rocsparse_int* h_row_ptr = sbp.hl.row_ptr;  // Device pointer
    rocsparse_int* h_col_ind = sbp.hl.col_ind;  // Device pointer
    real_t* h_values = sbp.hl.values;           // Device pointer

    std::size_t n2 = sbp.n * sbp.n;
    std::size_t face_size = sbp.n * sbp.n_blocks_dim;

    for (std::size_t block = 0; block < sbp.rank_limit_u; ++block) {
        std::size_t relblock = sbp.rank_index_u[block];

        // Compute offsets
        real_t* gi_dev  = (*g) + (n2 * relblock);
        real_t* gti_dev = gtemp + (n2 * relblock);

        constexpr std::size_t faces = 4;
        for (std::size_t face = 0; face < faces; ++face) {
            auto b_type_index = boundary_type_map[relblock][face];
            if (b_type_index != 0) {
                auto ti = face + ((b_type_index - 1) * faces);
                auto boundary = boundaries[ti];

                auto di = boundary_data_map[relblock][face] - 1;
                real_t* solution_dev = solutions + (face * face_size) + (di * sbp.n);

                // Perform sparse matrix-vector multiplication: gi += boundary * solution
                real_t alpha = 1.0;
                real_t beta = 1.0;
                rocsparse_dcsrmv(
                    handle,
                    rocsparse_operation_none,
                    sbp.n, sbp.n, sbp.n,  // Assuming square blocks
                    &al
                    pha,
                    boundary,
                    solution_dev,
                    &beta,
                    gi_dev
                );
            }
        }

        // Negate source values using HIP kernel
        hipLaunchKernelGGL(negate_kernel, (n2 + 255) / 256, 256, 0, 0, gti_dev, sources + (n2 * relblock), n2);

        // Perform sparse matrix-vector multiplication: gi += H * gti
        real_t alpha = 1.0;
        real_t beta = 1.0;
        rocsparse_dcsrmv(
            handle,
            rocsparse_operation_none,
            sbp.n, sbp.n, sbp.n,  // Assuming square blocks
            &alpha,
            h_descr,
            gti_dev,
            &beta,
            gi_dev
        );
    }

    // Destroy descriptors and free temporary memory
    rocsparse_destroy_mat_descr(h_descr);
    rocsparse_destroy_handle(handle);
    hipFree(gtemp);
}

// HIP Kernel for negating an array
__global__ void negate_kernel(real_t* output, real_t* input, std::size_t size) {
    std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        output[idx] = -input[idx];
    }
}
