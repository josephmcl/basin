#include "compute_b.h"


void compute_b(
    std::vector<sparse_matrix_t>    &  B, 
    components                      &sbp) {

    B.resize(8); // 4 for first-order 4 for second-order. Some may not be
                 // needed but it's easier to just compute them all for 
                 // now.

    compute_b1(B[0], sbp.hy, sbp.TAU_VALUE, sbp.lw, sbp.BETA_VALUE, sbp.bsx, sbp.n);
    compute_b1(B[1], sbp.hy, sbp.TAU_VALUE, sbp.le, sbp.BETA_VALUE, sbp.bsx, sbp.n);
    compute_b1(B[2], sbp.hx, sbp.TAU_VALUE, sbp.ls, sbp.BETA_VALUE, sbp.bsy, sbp.n);
    compute_b1(B[3], sbp.hx, sbp.TAU_VALUE, sbp.ln, sbp.BETA_VALUE, sbp.bsy, sbp.n);

    compute_b2(B[4], sbp.hy, sbp.TAU_VALUE, sbp.lw, sbp.BETA_VALUE, sbp.bsx, sbp.n);
    compute_b2(B[5], sbp.hy, sbp.TAU_VALUE, sbp.le, sbp.BETA_VALUE, sbp.bsx, sbp.n);
    compute_b2(B[6], sbp.hx, sbp.TAU_VALUE, sbp.ls, sbp.BETA_VALUE, sbp.bsy, sbp.n);
    compute_b2(B[7], sbp.hx, sbp.TAU_VALUE, sbp.ln, sbp.BETA_VALUE, sbp.bsy, sbp.n);
}

/* Compute H(a) * (TAU_VALUE * L(d)' - BETA_VALUE * BS(a)' * L(d)') for particular 
   directional orientations of L and axes orientations of BS. For 2-D 
   EWNS this includes 

    East:  H_y * (TAU_VALUE * LE' - BETA_VALUE * BS_x' * LE')
    West:  H_y * (TAU_VALUE * LW' - BETA_VALUE * BS_x' * LW')
    South: H_x * (TAU_VALUE * LS' - BETA_VALUE * BS_y' * LS')
    North: H_x * (TAU_VALUE * LN' - BETA_VALUE * BS_y' * LN')  
*/
void compute_b1(
    sparse_matrix_t    & B, 
    csr<real_t>        & H, 
    real_t       const   TAU_VALUE, 
    csr<real_t>        & L, 
    real_t       const   BETA_VALUE, 
    csr<real_t>        &BS,
    std::size_t  const   n) {
    
    sparse_matrix_t bs, l, temp1, temp2, lb, h;

    auto status = BS.mkl(&bs, BETA_VALUE);
    mkl_sparse_status(status);

    status = L.mkl(&l);
    mkl_sparse_status(status);

    matrix_descr da, db;

    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_TRANSPOSE, da, bs,
        SPARSE_OPERATION_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp1);
    mkl_sparse_status(status);

    status = L.mkl(&lb, TAU_VALUE);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
        SPARSE_OPERATION_TRANSPOSE, lb, -1., 
        temp1, &temp2);
    mkl_sparse_status(status);

    status = H.mkl(&h, -1);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, h,
        SPARSE_OPERATION_NON_TRANSPOSE, db, temp2,
        SPARSE_STAGE_FULL_MULT, &B);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(bs);
    mkl_sparse_status(status);
    
    status = mkl_sparse_destroy(l);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp2);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(lb);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(h);
    mkl_sparse_status(status);
    

    /*
    petsc_matrix LT, BST, temp1, temp2, temp3;
    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // temp1 := BETA_VALUE * BS_a^T
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n * n, n, nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, BETA_VALUE, BST, UNKNOWN_NONZERO_PATTERN);    

    // temp2 := temp1 * L_d^T 
    MatMatMult(temp1, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);   

    // temp3 := TAU_VALUE * L_d^T 
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n, n, nullptr, &temp3);
    finalize<fw>(temp3);    
    MatAXPY(temp3, TAU_VALUE, LT, UNKNOWN_NONZERO_PATTERN);

    // temp3 := temp3 + (-1) * temp2 
    MatAXPY(temp3, -1., temp2, UNKNOWN_NONZERO_PATTERN);

    // B := H * temp3
    MatMatMult(H, temp3, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);
    finalize<fw>(B);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
    */
}


/* Compute H(a) * (BETA_VALUE * L(d)' - 1/TAU_VALUE * BS(a')' * L(d)') for particular 
   directional of L and axes orientations of H and BS. For 2-D EWNS this 
   includes 

    South: H_x * (BETA_VALUE * LS' - 1/TAU_VALUE * BS_y' * LS')
    North: H_x * (BETA_VALUE * LN' - 1/TAU_VALUE * BS_y' * LN')  */
void compute_b2(
    sparse_matrix_t    & B, 
    csr<real_t>        & H, 
    real_t       const   TAU_VALUE, 
    csr<real_t>        & L, 
    real_t       const   BETA_VALUE, 
    csr<real_t>        &BS,
    std::size_t  const   n) {
    
    sparse_matrix_t bs, l, temp1, temp2, lb, h;

    auto status = BS.mkl(&bs, 1. / TAU_VALUE);
    mkl_sparse_status(status);

    status = L.mkl(&l);
    mkl_sparse_status(status);

    matrix_descr da, db;

    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_TRANSPOSE, da, bs,
        SPARSE_OPERATION_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp1);
    mkl_sparse_status(status);

    status = L.mkl(&lb, BETA_VALUE);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
        SPARSE_OPERATION_TRANSPOSE, lb, -1., 
        temp1, &temp2);
    mkl_sparse_status(status);

    status = H.mkl(&h, -1);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, h,
        SPARSE_OPERATION_NON_TRANSPOSE, db, temp2,
        SPARSE_STAGE_FULL_MULT, &B);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(bs);
    mkl_sparse_status(status);
    
    status = mkl_sparse_destroy(l);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp2);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(lb);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(h);
    mkl_sparse_status(status);

    /*
    sparse_index_base_t i;
    MKL_INT r, c, *rs, *re, *ci;
    real_t *v;
    rv = mkl_sparse_d_export_csr(bs, &i, &r, &c, &rs, &re, &ci, &v);
    std::cout << ((rv == SPARSE_STATUS_SUCCESS)? ":)": ">:(") << std::endl;

    std::cout << r << ", " << c << std::endl; 
    int o = 0;
    MKL_INT *rt = rs;
    while (rt != re) {
        std::cout << *rt <<  " ";
        rt += 1;
        o += 1;
    }
    std::cout << std::endl;
    
    std::cout << rs[r] << std::endl;
        for (int i = 0; i < rs[r - 1]; ++i) {
        std::cout << *(ci + i) <<  " ";
    }
    std::cout << std::endl;
    
    for (int i = 0; i < rs[r - 1]; ++i) {
        std::cout << *(v + i) <<  " ";
    }
    std::cout << std::endl;*/

    /*
    petsc_matrix LT, BST, temp1, temp2, temp3;
    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // temp1 := 1. / TAU_VALUE * BS(a')^T
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n * n, n, nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, 1. / TAU_VALUE, BST, UNKNOWN_NONZERO_PATTERN);    

    // temp2 := temp1 * L_d^T 
    MatMatMult(temp1, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);   

    // temp3 := BETA_VALUE * L_d^T 
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n, n, nullptr, &temp3);
    finalize<fw>(temp3);    
    MatAXPY(temp3, BETA_VALUE, LT, UNKNOWN_NONZERO_PATTERN);

    // temp3 := temp3 + (-1) * temp2 
    MatAXPY(temp3, -1., temp2, UNKNOWN_NONZERO_PATTERN);

    // B := H * temp3
    MatMatMult(H, temp3, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);
    finalize<fw>(B);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
    */
}