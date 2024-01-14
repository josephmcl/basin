#include "compute_f.h"

void compute_f(
  std::vector<sparse_matrix_t> &F_sparse,
  std::vector<real_t *>        &F_dense,
  components                   &sbp) {

   for (std::size_t index = 0; index != F_dense.size(); ++index) {
       F_dense[index] = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);
   }

  // (-τ * LN + β * LN* BS_y) * H_x 
  fcompop(&F_sparse[3], F_dense[3], sbp.ln, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LS + β * LS* BS_x) * H_x 
  fcompop(&F_sparse[2], F_dense[2], sbp.ls, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LE + β * LE* BS_y) * H_y 
  fcompop(&F_sparse[1], F_dense[1], sbp.le, sbp.bsx, sbp.hy, sbp.τ, sbp.β);
  // (-τ * LW + β * LW* BS_y) * H_y 
  fcompop(&F_sparse[0], F_dense[0], sbp.lw, sbp.bsx, sbp.hy, sbp.τ, sbp.β);

}

void fcompop(
  sparse_matrix_t   *f, 
  real_t            *f_dense, 
  csr<real_t>       &L, 
  csr<real_t>       &B,
  csr<real_t>       &H,
  real_t       const τ, 
  real_t       const β) {

    sparse_matrix_t tl, bl, b, h, temp1, temp2;

    matrix_descr da, db;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    da.mode = SPARSE_FILL_MODE_UPPER;
	da.diag = SPARSE_DIAG_NON_UNIT;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.mode = SPARSE_FILL_MODE_UPPER;
	db.diag = SPARSE_DIAG_NON_UNIT;

    auto status = L.mkl(&tl, -τ);
    mkl_sparse_status(status);
    
    status = L.mkl(&bl, β);
    mkl_sparse_status(status);

    status = B.mkl(&b);
    mkl_sparse_status(status);

    status = H.mkl(&h);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, bl,
        SPARSE_OPERATION_NON_TRANSPOSE, db, b,
        SPARSE_STAGE_FULL_MULT, &temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
    SPARSE_OPERATION_NON_TRANSPOSE, tl, 1., 
        temp1, &temp2);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp2,
        SPARSE_OPERATION_NON_TRANSPOSE, db, h,
        SPARSE_STAGE_FULL_MULT, f);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp2,
        SPARSE_OPERATION_NON_TRANSPOSE, db, h,
        SPARSE_STAGE_FULL_MULT, f);
    mkl_sparse_status(status);

    // sparse_status_t mkl_sparse_d_spmmd (const sparse_operation_t operation, const sparse_matrix_t A, const sparse_matrix_t B, const sparse_layout_t layout, double *C, const MKL_INT ldc);
    status = mkl_sparse_d_spmmd(
        SPARSE_OPERATION_NON_TRANSPOSE, 
        temp2, h, SPARSE_LAYOUT_COLUMN_MAJOR, 
        f_dense, CblasRowMajor);
    /*status = mkl_sparse_sp2md(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp2,
        SPARSE_OPERATION_NON_TRANSPOSE, db, h,
        1., f_dense, CblasRowMajor);*/
    mkl_sparse_status(status);
    
    status = mkl_sparse_destroy(tl);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(bl);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(b);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(h);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_destroy(temp2);
    mkl_sparse_status(status);

    //status = mkl_sparse_destroy(temp3);
    //mkl_sparse_status(status);

    /*
    petsc_matrix t;
    MatMatMult(l, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &t);
    MatScale(t, β);
    MatAXPY(t, -τ, l, UNKNOWN_NONZERO_PATTERN);
    MatMatMult(t, h, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &f);
    finalize<fw>(f);
    destroy<fw>(t);
    */
}

