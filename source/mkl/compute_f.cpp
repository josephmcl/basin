#include "compute_f.h"

void compute_f(
  std::vector<sparse_matrix_t> &F_sparse,
  std::vector<real_t *>        &F_dense,
  components                   &sbp) {

   for (std::size_t index = 0; index != F_dense.size(); ++index) {
       F_dense[index] = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);
   }

  // (-τ * LN + β * LN* BS_y) * H_x 
  fcompop(&F_sparse[3], F_dense[3], sbp.ln, sbp.bsy, sbp.hx, sbp.τ, sbp.β, sbp.n);
  // (-τ * LS + β * LS* BS_x) * H_x 
  fcompop(&F_sparse[2], F_dense[2], sbp.ls, sbp.bsy, sbp.hx, sbp.τ, sbp.β, sbp.n);
  // (-τ * LE + β * LE* BS_y) * H_y 
  fcompop(&F_sparse[1], F_dense[1], sbp.le, sbp.bsx, sbp.hy, sbp.τ, sbp.β, sbp.n);
  // (-τ * LW + β * LW* BS_y) * H_y 
  fcompop(&F_sparse[0], F_dense[0], sbp.lw, sbp.bsx, sbp.hy, sbp.τ, sbp.β, sbp.n);

}

void fcompop(
  sparse_matrix_t   *f, 
  real_t            *f_dense, 
  csr<real_t>       &L, 
  csr<real_t>       &B,
  csr<real_t>       &H,
  real_t       const τ, 
  real_t       const β,
  std::size_t        n) {

    sparse_matrix_t tl, bl, b, h, temp1, temp2;

    matrix_descr da, db;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;

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
        temp2, h, SPARSE_LAYOUT_ROW_MAJOR, 
        f_dense, n * n);
    /*status = mkl_sparse_d_sp2md(
        SPARSE_OPERATION_TRANSPOSE, da, temp2,
        SPARSE_OPERATION_TRANSPOSE, db, h,
        1., 1., f_dense, SPARSE_LAYOUT_ROW_MAJOR, n);
    mkl_sparse_status(status);*/

    /*
    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vals;
    status = mkl_sparse_d_export_csr(
      *f, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vals);
    mkl_sparse_status(status);

    ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((16) + 1), 64);
    std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * 16);
    ia[16] = rowe[16 - 1];
    */
    
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

