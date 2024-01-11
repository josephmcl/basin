void sbp_sat::x2::make_m(
  sparse_matrix_t                  &M,
  components                       &sbp, 
  std::array<std::size_t, 4> const &boundary) {

  petsc_matrix temp1, temp2; 
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, nullptr, 
    &temp1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, nullptr, 
    &temp2);
  finalize<fw>(temp1);
  finalize<fw>(temp2);

  MatAXPY(temp1, 1., sbp.d2x, UNKNOWN_NONZERO_PATTERN);
  MatAXPY(temp1, 1., sbp.d2y, UNKNOWN_NONZERO_PATTERN);
  MatAXPY(temp2, -1., sbp.hl, UNKNOWN_NONZERO_PATTERN);

  MatMatMult(temp2, temp1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &M);

  make_M_boundary(M, sbp, 1, boundary[0]);
  make_M_boundary(M, sbp, 2, boundary[1]);
  make_M_boundary(M, sbp, 3, boundary[2]);
  make_M_boundary(M, sbp, 4, boundary[3]);

  finalize<fw>(M);
}

void sbp_sat::x2::make_M_boundary(
  petsc_matrix &M,
  components const &sbp, 
  std::size_t const direction,
  std::size_t const boundary) {

  auto L = direction == 1 ? sbp.lw 
         : direction == 2 ? sbp.le
         : direction == 3 ? sbp.ls 
         : sbp.ln; 

  auto H = direction < 3 ? sbp.hy : sbp.hx;
  auto BS = direction < 3 ? sbp.bsx : sbp.bsy;

  if (boundary == dirichlet) {

    sparse_matrix_t th, bh, l, bs, temp1, temp2, temp3, temp4, temp5;

    // MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    // MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    // finalize<fw>(LT);
    // finalize<fw>(BST);

    status = H.mkl(&th, sbp.τ);
    mkl_sparse_status(status);

    status = H.mkl(&bh, sbp.β);
    mkl_sparse_status(status);

    status = L.mkl(&l);
    mkl_sparse_status(status);

    status = BS.mkl(&bs);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, th,
        SPARSE_OPERATION_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp1,
        SPARSE_OPERATION_NON_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp2);
    mkl_sparse_status(status);

    // τ*H_y*LW'*LW 
    /*
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, sbp.τ, H, UNKNOWN_NONZERO_PATTERN);    
    MatMatMult(LT, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);    
    MatMatMult(temp1, temp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp3);
    finalize<fw>(temp3);
    MatAXPY(M, 1., temp3, UNKNOWN_NONZERO_PATTERN);
    */

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, bh,
        SPARSE_OPERATION_TRANSPOSE, db, bs,
        SPARSE_STAGE_FULL_MULT, &temp3);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp3,
        SPARSE_OPERATION_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp4);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, temp4,
        SPARSE_OPERATION_TRANSPOSE, db, l,
        SPARSE_STAGE_FULL_MULT, &temp5);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, temp2, 1., 
        temp5, &M);
    mkl_sparse_status(status);

    // -β*H_y*BS_x'*LW'*LW 
    /*
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp4);
    finalize<fw>(temp4);
    MatAXPY(temp4, sbp.β, H, UNKNOWN_NONZERO_PATTERN);
    MatMatMult(temp4, BST, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp5);
    finalize<fw>(temp5);
    MatMatMult(temp5, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp6);
    finalize<fw>(temp6);
    MatMatMult(temp6, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp7);
    finalize<fw>(temp7);
    MatAXPY(M, -1., temp7, UNKNOWN_NONZERO_PATTERN);    
    */  
    // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(temp7, PETSC_VIEWER_STDOUT_SELF);

    status = mkl_sparse_destroy(th);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(bh);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(l);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(bs);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(temp1);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(temp2);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(temp3);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(temp4);
    mkl_sparse_status(status);
    status = mkl_sparse_destroy(temp5);
    mkl_sparse_status(status);
  }
  else {

    sparse_matrix_t h, th, l, bs, temp1, temp2, temp3, temp4, temp5;

    status = H.mkl(&h);
    mkl_sparse_status(status);

    status = H.mkl(&th, 1 / sbp.τ);
    mkl_sparse_status(status);

    status = L.mkl(&l);
    mkl_sparse_status(status);

    status = BS.mkl(&bs);
    mkl_sparse_status(status);

    // TODO

    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // H_x * LN' * LN * BS_y
    MatMatMult(H, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp1);
    finalize<fw>(temp1);
    MatMatMult(temp1, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);
    MatMatMult(temp2, BS, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp3);
    finalize<fw>(temp3);
    MatAXPY(M, 1., temp3, UNKNOWN_NONZERO_PATTERN);

    // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(temp3, PETSC_VIEWER_STDOUT_SELF);

    // -1/τ * H_x * BS_y' * LN' * LN * BS_y 
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp4);
    finalize<fw>(temp4);
    MatAXPY(temp4, 1. / sbp.τ, H, UNKNOWN_NONZERO_PATTERN);
    MatMatMult(temp4, BST, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp5);
    finalize<fw>(temp5);
    MatMatMult(temp5, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp6);
    finalize<fw>(temp6);
    MatMatMult(temp6, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp7);
    finalize<fw>(temp7);
    MatMatMult(temp7, BS, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp8);
    finalize<fw>(temp8);
    MatAXPY(M, -1., temp8, UNKNOWN_NONZERO_PATTERN);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
    destroy<fw>(temp4);
    destroy<fw>(temp5);
    destroy<fw>(temp6);
    destroy<fw>(temp7);
    destroy<fw>(temp8);

  }
}