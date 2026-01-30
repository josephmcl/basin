#include "compute_m.h"

namespace {
constexpr std::size_t dirichlet = 1;
constexpr std::size_t neumann = 2;
}

void make_m(
    sparse_matrix_t *M,
    components &sbp,
    std::array<std::size_t, 4> const &boundary) {

    sparse_matrix_t Mb1, Mb2, Mb3, Mb4;
    sparse_matrix_t temp1, temp3, temp4, temp5, temp6;
    sparse_matrix_t dd2x, dd2y, hhl;

    matrix_descr da, db;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;

    sbp.d2x.mkl(&dd2x);
    sbp.d2y.mkl(&dd2y);
    sbp.hl.mkl(&hhl, -1);

    // M = -H^{-1}(D2x + D2y) + boundary terms
    auto status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, dd2x, 1., dd2y, &temp1);
    mkl_sparse_status(status);

    status = mkl_sparse_sp2m(
        SPARSE_OPERATION_NON_TRANSPOSE, da, hhl,
        SPARSE_OPERATION_NON_TRANSPOSE, db, temp1,
        SPARSE_STAGE_FULL_MULT, &temp3);
    mkl_sparse_status(status);

    // Compute boundary SAT terms for each direction
    make_m_boundary(&Mb1, sbp, 1, boundary[0]);
    make_m_boundary(&Mb2, sbp, 2, boundary[1]);
    make_m_boundary(&Mb3, sbp, 3, boundary[2]);
    make_m_boundary(&Mb4, sbp, 4, boundary[3]);

    // Combine boundary terms
    status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, Mb1, 1., Mb2, &temp4);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, Mb3, 1., Mb4, &temp5);
    mkl_sparse_status(status);

    status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, temp4, 1., temp5, &temp6);
    mkl_sparse_status(status);

    // Final M = interior + boundary
    status = mkl_sparse_d_add(
        SPARSE_OPERATION_NON_TRANSPOSE, temp3, 1., temp6, M);
    mkl_sparse_status(status);

    // Cleanup temporaries
    mkl_sparse_destroy(Mb1);
    mkl_sparse_destroy(Mb2);
    mkl_sparse_destroy(Mb3);
    mkl_sparse_destroy(Mb4);
    mkl_sparse_destroy(temp1);
    mkl_sparse_destroy(temp3);
    mkl_sparse_destroy(temp4);
    mkl_sparse_destroy(temp5);
    mkl_sparse_destroy(temp6);
    mkl_sparse_destroy(dd2x);
    mkl_sparse_destroy(dd2y);
    mkl_sparse_destroy(hhl);
}

void make_m_boundary(
    sparse_matrix_t *M,
    components const &sbp,
    std::size_t const direction,
    std::size_t const boundary) {

    // Select operators based on direction
    auto L = (direction == 1) ? sbp.lw
           : (direction == 2) ? sbp.le
           : (direction == 3) ? sbp.ls
           : sbp.ln;

    auto H = (direction < 3) ? sbp.hy : sbp.hx;
    auto BS = (direction < 3) ? sbp.bsx : sbp.bsy;

    matrix_descr da, db;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    db.type = SPARSE_MATRIX_TYPE_GENERAL;

    if (boundary == dirichlet) {
        // Dirichlet: tau*H*L'*L - beta*H*BS'*L'*L
        sparse_matrix_t th, bh, l, bs, temp1, temp2, temp3, temp4, temp5;

        H.mkl(&th, sbp.TAU_VALUE);
        H.mkl(&bh, -sbp.BETA_VALUE);
        L.mkl(&l);
        BS.mkl(&bs);

        // tau*H*L'*L
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, th,
            SPARSE_OPERATION_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp1);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp1,
            SPARSE_OPERATION_NON_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp2);

        // -beta*H*BS'*L'*L
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, bh,
            SPARSE_OPERATION_TRANSPOSE, db, bs,
            SPARSE_STAGE_FULL_MULT, &temp3);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp3,
            SPARSE_OPERATION_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp4);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp4,
            SPARSE_OPERATION_NON_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp5);

        mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, temp2, 1., temp5, M);

        mkl_sparse_destroy(th);
        mkl_sparse_destroy(bh);
        mkl_sparse_destroy(l);
        mkl_sparse_destroy(bs);
        mkl_sparse_destroy(temp1);
        mkl_sparse_destroy(temp2);
        mkl_sparse_destroy(temp3);
        mkl_sparse_destroy(temp4);
        mkl_sparse_destroy(temp5);
    } else {
        // Neumann: H*L'*L*BS - (1/tau)*H*BS'*L'*L*BS
        sparse_matrix_t h, th, l, bs, temp1, temp2, temp3, temp4, temp5, temp6, temp7;

        H.mkl(&h);
        H.mkl(&th, -1.0 / sbp.TAU_VALUE);
        L.mkl(&l);
        BS.mkl(&bs);

        // H*L'*L*BS
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, h,
            SPARSE_OPERATION_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp1);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp1,
            SPARSE_OPERATION_NON_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp2);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp2,
            SPARSE_OPERATION_NON_TRANSPOSE, db, bs,
            SPARSE_STAGE_FULL_MULT, &temp3);

        // -(1/tau)*H*BS'*L'*L*BS
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, th,
            SPARSE_OPERATION_TRANSPOSE, db, bs,
            SPARSE_STAGE_FULL_MULT, &temp4);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp4,
            SPARSE_OPERATION_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp5);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp5,
            SPARSE_OPERATION_NON_TRANSPOSE, db, l,
            SPARSE_STAGE_FULL_MULT, &temp6);
        mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, da, temp6,
            SPARSE_OPERATION_NON_TRANSPOSE, db, bs,
            SPARSE_STAGE_FULL_MULT, &temp7);

        mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, temp3, 1., temp7, M);

        mkl_sparse_destroy(h);
        mkl_sparse_destroy(th);
        mkl_sparse_destroy(l);
        mkl_sparse_destroy(bs);
        mkl_sparse_destroy(temp1);
        mkl_sparse_destroy(temp2);
        mkl_sparse_destroy(temp3);
        mkl_sparse_destroy(temp4);
        mkl_sparse_destroy(temp5);
        mkl_sparse_destroy(temp6);
        mkl_sparse_destroy(temp7);
    }
}