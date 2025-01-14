#pragma once

/* Local headers */
#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"
#include "cholesky.h"

/* External headers */
#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"
#include "omp.h"

/* Implementations */
void compute_mg_sqr(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp);

void compute_mg_dc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp);

void compute_mg_rfpc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp);

void compute_mg_mpi(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp);

template <solver s>
struct mg {
    template <typename v, typename m>
    static inline void solve(
        components &sbp, v Mg, m M, v g) {
        return;
    }
    static std::size_t bytes(components &sbp) {
        return sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static std::size_t size(components &sbp) {
        return sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static double *allocate(components &sbp) {
        return static_cast<double *>(
            mkl_malloc(bytes(sbp), 64));
    }
};

template<>
struct mg<solver::sparse_qr> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        vv<sparse_matrix_t> &M,
        real_t *g) {

        compute_mg_sqr(Mg, M, g, sbp);
        return;
    }
    static std::size_t bytes(components &sbp) {
        return sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static std::size_t size(components &sbp) {
        return sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static double *allocate(components &sbp) {
        return static_cast<double *>(
            mkl_malloc(bytes(sbp), 64));
    }
};

template<>
struct mg<solver::dense_cholesky> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        std::vector<double *> &M,
        real_t *g = nullptr) {
        compute_mg_dc(Mg, M, sbp);
        return;
    }
    static std::size_t bytes(components &sbp) {
        return sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static std::size_t size(components &sbp) {
        return sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static double *allocate(components &sbp) {
        return static_cast<double *>(
            mkl_malloc(bytes(sbp), 64));
    }
};

template<>
struct mg<solver::rfp_cholesky> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        std::vector<double *> &M,
        real_t *g = nullptr) {
        compute_mg_rfpc(Mg, M, sbp);
        return;
    }
    static std::size_t bytes(components &sbp) {
        return sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static std::size_t size(components &sbp) {
        return sbp.n * sbp.n * sbp.rank_limit_u;
    }

    static double *allocate(components &sbp) {
        return static_cast<double *>(
            mkl_malloc(bytes(sbp), 64));
    }
};

