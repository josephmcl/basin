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
void compute_u_sqr(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp);

void compute_u_dc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp);

void compute_u_rfpc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp);

/* Compile time wrapper */
template <solver s>
struct mrhs {
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
struct mrhs<solver::sparse_qr> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        vv<sparse_matrix_t> &M,
        real_t *g) {

        compute_u_sqr(Mg, M, g, sbp);
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
struct mrhs<solver::dense_cholesky> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        std::vector<double *> &M,
        real_t *g = nullptr) {
        compute_u_dc(Mg, M, sbp);
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
struct mrhs<solver::rfp_cholesky> {
    static inline void solve(
        components &sbp,
        real_t *Mg, 
        std::vector<double *> &M,
        real_t *g = nullptr) {
        compute_u_rfpc(Mg, M, sbp);
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