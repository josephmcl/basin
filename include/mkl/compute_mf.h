#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"
#include "cholesky.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"

#include <cstring>

#include "omp.h"

using real_t = type::real_t;

void compute_mf_sqr(
  std::vector<real_t *>        &x,
  vv<sparse_matrix_t>          &m,
  std::vector<real_t *>        &f,
  components             const &sbp);

void compute_mf_dc(
  std::vector<real_t *>        &x,
  std::vector<real_t *>        &m,
  components             &sbp);

void compute_mf_rfpc(
  std::vector<real_t *>        &x,
  std::vector<real_t *>        &m,
  components             &sbp);

template <solver s>
struct mf {
    template <typename v, typename m>
    static inline void solve(
        components &sbp, v Mg, m M, v g) {
        return;
    }
};

template<>
struct mf<solver::sparse_qr> {
  static inline void solve(
    components            &sbp,
    std::vector<real_t *> &Mf,
    vv<sparse_matrix_t>   &M,
    std::vector<real_t *> &f) {
    compute_mf_sqr(Mf, M, f, sbp);
    return;
  }
};

static std::vector<double *> _f_dummy;

template<>
struct mf<solver::dense_cholesky> {
  static inline void solve(
    components            &sbp,
    std::vector<real_t *> &Mf,
    std::vector<real_t *> &M,
    std::vector<real_t *> &f = _f_dummy) {
    compute_mf_dc(Mf, M, sbp);
    return;
  }
};

template<>
struct mf<solver::rfp_cholesky> {
  static inline void solve(
    components &sbp,
    std::vector<real_t *> &Mf,
    std::vector<real_t *> &M,
    std::vector<real_t *> &f = _f_dummy) {
      compute_mf_rfpc(Mf, M, sbp);
      return;
  }
};