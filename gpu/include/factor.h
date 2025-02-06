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


template <solver s>
struct factor {
    template <typename m>
    static inline real_t *ize(
        components &sbp, m M) {
        return nullptr;
    }
};

template<>
struct factor<solver::sparse_qr> {
    static inline real_t *ize(
        components &sbp,
        sparse_matrix_t &M) {

        matrix_descr dc;
        dc.type = SPARSE_MATRIX_TYPE_GENERAL;
        auto status = mkl_sparse_qr_reorder(M, dc);
        mkl_sparse_status(status);
        status = mkl_sparse_d_qr_factorize(M, nullptr);
        mkl_sparse_status(status);

        return nullptr;
    }
};

template<>
struct factor<solver::dense_cholesky> {
    static inline real_t *ize(
        components &sbp,
        sparse_matrix_t &M) {

        real_t *res;
        to_dense(res, M, sbp.n * sbp.n);

        auto error = cholesky::factor(res, sbp);
        if (error) {
            std::cout << "Cholesky factor failed with code " 
                << *error << std::endl; 
        }    
        
        return res;
    }
    static inline real_t *ize(
        components &sbp,
        real_t *M) {

        real_t *res = M;

        auto error = cholesky::factor(res, sbp);
        if (error) {
            std::cout << "Cholesky factor failed with code " 
                << *error << std::endl; 
        }    
        
        return res;
    }
};

template<>
struct factor<solver::rfp_cholesky> {
    static inline real_t *ize(
        components &sbp,
        sparse_matrix_t &M) {
        
        real_t *res;
        to_rfp(res, M, sbp.n * sbp.n);
    
        auto error = cholesky::factor_rfp(res, sbp);
        if (error) {
            std::cout << "Cholesky RFP factor failed with code " 
            << *error << std::endl; 
        }

        return res;
    }
};

