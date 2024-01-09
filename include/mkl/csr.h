#pragma once 

// Helper class for generating arrays in CSR format. 

#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <iostream>

#include "mkl.h"
#include "mkl_spblas.h"

template<typename T>
struct csr {
    std::size_t n, m;
    std::vector<T> v;
    std::vector<MKL_INT> r, c;
    T *vv; MKL_INT *rr, *cc, *ee;

    std::size_t cols() const {
        return m;
    }
    std::size_t rows() const {
        return n;
    }
    std::size_t nnz() const { 
        return v.size();
    }
    std::size_t col_index_size() const { 
        return c.size();
    }
    std::size_t row_index_size() const {
        return r.size();
    }
    T *val_data() {
        return &v[0];
    }
    MKL_INT *col_index_data() const {
        return &c[0];
    }
    MKL_INT *row_index_data() const {
        return &r[0];
    }
    sparse_status_t mkl(
        sparse_matrix_t *mkls, 
        T alpha = T()) {

        for (auto &e: r)
            std::cout << e << " ";
        std::cout << std::endl; 
        std::cout << nnz() << std::endl; 

        
        T *data = val_data();
        std::vector<T> temp;
        if (alpha != T()) {
            temp = v;
            for (auto &e : temp) 
                e *= alpha;
            data = &temp[0];
        }

        vv = (T *) mkl_malloc(sizeof(T) * v.size(), 64);
        cc = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * c.size(), 64);
        rr = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * r.size(), 64);
        //ee = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * r.size(), 64);

        std::memcpy(vv, &data[0], sizeof(T) * v.size());
        std::memcpy(cc, &c[0], sizeof(MKL_INT) * c.size());
        std::memcpy(rr, &r[0], sizeof(MKL_INT) * r.size());
        //std::memcpy(ee, &r[1], sizeof(MKL_INT) * r.size());

        auto rv = mkl_sparse_d_create_csr(
            mkls, 
            SPARSE_INDEX_BASE_ZERO,
            n,
            m,
            &rr[0],
            &rr[1],
            &cc[0],
            &vv[0]);
        return rv;
    }

    csr &operator()(T value, std::size_t row, std::size_t column) {
        if (row >= n || row < 0) {
            throw std::invalid_argument("row out bounds");
        }
        else if (column >= m || column < 0) {
            throw std::invalid_argument("column out bounds");
        }
        else if (value == static_cast<T>(0)) {
            return *this;
        }
        MKL_INT rowi = r[row];
        MKL_INT coli = 0;
        for (; rowi < r[row + 1]; rowi++) {
			coli = c[rowi];
			if (static_cast<std::size_t>(coli) >= column) {
				break;
			}
		}

        // overwrite existing value
        if (static_cast<std::size_t>(coli) == column && nnz() > 0) { 
            v[rowi] = value;
        }
        else {
            v.insert(v.begin() + rowi, value);
            c.insert(c.begin() + rowi, column);
            for (std::size_t i = row + 1; i <= n; i++) {
                r[i] += 1;
            }
        }
        return *this;
    }
    csr(std::size_t n, std::size_t m) : n(n), m(m) {
        this->v = std::vector<T>();
		this->c = std::vector<MKL_INT>();
		this->r = std::vector<MKL_INT>(n + 1, 0);
    }
    csr<T> &operator=(const csr<T> &that) {
        this->n = that.n;
        this->m = that.m;
        this->v = that.v;
		this->c = that.c;
		this->r = that.r;
        return *this;
    }
    csr(): n(0), m(0) { }
    csr(csr const &that) { 
        *this = that;
    }
};
