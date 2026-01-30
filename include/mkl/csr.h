#pragma once

// Helper class for generating arrays in CSR format.

#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <type_traits>

#include "mkl.h"
#include "mkl_spblas.h"

template<typename T>
struct csr {
    std::size_t n, m;
    std::vector<T> v;
    std::vector<MKL_INT> r, c;

    std::vector<T *> _a;
    std::vector<MKL_INT *> _b;

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
        return v.empty() ? nullptr : v.data();
    }
    MKL_INT *col_index_data() {
        return c.empty() ? nullptr : c.data();
    }
    MKL_INT *row_index_data() {
        return r.empty() ? nullptr : r.data();
    }
    sparse_status_t mkl(
        sparse_matrix_t *mkls,
        T alpha = T()) {

        T *data = val_data();
        std::vector<T> temp;
        if (alpha != T()) {
            temp = v;
            for (auto &e : temp)
                e *= alpha;
            data = temp.data();
        }

        T *vv; MKL_INT *rr, *rrr, *cc;
        vv = (T *) mkl_malloc(sizeof(T) * v.size(), 64);
        cc = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * c.size(), 64);
        rr = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * (r.size() + 1), 64);
        rrr = (MKL_INT *) mkl_malloc(sizeof(MKL_INT) * (r.size() + 1), 64);
        memset(vv, 0, sizeof(T) * v.size());
        memset(cc, 0, sizeof(MKL_INT) * c.size());
        memset(rr, 0, sizeof(MKL_INT) * (r.size() + 1));
        memset(rrr, 0, sizeof(MKL_INT) * (r.size() + 1));

        std::memcpy(vv, data, sizeof(T) * v.size());
        std::memcpy(cc, c.data(), sizeof(MKL_INT) * c.size());
        std::memcpy(rr, r.data(), sizeof(MKL_INT) * (r.size() + 1));
        std::memcpy(rrr, r.data(), sizeof(MKL_INT) * (r.size() + 1));

        _a.push_back(vv);
        _b.push_back(rr);
        _b.push_back(rrr);
        _b.push_back(cc);

        sparse_status_t rv;
        if constexpr (std::is_same_v<T, float>) {
            rv = mkl_sparse_s_create_csr(
                mkls,
                SPARSE_INDEX_BASE_ZERO,
                n,
                m,
                rr,
                rrr + 1,
                cc,
                vv);
        } else {
            rv = mkl_sparse_d_create_csr(
                mkls,
                SPARSE_INDEX_BASE_ZERO,
                n,
                m,
                rr,
                rrr + 1,
                cc,
                vv);
        }

        return rv;
    }

    csr &operator()(T value, std::size_t row, std::size_t column) {

        if (row >= n) {
            throw std::invalid_argument("row out bounds");
        }
        else if (column >= m) {
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

        v.insert(v.begin() + rowi, value);
        c.insert(c.begin() + rowi, column);
        for (std::size_t i = row + 1; i <= n; i++) {
            r[i] += 1;
        }
        return *this;
    }

    csr() : n(0), m(0) { }

    csr(std::size_t n, std::size_t m) : n(n), m(m) {
        this->v = std::vector<T>();
        this->c = std::vector<MKL_INT>();
        this->r = std::vector<MKL_INT>(n + 1, 0);
    }

    // Copy constructor - copies data but not MKL allocations
    csr(const csr<T> &that) : n(that.n), m(that.m), v(that.v), r(that.r), c(that.c) {
        // _a and _b intentionally not copied - new object manages its own MKL allocations
    }

    // Copy assignment - copies data but not MKL allocations
    csr<T> &operator=(const csr<T> &that) {
        if (this != &that) {
            // Free existing MKL allocations
            for (auto &e : _a) MKL_free(e);
            for (auto &e : _b) MKL_free(e);
            _a.clear();
            _b.clear();

            n = that.n;
            m = that.m;
            v = that.v;
            c = that.c;
            r = that.r;
        }
        return *this;
    }

    // Move constructor
    csr(csr<T> &&that) noexcept
        : n(that.n), m(that.m),
          v(std::move(that.v)), r(std::move(that.r)), c(std::move(that.c)),
          _a(std::move(that._a)), _b(std::move(that._b)) {
        that.n = 0;
        that.m = 0;
        that._a.clear();
        that._b.clear();
    }

    // Move assignment
    csr<T> &operator=(csr<T> &&that) noexcept {
        if (this != &that) {
            // Free existing MKL allocations
            for (auto &e : _a) MKL_free(e);
            for (auto &e : _b) MKL_free(e);

            n = that.n;
            m = that.m;
            v = std::move(that.v);
            r = std::move(that.r);
            c = std::move(that.c);
            _a = std::move(that._a);
            _b = std::move(that._b);

            that.n = 0;
            that.m = 0;
            that._a.clear();
            that._b.clear();
        }
        return *this;
    }

    ~csr() {
        for (auto &e: _a) {
            MKL_free(e);
        }
        for (auto &e: _b) {
            MKL_free(e);
        }
    }
};
