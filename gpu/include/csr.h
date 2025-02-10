#pragma once 

// Helper class for generating arrays in CSR format. 

#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <fstream>


//#include "mkl.h"
//#include "mkl_spblas.h"



template<typename T>
struct csr {
    std::size_t n, m;
    std::vector<T> v;
    std::vector<int> r, c;

    std::vector<T *> _a;
    std::vector<std::size_t *> _b;

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
    std::size_t nnz(int sz) { 
        v.resize(static_cast<std::size_t>(sz));
        return v.size();
    }
    std::size_t col_index_size(int sz) { 
        c.resize(static_cast<std::size_t>(sz));
        return c.size();
    }
    std::size_t row_index_size(int sz) { 
        r.resize(static_cast<std::size_t>(sz));
        return r.size();
    }
    T *val_data() {
        return &v[0];
    }
    int *col_index_data() {
        return &c[0];
    }
    int *row_index_data() {
        return &r[0];
    }

    /*
    sparse_status_t mkl(
        sparse_matrix_t *mkls, 
        T alpha = T()) {
        
        T *data = val_data();
        std::vector<T> temp;
        if (alpha != T()) {
            temp = v;
            for (auto &e : temp) 
                e *= alpha;
            data = &temp[0];
        }

        T *vv; std::size_t *rr, *cc;
        vv = (T *) mkl_malloc(sizeof(T) * v.size(), 64);
        cc = (std::size_t *) mkl_malloc(sizeof(std::size_t) * c.size(), 64);
        rr = (std::size_t *) mkl_malloc(sizeof(std::size_t) * r.size(), 64);
        memset(vv, 0, sizeof(T) * v.size());
        memset(cc, 0, sizeof(std::size_t) * c.size());
        memset(rr, 0, sizeof(std::size_t) * r.size());

        std::memcpy(vv, &data[0], sizeof(T) * v.size());
        std::memcpy(cc, &c[0], sizeof(std::size_t) * c.size());
        std::memcpy(rr, &r[0], sizeof(std::size_t) * r.size());

        _a.push_back(vv);
        _b.push_back(rr);
        _b.push_back(cc);

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
    */

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
        std::size_t rowi = r[row];
        std::size_t coli = 0;
        for (; rowi < r[row + 1]; rowi++) {
			coli = c[rowi];
			if (static_cast<std::size_t>(coli) >= column) {
				break;
			}
		}

        // overwrite existing value
        //if (static_cast<std::size_t>(coli) == column && nnz() > 0) { 
        //    v[rowi] = value;
        //}
        //else {
            v.insert(v.begin() + rowi, value);
            c.insert(c.begin() + rowi, column);
            for (std::size_t i = row + 1; i <= n; i++) {
                r[i] += 1;
            }
        //}
        return *this;
    }
    csr(std::size_t n, std::size_t m) : n(n), m(m) {
        this->v = std::vector<T>();
		this->c = std::vector<std::size_t>();
		this->r = std::vector<std::size_t>(n + 1, 0);
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
    /*
    ~csr() {
        for (auto &e: _a) {
            free(e);
        }
        for (auto &e: _b) {
            free(e);
        }
    }
    */
};

template<typename T>
void load_operator(csr<T> &mat, std::string key, 
  std::size_t block_size, std::size_t block_count) {
  
  const std::size_t total = 
    block_size * block_size * block_count * block_count;
  std::string head = "../operator/V_" + std::to_string(total) + 
    "_N_" + std::to_string(block_size) + "_L_" + 
    std::to_string(block_count) + "/";

  int rsize, csize;
  std::ifstream ifs(head + key + ".rsize", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&rsize), sizeof(int));
  ifs.close();
  ifs = std::ifstream (head + key + ".csize", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&csize), sizeof(int));
  ifs.close();
  
  mat.n = rsize;
  mat.m = rsize;
  mat.row_index_size(rsize + 1);
  mat.col_index_size(csize);
  mat.nnz(csize);
  
  ifs = std::ifstream (head + key + ".row", std::ios::binary);
  ifs.read(
    reinterpret_cast<char*>(mat.row_index_data()), 
    sizeof(int) * (rsize + 1));
  ifs.close();

  ifs = std::ifstream (head + key + ".col", std::ios::binary);
  ifs.read(
    reinterpret_cast<char*>(mat.col_index_data()), 
    sizeof(int) * csize);
  ifs.close();
  
  ifs = std::ifstream (head + key + ".val", std::ios::binary);
  ifs.read(
    reinterpret_cast<char*>(mat.val_data()), 
    sizeof(T) * csize);
  ifs.close();

  double gb = ((sizeof(int) * (rsize + 1)) + (sizeof(int) * csize) 
    + (sizeof(T) * csize)) / 1e9;
  std::cout << "Loaded csr data \"" << head + key << "*\" (" 
    << std::fixed << gb << " GB)" << std::endl;
}


