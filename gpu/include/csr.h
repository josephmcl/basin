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
    std::vector<int> r, c, p;

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

    T *dense() {
        T *rv = (T *) malloc(sizeof(T) * n * m);
        int i,j;
        #pragma omp parallel for target teams distribute
        for(i=0; i< m * n; i++)
            rv[i] = 0.0;
        #pragma omp parallel for collapse(2) target teams distribute
        for(i = 0; i < n; i++) {
            for(j = r[i]; j < r[i+1]; j++) {
                int col_ind = c[j];
                rv[i + n * col_ind] = v[j];    /*column major*/
            }
        }
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
		this->c = std::vector<int>();
		this->r = std::vector<int>(n + 1, 0);
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
    

    csr(csr const &that, bool transpose = false) { 
        if (!transpose) {
            *this = that;
            return;
        }
        else {

            this->n = that.m;
            this->m = that.n;
            this->v.resize(that.nnz());
            this->r.resize(that.m + 2);
            this->c.resize(that.nnz());
            
            int count; 
            for (std::size_t i = 0; i < that.nnz(); ++i) {
                ++this->r[that.c[i] + 2];

            }

            for (std::size_t i = 2; i < this->r.size(); ++i) {
                // create incremental sum
                this->r[i] += this->r[i - 1];
            }

            for (int i = 0; i < that.n; ++i) {
                for (int j = that.r[i]; j < that.r[i + 1]; ++j) {
                    const int new_index = this->r[that.c[j] + 1]++;
                    this->r[new_index] = that.r[j];
                    this->c[new_index] = i;
                }
            }
            this->r.pop_back(); 
        }
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

  int ndim, mdim;
  std::ifstream ifs(head + key + ".ndim", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&ndim), sizeof(int));
  ifs.close();
  ifs = std::ifstream (head + key + ".mdim", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&mdim), sizeof(int));
  ifs.close();

  int rsize, csize;
  ifs = std::ifstream(head + key + ".rsize", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&rsize), sizeof(int));
  ifs.close();
  ifs = std::ifstream (head + key + ".csize", std::ios::binary);
  ifs.read(reinterpret_cast<char*>(&csize), sizeof(int));
  ifs.close();
  
  mat.n = ndim;
  mat.m = mdim;
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
  std::cout << "Loaded CSR data \"" << head + key << "*\" (" 
    << std::fixed << gb << " GB)" << std::endl;

    //for (int i = 0; i < csize; ++i) {
    //    std::cout << mat.val_data()[i] << " ";
    //}
    //std::cout << std::endl;
}


