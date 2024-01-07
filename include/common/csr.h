#pragma once 

// Helper class for generating arrays in CSR format. 

#include <cstdlib>
#include <vector>
#include <stdexcept>

template<typename T>
struct csr {
    std::size_t n, m;
    std::vector<T> v;
    std::vector<std::size_t> r, c;

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
    T *val_data() const {
        return &v[0];
    }
    T *col_index_data() const {
        return &c[0];
    }
    T *row_index_data() const {
        return &r[0];
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
			if (coli >= column) {
				break;
			}
		}

        // overwrite existing value
        if (coli == column && nnz() > 0) { 
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
    csr(csr &that): n(that.n), m(that.m) { }
};
