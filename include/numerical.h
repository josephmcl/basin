#pragma once 
#include <vector>
#include <cmath> 
#include <tgmath.h>
#include <functional>
#include <iostream>
#include <array>
#include <algorithm>
#include <tuple>

#include "definitions.h"
#include "ranges.h"

#include <petscsys.h>
#include <petscmat.h> 

namespace numerical {

using namespace type;

struct operators {

    using ℤ = std::size_t;
    using ℝ = type::real_t; 

    operators(){};
    
    /* Returns finite difference H matrix. */
    std::vector<ℝ> H(
        ℤ nodes, ℤ order=2, ℝ left=-1., ℝ right=1.);

    std::vector<ℝ> H_inverse(
        ℤ nodes, ℤ order=2, ℝ left=-1., ℝ right=1.);

    struct sbp {
        using row_p = std::tuple<ℤ, std::vector<ℝ> const *>;
        using row_t = std::tuple<ℤ, std::vector<ℝ> const>;
        ℤ const size, order; 
        ℝ const left, right, grid_size;
        std::vector<ℝ> d;
        std::vector<ℝ> h, hi;
        std::vector<std::vector<ℝ>> top, bottom;
        Mat _petsc_repr;
        sbp(ℤ const size, ℤ const order, ℝ const left, ℝ const right);
        /*  Given a const reference to an index, return a reference to 
            a vector of tuples of (size_t, long double) representing 
            the index and the value of rows. */
        row_p row(ℤ const index) const;

        row_t rowf(ℤ const index) const;

        auto product() {
            return [size=size, d=d, h=h, top=top, bottom=bottom]
            (std::vector<ℝ> const &rhs, std::vector<ℝ> &lhs){

                for (std::size_t i = 0; i != lhs.size(); ++i)
                    lhs[i] = 0.;

                for (std::size_t i = 0; i != lhs.size(); ++i) {

                    std::size_t j;
                    auto res = std::vector<real_t>();

                    if (i < top.size()) {
                        j = 0;
                        res = top[i];
                    }  
                    else if (i < size - bottom.size()) {
                        // Adjust column index for kernel size. 
                        j = i - ((d.size() - 1) / 2);
                        res = d;
                    }
                    else if (i < size) {
                        // Find logical index from explicitly stored kernel.
                        auto ii = i - (size - bottom.size());
                        // Adjust column index for kernel size.
                        j = size - bottom[ii].size();
                        res = bottom[ii];
                    }
                    else { throw; }

                    // // Apply the grid spacing "H" matrix. 
                    // for (std::size_t ii = 0; ii < res.size(); ++ii) {
                    //     res[ii] *= h[j + ii];
                    // } 

                    for (std::size_t jj = 0; jj != res.size(); ++jj) {
                        lhs[i] += res[jj] * rhs[j + jj];
                    }
                }
            };
        }

    };
    struct d1: sbp {
        d1(ℤ const size, ℤ const order=2, ℝ const left=-1., 
           ℝ const right=1.) : sbp(size, order, left, right) {
            load_operator(); 
        };
        void load_operator();
    };

    struct d2: sbp {

        
        std::vector<ℝ> d1_interior;
        std::vector<ℝ> h, hi;
        std::vector<std::vector<ℝ>> d1_top, d1_bottom;

        d2(ℤ const size, ℤ const order=2, ℝ const left=-1., 
           ℝ const right=1.) : sbp(size, order, left, right) {
            load_operator(); 
            h = std::vector<ℝ>(size, 1.); 
            hi = std::vector<ℝ>(size, 1.); 
        };
        void load_operator();
        void fuse(std::vector<ℝ> const &diag);
        void fuse_hi(std::vector<ℝ> const &diag);
        void left_boundary(ℝ value, ℤ degree=1);
        void right_boundary(ℝ value, ℤ degree=1);
    };


}; /* struct operators */
    
}; /* numerical:: */

namespace linalg {
    
    // In-place vector vector difference
    template <typename VectorType>
    void ip_vv_diff(VectorType &a, VectorType const &b) {
        for (std::size_t i = 0; i != a.size() || i != b.size(); ++i) 
            a[i] -= b[i];
    }

    // vT * v multiplication
    template <typename VectorType>
    auto vtv_prod(VectorType const &a, VectorType const &b) {
        auto res = a[0]; res = 0;
        for (std::size_t i = 0; i != a.size() || i != b.size(); ++i) 
            res += a[i] * b[i];
        return res;
    }

    template <typename VectorType, typename DataType>
    auto vs_prod(VectorType const &a, DataType const b, VectorType &c) {
        for (std::size_t i = 0; i != a.size() || i != c.size(); ++i) 
            c[i] = a[i] * b; 
    }

    template <typename VectorType>
    auto ip_vv_sum(VectorType &a, VectorType const &b) {
        for (std::size_t i = 0; i != a.size() || i != b.size(); ++i) 
            a[i] += b[i];
    }
}; /* linalg:: */



