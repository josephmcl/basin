#pragma once

#include <cmath>
#include <numbers>
#include <vector>

#include "definitions.h"
#include "ranges.h"

/* setup the Poisson Equation in 2D */
namespace poisson1d {

    using ℤ = std::size_t;
    using ℝ = type::real_t;  

    void analytical_solution(ℤ const size, ℝ const left=0., 
                             ℝ const right=1., std::vector<ℝ> &res) {
        auto h = (right - left) / static_cast<ℤ>(size - 1);
        for (std::size_t i = 0; i < size; ++i) 
            res.push_back(std::sin(h * std::numbers::pi_v<ℝ> * i));
        return;
    }

    void problem(ℤ size, ℝ left=0., ℝ right=1.) {

        auto domain = numerics::linrange<ℝ>(left, right, size); 

        ℝ left_boundary = 0.;
        ℝ right_boundary = -std::numbers::pi_v<ℝ>;

    }
};
