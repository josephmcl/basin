#pragma once

#include <cmath>
#include <numbers>

#include "definitions.h"


/* setup the Poisson Equation in 2D */
namespace poisson2d {

    using ℤ = std::size_t;
    using ℝ = type::real_t;  

    /* */
    template <typename T>
    T analytical_solution(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uxx(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return -(cx * cx) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uyy(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return -(cy * cy) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T ux(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return cx * pi * std::cos(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uy(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return cy * pi * std::cos(cx * pi * x + cy * pi * y); 
    }

};
