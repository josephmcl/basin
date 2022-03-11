#pragma once
#include <cmath>

namespace physical {

template<double in, double out, double c, double r̄, double rw> 
auto constexpr μ() noexcept {
    return [](double x, double y){ return (out - in) / 2. * (tanh(( x * 
        x + c * c * y * y - r̄) / rw) + 1.) + in; }; }


template<double in, double out, double c, double r̄, double rw> 
auto constexpr ρ() noexcept {
    return [](double x, double y){ return (out - in) / 2. * (tanh(( x * 
        x + c * c * y * y - r̄) / rw) + 1.) + in; }; }

} /* physical:: */
