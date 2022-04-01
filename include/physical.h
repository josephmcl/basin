#pragma once
#include <cmath>
#include <tuple>
#include <vector>
#include <iomanip>

#include "ranges.h"
#include "domain.h"

namespace physical {

template<double in, double out, double c, double r̄, double rw> 
auto constexpr μ() noexcept {
    return [](double x, double y){ return (out - in) / 2. * (tanh(( x * 
        x + c * c * y * y - r̄) / rw) + 1.) + in; }; }


template<double in, double out, double c, double r̄, double rw> 
auto constexpr ρ() noexcept {
    return [](double x, double y){ return (out - in) / 2. * (tanh(( x * 
        x + c * c * y * y - r̄) / rw) + 1.) + in; }; }


template <typename T>
using fault_param = std::tuple<std::size_t, std::size_t, std::size_t, 
    std::tuple<double, double, numerics::linrange<T>, double, double,
    double, double, double>>;

template<typename T>
fault_param<T> fault_params(numerics::linrange<T> const &face) {

    using namespace numerics;

    double constexpr Wf   = 24.;
    double constexpr Hvw  = 12.;
    double constexpr Ht   =  6.;
    double constexpr b0   =  0.02;
    double constexpr bmin =  0.;

    linrange<T> faceδ = face; 
    faceδ |= [](T x) {return std::abs(Wf - x);};
    auto δ = std::min_element(faceδ.begin(), faceδ.end());

    linrange<T> faceg = face; 
    faceg |= [](T x) {return std::abs(16 - x);};
    auto g = std::min_element(faceg.begin(), faceg.end());

    linrange<T> facev = face; 
    facev |= [](T x) {return std::abs((Hvw + Ht) - x);};
    auto v = std::min_element(facev.begin(), facev.end());

    auto faceb = face; face.clip(0, δ.index); // linrange(*face[0], *face[δ.index], δ.index + 1);
    faceb |= [](T x) -> T {
        if (0 <= x && x < Hvw) 
            return b0;
        else if (Hvw <= x && x < Hvw + Ht) 
            return b0 + (bmin - b0) * (x - Hvw) / Ht;
        else if (Hvw + Ht <= x && x < Wf)
            return bmin;
        else
            return bmin; };

    double constexpr σn    = 50.;
    double constexpr a     =  0.015;
    double constexpr Dc    =  0.008;
    double constexpr f0    =  0.6;
    double constexpr V0    =  1e-6;
    double constexpr τ_inf = 24.82;
    double constexpr Vp    =  1e-9;

    return std::make_tuple(δ.index, g.index, v.index, std::make_tuple(
        σn, a, faceb, Dc, f0, V0, τ_inf, Vp)); }

} /* physical:: */
