#pragma once 
#include <vector>
#include <cmath> 
#include <tgmath.h>
#include <functional>
#include <iostream>
#include <iomanip>

#include "physical.h"
#include "ranges.h"

namespace domain {

template<double length, double r̂, double l> 
auto constexpr transform_e1() {
    double constexpr a = (length - length * r̂ - length) 
        / (2 * tanh((r̂ - 1) / l) + tanh(-2 / l) * (r̂ - 1));
    double constexpr b = (a * tanh(-2 / l) + length) / 2;
    double constexpr c = length - b;
    return [] (double x) {return a * tanh((x - 1) / l) + b * x + c;}; }

template<double length, double r̂, double l> 
auto constexpr transform_e2() {
    double constexpr a = (length - length * r̂ - length) 
        / (2 * tanh((r̂ - 1) / l) + tanh(-2 / l) * (r̂ - 1));
    double constexpr b = (a * tanh(-2 / l) + length) / 2;
    return [] (double x) {return ((a * pow(1 / (cosh((x - 1) / l)), 
        2.)) / l) + b;}; }



struct data {
    double length;
    double basin_depth;
    /* NOTE: This is r HAT */
    double r̂; 
    double l;
    double μ_out;
    double μ_in;
    double ρ_out;
    double ρ_in;
    double c(){ return (length / 2) / basin_depth; }
    /* NOTE: This is r BAR */
    double r̄(){ return (length / 2) * (length / 2); } 
    double rw(){ return 1 + (length / 2) / basin_depth; }
};


    /* Here we decompose a disk from its physical space into a 
    reference block space.

    y                          S
    ╷                          ╷     face 4
    │    o ┄┄┄┄┄ o             ├───┬───┬───┬─┄─┐
    │   o ╲   4 ╱╱╲          f │   ┆   ┆   ┆   │ f     
    │  o   o┄┄┄┼┼┼┼┤         a │┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ a
    │  ┊   ┊ 1 ├┼┼┼┤ 2       c │   ┆   ┆   ┆   │ c
    │  o   o┄┄┄┼┼┼┼┤    -->  e ├┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ e
    │   o ╱   3 ╲╲╱            │   ┆   ┆   ┆   │ 
    │    o ┄┄┄┄┄ o           1 ├┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ 2
    │                          │   ┆   ┆   ┆   │
    └───────────────╴x         └───┴───┴───┴───┴──╴r  
                                     face 3 */

template<data D>
struct metrics {
    metrics(double r, double s) : r_size(r), s_size(s) {}

    double r_size;
    double s_size;

    /*  In the original problem, r and s are given as a 2D matrices,

            ┌ -1 .. -1 ┐       ┌ -1 .. 1 ┐
        r = │  :     : │,  s = │  :    : │
            └  1 ..  1 ┘       └ -1 .. 1 ┘ 
    
    or otherwise written in a Kronecker product form,

        r = ones(N + 1) ⊗ [-1:1:N + 1]
        s = [-1:1:N + 1]' ⊗ ones(N + 1).

    Instead of storing these outright we just store a single row and
    column of each in an iterator. */

    auto x() { return numerics::linrange(-1., 1., r_size + 1) 
                      | transform_e1<D.length, D.r̂, D.l>(); }

    auto y() { return numerics::linrange(-1., 1., s_size + 1) 
                      | transform_e1<D.length, D.r̂, D.l>(); }

    auto xr() { return numerics::linrange(-1., 1., r_size + 1) 
                      | transform_e2<D.length, D.r̂, D.l>(); }

    auto ys() { return numerics::linrange(-1., 1., s_size + 1) 
                      | transform_e2<D.length, D.r̂, D.l>(); }


    auto j() { 
        auto _xr = xr(); 
        auto _ys = ys(); 
        return [_xr, _ys](std::size_t ixr, std::size_t iys) {
        return *_xr[ixr] * *_ys[iys]; }; }

    auto ji() { 
        auto _xr = xr(); 
        auto _ys = ys(); 
        return [_xr, _ys](std::size_t ixr, std::size_t iys) {
        return 1. / (*_xr[ixr] * *_ys[iys]); }; }

    auto rx() {
        auto _ys = ys();
        auto _j  =  j();
        return [_ys, _j](std::size_t i, std::size_t j) {
        return *_ys[j] / _j(i, j); }; }

    auto sy() {
        auto _xr = xr();
        auto _j  =  j();
        return [_xr, _j](std::size_t i, std::size_t j) {
        return *_xr[i] / _j(i, j); }; }

    auto crr() {
        auto _ys = ys();
        auto _j  =  j();
        auto _x  =  x();
        auto _y  =  y();
        double constexpr c  = (D.length / 2) / D.basin_depth;
        double constexpr r̄  = (D.length / 2) * (D.length / 2); 
        double constexpr rw = 1 + (D.length / 2) / D.basin_depth; 
        auto _μ  =  physical::μ<D.μ_in, D.μ_out, c, r̄, rw>();
        return [_ys, _j, _x, _y, _μ](std::size_t i, std::size_t j){
        return (*_ys[j] * *_ys[j] * _μ(*_x[i], *_y[j])) / _j(i, j);};}

    auto css() {
        auto _xr = xr();
        auto _j  =  j();
        auto _x  =  x();
        auto _y  =  y();
        double constexpr c  = (D.length / 2) / D.basin_depth;
        double constexpr r̄  = (D.length / 2) * (D.length / 2); 
        double constexpr rw = 1 + (D.length / 2) / D.basin_depth; 
        auto constexpr   _μ = physical::μ<D.μ_in, D.μ_out, c, r̄, rw>();
        return [_xr, _j, _x, _y, _μ](std::size_t i, std::size_t j){
        return (*_xr[i] * *_xr[i] * _μ(*_x[i], *_y[j])) / _j(i, j);};}

    auto η() {
        auto _x = *x()[0];
        auto _y = y();
        double constexpr c  = (D.length / 2) / D.basin_depth;
        double constexpr r̄  = (D.length / 2) * (D.length / 2); 
        double constexpr rw = 1 + (D.length / 2) / D.basin_depth; 
        auto constexpr _μ = physical::μ<D.μ_in, D.μ_out, c, r̄, rw>();
        auto constexpr _ρ = physical::ρ<D.ρ_in, D.ρ_out, c, r̄, rw>();
        return [_x, _y, _μ, _ρ](std::size_t i){
                double μ = _μ(_x, *_y[i]);
        return μ / (2 * std::sqrt(μ / _ρ(_x, *_y[i])));};}

    auto μf2() {
        auto _x = *x()[r_size + 1];
        auto _y = y();
        double constexpr c  = (D.length / 2) / D.basin_depth;
        double constexpr r̄  = (D.length / 2) * (D.length / 2); 
        double constexpr rw = 1 + (D.length / 2) / D.basin_depth;
        auto constexpr _μ = physical::μ<D.μ_in, D.μ_out, c, r̄, rw>();
        return [_x, _y, _μ](std::size_t i){
        return _μ(_x, *_y[i]);};}

    auto sj1() {
        auto _ys = ys();
        return [_ys](std::size_t i){ return *_ys[i]; }; }

    auto sj2() {
        auto _ys = ys();
        return [_ys](std::size_t i){ return *_ys[i]; }; }

    auto sj3() {
        auto _ys = ys();
        return [_ys](std::size_t i){ return *_ys[i]; }; }

    auto sj4() {
        auto _ys = ys();
        return [_ys](std::size_t i){ return *_ys[i]; }; }

    };


template<typename T, data D> 
void intitial_conditions(metrics<D> &m, physical::fault_param<T> 
    const &fp, std::vector<T> &res) {


    auto [δ, g, v, rp] = fp;
    auto [σn, a, faceb, Dc, f0, V0, τ_inf, Vp] = rp;
    auto η = m.η();

    auto size = (m.r_size + 1) * 2;
    res = std::vector(size, 0.);

    for (auto i = 0; i < m.r_size + 1; ++i)
        res[i] = a * std::log(2 * (V0 / Vp) * std::sinh((τ_inf - η(i) 
                 * Vp) / (σn * a)));

    for (auto i = δ + 1; i < size; ++i)
        res[i] = 0.;

    return;
}

// ψδ = zeros(2nn)
//      for n in 1:nn
//             ψδ[n] = RS.a * log(2*(RS.V0/RS.Vp) * sinh((RS.τ_inf - η[n]*RS.Vp)/(RS.σn*RS.a)))
//       end
//         ψδ[δNp .+ (1:nn)] .= 0


} /* namespace domain */
