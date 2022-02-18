#include "ranges.h"
#include "domain.h"
#include "physical.h"

#include <ranges>
#include <vector>
#include <algorithm>

int main (int argc, char *argv[]) {

    using namespace numerics;

    auto square = [](double i) { return i * i; };
    auto plusone = [](double i) { return i + 1; };

    double constexpr N           = 400.;
    double constexpr length      =  24.;
    double constexpr basin_depth =   4.;
    /* NOTE: This is r HAT */
    double constexpr r̂           =   0.75; 
    double constexpr l           =   0.05;
    double constexpr μ_out       =  36.;
    double constexpr μ_in        =  8.;
    double constexpr c           = (length / 2) / basin_depth;
    /* NOTE: This is r BAR */
    double constexpr r̄           = (length / 2) * (length / 2); 
    double constexpr rw          = 1 + (length / 2) / basin_depth;


    /* Here we decompose a disk from its physical space into a 
    reference block space.

    y                          S
    ╷                          ╷     face 4
    │    o ┄┄┄┄┄ o             ├───┬───┬───┬─┄─┐
    │   o ╲     ╱╱╲          f │   ┆   ┆   ┆   │ f     
    │  o   o┄┄┄┼┼┼┼┤         a │┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ a
    │  ┊   ┊   ├┼┼┼┤         c │   ┆   ┆   ┆   │ c
    │  o   o┄┄┄┼┼┼┼┤    -->  e ├┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ e
    │   o ╱     ╲╲╱            │   ┆   ┆   ┆   │ 
    │    o ┄┄┄┄┄ o           1 ├┄┄┄┼┄┄┄┼┄┄┄┼┄┄┄┤ 2
    │                          │   ┆   ┆   ┆   │
    └───────────────╴x         └───┴───┴───┴───┴──╴r  
                                     face 3 */

    /*  In the original problem, r and s are given as a 2D matrices,

            ┌ -1 .. -1 ┐       ┌ -1 .. 1 ┐
        r = │  :     : │,  s = │  :    : │
            └  1 ..  1 ┘       └ -1 .. 1 ┘ 
    
    or otherwise written in a Kronecker product form,

        r = ones(N + 1) ⊗ [-1:1:N + 1]
        s = [-1:1:N + 1]' ⊗ ones(N + 1).

    Instead of storing these outright we just store a single row and
    column of each in an iterator. */

    auto r = linrange(-1., 1., N + 1);
    auto s = linrange(-1., 1., N + 1);

    /*  Subsequent x, xr, y, and ys, matrices are element-wise 
    compositions using the transform_e1 and transform_e2 functions. 
    These matrices still repeat over one dimension so we can apply the
    transforms as functions. 

    xs and yr operators also exist in the original code, but are purely
    zero matrices, so we can ignore them for now. */
    auto constexpr te1 = domain::transform_e1<length, r̂, l>();
    auto constexpr te2 = domain::transform_e2<length, r̂, l>();
    auto x  = linrange(-1., 1., N + 1) | te1;
    auto xr = linrange(-1., 1., N + 1) | te2;
    auto y  = linrange(-1., 1., N + 1) | te1;
    auto ys = linrange(-1., 1., N + 1) | te2;  

    /* The J matrix from the original code is computed as the element-
    wise product of xr * ys. Here, the expansion of the matrices 
    matters, so J is represented as function both. */
    auto J = [&xr, &ys](std::size_t ixr, std::size_t iys) {
        return *xr[ixr] * *ys[iys]; };
    auto JI = [&xr, &ys](std::size_t ixr, std::size_t iys) {
        return 1. / (*xr[ixr] * *ys[iys]); };

    /* Since J varies on both dimensions, so do subsequent matrices. */
    auto rx = [&ys, &J](std::size_t i, std::size_t j) {
        return *ys[j] / J(i, j); }; /* NOTE: ys varies on the j dim. */
    auto sy = [&xr, &J](std::size_t i, std::size_t j) {
        return *xr[i] / J(i, j); }; 

    /* The same constraints on J apply to μxy. */
    auto constexpr μ = physical::μ< μ_in, μ_out, c, r̄, rw>();
    auto μxy = [&x, &y, &μ](std::size_t ix, std::size_t iy) {
        return μ(*x[ix], *y[iy]); };    

    auto crr = [&J, &rx, &μxy](std::size_t i, std::size_t j) {
        return J(i, j) * (rx(i, j) * μxy(i, j) * rx(i, j)); };\
    auto css = [&J, &sy, &μxy](std::size_t i, std::size_t j) {
        return J(i, j) * (sy(i, j) * μxy(i, j) * sy(i, j)); };


    std::cout << "crr 0, 0: " << crr(0, 0) << std::endl;
    std::cout << "crr 200, 0: " << crr(200, 0) << std::endl;
    std::cout << "crr 0, 200: " << crr(0, 200) << std::endl;
    std::cout << "crr 200, 200: " << crr(200, 200) << std::endl;

    std::cout << "css 0, 0: " << css(0, 0) << std::endl;
    std::cout << "css 200, 0: " << css(200, 0) << std::endl;
    std::cout << "css 0, 200: " << css(0, 200) << std::endl;
    std::cout << "css 200, 200: " << css(200, 200) << std::endl;

    
    return 0;
}