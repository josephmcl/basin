#pragma once 
#include <vector>
#include <cmath>

#include <iostream>

namespace domain {

    //struct _domain {

    //};

    //struct _domain TheDomain;



template<int length, double r̂, double l> auto 
transform(std::vector<double> &r, std::vector<double> &s) -> void {

    double constexpr a = (length - length * r̂ - length) 
        / (2 * tanh((r̂ - 1) / l) + tanh(-2 / l) * (r̂ - 1));
    double constexpr b = (a * tanh(-2 / l) + length) / 2;
    double constexpr c = length - b;

    // (A .* tanh.((r .- 1) ./ l) .+ b .* r .+ c,

    std::cout << 
        "a: " << a << std::endl <<
        "b: " << b << std::endl <<
        "c: " << c << std::endl;}


} /* namespace domain */