#pragma once 
#include <vector>
#include <cmath>

#include <iostream>

namespace domain {

    //struct _domain {

    //};

    //struct _domain TheDomain;

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


} /* namespace domain */