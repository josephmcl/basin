#include "poisson_2d.h"

#include "csr.h"
#include <iostream>

#include <tuple>
#include <vector>
#include <string>

int main(int argc, char **argv) {

    /* 
    std::vector<std::tuple<std::size_t, std::size_t>> params = {
        {280, 3}, {210, 4}, {168, 5},  {140, 6}, {120, 7}, {105, 8}, {84, 10},
        {70, 12}, {60, 14}, {56, 15}, {42, 20}, {40, 21},  {35, 24}, 
        {30, 28}, {28, 30}, {24, 35}, {21, 40}
    };
    */
    
   if (argc < 3) exit(-1);

    //for (auto &e : params)
    poisson_2d::problem(std::stoi(argv[1]), std::stoi(argv[2]));
    //poisson_2d::problem(std::get<0>(e), std::get<1>(e));
}
