#include "ranges.h"
#include "domain.h"
#include "physical.h"

#include <ranges>
#include <vector>
#include <algorithm>
#include <tgmath.h>
#include <iomanip>

int main (int argc, char *argv[]) {
 

    domain::data constexpr D = {
        .length      = 24.   ,
        .basin_depth =   4.  ,
        .r̂           =   0.75, 
        .l           =   0.05,
        .μ_out       =  36.  ,
        .μ_in        =   8.  ,
        .ρ_out       =   2.8 ,
        .ρ_in        =   2.  
    };

    auto metrics = domain::metrics<D>(400, 400);


    auto j = metrics.j();

    std::cout << j(100, 100) << std::endl;

    return 0;
}