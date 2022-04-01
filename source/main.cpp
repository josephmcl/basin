#include <ranges>
#include <vector>
#include <algorithm>
#include <tgmath.h>
#include <iomanip>

#include "infrastructure.h" 
#include "ranges.h"
#include "domain.h"
#include "physical.h"
#include "numerical.h"

int main (int argc, char *argv[]) {
 
    if (!infrastructure::initialize()) exit(-1);; {

    domain::data constexpr D = {
        .length      =  24.  ,
        .basin_depth =   4.  ,
        .r̂           =   0.75, 
        .l           =   0.05,
        .μ_out       =  36.  ,
        .μ_in        =   8.  ,
        .ρ_out       =   2.8 ,
        .ρ_in        =   2.  };

    /* Construct the grid. */
    auto metrics = domain::metrics<D>(400, 400);

    auto fault_params = physical::fault_params(metrics.y());

    numerical::operators operators = {};

    std::vector<double> ψδ;
    domain::intitial_conditions(metrics, fault_params, ψδ);

    std::cout << std::setprecision(15);
    for (auto &i: ψδ) 
        std::cout << i << ", "; 
    std::cout << std::endl;

    } infrastructure::cleanup(); return 0;

}