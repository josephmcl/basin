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
#include "definitions.h"
#include "poisson2d.h"
// #include "plot.h"

#include "solve.h"


#include "poisson.h"

int main (int argc, char *argv[]) {
    
    if (!infrastructure::initialize()) exit(-1);; { 


    /* We are solving the Poisson Equation in 1D

        ⎧ u_xx π^2 sin(πx) = 0 , 0 < x < 1
        ⎨ u = 0    (Dirichlet) , x = 0 
        ⎩ u_x = -π (Neumann)   , x = 1 
    */


    std::vector<type::real_t> exact;
    poisson1d::analytical_blocked(exact, 100, 2);

    std::vector<type::real_t> result;
    poisson1d::petsc_problem(result, 100, 2);


    std::cout << std::setprecision(14);

    long double error = 0;
    long double temp;
    for (std::size_t i = 0; i < result.size(); ++i) {
        temp = result[i] - exact[i];
        error += temp * temp;
    }
   
    error = std::sqrt(error  / (result.size() / 2)) ;

    std::cout << exact.size() << " : " << error << std::endl;

    // TODO: Implement convergence tests, normalized with H

    } infrastructure::cleanup(); return 0;

}