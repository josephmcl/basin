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

int main (int argc, char *argv[]) {
    
    if (!infrastructure::initialize()) exit(-1);; {

    std::size_t N = 2187 * 3 ;

    auto grid = numerics::linrange(0., 1., N);

    std::vector<type::real_t> gridm;
    std::vector<type::real_t> soln;
    
    /* 
    for (auto &i : grid) {
        gridm.push_back(i);
        for (auto &j : grid) {
            soln.push_back(poisson2d::analytical_solution(i, j));
        }
    }
    */ 

    /* We are solving the Poisson Equation in 1D

        ⎧ u_xx π^2 sin(πx) = 0 , 0 < x < 1
        ⎨ u = 0    (Dirichlet) , x = 0 
        ⎩ u_x = -π (Neumann)   , x = 1 
    */

    auto operators = numerical::operators();

    std::size_t size = 100;
    std::size_t order = 2;
        
    auto h = operators.H(size, order, 0., 1.);
    auto hi = operators.H_inverse(size, order, 0., 1.);
    auto d1 = numerical::operators::d1(size, order, 0., 1.);
    auto d2 = numerical::operators::d2(size, order, 0., 1.);
    d2.fuse(h);
    d2.fuse_hi(hi);
    d2.left_boundary(-40);
    d2.right_boundary(1, 2);

    std::vector<type::real_t> lhs(size, 0.);
    std::vector<type::real_t> rhs(size, 1.);

    for (std::size_t i = 0; i != size; ++i) {
        rhs[i] = std::sin(2 * 3.14159 * i / size);
    }

    // d2.product(rhs, lhs);

    using vec = std::vector<type::real_t>;
    using dat = type::real_t;

    std::vector<type::real_t> x(size, 0.);
    std::vector<type::real_t> b(size, 1.);

    for (std::size_t i = 0; i != size; ++i) {
        b[i] = -3.14159 * 3.14159 * std::sin(3.14159 * (static_cast<
            type::real_t>(i) / static_cast<type::real_t>(size)));
    }
    // b[size - 1] = hi[size - 1] - b[size - 1];

    auto cg = solve::conj_grad {
        .size            = size,
        .matrix_action   = d2.product(),
        .inplace_vv_diff = linalg::ip_vv_diff<vec>,
        .vtv_prod        = linalg::vtv_prod<vec>,
        .vs_prod         = linalg::vs_prod<vec, dat>,
        .inplace_vv_sum  = linalg::ip_vv_sum<vec>
    };

    type::real_t epsilon = 10e-10;
    cg(x, b, epsilon, size * size);

    (cg.converged)
        ? std::cout << "Converged." << std::endl
        : std::cout << "Failed to converge." << std::endl;
    
    for (auto v = x.begin(); v != x.end(); ++v) 
        std::cout << *v << ", ";
    std::cout << std::endl;
    

    /*
    for (std::size_t i = 0; i < size; ++i) {
        std::cout << h[i] << std::endl;
    }

    for (std::size_t i = 0; i < size; ++i) {
        
        auto [offset, values] = d2.rowf(i);
        std::cout << offset << " : " ;
        for (auto v = values.begin(); v != values.end(); ++v) 
            std::cout << *v << ", ";
        std::cout << std::endl;
        
    }
    */
    
    // plot::matrix(gridm, gridm, soln);

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

    

    std::vector<double> ψδ;
    domain::intitial_conditions(metrics, fault_params, ψδ);

    /*
        -∇ * (b∇u) = f
        u = g_D
        n * b∇u = gN

    

        [M     F][u]   [g    ]
        [       ][ ] = [     ] 
        [F^(t) D][λ]   [g_(δ)]

        H σ_(k) + G_(k) u - H τ_(k) (L_(k) u - λ_(k))

        M u = q 

        q = kron(H, H) J f - sum(1...4, F_k, λ_k) 

        (D - F^(t) M^(-1) F) λ = g_(δ) - F^(t) M^(-1) g
    

    std::cout << std::setprecision(15);
    for (auto &i: ψδ) 
        std::cout << i << ", "; 
    std::cout << std::endl;
    */

    } infrastructure::cleanup(); return 0;

}