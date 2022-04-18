#pragma once 

#include <cmath>

/* Given a positive definite matrix A, an guess vector x, a vector b, a 
tolerance ϵ, and an maximum iteration count itermax, solve the form 
Ax=b using the conjugate gradient method to the minimized error 
specified by ϵ or for itermax iterations. */
namespace solve {

template <
    typename MatrixAction,
    typename InplaceVVDiff,
    typename VTVProd,
    typename VSProd, 
    typename InplaceVVPSum>
struct conjgrad {

    MatrixAction matrix_action;
    InplaceVVDiff inplace_vv_diff;
    VTVProd vtv_prod;
    VSProd vs_prod;
    InplaceVVPSum inplace_vv_sum;

    template <
        typename DataType,
        typename VectorType>
    void operator() (
        VectorType &x, 
        VectorType const &b, 
        DataType ϵ, 
        std::size_t size,
        std::size_t itermax=100) {

        VectorType temp(size), r(size), r_prior(size), p(size), 
            q(size), p_prior(size);
        DataType ρ_prior, ρ, δ, error;

        /* r_prior = A * x - b */
        matrix_action(x, r_prior);
        inplace_vv_diff(r_prior, b);
        r = r_prior;

        ρ_prior = static_cast<DataType>(0.);

        /* ρ = r' * r */
        ρ = vtv_prod(r, r);
        p = r;

        for (std::size_t i = 0; i != itermax; ++i) {
        
            /* q = A * p */
            matrix_action(p, q);

            /* δ = ρ / (p' * q) */
            δ = ρ / vtv_prod(p, q);

            /* x = x - δ * p */
            vs_prod(p, δ, temp);
            inplace_vv_diff(x, temp);

            ρ_prior = ρ;

            /* r_prior = r_prior - δ * q */
            vs_prod(q, δ, temp);
            inplace_vv_diff(r_prior, temp);

            p_prior = p;

            /* Compute error. */ 

            // num = A * x - b
            // num = sqrt(num' * num) 
            // den = sqrt(x' * x)
            // return (num / den)[1]

            /* r = A * x - b */
            matrix_action(x, r);
            inplace_vv_diff(r, b);

            error = std::sqrt(vtv_prod(r, r)) 
                  / std::sqrt(vtv_prod(x, x));

            if (error <= ϵ) break;

            // ρ = r' * r
            ρ = vtv_prod(r, r);

            // β = ρ / ρ_prior
            // p = r + β * p_prior 
            vs_prod(p_prior, ρ / ρ_prior, p);
            inplace_vv_sum(p, r);
        }
    }
}; /* struct conjgrad */
}; /* solve:: */