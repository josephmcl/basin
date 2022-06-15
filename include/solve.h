#pragma once 

#include <cmath>
#include <iostream>

/* Given a positive definite matrix A, an guess vector x, a vector b, a 
tolerance ϵ, and an maximum iteration count itermax, solve the form 
Ax=b using the conjugate gradient method to the minimized error 
specified by ϵ or for itermax iterations. */



namespace solve {

    // TODO: This template could be set up better so that the DataType
    //       VectorType template parameters are known at the 
    //       instantiation of the struct. The wrong template parameters
    //       will still prevent compilation, but not until the functor
    //       is called. Still you have to know the form of the function 
    //       signatures of each template type. This can possibly be 
    //       fixed in the future with concepts (?) 

    template <
        typename MatrixAction,
        typename InplaceVVDiff,
        typename VTVProd,
        typename VSProd, 
        typename InplaceVVPSum>

    struct conj_grad {

        // The problem size. The length of x and b. The length of any 
        // side of A. 
        std::size_t const size;
        bool converged = false;

        // matrix_action(T const &b, T &a) ie. a = Ab
        MatrixAction matrix_action; 

        // inplace_vv_diff(T &a, T const &b) ie. a -= b
        InplaceVVDiff inplace_vv_diff;

        // D c = vtv_prod(T const &a, T const &b) ie. c = aTb
        VTVProd vtv_prod;

        // vs_prod(T const &a, D b, T &c) ie. c = a * b 
        VSProd vs_prod;

        // inplace_vv_sum(T &a, T const &b) ie. a += b
        InplaceVVPSum inplace_vv_sum;

        /* Functor for conj_grad. Run conjugate gradient for Ax=b,
           where A is defined by maxtrix_action. */
        template <
            typename DataType,
            typename VectorType>

        void operator() (
            VectorType &x, 
            VectorType const &b, 
            DataType ϵ, 
            std::size_t itermax=100) {

            VectorType temp(size), r(size), r_prior(size), p(size), 
                q(size), p_prior(size);
            DataType ρ_prior, ρ, δ, error;

            converged = false;

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

                /* r = A * x - b */
                matrix_action(x, r);
                inplace_vv_diff(r, b);

                /* Compute error. */ 
                error = std::sqrt(vtv_prod(r, r)) 
                      / std::sqrt(vtv_prod(x, x));

                if (error <= ϵ) {
                    converged = true;
                    break;
                }

                // ρ = r' * r
                ρ = vtv_prod(r, r);

                // β = ρ / ρ_prior
                // p = r + β * p_prior 
                vs_prod(p_prior, ρ / ρ_prior, p);
                inplace_vv_sum(p, r);
            }
        }
    }; /* struct conj_grad */
}; /* solve:: */