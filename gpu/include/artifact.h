#include <string>
#include <filesystem>

#include "csr.h"
#include "io.h"
#include "components.h"
#include "definitions.h"

#include "rocsolver.h"
#include "rocblas.h"

using namespace io;

template <bool _>
struct artifact {
    artifact(std::string n, std::size_t s, std::size_t c) : 
        name(n),
        block_size_dim(s), 
        block_count_dim(c) {

            locate_data();

        }
    const std::string name;
    std::size_t block_size_dim  = 0;
    std::size_t block_count_dim = 0;
    struct {
        bool unfactored_sparse       = false;
        bool cholesky_factors_dense  = false;
        bool cholesky_factors_sparse = false;
    } loaded;
    struct {
        bool unfactored_sparse       = false;
        bool cholesky_factors_dense  = false;
        bool cholesky_factors_sparse = false;
    } exists;
    csr_t   unfactored_sparse;
    real_t *cholesky_factors_dense;
    csr_t   cholesky_factors_sparse;
      
    void locate_data() {
        std::size_t total = block_size_dim * block_size_dim * 
        block_count_dim * block_count_dim;
        std::string dir = "../operator/V_" + std::to_string(total) + 
            "_N_" + std::to_string(block_size_dim) + "_L_" + 
            std::to_string(block_count_dim) + "/";
        std::string ufs =  dir + name + ".val";
        std::string cfs =  dir + "sparse_" + name + ".val";
        std::string cfd =  dir + name + ".cholesky_factors";

        namespace fs = std::filesystem;
        
        if (fs::exists(ufs) && fs::is_regular_file(ufs)) {
            exists.unfactored_sparse = true;
        }
        if (fs::exists(cfd) && fs::is_regular_file(cfd)) {
            exists.cholesky_factors_dense = true;
        }
        if (fs::exists(cfs) && fs::is_regular_file(cfs)) {
            exists.cholesky_factors_sparse = true;
        }
    } 

    void load_unfactored_sparse() {
        if (!loaded.unfactored_sparse && exists.unfactored_sparse) {

            load_operator(unfactored_sparse, name, block_size_dim, 
                block_count_dim);

            loaded.unfactored_sparse = true;
        }
    }

    void load_cholesky_factors_dense(
        rocblas_handle &rb_handle, components &sbp) {

        if (!loaded.cholesky_factors_dense) {

            load_unfactored_sparse();

            int info;
            rocblas_status rbstatus;

            cholesky_factors_dense = new real_t[unfactored_sparse.n * unfactored_sparse.m];

            // Generate and store
            if (!exists.cholesky_factors_dense) {
                rbstatus = rocsolver_dpotf2(
                    rb_handle, 
                    rocblas_fill_lower, 
                    unfactored_sparse.n, 
                    cholesky_factors_dense, 
                    unfactored_sparse.n, 
                    &info);

                hipDeviceSynchronize();
                if (rbstatus != 0) {
                    std::cout << "rocsolver_dpotf2 (runtime error)" 
                        << std::endl;
                    return;
                }
                else if (info != 0) {
                    std::cout << "rocsolver_dpotf2 (factorization "
                        << "error)" << std::endl;
                    return;
                }
                else {
                    std::cout << "Factorized " << name << std::endl;
                    write_chol(
                        cholesky_factors_dense, 
                        unfactored_sparse.n, 
                        sbp, 
                        name);
                    exists.cholesky_factors_dense = true;
                }
                    
            }

            // Load existing data
            else {
                load_factors(
                    cholesky_factors_dense, 
                    unfactored_sparse.n,
                    sbp, 
                    name);
            }

            loaded.cholesky_factors_dense = true;
        }
    }

    void load_cholesky_factors_sparse(
        rocblas_handle &rb_handle, components &sbp) {
    
        if (!loaded.cholesky_factors_sparse) {

            load_unfactored_sparse();
            
            // Generate and store
            if (!exists.cholesky_factors_sparse) {

                load_cholesky_factors_dense(rb_handle, sbp);
                
                const std::string n = "sparse_" + name;

                print_csr(
                    cholesky_factors_dense, 
                    unfactored_sparse.n, 
                    sbp, 
                    n);

                exists.cholesky_factors_sparse = true;
            }

            // Load existing data
            else {
                load_operator(
                    cholesky_factors_sparse, 
                    "sparse_" + name, 
                    block_size_dim, 
                    block_count_dim);
            }

            loaded.cholesky_factors_sparse = true;
        }
    }

    ~artifact() {
        if (loaded.cholesky_factors_dense) {
            delete[] cholesky_factors_dense;
        }
    }
}; 