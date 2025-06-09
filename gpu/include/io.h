#pragma once 
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "definitions.h"
#include "components.h"

namespace io {

    const std::string operator_directory = "../operator/";
    const std::string dense_chol_suffix = ".cholesky_factors";

    std::string make_operator_path(
        std::size_t  quant, 
        std::size_t  size);

    void write_chol(
        real_t *m, 
        std::size_t sz, 
        components &sbp, 
        std::string key);
  
    bool load_factors(
        real_t *m, 
        std::size_t sz, 
        components &sbp, 
        std::string key);
        
        void write_chol_sparse(
            real_t *m, 
            std::size_t sz, 
            components &sbp, 
            std::string key);

    void print_csr(real_t *m, std::size_t sz, components &sbp, std::string key, std::size_t szb = 0);
}