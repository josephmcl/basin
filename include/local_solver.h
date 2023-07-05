#pragma once
// #include "hybrid_sbp_sat_2d.h"
#include "components.h"

namespace sbp_sat { namespace x2 {

void make_local_solver(
    petsc_solver       &sol,  
    petsc_matrix const &mat, 
    std::size_t  const  n);

};};