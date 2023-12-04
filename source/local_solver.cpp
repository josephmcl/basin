#include "local_solver.h"

void sbp_sat::x2::make_local_solver(
    petsc_solver       &sol,  
    petsc_matrix const &mat, 
    std::size_t  const  n) {

    sol = KSP();
    KSPCreate(PETSC_COMM_WORLD, &sol); 
    KSPSetOperators(sol, mat, mat);
    //KSPSetFromOptions(sol);
    KSPSetType(sol, KSPPREONLY);
    PC pc;
    KSPGetPC(sol, &pc);
    PCSetType(pc, PCLU);
    KSPSetPC(sol, pc);
    KSPSetUp(sol);

    // Call once to factorize
    petsc_vector x, y;
    VecCreateSeq(PETSC_COMM_WORLD, n * n, &x); 
    VecCreateSeq(PETSC_COMM_WORLD, n * n, &y); 

    KSPSolve(sol, x, y);

    destroy<fw>(x);
    destroy<fw>(y);
}