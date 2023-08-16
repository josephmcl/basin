#include "verify.h"

bool verify::
petsc_matrix(
    Mat const &A, 
    Mat const &B) {

    PetscBool equiv;
    MatEqual(A, B, &equiv);
    if (equiv == PETSC_TRUE) {
        std::cout << "PETSc matrices match." << std::endl;
        return true;
    }
    else {
        std::cout << "PETSc matrices differ." << std::endl;
        return false;
    }
} 