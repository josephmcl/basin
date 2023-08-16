#pragma once
#include <iostream>
#include "petscmat.h"

namespace verify {

    bool petsc_matrix(
        Mat const &A, 
        Mat const &b);

} /* namespace verify */