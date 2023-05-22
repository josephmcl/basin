#pragma once

#include <petscsys.h>

namespace infrastructure {
    enum class error : std::size_t { 
        nil                    = 0,
        petsc_init_failure     = 1,
        openmp_install_failure = 2};
    bool operator !(error e);
    error initialize(int c, char **v);
    void cleanup();
}