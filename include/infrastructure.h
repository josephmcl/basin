#pragma once

#include <petscsys.h>

namespace infrastructure {
    enum class error : std::size_t { 
        nil                 = 0,
        petsc_init_failure  = 1};
    bool operator !(error e);
    error initialize();
    void cleanup();
}