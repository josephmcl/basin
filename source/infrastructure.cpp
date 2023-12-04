#include "infrastructure.h"

#include <papi.h>
#include <iostream>

bool infrastructure::operator !(infrastructure::error e) {
    return static_cast<std::size_t>(e) == 0? false: true;
}

infrastructure::error
infrastructure::initialize(int c, char **v) {
    using namespace infrastructure;

    bool clean = true;

    // Initialize PETSc.
    PetscBool petsc_success;
    PetscInitialize(&c, &v, nullptr, nullptr);
    PetscInitialized(&petsc_success); 
    clean = static_cast<bool>(petsc_success);

    // PetscLogNestedBegin();
    // PetscMemorySetGetMaximumUsage();
    
    if (!clean) return error::petsc_init_failure;

    // Check that OPENMP is installed. 
    #ifdef _OPENMP 
        clean = true;
    #else
        clean = false;
    #endif 

    int retval;

    
    retval=PAPI_library_init(PAPI_VER_CURRENT);
    if (retval!=PAPI_VER_CURRENT) {
            fprintf(stderr,"Error initializing PAPI! %s\n",
                    PAPI_strerror(retval));
    }
    else {
        std::cout << "Loaded papi." << std::endl;
    }
    
    return !clean ?error::openmp_install_failure :error::nil; 
}


void infrastructure::cleanup() {
    using namespace infrastructure;
    // PetscLogDump("petsc.log");
    // PetscLogView(PETSCVIEWERASCII);
    PetscFinalize();
    return;
}
