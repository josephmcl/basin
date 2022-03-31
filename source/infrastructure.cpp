#include "infrastructure.h"

bool infrastructure::operator !(infrastructure::error e) {
    return static_cast<std::size_t>(e) == 0? false: true;
}

infrastructure::error
infrastructure::initialize() {
    using namespace infrastructure;
    int c = 0;
    char **v = nullptr;
    PetscErrorCode ierr;
    PetscBool      success;

    ierr = PetscInitialize(&c, &v, nullptr, nullptr);
    ierr = PetscInitialized(&success); 


    return !success ?error::petsc_init_failure :error::nil; 
}


void infrastructure::cleanup() {
    using namespace infrastructure;
    PetscErrorCode ierr;
    ierr = PetscFinalize();
    return;
}
