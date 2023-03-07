#include "infrastructure.h"

bool infrastructure::operator !(infrastructure::error e) {
    return static_cast<std::size_t>(e) == 0? false: true;
}

infrastructure::error
infrastructure::initialize(int c, char **v) {
    using namespace infrastructure;
    PetscBool      success;

    PetscInitialize(&c, &v, nullptr, nullptr);
    PetscInitialized(&success); 


    return !success ?error::petsc_init_failure :error::nil; 
}


void infrastructure::cleanup() {
    using namespace infrastructure;
    PetscFinalize();
    return;
}
