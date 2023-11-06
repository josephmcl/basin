#include "infrastructure.h"

bool infrastructure::operator !(infrastructure::error e) {
    return static_cast<std::size_t>(e) == 0? false: true;
}

infrastructure::error
infrastructure::initialize(int c, char **v) {
    using namespace infrastructure;
    bool clean = true;
    // Check that OPENMP is installed. 
    #ifdef _OPENMP 
        clean = true;
    #else
        clean = false;
    #endif 
    return !clean ?error::openmp_install_failure :error::nil; 
}


void infrastructure::cleanup() {
    using namespace infrastructure;
    return;
}
