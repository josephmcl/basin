#include "definitions.h"

bool solver::overwrite(solver_t s) {
    return solver_info_key[
        static_cast<std::underlying_type<solver_t>::type>(s)].overwrite_b;
}

storage_type solver::storage(solver_t s) {
    return solver_info_key[
        static_cast<std::underlying_type<solver_t>::type>(s)].storage;
}
