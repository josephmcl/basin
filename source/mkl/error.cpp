#include "error.h"
void mkl_sparse_status(sparse_status_t status) {

    if (status == SPARSE_STATUS_SUCCESS) {
        return;
    }
    else if (status == SPARSE_STATUS_NOT_INITIALIZED) {
        std::cout << "MKL sparse status: not initialized." << std::endl;
    }
    else if (status == SPARSE_STATUS_ALLOC_FAILED) {
        std::cout << "MKL sparse status: allocation failure." << std::endl;
    }
    else if (status == SPARSE_STATUS_INVALID_VALUE) {
        std::cout << "MKL sparse status: invalid value." << std::endl;
    }
    else if (status == SPARSE_STATUS_EXECUTION_FAILED) {
        std::cout << "MKL sparse status: execution failure." << std::endl;
    }
    else if (status == SPARSE_STATUS_INTERNAL_ERROR) {
        std::cout << "MKL sparse status: internal error." << std::endl;
    }
    else if (status == SPARSE_STATUS_NOT_SUPPORTED) {
        std::cout << "MKL sparse status: not supported." << std::endl;
    }
}