#pragma once

#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

// Abort with file/line and the failing call if a CUDA runtime call fails.
// Kernel launches don't return an error: follow them with
// CUDA_CHECK(cudaGetLastError()).
#define CUDA_CHECK(call)                                                      \
  do {                                                                        \
    cudaError_t err_ = (call);                                                \
    if (err_ != cudaSuccess) {                                                \
      std::fprintf(stderr, "CUDA error \"%s\" at %s:%d\n  in: %s\n",          \
                   cudaGetErrorString(err_), __FILE__, __LINE__, #call);      \
      std::exit(EXIT_FAILURE);                                                \
    }                                                                         \
  } while (0)
