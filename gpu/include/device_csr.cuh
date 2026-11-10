#pragma once

#include "csr.h"
#include "definitions.h"

// A CSR matrix whose arrays live in GPU memory (cudaMalloc).
//
// This is a plain handle: copying it copies the pointers, not the data.
// Create with to_device(), release with free_device().
struct device_csr {
  int n   = 0;  // rows
  int m   = 0;  // columns
  int nnz = 0;
  int    *row = nullptr;  // n + 1 row offsets
  int    *col = nullptr;  // nnz column indices
  real_t *val = nullptr;  // nnz values
};

// Allocate device arrays and copy a host CSR matrix into them.
device_csr to_device(csr_t const &h);

// Free the device arrays and reset the handle.
void free_device(device_csr &d);
