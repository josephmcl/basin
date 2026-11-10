#include "device_csr.cuh"
#include "cuda_check.cuh"

#include <stdexcept>
#include <string>

device_csr to_device(csr_t const &h) {
  if (h.r.size() < h.n + 1) {
    throw std::invalid_argument(
        "to_device: row offsets have " + std::to_string(h.r.size()) +
        " entries, expected n + 1 = " + std::to_string(h.n + 1));
  }

  device_csr d;
  d.n   = static_cast<int>(h.n);
  d.m   = static_cast<int>(h.m);
  d.nnz = static_cast<int>(h.v.size());

  CUDA_CHECK(cudaMalloc(&d.row, sizeof(int) * (h.n + 1)));
  CUDA_CHECK(cudaMemcpy(d.row, h.r.data(), sizeof(int) * (h.n + 1),
                        cudaMemcpyHostToDevice));

  // cudaMalloc(0) is legal, but skip it so empty matrices hold nullptr.
  if (d.nnz > 0) {
    CUDA_CHECK(cudaMalloc(&d.col, sizeof(int) * d.nnz));
    CUDA_CHECK(cudaMalloc(&d.val, sizeof(real_t) * d.nnz));
    CUDA_CHECK(cudaMemcpy(d.col, h.c.data(), sizeof(int) * d.nnz,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.val, h.v.data(), sizeof(real_t) * d.nnz,
                          cudaMemcpyHostToDevice));
  }
  return d;
}

void free_device(device_csr &d) {
  CUDA_CHECK(cudaFree(d.row));
  CUDA_CHECK(cudaFree(d.col));
  CUDA_CHECK(cudaFree(d.val));
  d = device_csr{};
}
