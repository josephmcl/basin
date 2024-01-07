#pragma once 
#include <cstddef>
namespace infrastructure {
  enum class error : std::size_t { 
    nil                      = 0,
    mkl_failure              = 1,
    openmp_failure           = 2};
  bool operator !(error e);
  error initialize(int c, char **v);
  void cleanup();
} /* namespace infrastructure */
