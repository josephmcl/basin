export module sbp_sat;

/* stdlib imports */
export import <cstddef>;
export import <vector>;

#include "petscsys.h"
#include "petscmat.h" 

export struct d2_1d {

  using dat_t = long double;
  using vec_t = std::vector<dat_t>;

  using petsc_matrix = Mat;

  d2_1d() = default;
  /*  Construct the 1-D operator given a domain size, an order of 
      accuracy, and left and right domain boundaries. */
  d2_1d(
    std::size_t const size, 
    std::size_t const order, 
    dat_t const left, 
    dat_t const right)
  :
    size(size), 
    order(order), 
    left(left), 
    right(right), 
    spacing((right - left) / static_cast<dat_t>(size - 1)) { }

  std::size_t const size, order;
  dat_t const left, right, spacing;

  struct iterator { 
    std::size_t const size, current;
    std::size_t const column_start = 0;
    vec_t const row = {};
  };

  iterator begin() const { return iterator{size, 0};    }
  iterator end()   const { return iterator{size, size}; }

  void write_petsc(petsc_matrix *pm) {
    for (auto [column_start, row_data] : rows) {

    }
  }  

}; /* struct d2_1d */  
