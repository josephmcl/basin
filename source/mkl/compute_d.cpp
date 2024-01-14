#include "compute_d.h"


void compute_d(
  sparse_matrix_t          *D, 
  components               &sbp, 
  vv<std::size_t> const &interfaces) {

  /* D = [H1y * 2 * τ,                              ]
         [             H1y * 2 * τ,                 ]
         [                        , ...,            ]
         [                               H1y * 2 * τ] 

  NOTE: either H1x or H1y for corresponding interface 
        orientations (NS or EW). For now we use h1v because 
        h1x and h1y are identical in this current config. */

  std::size_t n_interfaces = 0;
  for (std::size_t row = 0; row != interfaces.size(); ++row) {
    for (std::size_t col = 0; col != interfaces.size(); ++col) {
      std::size_t interface = interfaces[row][col];
      if (interface != 0) n_interfaces += 1;
    }
  }

  std::size_t size = sbp.n * n_interfaces;
  // MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 1, nullptr, &D);
  csr<double> d{size, size};
  std::size_t interface, index;
  real_t value;
  for (std::size_t row = 0; row != interfaces.size(); ++row) {
    for (std::size_t col = 0; col != interfaces.size(); ++col) {
      interface = interfaces[row][col];
      if (interface != 0 and row == col - 1) {  // NS 
        index = (interface - 1) * sbp.n; 
        for (std::size_t k = 0; k != sbp.n; ++k) {
          value = sbp.h1v[k] * 2. * sbp.τ;
          // MatSetValue(D, index + k,  index + k, value, ADD_VALUES);
          // std::cout << value << " " << index + k << std::endl;
          d(value, index + k,  index + k);
        }
      }
      else if (interface != 0) {  // EW
        index = (interface - 1) * sbp.n; 
        for (std::size_t k = 0; k != sbp.n; ++k) {
          value = sbp.h1v[k] * 2. * sbp.τ;
          // std::cout << value << " " << index + k << std::endl;
          // MatSetValue(D, index + k,  index + k, value, ADD_VALUES);
          d(value, index + k,  index + k);
        }
      }
    }
  }
  d.mkl(D);
  // finalize<fw>(D);
}