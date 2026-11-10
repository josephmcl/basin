#include "compute_f_symbols.h"

void compute_f_symbols(
  vv<std::size_t>       &F_symbols,
  vv<std::size_t>       &FT_symbols, 
  vv<std::size_t> const &interfaces,
  components      const &sbp) {

  constexpr std::size_t w = 1; constexpr std::size_t e = 2; 
  constexpr std::size_t s = 3; constexpr std::size_t n = 4; 

  for (std::size_t row = 0; row != sbp.n_blocks; ++row) {
    for (std::size_t col = 0; col != sbp.n_blocks; ++col) {
      std::size_t interface = interfaces[row][col];
      if (interface != 0 and row == col - 1) {  
        F_symbols[row][interface - 1] = n;
        F_symbols[col][interface - 1] = s;
        FT_symbols[interface - 1][row] = n;
        FT_symbols[interface - 1][col] = s;
      }
      else if (interface != 0 and row == col - sbp.n_blocks_dim) {  // BIG OOPS
        F_symbols[row][interface - 1] = e;
        F_symbols[col][interface - 1] = w;
        FT_symbols[interface - 1][row] = e;
        FT_symbols[interface - 1][col] = w;
      }
    }
  }
}