
#include<vector>

#include "definitions.h"
#include "components.h"

void compute_f_symbols(
  vv<std::size_t>       &F_symbols,
  vv<std::size_t>       &FT_symbols, 
  vv<std::size_t> const &interfaces,
  components      const &sbp);