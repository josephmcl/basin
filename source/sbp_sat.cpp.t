export module sbp_sat;

/* stdlib imports */
export import <cstddef>;

export struct d2_1d {

  d2_1d() = default;
  /*  Construct the 1-D operator given a domain size, an order of 
      accuracy, and left and right domain boundaries. */
  d2_1d(
    std::size_t const size, 
    std::size_t const order, 
    long double const left, 
    long double const right)
  :
    size(size), 
    order(order), 
    left(left), 
    right(right), 
    spacing((right - left) / static_cast<long double>(size - 1)) { }

  std::size_t const size, order;
  long double const left, right, spacing;


}; /* struct d2_1d */  
