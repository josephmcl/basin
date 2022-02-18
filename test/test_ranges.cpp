#include "test_ranges.h"

/*  TEST_STATIC_ASSERT_LINRANGE 
    Run static asserts for the numerics::linrange object to ensure it 
    complies with several range concepts. */
DEFINE_PROFILE(TEST_STATIC_ASSERT_LINRANGE)
using namespace numerics; 
static_assert(std::input_iterator<linrange<double>::iterator>);
static_assert(std::ranges::input_range<linrange<double>>);
static_assert(std::forward_iterator<linrange<double>::iterator>);
static_assert(std::ranges::forward_range<linrange<double>>);
static_assert(std::ranges::sized_range<linrange<double>>);
auto lr = linrange(0., 1., 13);
static_assert(std::is_same_v<decltype(lr.begin() == lr.end()), bool>);
static_assert(std::is_same_v<decltype(lr.end() == lr.begin()), bool>);
static_assert(std::ranges::range<linrange<double> &>);
ENDDEF_PROFILE

