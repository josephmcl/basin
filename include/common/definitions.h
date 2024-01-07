#pragma once 

#include <tuple>

namespace type {
    using real_t = long double;
};

template <typename T>
using type_p = std::tuple<std::size_t, T*>;
using real_p = type_p<type::real_t>;