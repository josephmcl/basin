#pragma once 

#include <vector>
#include <tuple>

namespace type {
    using real_t = double;
};

template <typename T>
using type_p = std::tuple<std::size_t, T*>;
using real_p = type_p<type::real_t>;

template <typename T>
using vv = std::vector<std::vector<T>>;