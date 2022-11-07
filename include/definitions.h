#pragma once

#include <vector>

namespace type {

    using real_t = long double;

    template <typename t> 
    using vv = std::vector<std::vector<t>>;

    template <typename t> 
    auto make_vv = [](std::size_t n, std::size_t m) {
        return vv<t>(n, std::vector<t>(m));};

}