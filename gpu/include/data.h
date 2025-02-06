#include "definitions.h"


template<typename t>
struct data {
    static constexpr t π = std::numbers::pi_v<t>;
    static constexpr auto w = [](t x, t y){
        return std::sin(π * x + π * y);};
    static constexpr auto e = [](t x, t y){
        return std::sin(π * x + π * y);};
    static constexpr auto s = [](t x, t y){
        return -π * std::cos(π * x + π * y);};
    static constexpr auto n = [](t x, t y){
        return  π * std::cos(π * x + π * y);};
    static constexpr auto source_function = [](t x, t y) {
        return 2. * (-π * π * std::sin(π * x + π * y));};
};