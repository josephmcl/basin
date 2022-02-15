#include "domain.h"

#include "ranges.h"
#include <ranges>


#include <vector>
#include <algorithm>

int main (int argc, char *argv[]) {

    // Domain.transforms_e()

    /*
    int constexpr length = 24;
    double constexpr r̂ = 0.75;
    double constexpr l = 0.05;

    std::vector<double> r = {};
    std::vector<double> s = {};
    */

    // domain::transform<length, r̂, l>(r, s); 

    using namespace numerics;

    //for (auto i: range(5., 0.1, 9.))
    //    std::cout << i << ' ';    // prints 5 6 7 8 

    // for (auto i: linrange(-1., 1., 10)) {
    //    std::cout << i << ' ';    // prints 5 6 7 8 
    // }

    // auto even = [](double d) { return 0 == d % 2; };
    
    
    static_assert(std::input_iterator<linrange<double>::iterator>);
    static_assert(std::ranges::input_range<linrange<double>>);

    static_assert(std::forward_iterator<linrange<double>::iterator>);
    static_assert(std::ranges::forward_range<linrange<double>>);
    static_assert(std::ranges::sized_range<linrange<double>>);

    auto lr = linrange(0., 1., 13);

    static_assert(std::is_same_v<decltype(lr.begin() == lr.end()), bool>);
    static_assert(std::is_same_v<decltype(lr.end() == lr.begin()), bool>);
    static_assert(std::ranges::range<linrange<double> &>, "It is not a range");
    // static_assert(std::ranges::input_range<linrange<double>>);
    // static_assert(std::ranges::forward_range<linrange<double>>);

    // "pipe" syntax of composing the views:

    auto square = [](double i) { return i * i; };
    auto plusone = [](double i) { return i + 1; };

    for (auto i : linrange(1., 2., 15) | plusone | square ) {
        std::cout << " " << i << std::endl;
    }
    std::cout << std::endl;
    


    return 0;
}