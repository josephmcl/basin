/* profile.h 
Core profiling and testing infrastructure for small C++ projects. */
#pragma once

#include <iomanip>
#include <iostream>
#include <source_location>
#include <string>
#include <sstream>
#include <tuple>
#include <cmath>
#include <cstdlib>

/* High-level usage: 

test.h:

    DECLARE_PROFILE(MY_TEST);

test.cpp:

    #include "test.h"
    DEFINE_PROFILE(MY_TEST)
        ASSERT(1 != 1, "Inanity.");
    ENDDEF_PROFILE

test_main.cpp:

    #include "test.h"
    int main() { test::test(MY_TEST); return 0; }

*/
#define DECLARE_PROFILE(x) std::tuple<std::source_location, bool> \
    const x(void);
#define DEFINE_PROFILE(x) std::tuple<std::source_location, bool>  \
    const x(void) { auto _res_ = std::source_location::current(); \
    bool _passed_ = true;
#define ENDDEF_PROFILE return std::make_tuple(_res_, _passed_); }

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { if (! (condition)) { std::cerr << std::setprecision(14) \
    << std::scientific << "\033[31mAssertion failed\033[0m in "  \
    << __FILE__ <<  " line " << __LINE__ << std::endl << message \
    << std::endl;        \
    return std::make_tuple(_res_, false); } } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace test {

    static bool picky = false;

    template<typename T>
    std::tuple<bool, std::string> approx(T a, T b, T ϵ=1e-9) {
        std::stringstream ss;
        ss << std::setprecision(14) << std::scientific 
           << a << " != " << b << " (ϵ=" << std::setprecision(0) 
           << ϵ << ")";
        return make_tuple(std::abs(a - b) < ϵ, ss.str()); }

    template<typename T>
    void test(T f) {
        // TODO: Add TSC timer. 
        auto [sl, success] = f(); 

        std::stringstream ss(sl.function_name());
        std::string s; ss >> s; ss >> s; ss >> s; ss >> s; ss >> s; 
        if (success)
            std::cout << "[\033[92m ✓ \033[0m] " << s << " passed." 
                      << std::endl;
        else
            std::cout << "[\033[31m ✗ \033[0m] " << s << " failed." 
                      << std::endl;

        if (picky && !success) exit(-1);
    }

}