#include "timing.h"

using namespace timing;

#if defined(__aarch64__)
    timing::tsc_count timing::read(){return tsc_count(read_<architecture::arm>());};
#elif defined(__x86_64__)
    timing::tsc_count timing::read(){return tsc_count(read_<architecture::x86>());};
#endif

double tsc_timer_ticks_per_second = 0.0;

auto timing::operator -(tsc_count lhs, tsc_count rhs) -> double {
    return static_cast<double>(lhs.value - rhs.value) 
        / tsc_timer_ticks_per_second;
}

void timing::init(void) {
    
    auto start = read();
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    auto stop = read();
    tsc_timer_ticks_per_second = static_cast<double>(
        stop.value - start.value);
    std::cout << "tsc_timer_ticks_per_second: " 
        << tsc_timer_ticks_per_second << std::endl;
}