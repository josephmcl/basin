#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>

namespace timing {
	
	enum class architecture:std::size_t {

		i386 = 0, // __i386__",
		x86 = 1, // "__x86_64__",
		powerpc = 2, //"__powerpc__",
		arm = 3 // "__aarch64__"

	};

	static double tsc_timer_ticks_per_second = 0.0;
	
    class tsc_count {
    public:
        uint64_t value;
        tsc_count(uint64_t value): value{value} {};
    };

    auto operator -(tsc_count lhs, tsc_count rhs) -> double {
        return static_cast<double>(lhs.value - rhs.value) 
            / tsc_timer_ticks_per_second;
    }

	template<architecture arch> uint64_t read_() {

		if constexpr (arch == architecture::arm) {
			uint64_t res;
			__asm__ __volatile__("mrs %0, cntvct_el0" : "=r"(res));
			return res;
		}
		else if constexpr (arch == architecture::x86) {
			uint32_t hi, lo;
		    __asm__ __volatile__("rdtsc":"=a"(lo), "=d"(hi));
		    return ((uint64_t) lo) | (((uint64_t) hi) << 32);
		}
		else if constexpr (arch == architecture::x86) {
			uint32_t hi, lo;
		    __asm__ __volatile__("rdtsc":"=a"(lo), "=d"(hi));
		    return ((uint64_t) lo) | (((uint64_t) hi) << 32);
		}
		else if constexpr (arch == architecture::powerpc) {
			uint64_t result = 0;
		    uint64_t upper, lower, tmp;
		    __asm__ __volatile__("0:                  \n\tmftbu   %0     "
		    	"      \n\tmftb    %1           \n\tmftbu   %2           "
		    	"\n\tcmpw    %2,%0        \n\tbne     0b         \n" : 
		    	"=r"(upper), "=r"(lower), "=r"(tmp));
		    result = upper;
		    result = result << 32;
		    result = result | lower;
		    return result;
		}
	}	

	#if defined(__aarch64__)
		auto read = [](){return tsc_count(read_<architecture::arm>());};
	#elif defined(__x86_64__)
		auto read = [](){return tsc_count(read_<architecture::x86>());};
	#endif

    void init(void) {
        auto start = read();
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        auto stop = read();
        tsc_timer_ticks_per_second = static_cast<double>(
            stop.value - start.value);
        std::cout << "tsc_timer_ticks_per_second: " 
            << tsc_timer_ticks_per_second << std::endl;
    }
}

