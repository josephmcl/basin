#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <functional>
#include <optional>

namespace timing {
	
	enum class architecture:std::size_t {

		i386 = 0, // __i386__",
		x86 = 1, // "__x86_64__",
		powerpc = 2, //"__powerpc__",
		arm = 3 // "__aarch64__"

	};
	
    class tsc_count {
    public:
        uint64_t value;
        tsc_count(uint64_t value): value{value} {};
    };

    auto operator -(tsc_count lhs, tsc_count rhs) -> double;

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

	tsc_count read();

    void init(void);


	/* This struct instruments functions, pass the parameter types 
	   into the template and plug in the pre, call, and post 
	   functions. Intended to clean up the mess of code at the top
	   level. 
	   
	   NOTE: We will have to extend this to measure things other than 
	         time.  */
	template <typename... arguments>
	struct instrument {

		using f = std::function<std::size_t(arguments&&... a)>;
		
		instrument(std::optional<f> pre, f call, std::optional<f> post){
			this->pre = pre;
			this->call = call;
			this->post = post;
		}

		std::optional<f> pre;
		f call;
		std::optional<f> post;
		
		tsc_count begin = read();   // TODO: default constructor for 
		tsc_count end = read();     //       this type.
		double elapsed;
		
		instrument &operator()(arguments&&... a){
			if (pre) (*pre)(a...);
			begin = read(); 
			call(a...);
			end = read();
			if (post) (*post)(a...);
			elapsed = end - begin;
		}
	};
	
}

