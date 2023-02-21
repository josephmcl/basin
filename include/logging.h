#pragma once 
#include <iostream>
#include <fstream>
namespace logging {
	// bool constexpr on = true; 
	extern std::ostream &out;
};

#define logging_on true