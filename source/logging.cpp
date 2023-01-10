#include "logging.h"

#if logging_on
	std::ostream &logging::out = std::cout;
#else 
	std::ofstream dev_null("/dev/null");
	std::ostream &logging::out = dev_null;
#endif
