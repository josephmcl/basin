#include "test_operators.h"

static std::string test_data_path = "./test_data/";

DEFINE_PROFILE(TEST_HI_OPERATOR)

    using namespace type;
    using namespace numerical;

    std::ifstream f;
    f.open(test_data_path + "HI.csv", std::ifstream::in);
    auto [rows, columns, test_data] = test::read_csv<real_t>(f);

    ASSERT(rows % 2 == 0, "H file not loaded correctly. Corrupted?");

    operators o = {};

    std::size_t offset = 0;
    for (std::size_t row = 0; row < rows ; row += 2) {

        // Extract test parameters and data.
        auto size  = static_cast<std::size_t>(test_data->at(offset));
        auto order = static_cast<std::size_t>(test_data->at(offset + 1));
        auto left  = test_data->at(offset + 2);
        auto right = test_data->at(offset + 3);
        auto iH_test = test_data->begin() + offset + 4;
        offset += 5 + size;
        
        // Generate H Inverse matrix.
        auto H = o.H_inverse(size + 1, order, left, right);

        // Test corresponding values of the matrix entries. 
        auto iH = H.begin();
        while(iH != H.end()) {
            auto[passed, message] = test::approx(*iH, *iH_test);
            ASSERT(passed, message);
            *iH++; iH_test++;
        }
    }
    delete test_data;

ENDDEF_PROFILE
