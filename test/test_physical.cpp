#include "test_physical.h"

static std::string test_input_path = "/Users/josephmcl/phd/thrase/Basin_Simulations/";

static domain::data constexpr D = {
    .length      = 24.   ,
    .basin_depth =   4.  ,
    .r̂           =   0.75, 
    .l           =   0.05,
    .μ_out       =  36.  ,
    .μ_in        =   8.  ,
    .ρ_out       =   2.8 ,
    .ρ_in        =   2.  
};

DEFINE_PROFILE(TEST_PHYSICAL_PARAMS_BFACE)
    
    std::ifstream f;
    f.open(test_input_path + "test/physical_params.b.csv", 
        std::ifstream::in);
    auto [columns, rows, test_data] = test::read_csv<double>(f);

    ASSERT(columns * rows == test_data->size(), 
        "Sanity check failed, input matrix x dims N × M ≠ |x|.");

    auto met = domain::metrics<D>(400., 400.);
    auto y = met.y();

    auto [δ, g, v, RP] = physical::fault_params(y);

    auto [σn, a, faceb, Dc, f0, V0, τ_inf, Vp] = RP;

    ASSERT(rows == faceb.size(), "faceb linrange is a "
        "different size than test input b.");

    for (std::size_t i = 0; i < rows; ++i) {
        auto facebi = *faceb[i];
        auto test = test_data->at(i);
        auto[passed, message] = test::approx(facebi, test);
        ASSERT(passed, message);
    }

    delete test_data;

ENDDEF_PROFILE

