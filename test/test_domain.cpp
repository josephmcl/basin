#include "test_domain.h"


static std::string test_input_path = "/Users/josephmcl/phd/thrase/Basin_Simulations/";

double constexpr test_domain_N  = 400.;

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

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_COORD_X)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.coord.x.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(), 
        "Sanity check failed, input matrix x dims N × M ≠ |x|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    auto x = met.x();

    ASSERT(std::get<1>(m) == x.size(), "x linrange is a different "
        "size than test input x.");

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto xi = *x[i];
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(xi, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE
/*                                           */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_COORD_Y)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.coord.y.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix y dims N × M ≠ |y|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    auto y = met.y();
    
    ASSERT(std::get<1>(m) == y.size(), "y linrange is a different "
        "size than test input y.");

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto yj = *y[j];
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(yj, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE
/*                                     */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_J)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.j.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix j dims N × M ≠ |j|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);


    auto J = met.j();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto jij = J(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(jij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_JI)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.ji.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix j dims N × M ≠ |j|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);

    
    auto JI = met.ji();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto jiij = JI(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(jiij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_CRR)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.crr.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);

    
    auto crr = met.crr();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto crrij = crr(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(crrij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE  

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_CSS)
    
    std::ifstream f;
    f.open(test_input_path + "test/metrics.css.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);

    auto css = met.css();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto crrij = css(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(crrij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE  

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_RX)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.rx.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix ji dimsz N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto rx = met.rx();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto rxij = rx(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(rxij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE  

DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_SY)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.sy.csv", 
        std::ifstream::in);
    auto m = test::read_csv<double>(f);

    ASSERT(std::get<0>(m) * std::get<1>(m) == std::get<2>(m)->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto sy = met.sy();

    for (std::size_t i = 0; i < std::get<0>(m); ++i) {
        for (std::size_t j = 0; j < std::get<1>(m); ++j) {
            auto syij = sy(i,j);
            auto test = std::get<2>(m)->at((std::get<0>(m) * i) + j);
            auto[passed, message] = test::approx(syij, test);
            ASSERT(passed, message);
        }    
    }

    delete std::get<2>(m);

ENDDEF_PROFILE  


DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_ETA)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.eta.csv", 
        std::ifstream::in);
    auto [N, M, test_eta] = test::read_csv<double>(f);

    ASSERT(N * M == test_eta->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto eta = met.η();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto etaij = eta(i);
            auto test = test_eta->at((M * i) + j);
            auto[passed, message] = test::approx(etaij, test);
            ASSERT(passed, message);
        }    
    }

    delete test_eta;

ENDDEF_PROFILE  


DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_MU_FACE_2)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.muf2.csv", 
        std::ifstream::in);
    auto [N, M, test_muf2] = test::read_csv<double>(f);

    ASSERT(N * M == test_muf2->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto muf2 = met.μf2();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto muf2i = muf2(i);
            auto test = test_muf2->at((M * i) + j);
            auto[passed, message] = test::approx(muf2i, test);
            ASSERT(passed, message);
        }    
    }

    delete test_muf2;

ENDDEF_PROFILE
/* */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_1)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.sj1.csv", 
        std::ifstream::in);
    auto [N, M, test_sj1] = test::read_csv<double>(f);

    ASSERT(N * M == test_sj1->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto sj1 = met.sj1();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto sj1i = sj1(i);
            auto test = test_sj1->at((M * i) + j);
            auto[passed, message] = test::approx(sj1i, test);
            ASSERT(passed, message);
        }    
    }

    delete test_sj1;

ENDDEF_PROFILE
/* */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_2)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.sj2.csv", 
        std::ifstream::in);
    auto [N, M, test_sj2] = test::read_csv<double>(f);

    ASSERT(N * M == test_sj2->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto sj2 = met.sj2();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto sj2i = sj2(i);
            auto test = test_sj2->at((M * i) + j);
            auto[passed, message] = test::approx(sj2i, test);
            ASSERT(passed, message);
        }    
    }

    delete test_sj2;

ENDDEF_PROFILE
/* */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_3)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.sj2.csv", 
        std::ifstream::in);
    auto [N, M, test_sj3] = test::read_csv<double>(f);

    ASSERT(N * M == test_sj3->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto sj3 = met.sj3();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto sj3i = sj3(i);
            auto test = test_sj3->at((M * i) + j);
            auto[passed, message] = test::approx(sj3i, test);
            ASSERT(passed, message);
        }    
    }

    delete test_sj3;

ENDDEF_PROFILE
/* */
DEFINE_PROFILE(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_4)

    std::ifstream f;
    f.open(test_input_path + "test/metrics.sj4.csv", 
        std::ifstream::in);
    auto [N, M, test_sj4] = test::read_csv<double>(f);

    ASSERT(N * M == test_sj4->size(),
        "Sanity check failed, input matrix ji dims N × M ≠ |ji|.");

    auto met = domain::metrics<D>(test_domain_N, test_domain_N);
    
    auto sj4 = met.sj4();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            auto sj4i = sj4(i);
            auto test = test_sj4->at((M * i) + j);
            auto[passed, message] = test::approx(sj4i, test);
            ASSERT(passed, message);
        }    
    }

    delete test_sj4;

ENDDEF_PROFILE