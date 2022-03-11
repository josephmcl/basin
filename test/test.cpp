#include "test_ranges.h"
    
#include "test_domain.h"

int main(int argc, char **argv) {

    test::test(TEST_STATIC_ASSERT_LINRANGE);


    test::test(TEST_TEST_INPUT_METRICS_COORD_X);
    test::test(TEST_TEST_INPUT_METRICS_COORD_Y);
    test::test(TEST_TEST_INPUT_METRICS_J);
    test::test(TEST_TEST_INPUT_METRICS_JI);
    test::test(TEST_TEST_INPUT_METRICS_RX);
    test::test(TEST_TEST_INPUT_METRICS_SY);
    test::test(TEST_TEST_INPUT_METRICS_CRR);
    test::test(TEST_TEST_INPUT_METRICS_CSS);
}