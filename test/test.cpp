#include "test_ranges.h"
#include "test_domain.h"
#include "test_physical.h"
#include "test_operators.h"

int main(int argc, char **argv) {

    test::test(TEST_HI_OPERATOR);

    test::test(TEST_PHYSICAL_PARAMS_BFACE);

    test::test(TEST_STATIC_ASSERT_LINRANGE);
    test::test(TEST_MIN_ELEMENT);

    test::test(TEST_TEST_INPUT_METRICS_COORD_X);
    test::test(TEST_TEST_INPUT_METRICS_COORD_Y);
    test::test(TEST_TEST_INPUT_METRICS_J);
    test::test(TEST_TEST_INPUT_METRICS_JI);
    test::test(TEST_TEST_INPUT_METRICS_RX);
    test::test(TEST_TEST_INPUT_METRICS_SY);
    test::test(TEST_TEST_INPUT_METRICS_CRR);
    test::test(TEST_TEST_INPUT_METRICS_CSS);      
    test::test(TEST_TEST_INPUT_METRICS_ETA);
    test::test(TEST_TEST_INPUT_METRICS_MU_FACE_2);
    test::test(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_1);
    test::test(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_2);
    test::test(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_3);
    test::test(TEST_TEST_INPUT_METRICS_SURFACE_JACOBIAN_FACE_4);

}