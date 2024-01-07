
#include "poisson_2d.h"


void main_2d();

int main (int argc, char *argv[]) {

  /* Initialize libraries. */
  if (!infrastructure::initialize(argc, argv)) exit(-1);; { 

    // timing::init();
    // poisson_1d_hybrid_convergence_test(199, 3);

    main_2d();

  /* Cleanup libraries */
  } infrastructure::cleanup(); return 0;
}

void main_2d() {
  sbp_sat::x2::hybridized_poisson(16, 3);
}
