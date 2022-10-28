#include "components.h"

sbp_sat::x2::components::
components(
  std::size_t const points,
  real_t      const span,
  std::size_t const accuracy)
: n(points), span(span), accuracy(accuracy), spacing(span/(points - 1)) {

  make_bs();
  make_l();
  make_h1();
  make_h();
  make_hl();
  make_d2();

}

sbp_sat::x2::components::
~components() {

  destroy<fw>(hl);
  destroy<fw>(hx);
  destroy<fw>(hy);
  destroy<fw>(h1x);
  destroy<fw>(h1y);
  destroy<fw>(bsx);
  destroy<fw>(bsy);
  destroy<fw>(ln);
  destroy<fw>(ls);
  destroy<fw>(le);
  destroy<fw>(lw);
  destroy<fw>(d2x);
  destroy<fw>(d2y);
}

void sbp_sat::x2::components::
make_bs() {

  std::size_t const n2 = n * n;

  vector_t bs = { // hard code p = 2 for now.
    (3./2.) / span * (n - 1), 
        -2. / span * (n - 1),   
    (1./2.) / span * (n - 1)};
  
  MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 3, nullptr, &bsx);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 3, nullptr, &bsy);
  
  for (std::size_t i = 0; i < bs.size(); ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      set_matrix_value<fw>(bsx, j, (i * n) + j, bs[i]);
      set_matrix_value<fw>(bsx, n2 - 1 - j, n2 - 1 - (i * n) - j, bs[i]);

      set_matrix_value<fw>(bsy, j * n, j * n + i, bs[i]);
      set_matrix_value<fw>(bsy, (j * n) + n - 1, (j * n) + n - i - 1, bs[i]);
    }
  }  
  
  finalize<fw>(bsx);
  finalize<fw>(bsy);
}


void sbp_sat::x2::components::
make_l() {

  std::size_t const n2 = n * n;

  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n2, 1, nullptr, &ln);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n2, 1, nullptr, &ls);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n2, 1, nullptr, &le);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n2, 1, nullptr, &lw);

  for (std::size_t i = 0; i != n; ++i) {
    MatSetValue(ln, i, n * i,  1., ADD_VALUES);
    MatSetValue(ls, i, n * i + (n - 1),  1., ADD_VALUES);
    MatSetValue(le, i, n * (n - 1) + i, 1., ADD_VALUES);
    MatSetValue(lw, i, i, 1., ADD_VALUES);
  }

  finalize<fw>(ln);
  finalize<fw>(ls);
  finalize<fw>(le);
  finalize<fw>(lw);
}

void sbp_sat::x2::components::
 make_h1() { 

  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &h1x);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &h1y);

  real_t const v = (1. / 2.) * spacing;
  MatSetValue(h1x, 0, 0, v, ADD_VALUES);
  MatSetValue(h1y, 0, 0, v, ADD_VALUES);
  h1v.push_back(v);
  for (std::size_t i = 1; i != n - 1; ++i) {
    MatSetValue(h1x, i, i, spacing, ADD_VALUES);
    MatSetValue(h1y, i, i, spacing, ADD_VALUES);
    h1v.push_back(spacing);
  }
  MatSetValue(h1x, n - 1, n - 1, v, ADD_VALUES);
  MatSetValue(h1y, n - 1, n - 1, v, ADD_VALUES);
  h1v.push_back(v);

  finalize<fw>(h1x);
  finalize<fw>(h1y);
}

void sbp_sat::x2::components::
make_h() {

  std::size_t const n2 = n * n;
  MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hx);
  MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hy);

  real_t const v = (1. / 2.) * spacing;
  for (std::size_t i = 0; i != n; ++i) {
    MatSetValue(hx, i, i, v, ADD_VALUES);
    MatSetValue(hy, n * i, n * i, v, ADD_VALUES);
    for (std::size_t j = 1; j != n - 1; ++j) {
      MatSetValue(hx, n * j + i, n * j + i, spacing, ADD_VALUES);
      MatSetValue(hy, n * i + j, n * i + j, spacing, ADD_VALUES);
    }
    MatSetValue(hx, n2 - n + i, n2 - n + i, v, ADD_VALUES);
    MatSetValue(hy, n * (i + 1) - 1, n * (i + 1) - 1, v, ADD_VALUES);
  }
  
  finalize<fw>(hx);
  finalize<fw>(hy);
}

void sbp_sat::x2::components::
make_hl() {

  std::size_t const n2 = n * n;
  MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hl);

  real_t const v = (1. / 2.) * spacing;
  MatSetValue(hl, 0, 0, v * v, ADD_VALUES);
  for (std::size_t j = 1; j != n - 1; ++j) {
    MatSetValue(hl, j, j, v * spacing, ADD_VALUES);
  }
  MatSetValue(hl, n - 1, n - 1, v * v, ADD_VALUES);
    for (std::size_t i = 1; i != n - 1; ++i) {
    MatSetValue(hl, n * i, n * i, v * spacing, ADD_VALUES);
    for (std::size_t j = 1; j != n - 1; ++j) {
      MatSetValue(hl, n * i + j, n * i + j, spacing * spacing, ADD_VALUES);
    }
    MatSetValue(hl, n * i + n - 1, n * i + n - 1, v * spacing, ADD_VALUES);
  }
  MatSetValue(hl, n * (n - 1), n * (n - 1), v * v, ADD_VALUES);
  for (std::size_t j = 1; j != n - 1; ++j) {
    MatSetValue(hl, n * (n - 1) + j, n * (n - 1) + j, v * spacing, ADD_VALUES);
  }
  MatSetValue(hl, n * n - 1, n * n - 1, v * v, ADD_VALUES);
  
  finalize<fw>(hl);
}

void sbp_sat::x2::components:: 
make_d2() {

  /*
  auto h = numerical::operators::H(n, 2, 0, span);

  for (auto e :h) std::cout << e << " ";
  std::cout << std::endl;

  petsc_matrix d2;
  MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 3, nullptr, &d2);

  real_t spacing2 = (spacing * spacing) / ((n - 1) * (n - 1));
  real_t coeff = -2.;

  // Initialize the first local skew row. 
  MatSetValue(d2, 0, 0,  1. / spacing2 * h[0] * coeff, ADD_VALUES);
  MatSetValue(d2, 0, 1, -2. / spacing2 * h[0] * coeff, ADD_VALUES);
  MatSetValue(d2, 0, 2,  1. / spacing2 * h[0] * coeff, ADD_VALUES); 

  // Initialize the final local skew row. 
  MatSetValue(d2, n - 1, n - 3,  1. / spacing2 * h[n - 1] * coeff, ADD_VALUES);
  MatSetValue(d2, n - 1, n - 2, -2. / spacing2 * h[n - 1] * coeff, ADD_VALUES);
  MatSetValue(d2, n - 1, n - 1,  1. / spacing2 * h[n - 1] * coeff, ADD_VALUES); 

  // Initialize the interior local diagonal rows. 
  for (auto i = 1; i != n - 1; ++i) {  
    MatSetValue(d2, i, i - 1,  1. / spacing2 * h[i] * coeff, ADD_VALUES);
    MatSetValue(d2, i, i,     -2. / spacing2 * h[i] * coeff, ADD_VALUES);
    MatSetValue(d2, i, i + 1,  1. / spacing2 * h[i] * coeff, ADD_VALUES);  
  }

  finalize<fw>(d2);

  MatView(d2, PETSC_VIEWER_STDOUT_SELF);

  destroy<fw>(d2);
  */

}
