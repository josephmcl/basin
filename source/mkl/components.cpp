#include "components.h"


components::components(
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


void 
components::make_bs() {

  std::size_t const n2 = n * n;

  vector_t bs = { // hard code p = 2 for now.
    (3./2.) / span * (n - 1), 
        -2. / span * (n - 1),   
    (1./2.) / span * (n - 1)};
  
  /*
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
  */

  
  bsx = csr<double>{n2, n2};
  bsy = csr<double>{n2, n2};
  for (std::size_t i = 0; i < bs.size(); ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      bsx(bs[i], j, (i * n) + j);
      bsx(bs[i], n2 - 1 - j, n2 - 1 - (i * n) - j);

      bsy(bs[i], j * n, j * n + i);
      bsy(bs[i], (j * n) + n - 1, (j * n) + n - i - 1);
    }
  }  
}


void 
components::make_l() {

  std::size_t const n2 = n * n;

  ln = csr<double>{n, n2};
  ls = csr<double>{n, n2};
  le = csr<double>{n, n2};
  lw = csr<double>{n, n2};

  for (std::size_t i = 0; i != n; ++i) {
    ls(1., i, n * i);
    ln(1., i, n * i + (n - 1));
    le(1., i, n * (n - 1) + i);
    lw(1., i, i);
    //MatSetValue(ls, i, n * i,  1., ADD_VALUES);
    //MatSetValue(ln, i, n * i + (n - 1),  1., ADD_VALUES);
    //MatSetValue(le, i, n * (n - 1) + i, 1., ADD_VALUES);
    //MatSetValue(lw, i, i, 1., ADD_VALUES);
  }

  /*
  finalize<fw>(ln);
  finalize<fw>(ls);
  finalize<fw>(le);
  finalize<fw>(lw);*/
}

void 
components::make_h1() { 

  //MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &h1x);
  //MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &h1y);

  h1x = csr<double>{n, n};
  h1y = csr<double>{n, n};

  real_t const v = (1. / 2.) * spacing;
  // MatSetValue(h1x, 0, 0, v, ADD_VALUES);
  // MatSetValue(h1y, 0, 0, v, ADD_VALUES);
  h1x(v, 0, 0);
  h1y(v, 0, 0);
  h1v.push_back(v);
  for (std::size_t i = 1; i != n - 1; ++i) {
    //MatSetValue(h1x, i, i, spacing, ADD_VALUES);
    //MatSetValue(h1y, i, i, spacing, ADD_VALUES);
    h1x(v, i, i);
    h1y(v, i, i);
    h1v.push_back(spacing);
  }
  //MatSetValue(h1x, n - 1, n - 1, v, ADD_VALUES);
  //MatSetValue(h1y, n - 1, n - 1, v, ADD_VALUES);
  h1x(v, n - 1, n - 1);
  h1y(v, n - 1, n - 1);
  h1v.push_back(v);
  //finalize<fw>(h1x);
  //finalize<fw>(h1y);
}

void 
components::make_h() {

  std::size_t const n2 = n * n;
  //MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hx);
  //MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hy);
  hx = csr<double>{n2, n2};
  hy = csr<double>{n2, n2};

  real_t const v = (1. / 2.) * spacing;
  for (std::size_t i = 0; i != n; ++i) {
    //MatSetValue(hx, i, i, v, ADD_VALUES);
    //MatSetValue(hy, n * i, n * i, v, ADD_VALUES);
    hx(v, i, i);
    hy(v, n * i, n * i);
    for (std::size_t j = 1; j != n - 1; ++j) {
      //MatSetValue(hx, n * j + i, n * j + i, spacing, ADD_VALUES);
      //MatSetValue(hy, n * i + j, n * i + j, spacing, ADD_VALUES);
      hx(spacing, n * j + i, n * j + i);
      hy(spacing, n * i + j, n * i + j);
    }
    //MatSetValue(hx, n2 - n + i, n2 - n + i, v, ADD_VALUES);
    //MatSetValue(hy, n * (i + 1) - 1, n * (i + 1) - 1, v, ADD_VALUES);
    hx(v, n2 - n + i, n2 - n + i);
    hy(v, n * (i + 1) - 1, n * (i + 1) - 1);
  }
  //finalize<fw>(hx);
  //finalize<fw>(hy);
}

void  
components::make_hl() {

  std::size_t const n2 = n * n;
  // MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 1, nullptr, &hl);
  hl = csr<double>{n2, n2};

  real_t const v = (1. / 2.) * spacing;
  // MatSetValue(hl, 0, 0, v * v, ADD_VALUES);
  hl(v * v, 0, 0);
  for (std::size_t j = 1; j != n - 1; ++j) {
    //MatSetValue(hl, j, j, v * spacing, ADD_VALUES);
    hl(v * spacing, j, j);
  }
  // MatSetValue(hl, n - 1, n - 1, v * v, ADD_VALUES);
  hl(v * v, n - 1, n - 1);
  for (std::size_t i = 1; i != n - 1; ++i) {
    // MatSetValue(hl, n * i, n * i, v * spacing, ADD_VALUES);
    hl(v * spacing, n * i, n * i);
    for (std::size_t j = 1; j != n - 1; ++j) {
      // MatSetValue(hl, n * i + j, n * i + j, spacing * spacing, ADD_VALUES);
      hl(spacing * spacing, n * i + j, n * i + j);
    }
    // MatSetValue(hl, n * i + n - 1, n * i + n - 1, v * spacing, ADD_VALUES);
    hl(v * spacing, n * i + n - 1, n * i + n - 1);
  }
  // MatSetValue(hl, n * (n - 1), n * (n - 1), v * v, ADD_VALUES);
  hl(v * v, n * (n - 1), n * (n - 1));
  for (std::size_t j = 1; j != n - 1; ++j) {
    // MatSetValue(hl, n * (n - 1) + j, n * (n - 1) + j, v * spacing, ADD_VALUES);
    hl(v * spacing, n * (n - 1) + j, n * (n - 1) + j);
  }
  // MatSetValue(hl, n * n - 1, n * n - 1, v * v, ADD_VALUES);
  hl(v * v, n * n - 1, n * n - 1);
  
  //finalize<fw>(hl);

}

void 
components::make_d2() {

  std::size_t const n2 = n * n;
  //MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 3, nullptr, &d2x);
  //MatCreateSeqAIJ(PETSC_COMM_SELF, n2, n2, 3, nullptr, &d2y);
  d2x = csr<double>{n2, n2};
  d2y = csr<double>{n2, n2};

  real_t spacing2 = (spacing * spacing) / ((n - 1) * (n - 1));
  real_t coeff = 2.;

  std::size_t jn, in;
  for (std::size_t j = 0; j != n; ++j) {  
    jn = j * n;

    //MatSetValue(d2x, j + 0, j + 0,  1. / spacing2 * h1v[0] * coeff, ADD_VALUES);
    //MatSetValue(d2x, j + 0, j + 1 * n, -2. / spacing2 * h1v[0] * coeff, ADD_VALUES);
    //MatSetValue(d2x, j + 0, j + 2 * n,  1. / spacing2 * h1v[0] * coeff, ADD_VALUES); 
    //MatSetValue(d2x, j + n2 - n, j + n2 - 3 * n,  1. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES);
    //MatSetValue(d2x, j + n2 - n, j + n2 - 2 * n, -2. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES);
    //MatSetValue(d2x, j + n2 - n, j + n2 - 1 * n,  1. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES); 
    d2x(1. / spacing2 * h1v[0] * coeff, j + 0, j + 0);
    d2x(-2. / spacing2 * h1v[0] * coeff, j + 0, j + 1 * n);
    d2x(1. / spacing2 * h1v[0] * coeff, j + 0, j + 2 * n);
    d2x(1. / spacing2 * h1v[n - 1] * coeff, j + n2 - n, j + n2 - 3 * n);
    d2x(-2. / spacing2 * h1v[n - 1] * coeff, j + n2 - n, j + n2 - 2 * n);
    d2x(1. / spacing2 * h1v[n - 1] * coeff, j + n2 - n, j + n2 - 1 * n);

    //MatSetValue(d2y, jn + 0, jn + 0,  1. / spacing2 * h1v[0] * coeff, ADD_VALUES);
    //MatSetValue(d2y, jn + 0, jn + 1, -2. / spacing2 * h1v[0] * coeff, ADD_VALUES);
    //MatSetValue(d2y, jn + 0, jn + 2,  1. / spacing2 * h1v[0] * coeff, ADD_VALUES); 
    //MatSetValue(d2y, jn + n - 1, jn + n - 3,  1. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES);
    //MatSetValue(d2y, jn + n - 1, jn + n - 2, -2. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES);
    //MatSetValue(d2y, jn + n - 1, jn + n - 1,  1. / spacing2 * h1v[n - 1] * coeff, ADD_VALUES); 
    d2y(1. / spacing2 * h1v[0] * coeff, jn + 0, jn + 0);
    d2y(-2. / spacing2 * h1v[0] * coeff, jn + 0, jn + 1);
    d2y(1. / spacing2 * h1v[0] * coeff, jn + 0, jn + 2);
    d2y(1. / spacing2 * h1v[n - 1] * coeff, jn + n - 1, jn + n - 3);
    d2y(-2. / spacing2 * h1v[n - 1] * coeff, jn + n - 1, jn + n - 2);
    d2y(1. / spacing2 * h1v[n - 1] * coeff, jn + n - 1, jn + n - 1);

    for (std::size_t i = 1; i != n - 1; ++i) {  
        in = i * n;

        //MatSetValue(d2x, in + j, in + j - n,  1. / spacing2 * h1v[i], ADD_VALUES);
        //MatSetValue(d2x, in + j, in + j,     -2. / spacing2 * h1v[i], ADD_VALUES);
        //MatSetValue(d2x, in + j, in + j + n,  1. / spacing2 * h1v[i], ADD_VALUES);  
        d2x(1. / spacing2 * h1v[i],  in + j, in + j - n);
        d2x(-2. / spacing2 * h1v[i], in + j, in + j);
        d2x(1. / spacing2 * h1v[i],  in + j, in + j + n);

        //MatSetValue(d2y, jn + i, jn + i - 1,  1. / spacing2 * h1v[i], ADD_VALUES);
        //MatSetValue(d2y, jn + i, jn + i,     -2. / spacing2 * h1v[i], ADD_VALUES);
        //MatSetValue(d2y, jn + i, jn + i + 1,  1. / spacing2 * h1v[i], ADD_VALUES);  
        d2y(1. / spacing2 * h1v[i],  jn + i, jn + i - 1);
        d2y(-2. / spacing2 * h1v[i], jn + i, jn + i);
        d2y(1. / spacing2 * h1v[i],  jn + i, jn + i + 1);
    }
  }

  //finalize<fw>(d2x);
  //finalize<fw>(d2y);
}
