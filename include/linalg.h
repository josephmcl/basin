#pragma once
/* linalg.h is intended to encapsulate capabilities common to several 
   different linear algebra libraries, allowing different libraries to 
   to arbitrarily swapped out with ease. */


/* stl headers */
#include <type_traits>

/* petsc headers */
#include "petscsys.h"
#include "petscmat.h" 
#include "petscviewer.h"
#include "petscvec.h"
#include "petscksp.h"



namespace linalg {

  /* framework enumerates the supported linear algebra paradigms. For 
     now this only supports PETSc. */
  enum class framework : std::size_t {
    petsc = 0
  };

  /* STL type aliases. */
  using st = std::size_t;

  /* PETSc type aliases. */
  using petsc_matrix = Mat;

  /* Evaluate the back-end linalg matrix type at compile time. Idea 
     robbed from stackoverflow.com/a/61177390 */
  template<framework f>
  constexpr auto pick_matrix_type() {
    if constexpr (f == framework::petsc) {
        return std::type_identity<petsc_matrix>{};
    } 
  }

  /* The compile time matrix type template. */
  template<framework f>
  using matrix = typename decltype(pick_matrix_type<f>())::type;

  /* Initialize a local sparse matrix. */
  template<framework f>
  inline void make_local_sparse_matrix(matrix<f> &m, st rows, st cols, st nz) {
    if constexpr (f == framework::petsc)
      MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, nz, nullptr, &m);
  }

  /* Finalize a matrix. */
  template<framework f>
  inline void finalize(matrix<f> &m) {
    if constexpr (f == framework::petsc) {
      MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(m,   MAT_FINAL_ASSEMBLY);
    }
  }

  /* Destroy a matrix. */
  template<framework f>
  inline void destroy(matrix<f> &m) {
    if constexpr (f == framework::petsc) MatDestroy(&m);
  }

  /* Set a matrix entry mat[row, col] = val for some matrix, mat and 
     value, val. */
  template<framework f, typename m, typename v>
  inline void set_matrix_value(m &mat, st row, st col, v val) {
    if constexpr (f == framework::petsc)
      MatSetValue(mat, row, col, val, ADD_VALUES);
  }

  template<framework f>
  inline void matmul(matrix<f> &a, matrix<f> &b, matrix<f> &c) {
    if constexpr (f == framework::petsc)
      MatMatMult(a, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c); 
  }
};