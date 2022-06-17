#pragma once

#include <type_traits>

/* petsc headers */
#include "petscsys.h"
#include "petscmat.h" 
#include "petscviewer.h"
#include "petscvec.h"
#include "petscksp.h"



namespace linalg {

  enum class framework : std::size_t {
    petsc = 0
  };

  /* petsc type aliases */
  using petsc_matrix = Mat;

  using st = std::size_t;

  /* NOTE: Idea for this was robbed from stackoverflow.com/a/61177390 */
  template<framework f>
  constexpr auto pick_matrix_type() {
    if constexpr (f == framework::petsc) {
        return std::type_identity<petsc_matrix>{};
    } 
  }

  template<framework f>
  using matrix = typename decltype(pick_matrix_type<f>())::type;

  /* set a matrix entry mat[row, col] = val for some matrix, mat and 
     value, val. */
  template<framework f, typename m, typename v>
  void set_matrix_value(m &mat, st row, st col, v val) {
    if constexpr (f == framework::petsc) {
      MatSetValue(mat, row, col, val, ADD_VALUES);
    }
  }


};