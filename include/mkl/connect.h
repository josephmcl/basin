#pragma once 

#include <vector>
#include <iostream>

#include "definitions.h"
#include "components.h"


// Given n, return the connectivity matrix for an [n x n] grid of 
// elements. Element ordering assumed column major, bottom left to 
// top right. 
std::size_t make_connectivity(
  std::vector<std::vector<std::size_t>> 
  &con, std::size_t const n);

void make_boundary_maps(
  std::vector<std::vector<std::size_t>> &data, 
  std::vector<std::vector<std::size_t>> &order,
  std::size_t const n);

void make_interface_list(
  vv<std::size_t>       &interface_list, 
  vv<std::size_t>          const &F_symbols,
  vv<std::size_t>          const &FT_symbols,
  components               const &sbp);