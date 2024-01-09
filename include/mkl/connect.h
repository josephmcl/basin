#pragma once 

#include <vector>

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