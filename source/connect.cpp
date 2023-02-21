#include "connect.h"

// Given n, return the connectivity matrix for an [n x n] grid of 
// elements. Element ordering assumed column major, bottom left to 
// top right. 
std::size_t x2::make_connectivity(
  std::vector<std::vector<std::size_t>> &con, std::size_t const n) {

  con.resize(n * n);
  for (auto &e: con) 
    e.resize(n * n);

  std::size_t edge = 1;
  for (std::size_t i = 0; i != n; ++i) {
    for (std::size_t j = 0; j != n - 1; ++j) {
      con[i * n + j][i * n + j + 1] = edge;
      edge += 1;
    }
  }

  for (std::size_t i = 0; i != n * (n - 1); ++i) {
    con[i][i + n] = edge;
    edge += 1;
  }

  return edge - 1;
}


void x2::make_boundary_maps(
  std::vector<std::vector<std::size_t>> &data, 
  std::vector<std::vector<std::size_t>> &order,
  std::size_t const n) { 

  data.resize(n * n);
  for (auto &e: data) 
    e.resize(4);
  order.resize(n * n);
  for (auto &e: order) 
    e.resize(4);

  for (std::size_t i = 0; i != n; ++i) {
    data[i][0] = i + 1;
    data[n * (n - 1) + i][1] = i + 1;
    data[n * i][2] = i + 1;
    data[n * i + n - 1][3] = i + 1;

    order[i][0] = 1;
    order[n * (n - 1) + i][1] = 1;
    order[n * i][2] = 2;
    order[n * i + n - 1][3] = 2;
  }
}