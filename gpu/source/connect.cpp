#include "connect.h"

// Given n, return the connectivity matrix for an [n x n] grid of 
// elements. Element ordering assumed column major, bottom left to 
// top right. 
std::size_t make_connectivity(
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


void make_boundary_maps(
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

void make_interface_list(
  vv<std::size_t>       &interface_list, 
  vv<std::size_t>          const &F_symbols,
  vv<std::size_t>          const &FT_symbols,
  components               const &sbp) {  


  interface_list.resize(sbp.n_interfaces);
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (F_symbols[k][i] > 0 && FT_symbols[j][k] > 0) {
          interface_list[i].push_back(j);

          //findex = FT_symbols[i][k] - 1;
          //auto r = k % sbp.n_blocks_dim == 0 ? 0 
          //      : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          //      : 1;
          //mindex = (r * 4) + F_symbols[k][j] - 1;

          // std::cout << findex << " " << mindex << std::endl;

        }
      }
    }
  }
}
