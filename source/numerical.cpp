#include "numerical.h"

using namespace type;

static std::vector<real_t> bhinv_p2 = {2.};
static std::vector<real_t> bhinv_p4 = {48./17., 48./59., 48./43., 
    48./49.};
static std::vector<real_t> bhinv_p6 = {43200./13649., 8640./12013.,  
    4320./2711., 4320./ 5359., 8640./7877., 43200./43801.};

std::vector<real_t> numerical::operators::H_inverse(std::size_t nodes, 
    std::size_t order, real_t left, real_t right) {

    real_t grid_size = (right - left) / static_cast<real_t>(nodes - 1);
    std::vector<real_t> bhinv = 
        (order == 2) ? bhinv_p2 :
        (order == 4) ? bhinv_p4 :
        (order == 6) ? bhinv_p6 : bhinv_p2;

    real_t t;
    auto H = std::vector<real_t>(nodes, 1.);
    for (std::size_t i = 0; i < bhinv.size(); ++i) {
        t = 1 / bhinv[i];
        H[i] = t;
        H[nodes - 1 - i] = t; }

    for (std::size_t i = 0; i < H.size(); ++i) 
        H[i] = 1 / (H[i] * grid_size);

    return H; }

std::vector<real_t> numerical::operators::H(std::size_t nodes, 
    std::size_t order, real_t left, real_t right) {

    real_t grid_size = (right - left) / static_cast<real_t>(nodes - 1);
    std::vector<real_t> bhinv = 
        (order == 2) ? bhinv_p2 :
        (order == 4) ? bhinv_p4 :
        (order == 6) ? bhinv_p6 : bhinv_p2;

    real_t t;
    auto H = std::vector<real_t>(nodes, 1.);
    for (std::size_t i = 0; i < bhinv.size(); ++i) {
        t = 1 / bhinv[i];
        H[i] = t;
        H[nodes - 1 - i] = t; }

    for (std::size_t i = 0; i < H.size(); ++i) 
        H[i] *= grid_size;

    return H; }
