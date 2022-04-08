#include "numerical.h"

using namespace type;

static std::vector<real_t> bhinv_p2 = {2.};
static std::vector<real_t> d_p2 = {-1./2., 0., 1./2.};
static std::vector<real_t> bd_p2 = {-1., 1.};

static std::vector<real_t> bhinv_p4 = {48./17., 48./59., 48./43., 
    48./49.};
static std::vector<real_t> d_p4 = {1./12., -2./3., 0., 2./3., -1./12.};
static std::vector<real_t> bd_p4 = {
    -24./17.,  59./34.,  -4./17., -3./34.,      0.,      0., 
     -1./ 2.,       0.,   1./ 2.,      0.,      0.,      0., 
      4./43., -59./86.,       0., 59./86., -4./43.,      0., 
      3./98.,       0., -59./98.,      0., 32./49., -4./49.};

static std::vector<real_t> bhinv_p6 = { 43200./13649., 8640./12013., 
    4320./2711., 4320./ 5359., 8640./7877., 43200./43801.};
static std::vector<real_t> d_p6 = {-1./60., 3./20., -3./4., 0., 3./4., 
    -3./20., 1./60.};
static std::vector<real_t> bd_p6 = { 
     -21600./ 13649., 104009./ 54596.,  30443./81894., -33311./ 27298.,   16863./ 27298., -15025./163788.,            0.,            0.,          0.,
    -104009./240260.,              0.,   -311./72078.,  20229./ 24026.,  -24337./ 48052.,  36661./360390.,            0.,            0.,          0., 
     -30443./162660.,    311./ 32532.,             0., -11155./ 16266.,   41287./ 32532., -21999./ 54220.,            0.,            0.,          0., 
      33311./107180., -20229./ 21436.,    485./ 1398.,              0.,    4147./ 21436.,  25427./321540.,    72./ 5359.,            0.,          0.,
     -16863./ 78770.,  24337./ 31508., -41287./47262.,  -4147./ 15754.,               0., 342523./472620., -1296./ 7877.,   144./ 7877.,          0.,
      15025./525612., -36661./262806.,  21999./87602., -25427./262806., -342523./525612.,              0., 32400./43801., -6480./43801., 720./43801.
};

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

void numerical::operators::D1(std::size_t const &size, 
    std::size_t const &order=2, real_t const left=-1., 
    real_t const right=1., result &std::vector<real_t>) {



}
