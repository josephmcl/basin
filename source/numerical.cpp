#include "numerical.h"

using namespace type;

const std::vector<real_t> bhinv_p2 = {2.};
const std::vector<real_t> d_p2 = {-1./2., 0., 1./2.};
const std::vector<std::vector<real_t>> bd_p2 = {{-1., 1.}};

const std::vector<real_t> bhinv_p4 = {48./17., 48./59., 48./43., 
    48./49.};
const std::vector<real_t> d_p4 = {1./12., -2./3., 0., 2./3., -1./12.};
const std::vector<std::vector<real_t>> bd_p4 = {
    {-24./17.,  59./34.,  -4./17., -3./34.,      0.,      0.}, 
    { -1./ 2.,       0.,   1./ 2.,      0.,      0.,      0.}, 
    {  4./43., -59./86.,       0., 59./86., -4./43.,      0.}, 
    {  3./98.,       0., -59./98.,      0., 32./49., -4./49.}};

const std::vector<real_t> bhinv_p6 = { 43200./13649., 8640./12013., 
    4320./2711., 4320./ 5359., 8640./7877., 43200./43801.};
const std::vector<real_t> d_p6 = {-1./60., 3./20., -3./4., 0., 3./4., 
    -3./20., 1./60.};
const std::vector<std::vector<real_t>> bd_p6 = { 
    { -21600./ 13649., 104009./ 54596.,  30443./81894., -33311./ 27298.,   16863./ 27298., -15025./163788.,            0.,            0.,          0.},
    {-104009./240260.,              0.,   -311./72078.,  20229./ 24026.,  -24337./ 48052.,  36661./360390.,            0.,            0.,          0.}, 
    { -30443./162660.,    311./ 32532.,             0., -11155./ 16266.,   41287./ 32532., -21999./ 54220.,            0.,            0.,          0.}, 
    {  33311./107180., -20229./ 21436.,    485./ 1398.,              0.,    4147./ 21436.,  25427./321540.,    72./ 5359.,            0.,          0.},
    { -16863./ 78770.,  24337./ 31508., -41287./47262.,  -4147./ 15754.,               0., 342523./472620., -1296./ 7877.,   144./ 7877.,          0.},
    {  15025./525612., -36661./262806.,  21999./87602., -25427./262806., -342523./525612.,              0., 32400./43801., -6480./43801., 720./43801.}};

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

numerical::operators::D1::
D1(ℤ const &size, ℤ const &order, ℝ const left, ℝ const right): 
    size(size), order(order), left(left), right(right), 
    grid_size((right - left) / static_cast<ℝ>(size - 1)) {

    /* Set the operator values in the struct. */
    if (order == 2) {
        d = d_p2;
        top = bd_p2;
        bottom = bd_p2; 
    }
    else if (order == 4) {
        d = d_p4;
        top = bd_p4;
        bottom = bd_p4; 
    }

    std::vector<std::vector<real_t>> temp(bottom.size());

    std::reverse_copy(std::begin(bottom), std::end(bottom), 
        std::begin(temp));
    bottom.clear();

    for (std::size_t i = 0; i != temp.size(); ++i) {
        bottom.push_back(temp[i]);
        std::reverse(bottom[i].begin(), bottom[i].end());
    } 

    for (std::size_t i = 0; i != d.size(); ++i) {
        d[i] /= grid_size;
    }

    for (std::size_t i = 0; i != top.size(); ++i) {
        for (std::size_t j = 0; j != top[i].size(); ++j) {
            top[i][j] /= grid_size;
            bottom[i][j] /= -grid_size;
        }
    }
}

std::tuple<std::size_t, std::vector<real_t> const *>
numerical::operators::D1::row(std::size_t const index) const {

    if (index < top.size()) {
        return std::make_tuple(0, &top[index]);
    }  
    else if (index < size - bottom.size()) {
        auto j = index - ((d.size() - 1) / 2);
        return std::make_tuple(j, &d);
    }
    else if (index < size) {
        auto i = index - (size - bottom.size());
        auto j = size - bottom[i].size();
        return std::make_tuple(j, &bottom[i]);
    }
    else { throw; }
}

