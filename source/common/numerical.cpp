#include "numerical.h"

using namespace type;

static const std::vector<real_t> bhinv_p2 = {2.};

/* second-order, first derivative matrix components. */ 
static const std::vector<real_t> d_p2 = {-1./2., 0., 1./2.};
static const std::vector<std::vector<real_t>> bd_p2 = {{-1., 1.}};

static const std::vector<real_t> bhinv_p4 = {48./17., 48./59., 48./43., 
    48./49.};

/* fourth-order, first derivative matrix components. */ 
static const std::vector<real_t> d_p4 = {1./12., -2./3., 0., 2./3., -1./12.};
static const std::vector<std::vector<real_t>> bd_p4 = {
    {-24./17.,  59./34.,  -4./17., -3./34.,      0.,      0.}, 
    { -1./ 2.,       0.,   1./ 2.,      0.,      0.,      0.}, 
    {  4./43., -59./86.,       0., 59./86., -4./43.,      0.}, 
    {  3./98.,       0., -59./98.,      0., 32./49., -4./49.}};

static const std::vector<real_t> bhinv_p6 = { 43200./13649., 8640./12013., 
    4320./2711., 4320./ 5359., 8640./7877., 43200./43801.};

/* sixth-order, first derivative matrix components. */ 
static const std::vector<real_t> d_p6 = {-1./60., 3./20., -3./4., 0., 3./4., 
    -3./20., 1./60.};
static const std::vector<std::vector<real_t>> bd_p6 = { 
    { -21600./ 13649., 104009./ 54596.,  30443./81894., -33311./ 27298.,   16863./ 27298., -15025./163788.,            0.,            0.,          0.},
    {-104009./240260.,              0.,   -311./72078.,  20229./ 24026.,  -24337./ 48052.,  36661./360390.,            0.,            0.,          0.}, 
    { -30443./162660.,    311./ 32532.,             0., -11155./ 16266.,   41287./ 32532., -21999./ 54220.,            0.,            0.,          0.}, 
    {  33311./107180., -20229./ 21436.,    485./ 1398.,              0.,    4147./ 21436.,  25427./321540.,    72./ 5359.,            0.,          0.},
    { -16863./ 78770.,  24337./ 31508., -41287./47262.,  -4147./ 15754.,               0., 342523./472620., -1296./ 7877.,   144./ 7877.,          0.},
    {  15025./525612., -36661./262806.,  21999./87602., -25427./262806., -342523./525612.,              0., 32400./43801., -6480./43801., 720./43801.}};


/* second-order, second derivative matrix components. */ 
static const std::vector<real_t> d2_p2 = {1., -2., 1.};
static const std::vector<std::vector<real_t>> bd2_p2 = {
    {1., -2., 1.},
    {0.,  0., 0.},
    {0.,  0., 0.}};

/* fourth-order, second derivative matrix components. */ 
static const std::vector<real_t> d2_p4 = {
    -1./12., 4./3., -5./2., 4./3., -1./12.};
static const std::vector<std::vector<real_t>> bd2_p4 = {
    {     2.,     -5.,        4.,       -1.,      0.,      0.},
    {     1.,     -2.,        1.,        0.,      0.,      0.},
    {-4./43., 59./43., -110./43.,   59./43., -4./43.,      0.},
    {-1./49.,      0.,   59./49., -118./49., 64./49., -4./49.}};

/* second-order, second derivative matrix SAT component. */ 
static const std::vector<real_t> d2_bs2 = {3./2., -2., 1./2.};
/* fourth-order, second derivative matrix SAT component.. */ 
static const std::vector<real_t> d2_bs4 = {11./6., -3., 3./2., -1./3.,};
    
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

numerical::operators::sbp::sbp(
    std::size_t const size, std::size_t const order, 
    real_t const left, real_t const right): 
    size(size), order(order), left(left), right(right), 
    grid_size((right - left) / static_cast<ℝ>(size - 1)) { }

std::tuple<std::size_t, std::vector<real_t> const *>
numerical::operators::sbp::row(std::size_t const index) const {

    if (index < top.size()) {
        return std::make_tuple(0, &top[index]);
    }  
    else if (index < size - bottom.size()) {
        // adjust column index for kernel size. 
        auto j = index - ((d.size() - 1) / 2);
        return std::make_tuple(j, &d);
    }
    else if (index < size) {
        // find logical index from explicitly stored kernel.
        auto i = index - (size - bottom.size());
        // adjust column index for kernel size.
        auto j = size - bottom[i].size();
        return std::make_tuple(j, &bottom[i]);
    }
    else { throw; }
}

std::tuple<std::size_t, std::vector<real_t> const >
numerical::operators::sbp::rowf(std::size_t const index) const {

    std::size_t row_start_index;
    auto res = std::vector<real_t>();

    // TODO: soft code this later
    if (index  == 0) { 
        row_start_index = 0;
        res = top[index];
    }  
    else if (index < size - 1) {
        row_start_index = index - ((d.size() - 1) / 2);
        res = d;
    }
    else if (index == size - 1) {
        row_start_index = size - bottom[2].size();
        res = bottom[2];
        
    }
    else { throw; }

    if (index == 1 or index == 2) {
        auto temp = std::vector<real_t>();
        temp = top[index];
        for (std::size_t i = 0; i != res.size(); ++i) {
            if (row_start_index + i < temp.size()) {
                temp[row_start_index + i] += res[i];
            }
            else {
                temp.push_back(res[i]);
            }
        }
        res = temp;
        row_start_index = 0;
    }

    if (index == size - 2 or index == size - 3) {
        // row_start_index = 0;
        auto f = index - (size - 3);
        std::cout << f << std::endl;
        auto temp_start_index = size - bottom[f].size();
        auto temp = std::vector<real_t>();
        // temp = bottom[index];
        auto ofs = bottom[f].size() - res.size();
        for (std::size_t i = row_start_index; i != size; ++i) {
            if (i < temp_start_index) {
                temp.push_back(res[i - temp_start_index]);
            }
            else {
                temp.push_back(res[i - temp_start_index] + bottom[f][i - temp_start_index - ofs]);
            }
        }
        res = temp;
        row_start_index = size - 3;
    }

    // Apply the grid spacing "H" matrix. 
    for (std::size_t i = 0; i < res.size(); ++i) {
        res[i] *= h[row_start_index + i];
    } 

    return std::make_tuple(row_start_index, res);
}

/*

void numerical::operators::sbp::
product(std::vector<ℝ> const &rhs, std::vector<ℝ> &lhs) {
    for (std::size_t i = 0; i != lhs.size(); ++i) {
        auto [column, r] = row(i);
        for (std::size_t j = 0; j != r->size(); ++j) {
            lhs[i] += r->at(j) * rhs[column + j];
        }
    }
}

*/

static void load_d1_operator(
    std::size_t order,
    real_t grid_size,
    std::vector<real_t> &d1_interior,
    std::vector<std::vector<real_t>> &d1_top,
    std::vector<std::vector<real_t>> &d1_bottom) {

    if (order == 2) {
        d1_interior = d_p2;
        d1_top      = bd_p2;
        d1_bottom   = bd_p2; 
    }
    else if (order == 4) {
        d1_interior = d_p4;
        d1_top      = bd_p4;
        d1_bottom   = bd_p4; 
    }


    std::vector<std::vector<real_t>> temp(d1_bottom.size());

    std::reverse_copy(std::begin(d1_bottom), std::end(d1_bottom), 
        std::begin(temp));
    d1_bottom.clear();

    for (std::size_t i = 0; i != temp.size(); ++i) {
        d1_bottom.push_back(temp[i]);
        std::reverse(d1_bottom[i].begin(), d1_bottom[i].end());
    } 

    for (std::size_t i = 0; i != d1_interior.size(); ++i) {
        d1_interior[i] /= grid_size;
    }

    for (std::size_t i = 0; i != d1_top.size(); ++i) {
        for (std::size_t j = 0; j != d1_top[i].size(); ++j) {
            d1_top[i][j] /= grid_size;
            d1_bottom[i][j] /= -grid_size;
        } 
    }
}

void numerical::operators::d1::load_operator() {
    /* Set the operator values in the struct. */

    load_d1_operator(order, grid_size, d, top, bottom);

}

void numerical::operators::d2::load_operator() {
    /* Set the operator values in the struct. */
    if (order == 2) {
        d = d2_p2;
        top = bd2_p2;
        bottom = bd2_p2; 
    }
    else if (order == 4) {
        d = d2_p4;
        top = bd2_p4;
        bottom = bd2_p4; 
    }

    std::vector<std::vector<real_t>> temp(bottom.size());

    std::reverse_copy(std::begin(bottom), std::end(bottom), 
        std::begin(temp));
    bottom.clear();

    for (std::size_t i = 0; i != temp.size(); ++i) {
        bottom.push_back(temp[i]);
        std::reverse(bottom[i].begin(), bottom[i].end());
    } 

    auto grid_size_square = grid_size * grid_size;

    for (std::size_t i = 0; i != d.size(); ++i) {
        d[i] /= grid_size_square;
    }

    for (std::size_t i = 0; i != top.size(); ++i) {
        for (std::size_t j = 0; j != top[i].size(); ++j) {
            top[i][j] /= grid_size_square;
            bottom[i][j] /= grid_size_square;
        }
    }

    // Load the first derivative finite difference operator needed for 
    // deriving boundary conditions.
    load_d1_operator(order, grid_size, d1_interior, d1_top, d1_bottom);
}

void numerical::operators::d2::load_h() {
    auto diag = numerical::operators::H(size, order, left, right);
    for (std::size_t i = 0; i < size; ++i) {
        h[i] = diag[i];
    }
}

void numerical::operators::d2::load_h_inverse() {
    auto diag = numerical::operators::H_inverse(
        size, order, left, right);
    for (std::size_t i = 0; i < size; ++i) {
        hi[i] = diag[i];
    }
}

void numerical::operators::d2::load_left_boundary_data() {
    
    std::vector<real_t> bs; 
    if (order == 2) {
        bs = d2_bs2;
    }
    else if (order == 4) {
        bs = d2_bs4;
    }

    for (auto &e: bs) {
        top_boundary_data.push_back(e / grid_size);
    }   
}

void numerical::operators::d2::load_right_boundary_data() {
    
    std::vector<real_t> bs; 
    if (order == 2) {
        bs = d2_bs2;
    }
    else if (order == 4) {
        bs = d2_bs4;
    }

    for (auto &e: bs) {
        bot_boundary_data.push_back(e / grid_size);
        std::reverse(bot_boundary_data.begin(),bot_boundary_data.end());
    }   
}

void numerical::operators::d2::fuse(std::vector<ℝ> const &diag) {
    for (std::size_t i = 0; i < size; ++i) h[i] = diag[i];
}

void numerical::operators::d2::fuse_hi(std::vector<ℝ> const &diag) {
    for (std::size_t i = 0; i < size; ++i) hi[i] = diag[i];
}

void numerical::operators::d2::
left_boundary(real_t value, std::size_t degree) {
    if (degree == 1) {
        top[0][0] += value * hi[0] * top_boundary_data[0];
        top[1][0] += value * hi[1] * top_boundary_data[1];
        top[2][0] += value * hi[2] * top_boundary_data[2];
    }
}

void numerical::operators::d2::
right_boundary(real_t value, std::size_t degree) {
    
    auto i = bottom.size() - 1;
    auto j = bottom[i].size() - 1;

    if (degree == 1) {
        
        bottom[i][j] = value * hi[size - 1];
    }
    if (degree == 2) {

        auto &bottom_row = bottom[bottom.size() - 1];
        auto offset = bottom_row.size() - 3;
        for (std::size_t i = 0; i != 3; ++i) {
            bottom_row[i + offset] += hi[size - 1] * bot_boundary_data[i] * value;
        }

    }
}

void numerical::operators::d2::left_sat(real_t value) {
    load_left_boundary_data();
    std::size_t i = 0;
    for (auto &e : top_boundary_data)
        e *= value * hi[i++];
}

void numerical::operators::d2::right_sat(real_t value) {
    load_right_boundary_data();
    std::size_t i = size;
    for (auto &e : bot_boundary_data)
        e *= value * hi[--i];
}



