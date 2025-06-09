#include "io.h"

std::string io::make_operator_path(
    std::size_t  quant, 
    std::size_t  size) {
    std::string res = "../operator/";
    res += "V_";
    res += std::to_string(size * size * quant * quant); 
    res += "_N_";
    res += std::to_string(size); 
    res += "_L_";
    res += std::to_string(quant); 
    return res;
}

void io::write_chol(real_t *m, std::size_t sz, components &sbp, std::string key) {
    std::string dir = make_operator_path(sbp.n_blocks_dim, sbp.n);
    namespace fs = std::filesystem;
    fs::create_directories(dir); // create dir if not already there. 
    std::ofstream out;
    out.open(dir + "/" + key +  dense_chol_suffix,
    std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&m[0]), sizeof(real_t) * sz * sz);
    out.close();
    std::cout << dir + dense_chol_suffix << std::endl;
    return;
}

void io::write_chol_sparse(real_t *m, std::size_t sz, components &sbp, std::string key) {
    std::string dir = make_operator_path(sbp.n_blocks_dim, sbp.n);
    namespace fs = std::filesystem;
    fs::create_directories(dir); // create dir if not already there. 
    std::ofstream out;
    out.open(dir + "/" + key +  dense_chol_suffix,
    std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&m[0]), sizeof(real_t) * sz * sz);
    out.close();
    std::cout << dir + dense_chol_suffix << std::endl;
    return;
}

bool io::load_factors(real_t *m, std::size_t sz, components &sbp, std::string key) {
    std::string dir = make_operator_path(sbp.n_blocks_dim, sbp.n);
    dir += "/";
    dir += key;
    dir += dense_chol_suffix;
    namespace fs = std::filesystem;
    const fs::path p{dir};
    if (fs::exists(p)) {
        auto ifs = std::ifstream(dir, std::ios::binary);
        ifs.read(
            reinterpret_cast<char*>(m), 
            sizeof(real_t) * sz * sz);
        ifs.close();
        double gb = (sizeof(real_t) * sz * sz) / 1e9;
        std::cout << "Loaded dense cholesky factors \"" << dir << "*\" (" 
        << std::fixed << gb << " GB)" << std::endl;
        return true;
    }
    else {
        return false;
    }
}

void io::print_csr(real_t *m, std::size_t sz, components &sbp, std::string key, std::size_t szb) {
  
    int ssz = static_cast<int>(sz);
  
    std::string dir = "../operator/";
    dir += "V_";
    dir += std::to_string(sbp.n * sbp.n * sbp.n_blocks_dim * sbp.n_blocks_dim); 
    dir += "_N_";
    dir += std::to_string(sbp.n); 
    dir += "_L_";
    dir += std::to_string(sbp.n_blocks_dim); 
  
    namespace fs = std::filesystem;
    
    fs::create_directories(dir);
  
    int *rs, *re, *ci; 
    double *v;
    int sszb = 0;
    if (szb != 0) {
      sszb = szb;
    }
    else {
      sszb = ssz;
    }


    csr_t mm = {m, sz, sz};
    //out.write( reinterpret_cast<const char*>( &f ), sizeof( float ));
    
    std::ofstream out;

  
    out.open( dir + "/" + key + ".ndim", std::ios::out | std::ios::binary);
    out.write( reinterpret_cast<const char*>(&ssz), sizeof(int));
    out.close();
    out.open( dir + "/" + key + ".mdim", std::ios::out | std::ios::binary);
    out.write( reinterpret_cast<const char*>(&sszb), sizeof(int));
    out.close();
  
    out.open( dir + "/" + key + ".row", std::ios::out | std::ios::binary);
    out.write( reinterpret_cast<const char*>(&mm.r[0]), sizeof(int) * (mm.row_index_size()));
    out.close();

    out.open(dir + "/" + key + ".col", std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&mm.c[0]), sizeof(int) * (mm.col_index_size()));
    out.close();
    
    out.open(dir + "/" + key + ".val", std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&mm.v[0]), sizeof(double) * (mm.nnz()));
    out.close();
  
    out.open(dir + "/" + key + ".rsize", std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&mm.n), sizeof(int));
    out.close();
  
    out.open(dir + "/" + key + ".csize", std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(&rs[sz]), sizeof(int));
    out.close();
  
    std::cout << dir + "/" + key + ".ndim " << ssz << std::endl;
    std::cout << dir + "/" + key + ".mdim " << sszb << std::endl;
    std::cout << dir + "/" + key + ".rsize " << mm.row_index_size() << std::endl;
    std::cout << dir + "/" + key + ".csize " << mm.col_index_size() << std::endl;

  }