#include <tuple>
#include <vector>
#include <fstream>
#include <sstream>

namespace test {


template<typename T> 
std::tuple<std::size_t, std::size_t, std::vector<T> const *> 
read_csv(std::istream &f) {
    std::size_t rows    = 0;
    std::size_t columns = 0;
    auto data = new std::vector<T>();
    std::string line;
    std::stringstream ss;
    while (std::getline(f, line)) {
        ss = std::stringstream(line);
        T tempt; std::size_t col = 0;
        while (ss >> tempt) {
            std::string temps;
            data->push_back(tempt);
            col++;
            if (!std::getline (ss, temps, ',')) break; 
        }
        columns = col > columns ? col : columns;
        rows++;
    }
    return std::make_tuple(rows, columns, data);
}

};
