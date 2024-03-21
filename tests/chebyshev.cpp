#include "../src/chebyshev.hpp"

int main(){
    std::vector<double> x0 = {1., 2., 3.};
    std::vector<double> v = {1., 2., 4., 4., 6.};
    std::vector<std::size_t> c = {0, 1, 1, 1, 2};
    std::vector<std::size_t> r = {0, 2, 3, 5};
    CSR_Matrix a{v, c, r};
    std::vector<double> b = {1., 1., 1.};
    std::vector<double> res = chebyshev(a, x0, b, 7, 0.000001, 1., 6.);
    for(size_t  i = 0; i < b.size(); i++){
        std::cout<<res[i]<<" ";
    }
}