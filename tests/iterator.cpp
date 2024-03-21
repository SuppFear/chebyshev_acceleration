#include "../src/chebyshev.hpp"
#include "../../simple_iterative_methods/src/iterations.hpp"
#include <fstream>

int main(){
    std::vector<double> x0 = {1., 2., 3.};
    std::vector<double> v = {1., 2., 4., 4., 6.};
    std::vector<std::size_t> c = {0, 1, 1, 1, 2};
    std::vector<std::size_t> r = {0, 2, 3, 5};
    CSR_Matrix a{v, c, r};
    std::vector<double> b = {1., 1., 1.};

    std::ofstream out;
    out.open("simple_iteration.txt");
    if (out.is_open())
    {
        for (size_t n = 1; n < 32; n++){
            std::cout << n << std::endl;
            std::vector<double> res = simple(a, x0, b, 0.25, n, 0.00000001);
            out << std::to_string(n) << " "  << std::to_string(std::sqrt(res*res)) << std::endl; 
        }   
    }
    out.close();

    out.open("jacobi_iteration.txt");
    if (out.is_open())
    {
        for (size_t n = 1; n < 32; n++){
            std::cout << n << std::endl;
            std::vector<double> res = jacobi(a, x0, b, 0.25, n, 0.00000001);
            out << std::to_string(n) << " "  << std::to_string(std::sqrt(res*res)) << std::endl; 
        }   
    }
    out.close();

    out.open("gauss_seidel_iteration.txt");
    if (out.is_open())
    {
        for (size_t n = 1; n < 32; n++){
            std::cout << n << std::endl;
            std::vector<double> res = gauss_seidel(a, x0, b, 0.25, n, 0.00000001);
            out << std::to_string(n) << " "  << std::to_string(std::sqrt(res*res)) << std::endl; 
        }   
    }
    out.close();

    out.open("simple_chebyshev_iteration.txt");
    if (out.is_open())
    {
        for (size_t n = 1; n < 5; n++){
            std::cout << n << std::endl;
            std::vector<double> res = chebyshev(a, x0, b, n, 0.00000001, 1, 6);
            out << std::to_string(n) << " "  << std::to_string(std::sqrt(res*res)) << std::endl; 
        }   
    }
    out.close();
}