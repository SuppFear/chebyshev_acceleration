#pragma once
#include "../../Dense-CSR/csr_matrix.h"
#include "../../Dense-CSR/vector_operations.h"
#include <cmath>

std::vector<double> chebyshev(const CSR_Matrix &A, std::vector<double> &x0, const std::vector<double> &b, const size_t &r_max, const double &tolerance, const double &lmin, const double &lmax){
    size_t n = pow(2, r_max);
    std::vector<size_t> ind = std::vector<size_t>(n);
    ind[0] = 0;
    ind[n/2] = 1;
    size_t tmp = n/2;
    for(size_t r = 2; r <= r_max; r++){
        
        for(size_t j = 0; j < n; j += tmp){
            ind[j + tmp/2] = pow(2, r) - 1 - ind[j];
        }
        tmp /= 2;
    }
    // for(size_t i = 0; i < n; i++){
    //     std::cout<<ind[i]<<" ";
    // }
    std::vector<double> t = std::vector<double>(n);
    for(size_t i = 0; i < n; i++){
        t[ind[i]] = 1 / (((lmin + lmax) / 2) + ((lmax - lmin) * cos(M_PI * (2 * i + 1) / (2 * n)) / 2));
    }
    // for(size_t i = 0; i < n; i++){
    //     std::cout<<t[i]<<" ";
    // }
    std::vector<double> x = x0;
    std::vector<double> error = std::vector<double>(b.size());
    double diff = 0.;
    for(size_t i = 0; i < n; i++){
        //std::cout<<diff<<std::endl;
        error = (A*x - b);
        x = x - t[i]*error;
        diff = 0;
        for(size_t j = 0; j < b.size(); j++){
            diff += (t[i] * error[j]) * (t[i] * error[j]);
        }
        //std::cout<<diff<<" ";
        if(std::sqrt(diff) < tolerance){
            //std::cout<<i<<" ";
            break;
        }
    }
    return x;
}