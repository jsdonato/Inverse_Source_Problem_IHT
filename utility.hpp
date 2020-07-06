//
//  utility.hpp
//  ISP3C
//
//  Created by Joseph Donato on 5/31/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef utility_hpp
#define utility_hpp

#include <stdio.h>
#include <cmath>
#include <complex>

std::complex<double> Ai(double rx, double ry, double rz, double theta, double phi){
    std::complex<double> num(0.0, freq * ((rx * std::cos(theta) * std::sin(phi)) + (ry * std::sin(theta) * std::sin(phi)) + (rz * std::cos(phi))));
    return exp(num);
}

std::complex<double> G(double rx, double ry, double rz, double rxs, double rys, double rzs){
    double mag = std::sqrt(std::pow(rx - rxs, 2) + std::pow(ry - rys, 2) + std::pow(rz - rzs, 2));
    double num = (1 / (4 * M_PI * mag));
    std::complex<double> numtop(0.0, freq * mag);
    return num * exp(numtop);
}

double coherence(arma::cx_mat A){
    std::vector<double> data(A.n_cols, 0);
    
    std::vector<double> cdots;
    for (int i = 0; i < A.n_cols; i++){
        cdots.push_back(real(sqrt(cdot(A.col(i), A.col(i)))));
    }
    
    for (int i = 0; i < A.n_cols; i++){
        double max = 0.0;
        std::complex<double> temp = cdots[i];
        for (int j = i; j < A.n_cols; j++){
            if (i != j){
                std::complex<double> num = abs(cdot(A.col(i), A.col(j))) / (cdots[j] * temp);
                if (real(num) > max){
                    max = real(num);
                }
            }
        }
        data[i] = max;
    }
    
    return *max_element(data.begin(), data.end());
}

#endif /* utility_hpp */
