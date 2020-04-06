//
//  lattice.cpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#include "lattice.hpp"

const double med_coeff = 1;

Lattice::Lattice(){
    size = 0;
    len = 0;
}

Lattice::Lattice(size_t len_n, double x_0, double y_0, double z_0, double radius){
    len = len_n;
    size = len * len *len;
    for (size_t i = 0; i < len; i++){
        for (size_t j = 0; j < len; j++){
            for (size_t k = 0; k < len; k++){
                make_cell(i, j, k, x_0, y_0, z_0, radius);
            }
        }
    }
}

void Lattice::make_cell(size_t i, size_t j, size_t k, double x_0, double y_0, double z_0, double radius){
    double x = (double)((2 * k) + 1)/(double)(2 * len);
    double y = (double)((2 * j) + 1)/(double)(2 * len);
    double z = (double)((2 * i) + 1)/(double)(2 * len);
    
    Point temp(x, y, z, 0);
    
    double sum = 0;
    for (size_t l = 0; l < i + 2; i++){
        for (size_t m = 0; m < j + 2; j++){
            for (size_t n = 0; n < k + 2; k++){
                double x = (double)n / (double)len;
                double y = (double)m / (double)len;
                double z = (double)l / (double)len;
                
                if (((x-x_0)*(x-x_0)) + ((y-y_0)*(y-y_0)) + ((z-z_0)*(z-z_0)) <= radius * radius){
                    sum += med_coeff;
                }
            }
        }
    }
    temp.set_Eta((double)sum / 8.0);
    lattice.push_back(temp);
}

Point Lattice::at(size_t i){
    return lattice[i];
}

size_t Lattice::Size(){
    return size;
    }
