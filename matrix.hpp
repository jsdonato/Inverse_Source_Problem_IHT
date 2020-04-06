//
//  matrix.hpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>
#include <complex>
#include <cmath>
#include <armadillo>
#include "detectors.hpp"
#include "lattice.hpp"

class Matrix{
public:
    Matrix(Lattice &l_n, Detectors &d_n);
    
    std::complex<double> Ai(double rx, double ry, double rz, double theta, double phi);
    
    std::complex<double> G(double rx, double ry, double rz, double rxs, double rys, double rzs);
    
    /*cx_mat Gv();*/
    
    arma::cx_mat V();
    
    /*cx_mat Gd();*/
    
    /*cx_mat Gs();*/
    
private:
    Lattice l = Lattice();
    Detectors d = Detectors();
};

#endif /* matrix_hpp */
