//
//  matrix.cpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#include "matrix.hpp"

const double freq = 1;

Matrix::Matrix(Lattice &l_n, Detectors &d_n){
    l = l_n;
    d = d_n;
}

std::complex<double> Matrix::Ai(double rx, double ry, double rz, double theta, double phi){
    std::complex<double> num(0.0, freq * ((rx * std::cos(theta) * std::sin(phi)) + (ry * std::sin(theta) * std::sin(phi)) + (rz * std::cos(phi))));
    return exp(num);
}

std::complex<double> Matrix::G(double rx, double ry, double rz, double rxs, double rys, double rzs){
    double mag = std::sqrt(std::pow(rx - rxs, 2) + std::pow(ry - rys, 2) + std::pow(rz - rzs, 2));
    double num = (1 / (4 * M_PI * mag));
    std::complex<double> numtop(0.0, freq * mag);
    return num * exp(numtop);
}

/*cx_mat Matrix::Gv(){
 //TODO: implement this
 }*/

arma::cx_mat Matrix::V(){
    arma::cx_mat V_N(l.Size(), l.Size(), arma::fill::zeros);
    for (int i = 0; i < l.Size(); i++){
        std::complex<double> num(freq * freq * l.at(i).Eta(), 0);
        V_N(i,i) = num;
    }
    return V_N;
}

/*cx_mat Matrix::Gd(){
 //TODO: implement this
 }*/

/*cx_mat Matrix::Gs(){
 //TODO: implement this
 }*/
