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

#endif /* utility_hpp */
