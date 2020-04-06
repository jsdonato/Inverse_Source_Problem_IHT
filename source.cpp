//
//  source.cpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#include "source.hpp"

Source::Source(double theta_n, double phi_n, double amplitude_n){
    theta = theta_n;
    phi = phi_n;
    amplitude = amplitude_n;
}

double Source::T(){
    return theta;
}

double Source::P(){
    return phi;
}

double Source::A(){
    return amplitude;
}
