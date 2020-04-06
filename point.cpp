//
//  point.cpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#include "point.hpp"

Point::Point(double x_n, double y_n, double z_n, double eta_n){
    x = x_n;
    y = y_n;
    z = z_n;
    eta = eta_n;
}

double  Point::X(){
    return x;
}

double Point::Y(){
    return y;
}

double Point::Z(){
    return z;
}

double Point::Eta(){
    return eta;
}

void Point::set_Eta(double eta_n){
    eta = eta_n;
}
