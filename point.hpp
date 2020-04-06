//
//  point.hpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef point_hpp
#define point_hpp

#include <stdio.h>

class Point{
public:
    Point(double x_n, double y_n, double z_n, double eta_n);
    
    double  X();
    
    double Y();
    
    double Z();
    
    double Eta();
    
    void set_Eta(double eta_n);
    
private:
    double x;
    double y;
    double z;
    double eta;
};

#endif /* point_hpp */
