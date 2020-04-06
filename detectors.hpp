//
//  detectors.hpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef detectors_hpp
#define detectors_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"

class Detectors{
public:
    Detectors();
    
    Detectors(size_t size_n, double x_0, double y_0, double z_0, double radius);
    
    size_t Size();
    
    
private:
    std::vector<Point> detectors;
    size_t size;
};

#endif /* detectors_hpp */
