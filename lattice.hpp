//
//  lattice.hpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef lattice_hpp
#define lattice_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"




class Lattice{
public:
    Lattice();
    
    Lattice(size_t len_n, double x_0, double y_0, double z_0, double radius);
    
    void make_cell(size_t i, size_t j, size_t k, double x_0, double y_0, double z_0, double radius);
    
    Point at(size_t i);
    
    size_t Size();
    
private:
    std::vector<Point> lattice;
    size_t len;
    size_t size;
};

#endif /* lattice_hpp */
