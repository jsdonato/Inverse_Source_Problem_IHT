//
//  source.hpp
//  ISP3C
//
//  Created by Joseph Donato on 4/5/20.
//  Copyright © 2020 Joseph Donato. All rights reserved.
//

#ifndef source_hpp
#define source_hpp

#include <stdio.h>

class Source{
public:
    Source(double theta_n, double phi_n, double amplitude_n);
    
    double T();
    
    double P();
    
    double A();
    
private:
    double theta;
    double phi;
    double amplitude;
};


#endif /* source_hpp */
