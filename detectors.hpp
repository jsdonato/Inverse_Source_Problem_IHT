#ifndef detectors_hpp
#define detectors_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"

class Detectors{
public:
    
    Detectors(double x_0, double y_0, double z_0, double radius){
        int num1 = 40;
        int num2 = 20;
        detectors.reserve(num1 * num2);
        for (int i = 0; i < num1; i++){
            for (int j = 0; j < num2; j++){
                double x = (radius * cos((i * 2 * M_PI) / num1) * sin((j * M_PI) / num2)) + x_0;
                double y = (radius * sin((i * 2 * M_PI) / num1) * sin((j * M_PI) / num2)) + y_0;
                double z = (radius * cos((j * M_PI) / num2)) + z_0;
                detectors.emplace_back(x, y, z, 0.0);
            }
        }
    }
    
    Point at(size_t i) const {
        return detectors[i];
    }
    
    size_t Size() const {
        return detectors.size();
    }
    
    
private:
    std::vector<Point> detectors;
};

#endif /* detectors_hpp */
