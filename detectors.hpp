#ifndef detectors_hpp
#define detectors_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"
#include "constants.hpp"

class Detectors{
public:
    //Default constructor
    Detectors() : detectors({}) {}
    
    //This constructor constructs a sphere of evenly spaced detectors around the medium centered at (x_0, y_0, z_0)
    Detectors(double x_0, double y_0, double z_0, double radius){
        detectors.reserve(2 * detect_div_pi * detect_div_pi);
        for (int i = 0; i < 2 * detect_div_pi; i++){
            for (int j = 0; j < detect_div_pi; j++){
                double x = (radius * cos((i * 2 * M_PI) / (2 * detect_div_pi)) * sin((j * M_PI) / detect_div_pi)) + x_0;
                double y = (radius * sin((i * 2 * M_PI) / (2 * detect_div_pi)) * sin((j * M_PI) / detect_div_pi)) + y_0;
                double z = (radius * cos((j * M_PI) / detect_div_pi)) + z_0;
                detectors.emplace_back(x, y, z, 0.0);
            }
        }
    }

    //This constructor constructs a sphere of evenly spaced detectors around the medium centered at (vec[0], vec[1], vec[2])
    //and radius vec[3]
    Detectors(std::vector<double> vec){
        detectors.reserve(2 * detect_div_pi *  detect_div_pi);
        for (int i = 0; i < 2 * detect_div_pi; i++){
            for (int j = 0; j < detect_div_pi; j++){
                double x = (vec[3] * cos((i * 2 * M_PI) / (2 * detect_div_pi)) * sin((j * M_PI) / detect_div_pi)) + vec[0];
                double y = (vec[3] * sin((i * 2 * M_PI) / (2 * detect_div_pi)) * sin((j * M_PI) / detect_div_pi)) + vec[1];
                double z = (vec[3] * cos((j * M_PI) / detect_div_pi)) + vec[2];
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
