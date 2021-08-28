#ifndef lattice_hpp
#define lattice_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"
#include "constants.hpp"




class Lattice{
public:
    //This constructor makes a lattice with uniform eta throughout.
    Lattice(){
        for (size_t i = 0; i < len; i++){
            for (size_t j = 0; j < len; j++){
                for (size_t k = 0; k < len; k++){
                    make_cell_uni(i, j, k);
                }
            }
        }
    }
    
    //This constructor makes a lattice with a spherical region filled set to eta and the rest set to zero.
    Lattice(double x_0, double y_0, double z_0, double radius){
        for (size_t i = 0; i < len; i++){
            for (size_t j = 0; j < len; j++){
                for (size_t k = 0; k < len; k++){
                    make_cell(i, j, k, x_0, y_0, z_0, radius);
                }
            }
        }
    }
    
    void make_cell(size_t i, size_t j, size_t k, double x_0, double y_0, double z_0, double radius){
        double x = (double)((2 * k) + 1)/(double)(2 * len);
        double y = (double)((2 * j) + 1)/(double)(2 * len);
        double z = (double)((2 * i) + 1)/(double)(2 * len);
        
        Point temp(x, y, z, 0);
        
        double sum = 0;
        for (size_t l = i; l < i + 2; l++){
            for (size_t m = j; m < j + 2; m++){
                for (size_t n = k; n < k + 2; n++){
                    double x = (double)n / (double)len;
                    double y = (double)m / (double)len;
                    double z = (double)l / (double)len;
                    
                    if (((x-x_0)*(x-x_0)) + ((y-y_0)*(y-y_0)) + ((z-z_0)*(z-z_0)) <= radius * radius){
                        sum += eta;
                    }
                }
            }
        }
        temp.set_Eta((double)sum / 8.0);
        lattice.push_back(temp);
    }
    
    void make_cell_uni(size_t i, size_t j, size_t k){
        double x = (double)((2 * k) + 1)/(double)(2 * len);
        double y = (double)((2 * j) + 1)/(double)(2 * len);
        double z = (double)((2 * i) + 1)/(double)(2 * len);
        
        Point temp(x, y, z, eta);
        lattice.push_back(temp);
    }
    
    Point at(size_t i) const {
        return lattice[i];
    }
    
    size_t Size() const {
        return lattice.size();
    }
    
private:
    std::vector<Point> lattice;
};

#endif /* lattice_hpp */
