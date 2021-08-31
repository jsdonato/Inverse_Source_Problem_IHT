#ifndef lattice_hpp
#define lattice_hpp

#include <stdio.h>
#include <vector>
#include "point.hpp"
#include "constants.hpp"




class Lattice{
public:
    Lattice() : lattice({}) {}
    
    virtual void make_cell(size_t i, size_t j, size_t k) = 0;
    
    
    Point at(size_t i) const {
        return lattice[i];
    }
    
    size_t Size() const {
        return lattice.size();
    }
    
    virtual ~Lattice() = default;
    
protected:
    std::vector<Point> lattice;
};

class UniformLattice : public Lattice {
public:
    //This constructor makes a lattice with uniform eta throughout.
    UniformLattice(){
        for (size_t i = 0; i < len; i++){
            for (size_t j = 0; j < len; j++){
                for (size_t k = 0; k < len; k++){
                    make_cell(i, j, k);
                }
            }
        }
    }

    void make_cell(size_t i, size_t j, size_t k){
        double x = (double)((2 * k) + 1)/(double)(2 * len);
        double y = (double)((2 * j) + 1)/(double)(2 * len);
        double z = (double)((2 * i) + 1)/(double)(2 * len);
        
        Point temp(x, y, z, eta);
        lattice.push_back(temp);
    }

};

class SphericalLattice : public Lattice {
private:
    double x_0, y_0, z_0, radius;

public:
    //This constructor makes a lattice with a spherical region filled set to eta and the rest set to zero.
    //This sphere is defined with radius <radius> and center at (x_0,y_0,z_0).
    SphericalLattice(double x0, double y0, double z0, double rad)
        : x_0(x0), y_0(y0), z_0(z0), radius(rad) {
        for (size_t i = 0; i < len; i++){
            for (size_t j = 0; j < len; j++){
                for (size_t k = 0; k < len; k++){
                    make_cell(i, j, k);
                }
            }
        }
    }


    //This constructor makes a lattice with a spherical region filled set to eta and the rest set to zero.
    //This sphere is defined with radius vec[3] and center at (vec[0],vec[1],vec[2]).
    SphericalLattice(std::vector<double> vec)
        : x_0(vec[0]), y_0(vec[1]), z_0(vec[2]), radius(vec[3]) {
        for (size_t i = 0; i < len; i++){
            for (size_t j = 0; j < len; j++){
                for (size_t k = 0; k < len; k++){
                    make_cell(i, j, k);
                }
            }
        }
    }

    void make_cell(size_t i, size_t j, size_t k){
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

};

#endif /* lattice_hpp */
