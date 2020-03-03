#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <armadillo>
#include <thread>
#include <fstream>
using namespace std;
using namespace arma;

class Point{
public:
    Point(double x_n, double y_n, double z_n, double eta_n){
        x = x_n;
        y = y_n;
        z = z_n;
        eta = eta_n;
    }
    
    double  X(){
        return x;
    }
    
    double Y(){
        return y;
    }
    
    double Z(){
        return z;
    }
    
    double Eta(){
        return eta;
    }
    
private:
    double x;
    double y;
    double z;
    double eta;
};


class Source{
public:
    Source(double theta_n, double phi_n, double amplitude_n){
        theta = theta_n;
        phi = phi_n;
        amplitude = amplitude_n;
    }
    
    double T(){
        return theta;
    }
    
    double P(){
        return phi;
    }
    
    double A(){
        return amplitude;
    }

private:
    double theta;
    double phi;
    double amplitude;
};


class Lattice{
public:
    Lattice(size_t size_n, double x_0, double y_0, double z_0, double radius){
        vector<Point>lattice_temp;
        for (size_t i = 0; i < (size_n + 1); i++){
            for (size_t j = 0; j < (size_n + 1); j++){
                for (size_t k = 0; k < (size_n + 1); k++){
                    lattice_temp.emplace_back(k / size_n, j / size_n ,i / size_n, 0);
                    /*if (){
                        
                    }
                    else{
                        
                    }*/
                }
            }
        }
        
    }
private:
    vector<Point> lattice;
};


class Matrix{
public:

private:
    
};


class ISP{
public:
    ISP(){
        
    }

private:
    
};

int main(){
    
    return 0;
}
