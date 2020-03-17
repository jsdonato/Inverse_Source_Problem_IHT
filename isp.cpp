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

const double freq = 1;
const double med_coeff = 1;

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
    
    void set_Eta(double eta_n){
        eta = eta_n;
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
                    
                    double x = (double)k / (double)size_n;
                    double y = (double)j / (double)size_n;
                    double z = (double)i / (double)size_n;
                    
                    lattice_temp.emplace_back(x, y , z, 0);
                    
                    if (((x-x_0)*(x-x_0)) + ((y-y_0)*(y-y_0)) + ((z-z_0)*(z-z_0)) <= radius * radius){
                        
                        lattice_temp[lattice_temp.size() - 1].set_Eta(med_coeff);
                    
                    }
                }
            }
        }
        
        for (size_t i = 0; i < size_n; i++){
            for (size_t j = 0; j < size_n; j++){
                for (size_t k = 0; k < size_n; k++){
                    //In this loop I will take the averages of eta we have in lattice_temp as we did in
                    //make_gridn in the previous routine
                    
                }
            }
        }
        
    }
private:
    vector<Point> lattice;
};


class Matrix{
public:
    Matrix(Lattice &l){
        
    }
    
    complex<double> Ai(double rx, double ry, double rz, double theta, double phi){
        complex<double> num(0.0, freq * ((rx * cos(theta) * sin(phi)) + (ry * sin(theta) * sin(phi)) + (rz * cos(phi))));
        return exp(num);
    }
    
    complex<double> G(double rx, double ry, double rz, double rxs, double rys, double rzs){
        double mag = sqrt(pow(rx - rxs, 2) + pow(ry - rys, 2) + pow(rz - rzs, 2));
        double num = (1 / (4 * M_PI * mag));
        complex<double> numtop(0.0, freq * mag);
        return num * exp(numtop);
    }
    
    cx_mat Gv(){
        
    }
    
    cx_mat V(){
        
    }
    
    cx_mat Gd(){
        
    }
    
    cx_mat Gs(){
        
    }

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
