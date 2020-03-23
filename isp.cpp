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

class Detectors{
public:
    Detectors(){
        size = 0;
    }
    
    Detectors(size_t size_n, double x_0, double y_0, double z_0, double radius){
        //TODO: implement this
    }
    
    size_t Size(){
        return size;
    }
    
    
private:
    vector<Point> detectors;
    size_t size;
};

class Lattice{
public:
    Lattice(){
        size = 0;
        len = 0;
    }
    
    Lattice(size_t len_n, double x_0, double y_0, double z_0, double radius){
        len = len_n;
        size = len * len *len;
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
        for (size_t l = 0; l < i + 2; i++){
            for (size_t m = 0; m < j + 2; j++){
                for (size_t n = 0; n < k + 2; k++){
                    double x = (double)n / (double)len;
                    double y = (double)m / (double)len;
                    double z = (double)l / (double)len;
                    
                    if (((x-x_0)*(x-x_0)) + ((y-y_0)*(y-y_0)) + ((z-z_0)*(z-z_0)) <= radius * radius){
                        sum += med_coeff;
                    }
                }
            }
        }
        temp.set_Eta((double)sum / 8.0);
        lattice.push_back(temp);
    }
    
    Point at(size_t i){
        return lattice[i];
    }
    
    size_t Size(){
        return size;
    }
    
private:
    vector<Point> lattice;
    size_t len;
    size_t size;
};



class Matrix{
public:
    Matrix(Lattice &l_n, Detectors &d_n){
        l = l_n;
        d = d_n;
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
        //TODO: implement this
    }
    
    cx_mat V(){
        cx_mat V_N(l.Size(), l.Size(), fill::zeros);
        for (int i = 0; i < l.Size(); i++){
            complex<double> num(freq * freq * l.at(i).Eta(), 0);
            V_N(i,i) = num;
        }
        return V_N;
    }
    
    cx_mat Gd(){
        //TODO: implement this
    }
    
    cx_mat Gs(){
        //TODO: implement this
    }

private:
    Lattice l = Lattice();
    Detectors d = Detectors();
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

