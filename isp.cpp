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

class Grid{
public:
    Grid(size_t size_n, double center, double radius){
        grid.resize(size_n);
        
        vector<Point> grid_temp;
        
    }
private:
    vector<Point> grid;
};

int main(){
    cout << M_PI << endl;
    return 0;
}
