#ifndef utility_hpp
#define utility_hpp

#include <stdio.h>
#include <cmath>
#include <complex>
#include <string>
#include <regex>
#include "source.hpp"
#include "point.hpp"

//Solution to wave equation with no scattering.
std::complex<double> Ai(const Point& p, double theta, double phi){
    std::complex<double> num(0.0, freq * ((p.X() * std::cos(theta) * std::sin(phi)) + (p.Y() * std::sin(theta) * std::sin(phi)) + (p.Z() * std::cos(phi))));
    return exp(num);
}

//Green's function
std::complex<double> G(const Point& p0, const Point& p1){
    double mag = std::sqrt(std::pow(p0.X() - p1.X(), 2) + std::pow(p0.Y() - p1.Y(), 2) + std::pow(p0.Z() - p1.Z(), 2));
    double num = (1 / (4 * M_PI * mag));
    std::complex<double> numtop(0.0, freq * mag);
    return num * exp(numtop);
}

//returns the Mutual Coherence of matrix.
//https://en.wikipedia.org/wiki/Mutual_coherence_(linear_algebra)
double coherence(arma::cx_mat A){
    std::vector<double> data(A.n_cols, 0);
    
    std::vector<double> cdots;
    for (int i = 0; i < A.n_cols; i++){
        cdots.push_back(real(sqrt(cdot(A.col(i), A.col(i)))));
    }
    
    for (int i = 0; i < A.n_cols; i++){
        double max = 0.0;
        std::complex<double> temp = cdots[i];
        for (int j = i; j < A.n_cols; j++){
            if (i != j){
                std::complex<double> num = abs(cdot(A.col(i), A.col(j))) / (cdots[j] * temp);
                if (real(num) > max){
                    max = real(num);
                }
            }
        }
        data[i] = max;
    }
    
    return *max_element(data.begin(), data.end());
}

//Returns the sum of the real parts of a complex vector.
double real_accumulate(arma::cx_vec x, int begin, int end){
    double sum = 0.0;
    for (int i = begin; i < end; i++){
        sum += real(x(i));
    }
    return sum;
}

//takes in a string of the form "(num1,num2,num3,...)" and returns the string converted to
//a std::vector object. 
std::vector<double> parse_point(std::string str) {
    std::regex re("[ (),\n]");
    std::sregex_token_iterator first{str.begin(), str.end(), re, -1}, last;
    std::vector<std::string> v{first, last};
    v.erase(remove(v.begin(), v.end(), "\0"), v.end());
    std::vector<double> vec;
    for (const auto& s : v) {
        vec.push_back(std::stod(s));
    }
    return vec;
}

#endif /* utility_hpp */
