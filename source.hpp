#ifndef source_hpp
#define source_hpp

#include <stdio.h>

class Source{
public:
    Source(double theta_n, double phi_n, double amplitude_n)
        : theta(theta_n), phi(phi_n), amplitude(amplitude_n) {}

    Source(std::vector<double> vec) 
        : theta(vec[0]), phi(vec[1]), amplitude(vec[2]) {}
    
    double Theta() const {
        return theta;
    }
    
    double Phi() const {
        return phi;
    }
    
    double Amp() const {
        return amplitude;
    }
    
    void set_Theta(double T_n){
        theta = T_n;
    }
    
    void set_Phi(double P_n){
        phi = P_n;
    }
    
    void set_Amp(double A_n){
        amplitude = A_n;
    }
    
private:
    double theta;
    double phi;
    double amplitude;
};


#endif /* source_hpp */
