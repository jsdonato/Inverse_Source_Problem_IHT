#include <iostream>
#include <complex>
#include "lattice.hpp"
#include "detectors.hpp"
#include "point.hpp"
#include "matrix.hpp"

class ISP{
public:
    ISP(){
        make_sources();
        make_init_matrix();
        solver();
        print_results();
    }
    
    void make_sources(){
        sources.emplace_back(1.0, 2.6, 2.0);
        sources.emplace_back(2.0, 1.5, 1.0);
    }
    
    void make_init_matrix(){
        arma::cx_mat G_d = Gd(lattice, detectors);
        arma::cx_mat G_s = Gs(lattice, sources);
        arma::cx_mat V_ = V(lattice);
        arma::cx_mat Ai_num = Ai(sources, detectors);
        
        
        A_num = arma::sum((G_d * V_ * G_s) + Ai_num, 1);
        A_s_part = G_d * V_;
    }
    
   void solver(){
        delta = (M_PI / (double)num_div);
       
        make_test_sources();
       
        arma::cx_mat G_s_test = Gs(lattice, test_sources);
        arma::cx_mat A_i_test = Ai(test_sources, detectors);
        arma::cx_mat A_test = (A_s_part * G_s_test) + A_i_test;
        
        x = arma::solve(A_test, A_num);
        
        for (int i = 0; i < 100; i++){
            std::cout << i << "\n";
            delta *= (7.0 / 8.0); //delta should be greater than 0.71
            std::vector<Source> new_test_sources;
            
            if (i % 2 == 0){
                refine(i, new_test_sources);
            }
            
            else if (i % 2 == 1){
                find_2_largest(new_test_sources);
            }
            
            test_sources = new_test_sources;
            
            G_s_test = Gs(lattice, test_sources);
            A_i_test = Ai(test_sources, detectors);
            A_test = (A_s_part * G_s_test) + A_i_test;
            
            x = solve(A_test, A_num);
            
        }
    }
    
    void make_test_sources(){
        for (int j = 0; j <= num_div; j++){
            for (int i = 0; i <= num_div; i++){
                test_sources.emplace_back((i * M_PI) / num_div, (j * M_PI) / num_div, 1.0);
            }
        }
    }
    
    void refine(const int &i, std::vector<Source> &new_test_sources){
        double threshold = 0.5 * real(arma::mean(x));
        
        for (int i = 0; i < test_sources.size(); i++){
            if (real(x(i)) > threshold){
                double theta = test_sources[i].Theta();
                double phi = test_sources[i].Phi();
                
                new_test_sources.emplace_back(theta, phi, 1.0);
                
                new_test_sources.emplace_back(theta + delta, phi, 1.0);
                new_test_sources.emplace_back(theta - delta, phi, 1.0);
                
                new_test_sources.emplace_back(theta, phi + delta, 1.0);
                new_test_sources.emplace_back(theta, phi - delta, 1.0);
                
                new_test_sources.emplace_back(theta + delta, phi + delta, 1.0);
                new_test_sources.emplace_back(theta + delta, phi - delta, 1.0);
                new_test_sources.emplace_back(theta - delta, phi + delta, 1.0);
                new_test_sources.emplace_back(theta - delta, phi - delta, 1.0);
            }
        }
    }
    
    void find_2_largest(std::vector<Source> &new_test_sources){
        for (int i = 0; i < test_sources.size(); i += 9){
            int end = i + 9;
            double max = 0.0;
            double max_2 = 0.0;
            int max_index = 0;
            int max_2_index = 0;
            for (int j = i; j < end; j++){
                if (real(x(j)) > max){
                    max_2 = max;
                    max = real(x(j));
                    max_2_index = max_index;
                    max_index = j;
                }
                else if (real(x(j)) > max_2){
                    max_2 = real(x(j));
                    max_2_index = j;
                }
            }
            new_test_sources.emplace_back(test_sources[max_index].Theta(), test_sources[max_index].Phi(), 1.0);
            new_test_sources.emplace_back(test_sources[max_2_index].Theta(), test_sources[max_2_index].Phi(), 1.0);
        }
    }
    
    void print_results(){
        std::vector<std::pair<std::pair<double, double>, double> > vals_counters;
        std::vector<double> temp;
        for (int i = 0; i < x.n_elem; i++){
            bool in = false;
            for (int j = 0; j < vals_counters.size(); j++){
                if (abs(test_sources[i].Theta() - vals_counters[j].first.first) < .2 && abs(test_sources[i].Phi() - vals_counters[j].first.second) < .2){
                    double CMA_theta = vals_counters[j].second * vals_counters[j].first.first;
                    double CMA_phi = vals_counters[j].second * vals_counters[j].first.second;
                    vals_counters[j].second++;
                    vals_counters[j].first.first =  (test_sources[i].Theta() + CMA_theta) / vals_counters[j].second;
                    vals_counters[j].first.second =  (test_sources[i].Phi() + CMA_phi) / vals_counters[j].second;
                    temp[j] += real(x(i));
                    in = true;
                }
            }
            if (!in){
                temp.push_back(real(x(i)));
                vals_counters.push_back({{test_sources[i].Theta(), test_sources[i].Phi()}, 1});
            }
        }
        for (int m = 0; m < vals_counters.size(); m++){
            std::cout << "(" << vals_counters[m].first.first << "," << vals_counters[m].first.second << "," << temp[m] << ")\n";
        }
    }

private:
    Lattice lattice = Lattice();
    Detectors detectors = Detectors(0.5, 0.5, 0.5, 4.0);
    
    std::vector<Source> test_sources;
    std::vector<Source> sources;
    
    arma::cx_vec A_num;
    arma::cx_mat A_s_part;
    arma::cx_vec x;
    
    double delta;
};

int main(){
    ISP();
    
    return 0;
}

