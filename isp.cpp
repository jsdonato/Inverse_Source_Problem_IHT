#include <iostream>
#include <complex>
#include "lattice.hpp"
#include "detectors.hpp"
#include "point.hpp"
#include "matrix.hpp"

class ISP{
public:
    ISP(){
        for (int i = 0; i <= 15; i++){
            std::cout << i << " & ";
            d_theta = (0.5 * i) / 15.0;
            
            make_sources();
            make_init_matrix();
            solver();
            print_results();
            
            std::cout << Coherence << " \\\\ \\hline" << std::endl;
            
        }
    }
    
    void make_sources(){
        if (sources.size() > 0){
            sources.clear();
        }
        
        sources.emplace_back(1.0 + d_theta, 1.5, 1.0);
        sources.emplace_back(2.0 - d_theta, 1.5, 1.0);
        //sources.emplace_back(1.5, 1.5, 1.0);
        
        std::cout << sources[0].Theta() << " & " << sources[1].Theta() << " & ";
    }
    
    void make_init_matrix(){
        arma::cx_mat G_d = Gd(lattice, detectors);
        arma::cx_mat G_s = Gs(lattice, sources);
        arma::cx_mat V_ = V(lattice);
        arma::cx_mat Ai_num = Ai(sources, detectors);
        
        
        
        A_num = arma::sum((G_d * V_ * G_s) + Ai_num, 1);
        arma::cx_mat Zero(lattice.Size(), lattice.Size(), arma::fill::zeros);
        
        A_s_part = G_d * V_;
        //A_s_part = G_d * Zero;
        
        Coherence = coherence((G_d * V_ * G_s) + Ai_num);
    }
    
   void solver(){
        delta = (M_PI / num_div);
       
        make_test_sources();
       
        arma::cx_mat G_s_test = Gs(lattice, test_sources);
        arma::cx_mat A_i_test = Ai(test_sources, detectors);
        arma::cx_mat A_test = (A_s_part * G_s_test) + A_i_test;
        x = arma::solve(A_test, A_num);
       
        Cond = arma::cond(A_test);
       
        //print_x(0);
       
        //std::cout << "\n";
        
        for (int i = 0; i < 100; i++){
            
            delta *= (7.0 / 8.0); //should be greater than 0.71
            std::vector<Source> new_test_sources;
            if (i % 2 == 0){
                
                refine(i, new_test_sources);
            }
            
            else if (i % 2 == 1){
                find_largest(i, new_test_sources);
            }
            
            test_sources = new_test_sources;
            
            //std::cout << i << "," << test_sources.size() << "\n";
            
            G_s_test = Gs(lattice, test_sources);
            A_i_test = Ai(test_sources, detectors);
            A_test = (A_s_part * G_s_test) + A_i_test;
            
            x = solve_regular(A_test);
            
            make_vec_x();
            
            //print_x(i + 1);
            
            
        }
    }
    
    void make_test_sources(){
        if (test_sources.size() > 0){
            test_sources.clear();
        }
        
        for (int j = 0; j <= num_div; j++){
            for (int i = 0; i <= num_div; i++){
                test_sources.emplace_back((i * M_PI) / num_div, (j * M_PI) / num_div, 1.0);
            }
        }
    }
    
    void refine(const int &k, std::vector<Source> &new_test_sources){
        double mean = real(arma::mean(x));
        double threshold = 0.5 * mean;
        
        /*if (k >= 50){
            threshold *= 2;
        }*/
        
        for (int i = 0; i < test_sources.size(); i++){
            if (real(x(i)) > threshold /*&& real(x(i)) > 0.0*/){
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
                
                new_test_sources.emplace_back(theta + (0.5 * delta), phi, 1.0);
                new_test_sources.emplace_back(theta - (0.5 * delta), phi, 1.0);
                
                new_test_sources.emplace_back(theta, phi + (0.5 * delta), 1.0);
                new_test_sources.emplace_back(theta, phi - (0.5 * delta), 1.0);
                
                new_test_sources.emplace_back(theta + (0.5 * delta), phi + (0.5 * delta), 1.0);
                new_test_sources.emplace_back(theta + (0.5 * delta), phi - (0.5 * delta), 1.0);
                new_test_sources.emplace_back(theta - (0.5 * delta), phi + (0.5 * delta), 1.0);
                new_test_sources.emplace_back(theta - (0.5 * delta), phi - (0.5 * delta), 1.0);
             
            }
        }
    }
    
    void find_largest(const int &k, std::vector<Source> &new_test_sources){
        for (int i = 0; i < test_sources.size(); i += 17){
            int end = i + 17;
            double max = 0.0;
            double max_2 = 0.0;
            double max_3 = 0.0;
            int max_index = 0;
            int max_2_index = 0;
            int max_3_index = 0;
            for (int j = i; j < end; j++){
                if (real(x(j)) > max){
                    max_3 = max_2;
                    max_2 = max;
                    max = real(x(j));
                    max_3_index = max_2_index;
                    max_2_index = max_index;
                    max_index = j;
                }
                else if (real(x(j)) > max_2){
                    max_3 = max_2;
                    max_2 = real(x(j));
                    max_3_index = max_2_index;
                    max_2_index = j;
                }
                else if (real(x(j)) > max_3){
                    max_3 = real(x(j));
                    max_3_index = j;
                }
            }
            new_test_sources.emplace_back(test_sources[max_index].Theta(), test_sources[max_index].Phi(), 1.0);
            //new_test_sources.emplace_back(test_sources[max_2_index].Theta(), test_sources[max_2_index].Phi(), 1.0);
        }
    }
    
    arma::cx_vec solve_regular(arma::cx_mat &A_test){
        double lambda_sq = 10e-10;
        arma::cx_mat I(A_test.n_cols, A_test.n_cols, arma::fill::eye);
        
        return solve((A_test.t() * A_test) + (lambda_sq * I), A_test.t() * A_num);
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
            std::cout << "(" << vals_counters[m].first.first << "," << vals_counters[m].first.second << "," << temp[m] << ")";
        }
        
        print_data(vals_counters, temp);
        
        std::cout << " & " << Cond << " & ";
        
    }
    
    void print_data(std::vector<std::pair<std::pair<double, double>, double> > vals_counters, std::vector<double> temp){
        if (vals_counters.size() != 2){
            std::cout << " & inf & inf & inf";
        }
        
        else{
            int index_1 = 0;
            int index_2 = 1;
            
            if (vals_counters[1].first.first < vals_counters[0].first.first){
                index_1 = 1;
                index_2 = 0;
            }
            
            double d_theta = std::max(abs(vals_counters[index_1].first.first - sources[0].Theta()), abs(vals_counters[index_2].first.first - sources[1].Theta()));
            double d_phi = std::max(abs(vals_counters[index_1].first.second - sources[0].Phi()), abs(vals_counters[index_2].first.second - sources[1].Phi()));
            double d_amp = std::max(abs(temp[index_1] - 1), abs(temp[index_2] - 1));
            
            std::cout << " & " << d_theta << " & " <<  d_phi << " & " << d_amp;
        }
    }
    
    void print_x(const int &i){
        std::ofstream out("output_" + std::to_string(i) + ".txt");
        
        for (int j = 0; j < x.n_elem; j++){
            out << test_sources[j].Theta() << " " << test_sources[j].Phi() << " " << real(x(j)) << "\n";
        }
        
        out.close();
    }
    
    void make_vec_x(){
        vec_x.clear();
        
        for (int i = 0; i < x.n_elem; i++){
            vec_x.push_back(real(x(i)));
        }
        
    }

private:
    Lattice lattice = Lattice();
    Detectors detectors = Detectors(0.5, 0.5, 0.5, 4.0);
    
    std::vector<Source> test_sources;
    std::vector<Source> sources;
    std::vector<double> vec_x;
    
    arma::cx_vec A_num;
    arma::cx_mat A_s_part;
    arma::cx_vec x;
    
    double delta;
    
    double d_theta;
    double Cond;
    double Coherence;
};

int main(){
    ISP();
    
    return 0;
}

