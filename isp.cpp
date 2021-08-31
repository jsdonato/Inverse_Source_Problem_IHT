#include <iostream>
#include <complex>
#include <getopt.h>
#include <fstream>
#include <memory>
#include <limits>
#include <sys/stat.h>
#include "lattice.hpp"
#include "detectors.hpp"
#include "point.hpp"
#include "matrix.hpp"
#include "utility.hpp"

class ISP{
public:
    ISP(int argc, char **argv){
        get_options(argc, argv);
        if (analysis) {
             mkdir("num_test_sources", 0777);
        }
        if (spherical){
            lattice = std::make_unique<SphericalLattice>(sphere_info);
        }
        else if (uniform) {
            lattice = std::make_unique<UniformLattice>();
        }
        detectors = Detectors(detectors_info);
        
        std::ifstream in;
        read_and_run(in);
    }

    void read_and_run(std::ifstream &in) {
        in.open(file_name);
        std::string entry;
        while (getline(in, entry)) {
            if (entry[0] == '#') {
                make_init_matrix();
                solver();
                std::cout << sim_num << " & ";
                print_sources();
                print_results();
            
                std::cout << Coherence << " \\\\ \\hline" << std::endl;
                sources.clear();
                sim_num++;
            }
            else {
                sources.emplace_back(parse_point(entry));
            }
        }
    }

    void get_options(int argc, char **argv){
        int option_index = 0, option = 0;
        struct option longOpts[] = {{"file", required_argument, nullptr, 'f'},
            {"detectors", required_argument, nullptr, 'd'},
            {"uniform_lattice", no_argument, nullptr, 'u'},
            {"spherical_lattice", required_argument, nullptr, 's'},
            {"analysis", no_argument, nullptr, 'a'},
            { nullptr, 0, nullptr, '\0'}};
        while ((option = getopt_long(argc, argv, "f:d:us:a", longOpts, &option_index)) != -1){
            switch (option){
                case 'f': {
                    file = true;
                    std::string str_f(optarg);
                    file_name = str_f;
                    break;
                }
                case 'd': {
                    detect = true;
                    std::string str_d(optarg);
                    detectors_info = parse_point(str_d);
                    break;
                }
                case 'u': {
                    uniform = true;
                    break;
                }
                case 's': {
                    spherical = true;
                    std::string str_s(optarg);
                    sphere_info = parse_point(str_s);
                    break;
                }
                case 'a': {
                    analysis = true;
                    break;
                }
                default:
                    exit(1);
            }
        }
        if (!file) {
            throw std::runtime_error("File name with sources is required.");
        }
        if (uniform == spherical) {
            throw std::runtime_error("Strictly one of uniform or spherical lattice must be chosen");
        }
        if (!detect) {
            throw std::runtime_error("Detector information is required");
        }

    }
    
    void print_sources() {
        for (const auto& s : sources) {
            std::cout << "(" <<s.Theta() << "," << s.Phi() << "," << s.Amp() << ")";
        }
        std::cout << " & ";
    }
    
    void make_init_matrix(){
        arma::cx_mat G_d = Gd(*lattice, detectors);
        arma::cx_mat G_s = Gs(*lattice, sources);
        arma::cx_mat V_ = V(*lattice);
        arma::cx_mat Ai_num = Ai(sources, detectors);
        
        
        
        A_num = arma::sum((G_d * V_ * G_s) + Ai_num, 1);
        
        A_s_part = G_d * V_;
        
        Coherence = coherence((G_d * V_ * G_s) + Ai_num);
    }
    
   void solver(){
        std::ofstream out;
        if (analysis) {
            out.open("num_test_sources/num_test_sources_" + std::to_string(sim_num) + ".txt");
        }
    
        delta = (M_PI / num_div);
       
        make_test_sources();
       
        arma::cx_mat G_s_test = Gs(*lattice, test_sources);
        arma::cx_mat A_i_test = Ai(test_sources, detectors);
        arma::cx_mat A_test = (A_s_part * G_s_test) + A_i_test;
        x = solve_regular(A_test);
       
        Cond = arma::cond(A_test);
       
        if (analysis) { 
            print_x(0);
        }
        
        for (int i = 0; i < num_iter; i++){
            
            delta *= (7.0 / 8.0); //should be greater than 0.71
            
            refine(i);
    
            G_s_test = Gs(*lattice, test_sources);
            A_i_test = Ai(test_sources, detectors);
            A_test = (A_s_part * G_s_test) + A_i_test;
            
            x = solve_regular(A_test);
            
            if (analysis){
                print_x(i + 1);
                out << i << "," << test_sources.size() << "\n";
            }
            
            threshold_max_of_groups(i);
            
            group_sources();
        
        }
       
       
        G_s_test = Gs(*lattice, test_sources);
        A_i_test = Ai(test_sources, detectors);
        A_test = (A_s_part * G_s_test) + A_i_test;
       
        x = solve_regular(A_test);
       
    }
    
    arma::cx_vec solve_regular(arma::cx_mat &A_test){
        arma::cx_mat I(A_test.n_cols, A_test.n_cols, arma::fill::eye);
        if (lambda == 0) {
            return solve(A_test, A_num);
        }
        return solve((A_test.t() * A_test) + (lambda * I), A_test.t() * A_num);
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
    
    void refine(const int &k){
        std::vector<Source> new_test_sources;
        
        for (int i = 0; i < test_sources.size(); i++){
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
        
        test_sources = new_test_sources;
    }
    
    void threshold_max_of_groups(const int &m){
        std::vector<Source> new_test_sources;
        
        std::vector<int> max_indicies;
        double sum = 0.0;
        for (int i = 0; i < test_sources.size(); i += 17){
            int end = i + 17;
            double max = 0.0;
            int max_index = 0;
            for (int j = i; j < end; j++){
                if (real(x(j)) > max){
                    max = real(x(j));
                    max_index = j;
                }
                
            }
            sum += max;
            max_indicies.push_back(max_index);
        }
        
        double threshold = threshold_weight * (sum / (double)max_indicies.size());
        for (int k = 0; k < max_indicies.size(); k++){
            if (real(x(max_indicies[k])) > threshold){
               new_test_sources.emplace_back(test_sources[max_indicies[k]].Theta(), test_sources[max_indicies[k]].Phi(), 1.0);
            }
        }
        
        test_sources = new_test_sources;
    }
    
    void group_sources(){
        std::vector<Source> new_test_sources;
        
        std::vector<double> vals_counters;
        for (int i = 0; i < test_sources.size(); i++){
            bool in = false;
            for (int j = 0; j < vals_counters.size(); j++){
                double d_T = test_sources[i].Theta() - new_test_sources[j].Theta();
                double d_P = test_sources[i].Phi() - new_test_sources[j].Phi();
                if (d_T < 0.001 && d_P < 0.001){
                    double CMA_theta = vals_counters[j] * new_test_sources[j].Theta();
                    double CMA_phi = vals_counters[j] * new_test_sources[j].Phi();
                    vals_counters[j]++;
                    new_test_sources[j].set_Theta((test_sources[i].Theta() + CMA_theta) / vals_counters[j]);
                    new_test_sources[j].set_Phi((test_sources[i].Phi() + CMA_phi) / vals_counters[j]);
                    in = true;
                }
            }
            if (!in){
                vals_counters.push_back(1);
                new_test_sources.emplace_back(test_sources[i].Theta(), test_sources[i].Phi(), 1.0);
            }
        }
        test_sources = new_test_sources;
        
    }
    
    void print_results(){
        std::vector<Source> result_sources;
        
        std::vector<double> vals_counters;
        for (int i = 0; i < test_sources.size(); i++){
            bool in = false;
            for (int j = 0; j < vals_counters.size(); j++){
                double d_T = test_sources[i].Theta() - result_sources[j].Theta();
                double d_P = test_sources[i].Phi() - result_sources[j].Phi();
                if (d_T < 0.1 && d_P < 0.1){
                    double CMA_theta = vals_counters[j] * result_sources[j].Theta();
                    double CMA_phi = vals_counters[j] * result_sources[j].Phi();
                    vals_counters[j]++;
                    result_sources[j].set_Theta((result_sources[i].Theta() + CMA_theta) / vals_counters[j]);
                    result_sources[j].set_Phi((result_sources[i].Phi() + CMA_phi) / vals_counters[j]);
                    result_sources[j].set_Amp(result_sources[j].Amp() + real(x(j)));
                    in = true;
                }
            }
            if (!in){
                vals_counters.push_back(1);
                result_sources.emplace_back(test_sources[i].Theta(), test_sources[i].Phi(), real(x(i)));
            }
        }
        for (const auto& s : result_sources) {
            std::cout << "(" << s.Theta() << "," << s.Phi() << "," << s.Amp() << ")";
        }
        
        print_data(result_sources);
        
        std::cout << " & " << Cond << " & ";
        
    }
    
    void print_data(std::vector<Source> result_sources){
        if (result_sources.size() != sources.size()){
            std::cout << " & inf & inf & inf";
        }
        
        else{
            std::vector<Source> match_sources;
            double d_theta = 0;
            double d_phi = 0;
            double d_amp = 0;
            for (const auto& s : sources) {
                size_t nearest_source_index = 0;
                double nearest_source_dist = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < result_sources.size(); i++) {
                    double num = (((result_sources[i].Theta() - s.Theta()) * (result_sources[i].Theta() - s.Theta())) +
                                  ((result_sources[i].Phi() - s.Phi()) * (result_sources[i].Phi() - s.Phi())));
                    if (num < nearest_source_dist) {
                        nearest_source_dist = num;
                        nearest_source_index = i;
                        
                    }
                }
                double num = abs(result_sources[nearest_source_index].Theta() - s.Theta());
                if (num > d_theta) {
                    d_theta = num;
                }
                num = abs(result_sources[nearest_source_index].Phi() - s.Phi());
                if (num > d_phi) {
                    d_phi = num;
                }
                num = abs(result_sources[nearest_source_index].Amp() - s.Amp());
                if (num > d_amp) {
                    d_amp = num;
                }
            }
            
            std::cout << " & " << d_theta << " & " <<  d_phi << " & " << d_amp;
        }
    }
    
    void print_x(const int &i){
        std::string dir = "IHT_output_" + std::to_string(sim_num);
        mkdir(dir.c_str(), 0777);
        std::ofstream out("IHT_output_" + std::to_string(sim_num) + "/output_" + std::to_string(sim_num) + "_" + std::to_string(i) + ".txt");
        
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
    std::unique_ptr<Lattice> lattice = nullptr;
    Detectors detectors = Detectors();
    
    std::vector<Source> test_sources;
    std::vector<Source> sources;
    std::vector<double> vec_x;
    
    arma::cx_vec A_num;
    arma::cx_mat A_s_part;
    arma::cx_vec x;

    std::string file_name;
    std::vector<double> sphere_info;
    std::vector<double> detectors_info;
    
    int sim_num=0;

    double delta;
    double d_theta;
    double Cond;
    double Coherence;

    bool file = false;
    bool detect = false;
    bool uniform = false;
    bool spherical = false;
    bool analysis = false;
};

int main(int argc, char **argv){
    ISP(argc, argv);
    return 0;
}
