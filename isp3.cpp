///////////// NOTE: Documents/armadillo-9.400.3/include/armadillo_bits/config.hpp  -->  ARMA_OPENMP_THREADS (default) = 10 
/////////////ARMA_OPENMP_THREADS (current) = 10000000000000000
/////////////compile with: g++ -std=c++11 isp3.cpp -o isp3 -O3 -I ~/Documents/armadillo-7.600.2/include -DARMA_DONT_USE_WRAPPER -L/usr/local/opt/openblas/lib -lopenblas -llapack -pthread
/////////////Mac OS Accelerate framework has blas built in but openblas is faster 
/////////////Make sure to install lapack, openblas and openmp.  This can be done using "brew install lapack", "brew install openblas", and "brew install libomp" respectively. 
/////////////Make sure the -O3 flag is included when compiling.  This will allow for a tail call optimization.  The stack from previous tail function calls will be deleted preventing a seg fault.
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

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

const int num_squaresx = 10;
const int num_squaresy = 10;
const int num_squaresz = 10;
const double dL = 1.0 / (num_squaresx);

const int num_div = 100;

const double circle_centerx = .50;
const double circle_centery = .50;
const double circle_centerz = .50;
const double circle_radius = .25;

const double med_coeff = 1;

const double freq = 1;

double coherence(cx_mat A){
    double max = 0.0;
    for (int i = 0; i < A.n_cols; i++){
        for (int j = 0; j < A.n_cols; j++){
            if (i != j){
                complex<double> num = abs(cdot(A.col(i), A.col(j))) / (sqrt(cdot(A.col(j), A.col(j))) * sqrt(cdot(A.col(i), A.col(i))));
                if (real(num) > max){
                    max = real(num);
                }
            }
        }
    }
    return max;
}

void do_join(std::thread& t){t.join();}
void join_all(std::vector<std::thread>& v){std::for_each(v.begin(),v.end(),do_join);}

void join_all(vector<vector<pair<double, double> > > &vectors){
    int total = 0;
    for (int i = 1; i < vectors.size(); i++){
        if (vectors[i].size() != 0){
            total += vectors[i].size();
            vectors[0].reserve(total);
            move(vectors[i].begin(), vectors[i].end(), back_inserter(vectors[0]));
        }
    }
}

void initialize_grid(double (&grid)[num_squaresy + 1][num_squaresx + 1][num_squaresz + 1][4], const int num_squaresx, const int num_squaresy, const int num_squaresz){
    for (int j = 0; j <= num_squaresy; j++){
        for (int i = 0; i <= num_squaresx; i++){
            for (int k = 0; k <= num_squaresz; k++){
                grid[j][i][k][0] = 0;
                grid[j][i][k][1] = (i / static_cast<double>(num_squaresx));
                grid[j][i][k][2] = ((num_squaresy - j) / static_cast<double>(num_squaresy));
                grid[j][i][k][3] = ((num_squaresz - k) / static_cast<double>(num_squaresz));
                
            }
        }
    }
}


void make_circle(double (&grid)[num_squaresy + 1][num_squaresx + 1][num_squaresz + 1][4], double centerx, double centery, double centerz, double radius, double medcoeff){
    for (int j = 0; j <= num_squaresy; j++){
        for (int i = 0; i <= num_squaresx; i++){
            for (int k = 0; k <= num_squaresz; k++){
                if ((pow(grid[j][i][k][1] - centerx, 2) + pow(grid[j][i][k][2] - centery, 2) + pow(grid[j][i][k][3] - centerz, 2)) <= pow(radius, 2)){
                    grid[j][i][k][0] = medcoeff;
                }
            }
        }
    }
}

void make_detectors(vector<vector<double> > &detectors){
    int num1 = 40;
    int num2 = 20;
    detectors.reserve(num1 * num2);
    for (int i = 0; i < num1; i++){
        for (int j = 0; j < num2; j++){
            double x = (4 * cos((i * 2 * pi) / num1) * sin((j * pi) / num2)) + 0.5;
            double y = ((4 * sin((i * 2 * pi) / num1)) * sin((j * pi) / num2)) + 0.5;
            double z = (4 * cos((j * pi) / num2)) + 0.5;
            detectors.push_back({x, y, z});
        }
    }
}

void make_thetas_test(vector<pair<double, double> > &thetas_test){
    for (int j = 0; j <= num_div; j++){
        for (int i = 0; i <= num_div; i++){
            thetas_test.push_back({(i * pi) / num_div, (j * pi) / num_div});
        }
    }
}

double average8(double num1, double num2, double num3, double num4, double num5, double num6, double num7, double num8){
    return (num1 + num2 + num3 + num4 + num5 + num6 + num7 + num8) / 8;
}

double average2(double num1, double num2){
    return ((num1 + num2) / 2);
}

void make_gridn(double (&grid)[num_squaresy + 1][num_squaresx + 1][num_squaresz + 1][4], double (&gridn)[num_squaresy][num_squaresx][num_squaresz][4]){
    for (int j = 0; j < num_squaresy; j++){
        for (int i = 0; i < num_squaresx; i++){
            for (int k = 0; k < num_squaresz; k++){
                gridn[j][i][k][0] =
                average8(grid[j][i][k][0], grid[j + 1][i][k][0], grid[j][i + 1][k][0], grid[j][i][k + 1][0], grid[j + 1][i + 1][k][0], grid[j + 1][i][k + 1][0],
                         grid[j][i + 1][k + 1][0], grid[j + 1][i + 1][k + 1][0]);
                
                gridn[j][i][k][1] = average2(grid[j][i][k][1], grid[j][i + 1][k][1]);
                gridn[j][i][k][2] = average2(grid[j][i][k][2], grid[j + 1][i][k][1]);
                gridn[j][i][k][3] = average2(grid[j][i][k][3], grid[j][i][k + 1][1]);
            }
        }
    }
}

void make_space(double (&grid)[num_squaresy + 1][num_squaresx + 1][num_squaresz + 1][4], double (&gridn)[num_squaresy][num_squaresx][num_squaresz][4]){
    initialize_grid(grid, num_squaresx, num_squaresy, num_squaresz);
    make_circle(grid, circle_centerx, circle_centery, circle_centerz, circle_radius, med_coeff);
    make_gridn(grid, gridn);
}

complex<double> Ai(double rx, double ry, double rz, double theta, double phi){
    complex<double> num(0.0, freq * ((rx * cos(theta) * sin(phi)) + (ry * sin(theta) * sin(phi)) + (rz * cos(phi))));
    return exp(num);
}

complex<double> G(double rx, double ry, double rz, double rxs, double rys, double rzs){
    double mag = sqrt(pow(rx - rxs, 2) + pow(ry - rys, 2) + pow(rz - rzs, 2));
    double num = (1 / (4 * pi * mag));
    complex<double> numtop(0.0, freq * mag);
    return num * exp(numtop);
}

void make_Gd(cx_mat& Gd, double gridn[num_squaresy][num_squaresx][num_squaresz][4], vector<vector<double> > detectors){
    for (int m = 0; m < detectors.size(); m++){
        int n = 0;
        for (int i = 0; i < num_squaresy; i++){
            for (int j = 0; j < num_squaresx; j++){
                for (int k = 0; k < num_squaresz; k++){
                    Gd(m, n) = G(gridn[i][j][k][1], gridn[i][j][k][2], gridn[i][j][k][3], detectors[m][0], detectors[m][1], detectors[m][2]) * pow(dL, 2);
                    n++;
                }
            }
        }
    }
}

void make_V(cx_mat& V, double gridn[num_squaresy][num_squaresx][num_squaresz][4]){
    int n = 0;
    for (int i = 0; i < num_squaresy; i++){
        for (int j = 0; j < num_squaresx; j++){
            for (int k = 0; k < num_squaresz; k++){
                complex<double> num((pow(freq, 2) * gridn[i][j][k][0]), 0);
                V(n, n) = num;
                n++;
            }
        }
    }
}

void make_Gs(cx_mat& Gs, double gridn[num_squaresy][num_squaresx][num_squaresz][4], vector<vector<double> > thetas){
    for (int m = 0; m < thetas.size(); m++){
        int n = 0;
        for (int i = 0; i < num_squaresy; i++){
            for (int j = 0; j < num_squaresx; j++){
                for (int k = 0; k < num_squaresz; k++){
                    Gs(n, m) = thetas[m][2] * Ai(gridn[i][j][k][1], gridn[i][j][k][2], gridn[i][j][k][3], thetas[m][0], thetas[m][1]);
                    n++;
                }
            }
        }
    }
}

void make_Gv(cx_mat& Gv, double gridn[num_squaresy][num_squaresx][num_squaresz][4]){
    mat temp(3, num_squaresy * num_squaresx * num_squaresz);
    int n = 0;
    for (int i = 0; i < num_squaresy; i++){
        for (int j = 0; j < num_squaresx; j++){
            for (int k = 0; k < num_squaresz; k++){
                temp(0, n) = gridn[i][j][k][1];
                temp(1, n) = gridn[i][j][k][2];
                temp(2, n) = gridn[i][j][k][3];
                n++;
            }
        }
    }
    
    for (int j = 0; j < num_squaresy * num_squaresx * num_squaresz; j++){
        for (int m = 0; m < num_squaresy * num_squaresx * num_squaresz; m++){
            if (j != m){
                Gv(j, m) = G(temp(0, j), temp(1, j), temp(2, j), temp(0, m), temp(1, m), temp(2, m)) * dL;
            }
            else{
                complex<double> num((2.38 / (4 * pi)) * pow(dL, 2), (freq / (4 * pi)) * pow(dL, 3));
                Gv(j, m) = num;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////FUNCTIONS FOR TREE SEARCH/////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void make_Ai_num(cx_mat& Ai_num, vector<vector<double> > thetas_num, vector<vector<double> > detectors){
    for (int i = 0; i < detectors.size(); i++){
        complex<double> sum(0, 0);
        for (int j = 0; j < thetas_num.size(); j++){
            sum += thetas_num[j][2] * Ai(detectors[i][0], detectors[i][1], detectors[i][2], thetas_num[j][0], thetas_num[j][1]);
        }
        Ai_num(i, 0) = sum;
    }
}

void Ai_test_helper(cx_mat& Ai_test, vector<pair<double, double> > thetas_test, vector<vector<double> > detectors, int i){
    for (int y = 0; y < thetas_test.size(); y++){
        Ai_test(i, y) = Ai(detectors[i][0], detectors[i][1], detectors[i][2], thetas_test[y].first, thetas_test[y].second);
    }
}

void make_Ai_test(cx_mat& Ai_test, vector<pair<double, double> > thetas_test, vector<vector<double> > detectors){
    vector<thread> threads;
    for (int i = 0; i < detectors.size(); i++){
        threads.push_back(thread(Ai_test_helper, ref(Ai_test), thetas_test, detectors, i));
    }
    join_all(threads);
}

void Gs_test_helper(cx_mat& Gs, double gridn[num_squaresy][num_squaresx][num_squaresz][4], vector<pair<double, double> > thetas_test, int n, int i, int j, int k){
    for (int y = 0; y < thetas_test.size(); y++){
        Gs(n, y) = Ai(gridn[i][j][k][1], gridn[i][j][k][2], gridn[i][j][k][3], thetas_test[y].first, thetas_test[y].second);
    }
}

void make_Gs_test(cx_mat& Gs, double gridn[num_squaresy][num_squaresx][num_squaresz][4], vector<pair<double, double> > thetas_test){
    vector<thread> threads;
    int n = 0;
    for (int i = 0; i < num_squaresy; i++){
        for (int j = 0; j < num_squaresx; j++){
            for (int k = 0; k < num_squaresz; k++){
                threads.push_back(thread(Gs_test_helper, ref(Gs), gridn, thetas_test, n, i, j, k));
                n++;
            }
        }
    }
    join_all(threads);
}

double f_xx(cx_vec x, int n){
    return (real(x(n + 1)) - (2.0 * real(x(n))) + real(x(n - 1))) / (pow(pi / num_div, 2));
}

double f_yy(cx_vec x, int n){
    return (real(x(n + num_div + 1)) - (2.0 * real(x(n))) + real(x(n - num_div - 1))) / (pow(pi / num_div, 2));
}

double f_xy(cx_vec x, int n){
    return (real(x(n + num_div + 2)) - real(x(n - num_div)) - real(x(n + num_div)) + real(x(n - num_div - 2))) / (4 * pow(pi / num_div, 2));
}

double D(cx_vec x, int n){
    return (f_xx(x, n) * f_yy(x, n) - pow(f_xy(x, n), 2));
}

void find_2_largest(int n, cx_vec x, vector<pair<double, double> > thetas_test, vector<pair<double, double> > &new_vals){
    double max = 0.0;
    double max_2 = 0.0;
    double max_3 = 0.0;
    int max_index = 0;
    int max_2_index = 0;
    int max_3_index = 0;
    int end = n + 9;
    for (int j = n; j < end; j++){l
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
    new_vals.push_back({thetas_test[max_index].first, thetas_test[max_index].second});
    new_vals.push_back({thetas_test[max_2_index].first, thetas_test[max_2_index].second});
    new_vals.push_back({thetas_test[max_3_index].first, thetas_test[max_3_index].second});
}

void first_groupings(int n, cx_vec x, double delta, double mean_vec, vector<pair<double, double> > thetas_test, vector<pair<double, double> > &new_vals){
    int end = n + num_div - 1;
    for (int j = n; j < end; j++){
        if (real(x(j)) > (30 * mean_vec) && D(x, j) > 0.5 && f_xx(x, j) < -0.5){
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second});
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second});
            
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second - delta});
            
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second - delta});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second - delta});
        }
    }
}

void second_groupings(int n, cx_vec x, double delta, double mean_vec, vector<pair<double, double> > thetas_test, vector<pair<double, double> > &new_vals){
    unsigned long int end = n + 20;
    if (thetas_test.size() - n - 1 < 19){end = thetas_test.size();}
    for (int j = n; j < end; j++){
        if (real(x(j)) > 3 * mean_vec){
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second});
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second});
            
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first, thetas_test[j].second - delta});
            
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second + delta});
            new_vals.push_back({thetas_test[j].first + delta, thetas_test[j].second - delta});
            new_vals.push_back({thetas_test[j].first - delta, thetas_test[j].second - delta});
        }
    }
}

void print_result(vector<pair<double, double> > &thetas_test, cx_mat x){
    vector<pair<pair<double, double>, double> > vals_counters;
    vector<double> temp;
    for (int i = 0; i < x.n_elem; i++){
        bool in = false;
        for (int j = 0; j < vals_counters.size(); j++){
            if (abs(thetas_test[i].first - vals_counters[j].first.first) < .2 && abs(thetas_test[i].second - vals_counters[j].first.second) < .2){
                double CMA_theta = vals_counters[j].second * vals_counters[j].first.first;
                double CMA_phi = vals_counters[j].second * vals_counters[j].first.second;
                vals_counters[j].second++;
                vals_counters[j].first.first =  (thetas_test[i].first + CMA_theta) / vals_counters[j].second;
                vals_counters[j].first.second =  (thetas_test[i].second + CMA_phi) / vals_counters[j].second;
                temp[j] += real(x(i));
                in = true;
            }
        }
        if (!in){
            temp.push_back(real(x(i)));
            vals_counters.push_back({{thetas_test[i].first, thetas_test[i].second}, 1});
        }
    }
    for (int m = 0; m < vals_counters.size(); m++){
        cout << vals_counters[m].first.first << " " << vals_counters[m].first.second << " " << temp[m] << endl;
        
    }
}

void algo(vector<pair<double, double> > thetas_test, vector<vector<double> > detectors, cx_mat As_part, cx_mat A_num, double gridn[num_squaresy][num_squaresx][num_squaresz][4], int i){
    cx_mat Gs_test (num_squaresy * num_squaresx * num_squaresz, thetas_test.size());
    cx_mat Ai_test (detectors.size(), thetas_test.size());
    thread t14 (make_Gs_test, ref(Gs_test), gridn, thetas_test);
    thread t15 (make_Ai_test, ref(Ai_test), thetas_test, detectors);
    t15.join();
    t14.join();
    cx_mat As_test = As_part * Gs_test;
    cx_mat A_test = As_test + Ai_test;
    cout << coherence(A_test) << " ";
    cx_vec x = solve(A_test, A_num);
    
    if (i == 1){
        ofstream out("output.txt");
        for (int j = 0; j < x.n_elem; j++){
            out << thetas_test[j].first << " " << thetas_test[j].second << " " << real(x(j)) << endl;
        }
        out.close(); 
    }
    
    
    else if (i == 10){
        cout << endl << endl;
        print_result(thetas_test, x);
        return;
    }
    
    double mean_vec = real(mean(x));
    vector<vector<pair<double, double> > > vectors;
    if (i == 1){
        vector<thread> threads;
        double delta = (pi / (3 * num_div));
        unsigned long long int end = x.n_elem - num_div;
        for (int m = num_div + 2; m < end; m += (num_div + 1)){
            vectors.reserve(num_div - 1);
            vectors.push_back(vector<pair<double, double> >());
            threads.push_back(thread(first_groupings, m, x, delta, mean_vec, thetas_test, ref(vectors[vectors.size() - 1])));
        }
        join_all(threads);
    }
    else if (i % 2 == 0){
        vector<thread> threads;
        for (int m = 0; m < thetas_test.size(); m += 9){
            vectors.reserve(thetas_test.size() / 9);
            vectors.push_back(vector<pair<double, double> >());
            threads.push_back(thread(find_2_largest, m, x, thetas_test, ref(vectors[vectors.size() - 1])));
        }
        join_all(threads);
    }
    else{
        int k = (i - 3) / 2;
        double delta = (pow(2, k + 1) * pi / (pow(3, k + 2) * num_div));
        vector<thread> threads;
        for (int m = 0; m < thetas_test.size(); m += 20){
            vectors.reserve(ceil(thetas_test.size() / 20.0));
            vectors.push_back(vector<pair<double, double> >());
            threads.push_back(thread(second_groupings, m, x, delta, mean_vec, thetas_test, ref(vectors[vectors.size() - 1])));
        }
        join_all(threads);
        k++;
    }
    join_all(vectors);
    return algo(vectors[0], detectors, As_part, A_num, gridn, ++i);
}


int main() {
    
    double grid[num_squaresy + 1][num_squaresx + 1][num_squaresz + 1][4];
    double gridn[num_squaresy][num_squaresx][num_squaresz][4];
    vector<vector<double> > detectors;
    vector<pair<double, double> > thetas_test;
    thread t0 (make_space, ref(grid), ref(gridn));
    thread t1 (make_detectors, ref(detectors));
    thread t2 (make_thetas_test, ref(thetas_test));
    t0.join();
    t1.join();
    t2.join();
    
    vector<vector<double> > thetas_num;
    thetas_num.push_back({1.5, 1.5, 1.0});
    thetas_num.push_back({2.0, 2.0, .8});
    thetas_num.push_back({3.0, 3.0, .5});
    thetas_num.push_back({1.0, 1.0, .5});
    thetas_num.push_back({1.0, 1.5, .5});
    thetas_num.push_back({3.0, 1.5, 1.0});
    
    cx_mat Gd(detectors.size(), num_squaresy * num_squaresx * num_squaresz);
    cx_mat V(num_squaresy * num_squaresx * num_squaresz, num_squaresy * num_squaresx * num_squaresz, fill::zeros);
    cx_mat Gv(num_squaresy * num_squaresx * num_squaresz, num_squaresy * num_squaresx * num_squaresz);
    cx_mat Gs(num_squaresy * num_squaresx * num_squaresz, thetas_num.size());
    cx_mat Ai_num (detectors.size(), 1);
    thread t3 (make_Gd, ref(Gd), gridn, detectors);
    thread t4 (make_V, ref(V), gridn);
    thread t5 (make_Gv, ref(Gv), gridn);
    thread t6 (make_Gs, ref(Gs), gridn, thetas_num);
    thread t7 (make_Ai_num, ref(Ai_num), thetas_num, detectors);
    t3.join();
    cout << "MADE GD" << endl;
    t4.join();
    cout << "MADE V" << endl;
    t5.join();
    cout << "MADE GV" << endl;
    t6.join();
    cout << "MADE GS" << endl;
    t7.join();
    cout << "MADE Ai_NUM (NONSCATTERED DATA INPUT)" << endl;
    
    
    cx_mat Gd_V;
    cx_mat Gv_V;
    cx_mat Id;
    cx_mat in;
    thread t8 ([&Gd_V, Gd, V](){Gd_V = Gd * V;});
    thread t9 ([&Gv_V, Gv, V](){Gv_V = Gv * V;});
    thread t10 ([&Id](){Id = eye<cx_mat>(num_squaresy * num_squaresx * num_squaresz, num_squaresy * num_squaresx * num_squaresz);});
    t9.join();
    t10.join();
    thread t11 ([&in, Id, Gv_V](){in = inv(Id - Gv_V);});
    t8.join();
    t11.join();
    
    cx_mat As_part = Gd_V * in;
    cx_mat As_num = sum(As_part * Gs, 1);
    cout << "MADE AS_NUM (SCATTERED DATA INPUT)" << endl;
    cx_mat A_num = As_num + Ai_num;
    cout << "MADE A_NUM (DATA INPUT)" << endl << endl;
    
    
    algo(thetas_test, detectors, As_part, A_num, gridn, 1);
    return 0;
}



