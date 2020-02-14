g++ -std=c++11 isp3.cpp -o isp3 -O3 -I ~/Documents/armadillo-7.600.2/include -DARMA_DONT_USE_WRAPPER -L/usr/local/opt/openblas/lib -lopenblas -llapack -pthread
set OMP_NUM_THREADS=10000000000000000
./isp3
