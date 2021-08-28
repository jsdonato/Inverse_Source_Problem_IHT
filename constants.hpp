#ifndef constants_hpp
#define constants_hpp

const double freq = 1.0; //Frequency of incoming sources.

const int num_div = 10; //The test sources in theta x phi are placed in a equally spaced num_div x num_div grid. 

//lattice

const int len = 10; //As the medium is discretized, len represents the length of the cubes making up the discretized mesh.

const double eta = 0; //eta (>0) is a measure of how much the medium scatters light with higher magnitude translating to more scattering. 

const double dL = 1.0 / (double)len;

#endif /* constants_hpp */
