#ifndef constants_hpp
#define constants_hpp

//sources

const double freq = 1.0; //Frequency of incoming sources.

const int num_div = 10; //The test sources in theta x phi are placed in a equally spaced num_div x num_div grid. 

//lattice

const int len = 10; //As the medium is discretized, the number of "cuts" or "discretizations" on a single side of the medium.
                    //Note that the medium is a unit cube in R^3 defined by [0,1]x[0,1]x[0,1].

const double eta = 0; //eta (>0) is a measure of how much the medium scatters light with higher magnitude translating to more scattering. 

const double dL = 1.0 / (double)len; //As the medium is discretized, dL represents the length of the cubes making up the discretized mesh.

//detectors

const int detect_div_pi = 20; //

//IHT

const int num_iter = 100; //Number of iterations for the Iterative Hard Thresholding Algorithm.

const double lambda = 10e-9; //The constant \lambda which is used in Tikhonov regularization.
                             //If lambda is set to zero then a regular solve is done instead.

const double threshold_weight = 0.15 //This number describes the threshold for the IHT algorithm.
                                     //In particular, this weight is multiplied by the average
                                     //of the maximums for each grouping.

#endif /* constants_hpp */
