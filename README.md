# ISP
## Table of contents
* [General info](#general-info)
* [Setup](#setup)

## General info
The ISP or the Inverse Source Problem is a problem which considers the following;  Suppose we have an unknown integer amount of sources penetrating a medium we are concerned with and we are able to retreive the information of the escaping plane wave and we know some information about the medium.  From here, how do we recover the number of sources that penetrated the medium and their location?  In the scope of this problem we will refer to the location of these sources as simply the angles theta and phi in spherical coordinates.

## Setup
In the trivial case (there is no scattering) we consider the following mathematical setup.  First, for the sake of simplicity lets suppose that 0<theta<pi and 0<phi.  From here, we consider a vector of size N_d which well tall v_t where N_d is the number of our detectors we have placed around our medium.  v_t contains the data at each detector of the out going plane wave.  From here we then consider a matrix of size N_d x M which we'll call A_test where M is large.  More specifically, we will discretely split up phi and theta into small intervals and we end up with a lattic in "theta-phi" space.  The number of nodes in this lattice will be M.  Considering what we have thus far if we have two sources say at (theta_1, phi_1) and (theta_2, phi_2) we would expect that when solving A_test * x = v_t and plot the resulting x in theta/phi space that there will be peaks near or at (theta_1, phi_1) and (theta_2, phi_2) which will then allow us to determine the location of these sources.  Furthermore, we take things a step further and remove nodes from our theta/phi lattice which we know arent viable solutions and from there refine our search over and over so that we are ultimately left with our sources locations and their amplitudes  
