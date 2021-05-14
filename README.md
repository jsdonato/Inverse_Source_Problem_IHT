# ISP

## General info
The ISP or the Inverse Source Problem is a problem which considers the following;  Suppose we have an unknown integer amount of sources penetrating a medium we are concerned with and we are able to retreive the information of the escaping plane wave and we know some information about the medium.  From here, how do we recover the number of sources that penetrated the medium and their location?  In the scope of this problem we will refer to the location of these sources as simply the angles theta and phi in spherical coordinates.

## Setup
In the trivial case (there is no scattering) we consider the following mathematical setup.  First, for the sake of simplicity lets suppose that 0<theta<pi and 0<phi<pi.  From here, we consider a vector of size N_d which well tall v_t where N_d is the number of our detectors we have placed around our medium.  v_t contains the data at each detector of the out going plane wave.  From here we then consider a matrix of size N_d x M which we'll call A_test where M is large.  More specifically, we will discretely split up phi and theta into small intervals and we end up with a lattic in "theta-phi" space.  The number of nodes in this lattice will be M.  Considering what we have thus far if we have two sources say at (theta_1, phi_1) and (theta_2, phi_2) we would expect that after solving A_test * x = v_t and plotting the resulting x in (theta, phi) space there will be peaks near or at (theta_1, phi_1) and (theta_2, phi_2) with height or magnitude equivalent to the amplitude of the sources.  In our algorithm, the process of solving A_test * x = v_t is repeated over and over as we refine the set of possible (theta, phi) by elemtns whose amplitudes dont exceed a certain threshold and adding elements closely surrounding elements which do exceed this threshold. 

## Results
The results of the research surrounding this codebase can be found in the following paper

https://arxiv.org/pdf/2012.12783.pdf
