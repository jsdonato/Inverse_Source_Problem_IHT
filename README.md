# ISP

## General info
The ISP or the Inverse Source Problem is a problem which considers the following;  Suppose we have an unknown integer amount of sources penetrating a medium we are concerned with and we are able to retreive the information of the escaping plane wave and we know some information about the medium.  From here, how do we recover the number of sources that penetrated the medium and their location?  In the scope of this problem we will refer to the location of these sources as simply the angles theta and phi in spherical coordinates.

## Setup
In the trivial case (there is no scattering) we consider the following mathematical setup.  First, for the sake of simplicity lets suppose that 0<theta<pi and 0<phi<pi (The code in this repository works under these assumptions too).  From here, we consider a vector of size N_d which we'll call A_num where N_d is the number of our detectors we have placed around our medium.  A_num contains the data at each detector of the out going plane wave.  From here we then consider a matrix of size N_d x M which we'll call A_test where M is large.  More specifically, we will discretely split up phi and theta into small intervals and we end up with a lattice of test cources in "theta-phi" space.  The number of nodes in this lattice will be M.  Considering what we have thus far if we have two sources say at (theta_1, phi_1) and (theta_2, phi_2) we would expect that after solving A_test * x = A_num and plotting the resulting x in (theta, phi) space there will be peaks near or at (theta_1, phi_1) and (theta_2, phi_2) with height or magnitude equivalent to the amplitude of the sources.  In our algorithm, the process of solving A_test * x = A_num is repeated over and over as we refine the set of possible (theta, phi) by elements whose amplitudes dont exceed a certain threshold and adding elements closely surrounding elements which do exceed this threshold.

## Usage
This code is meant as a means to test the algorithm mentioned in the paper listed under the Results section.  The `constants.hpp` contains fundamental mathematical and procedural constants with descriptions.  The file labelled `input_sources.txt` (but can have an arbitrary name as we'll soon see) allows the user to input which sources they would like to test the algorithm against and for which simulations.  For example, if the user wishes to run one simulation with input sources coming from directions `(theta,phi)=(1.2,2.3)` and `(theta,phi)=(0.5,1.5)` with amplitudes `1` and `2` respectively and a second simulation with sources coming from `(theta,phi)=(0.2,0.3)` and `(theta,phi)=(1.1,1.3)` with amplitudes `1` and `1` respectively then the user would input the following in `input_sources.txt`.
```
(1.2,2.3,1)
(0.5,1.5,2)
#
(0.2,0.3,1)
(1.1,1.3,1)
#
```
Note that the input file must be in this format.  That is, each line of the input file is either a single source of the form `(theta,phi,amplitude)` or a colon `#` with a colon signifying that the user wishes to run a simulation with the sources listed above/before it and after the previous colon if there is one (Note that a colon must appear on the last line of the input file in order to run the last simulation listed).
### Command Line Flags
| Flag | Description |
| --- | --- |
| `-f`,`--file` | This flag is required and the input file with the test sources is placed after it. |
| `-d`,`--detectors` | The code places a sphere of detectors around our medium, a string of the form `x,y,z,radius` is required after this flag to specify the parameters of this sphere.  This flag is required.|
| `-u`,`--uniform_lattice` | This flag specifies that the medium has a uniform scattering coefficient `eta` (defined in `constants.hpp`) throughout.  Strictly one of the flags `-u` or `-s` is required.|
| `-s`,`--spherical_lattice` | This flag specifies that the medium has a spherical region in which the scattering coefficient of `eta` throughout it and the rest of the lattice has a scattering coefficient of `0`.  A string of the form `x,y,z,radius` is required after this flag.  Strictly one of the flags `-u` or `-s` is required.|
| `-a`,`--analysis` | This flag results in a directory named `num_test_sources` being made and inside of it are files named `num_test_sources_i.txt` (where `i` corresponds to the simulation number starting from `0`) which hold two columns of data pertaining to the step in the IHT algorithm and the number of test sources in that step.  In addition, directories named `IHT_output_i` (where `i` corresponds to the simulation number starting from `0`) are mode and in each directory are a series of files which output three columns of data pertaining to the `theta  phi  amplitude` of each element in the set of test sources for each step of the IHT algorithm. |
### Command Line Output
For each simulation a line is printed of the following form 

`sim_num & input sources & sources recovered by algorithm & d_theta & d_phi & d_amp & cond & coherence \\ \hline`

The following is the key for the parameters listed above.
| Parameter | Description |
| --- | --- |
| `sim num` | The simulation number. |
| `input sources` | The sources inputted into the input file by the user to run the algorithm against. |
| `sources recovered by algorithm` | The sources that the algorithm recovers after it is ran against the scattered wave data generated by `input sources`. |
| `d_theta` | The error present in the parameter `theta`.  If the number of recovered sources matches the number of input sources then this is computed by matching the recovered sources with the input sources and determing the max difference in `theta` among those.  Otherwise if the number of recovered sources does not match the number of test sources then `inf` is printed. |
| `d_phi` | The error present in the parameter `phi`.  The procedure used to compute `d_theta` is used to compute `d_phi` but within the context of the `phi` parameter. |
| `d_amp` | The error present in the parameter `amp`.  The procedure used to compute `d_theta` is used to compute `d_amp` but within the context of the amplitude (`amp`) parameter. |
| `cond` | The conditional number of `A_test` before the algorithm begins.|
| `coherence` | The mutual coherence of `A_num` before each row is accumulated together.|

NOTE: The reason the parameters in the output are separated by an ampersand `&` and end with a `\\ \hline` is so that the user may easily tabulate the data in LaTeX


## Requirements
The Armadillo C++ library for linear algebra is required to run this code.  Information regarding the documentation and installation of this library can be found at http://arma.sourceforge.net/. 

## Results
The results of the research surrounding this codebase can be found in the following paper

https://arxiv.org/pdf/2012.12783.pdf
