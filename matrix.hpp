#include <armadillo>
#include <complex>
#include "source.hpp"
#include "detectors.hpp"
#include "lattice.hpp"
#include "utility.hpp"

arma::cx_mat Gv(const Lattice &lattice){
    arma::cx_mat G_v(lattice.Size(), lattice.Size());
    for (int i = 0; i < lattice.Size(); i++){
        for (int j = 0; j < lattice.Size(); j++){
            if (i != j){
                G_v(i, j) = G(lattice.at(i).X(), lattice.at(i).Y(), lattice.at(i).Z(), lattice.at(j).X(), lattice.at(j).Y(), lattice.at(j).Z()) * dL;
            }
            else{
                std::complex<double> num((2.38 / (4 * M_PI)) * pow(dL, 2), (freq / (4 * M_PI)) * dL * dL * dL);
                G_v(i, j) = num;
            }
        }
    }
    
    return G_v;
}

arma::cx_mat V(const Lattice &lattice){
    arma::cx_mat V_(lattice.Size(), lattice.Size(), arma::fill::zeros);
    for (int i = 0; i < lattice.Size(); i++){
        std::complex<double> num(freq * freq * lattice.at(i).Eta(), 0);
        V_(i,i) = num;
    }
    return V_;
}

arma::cx_mat Gd(const Lattice &lattice, const Detectors &detectors){
    arma::cx_mat G_d(detectors.Size(), lattice.Size());
    for (int i = 0; i < detectors.Size(); i++){
        for (int j = 0; j < lattice.Size(); j++){
            G_d(i, j) = G(lattice.at(j).X(), lattice.at(j).Y(), lattice.at(j).Z(), detectors.at(i).X(), detectors.at(i).Y(), detectors.at(i).Z()) * dL * dL;
        }
    }
    return G_d;
}

arma::cx_mat Gs(const Lattice &lattice, const std::vector<Source> &sources){
    arma::cx_mat G_s(lattice.Size(), sources.size());
    for (size_t i = 0; i < sources.size(); i++){
        for (size_t j = 0; j < lattice.Size(); j++){
            G_s(j, i) = sources[i].Amp() * Ai(lattice.at(j).X(), lattice.at(j).Y(), lattice.at(j).Z(), sources[i].Theta(), sources[i].Phi());
        }
    }
    return G_s;
}

arma::cx_mat Ai(const std::vector<Source> &sources, const Detectors &detectors){
    arma::cx_mat A_i(detectors.Size(), sources.size());
    for (int i = 0; i < detectors.Size(); i++){
        for (int j = 0; j < sources.size(); j++){
            A_i(i, j) = sources[j].Amp() * Ai(detectors.at(i).X(), detectors.at(i).Y(), detectors.at(i).Z(), sources[j].Theta(), sources[j].Phi());
        }
    }
    return A_i;
}
