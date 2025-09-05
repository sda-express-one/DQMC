#ifndef COMPUTATIONAL_METHODS_HPP
#define COMPUTATIONAL_METHODS_HPP
#include <iostream>
#include <iterator>
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"
#include "MC_data_structures.hpp"


// evaluates equality between two double precision values
inline bool isEqual(long double a, long double b, long double epsilon = 1e-9L) {return std::fabs(a - b) < epsilon;};

// returns eigenvalues and eigenvectors of LK matrix in correct order (band 0, 1 and 2)
inline void selectionRules(const double A_LK, const double B_LK, const double C_LK,
    Eigen::Matrix3d& eigenvector_matrix, Eigen::RowVector3d& eigenvalues){
    
    Eigen::Vector3i order;
    Eigen::Matrix3d eigenvectors_new;
    Eigen::RowVector3d eigenvalues_new;

    if(C_LK == 0 && A_LK == B_LK){
        return; // do nothing if C_LK_el == 0 and A_LK_el == B_LK_el
    }    
    else if(A_LK >= B_LK){
        order << 2, 1, 0;
        eigenvectors_new = eigenvector_matrix(Eigen::all, order);
        if(eigenvectors_new.determinant() > 0.01){
            eigenvector_matrix = eigenvectors_new;
            eigenvalues_new = eigenvalues(order);
            eigenvalues = eigenvalues_new;
            return;
        }
        else{
            order << 2, 0, 1;
            eigenvectors_new = eigenvector_matrix(Eigen::all, order);
            eigenvector_matrix = eigenvectors_new;
            eigenvalues_new = eigenvalues(order);
            eigenvalues = eigenvalues_new;
            return;
        }
    }
    else if(B_LK > A_LK){
        if(eigenvector_matrix.determinant() > 0.01){
            return;
        }
        else{
            order << 0, 2, 1;
            eigenvectors_new = eigenvector_matrix(Eigen::all, order);
            eigenvector_matrix = eigenvectors_new;
            eigenvalues_new = eigenvalues(order);
            eigenvalues = eigenvalues_new;
            return;
        }
    }
    else{
        return; // do nothing if error occurs
    }

    
};


// diagonalizes LK Hamiltonian and returns eigenvalues and eigenvectors in correct order
inline Eigen::Matrix<double, 4, 3> diagonalizeLKHamiltonian(const double kx, const double ky, const double kz,
    const double A_LK, const double B_LK, const double C_LK){
    
    Eigen::Matrix<double, 4, 3> result;
    
    // detects free propagators
    if(isEqual(kx,0) && isEqual(ky,0) && isEqual(kz,0)){
        result << -2, 1, 1,
                (1./3), 0, 0,
                0, (1./3), 0,
                0, 0, (1./3);
        return result; 
    }
    
    double k_modulus = std::sqrt(kx*kx + ky*ky + kz*kz);
    double kx_n = kx/k_modulus;
    double ky_n = ky/k_modulus;
    double kz_n = kz/k_modulus;

    // check if k vector has kx==ky or kx==kz or ky==kz since it may cause issues in eigenvector ordering
    if(isEqual(kx_n,ky_n) || isEqual(kx_n,kz_n) || isEqual(ky_n, kz_n)){
        result << -1, -1, -1,
                   0,  0,  0,
                   0,  0,  0,
                   0,  0,  0;
        return result; // return error values
    }

    // build LK Hamiltonian matrix (3x3)
    Eigen::Matrix3d LK_matrix;
    LK_matrix << A_LK*kx_n*kx_n + B_LK*(ky_n*ky_n + kz_n*kz_n), C_LK*kx_n*ky_n, C_LK*kx_n*kz_n,
                 C_LK*kx_n*ky_n, A_LK*ky_n*ky_n + B_LK*(kx_n*kx_n + kz_n*kz_n), C_LK*ky_n*kz_n,
                 C_LK*kx_n*kz_n, C_LK*ky_n*kz_n, A_LK*kz_n*kz_n + B_LK*(kx_n*kx_n + ky_n*ky_n);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver; 
    eigensolver.compute(LK_matrix); // compute eigenvalues and eigenvectors

    if(eigensolver.info() != 0){
        result << -1, -1, -1,
                   0,  0,  0,
                   0,  0,  0,
                   0,  0,  0;
        return result; // return error values
    }

    Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors(); // eigenvector matrix 
    Eigen::RowVector3d eigenvalues = eigensolver.eigenvalues().transpose(); // eigenvalues row vector (from smallest to largest

    selectionRules(A_LK, B_LK, C_LK, eigenvectors, eigenvalues); // reorder eigenvalues and eigenvectors according to band number

    result << eigenvalues, eigenvectors;
    return result;
};

// diagonalizes LK Hamiltonian and returns eigenvalues (no eigenvectors), for effective mass exact estimator
inline Eigen::RowVector3d diagonalizeLKHamiltonianEigenval(const double kx, const double ky, const double kz,
    const double A_LK, const double B_LK, const double C_LK){

    Eigen::RowVector3d eigenvalues;
        
    // LK Hamiltonian undefined for k=(0,0,0)
    if(isEqual(kx,0) && isEqual(ky,0) && isEqual(kz,0)){
        eigenvalues << -1, 1, 1;
        return eigenvalues; 
    }
    
    double k_modulus = std::sqrt(kx*kx + ky*ky + kz*kz);
    double kx_n = kx/k_modulus;
    double ky_n = ky/k_modulus;
    double kz_n = kz/k_modulus;


    // build LK Hamiltonian matrix (3x3)
    Eigen::Matrix3d LK_matrix;
    LK_matrix << A_LK*kx_n*kx_n + B_LK*(ky_n*ky_n + kz_n*kz_n), C_LK*kx_n*ky_n, C_LK*kx_n*kz_n,
                 C_LK*kx_n*ky_n, A_LK*ky_n*ky_n + B_LK*(kx_n*kx_n + kz_n*kz_n), C_LK*ky_n*kz_n,
                 C_LK*kx_n*kz_n, C_LK*ky_n*kz_n, A_LK*kz_n*kz_n + B_LK*(kx_n*kx_n + ky_n*ky_n);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver; 
    eigensolver.compute(LK_matrix, Eigen::EigenvaluesOnly); // compute eigenvalues and eigenvectors

    if(eigensolver.info() != 0){
        eigenvalues << -1, 1, 1;
        return eigenvalues; // return error values
    }

    Eigen::RowVector3d eigenvalues = eigensolver.eigenvalues().transpose(); // eigenvalues row vector (from smallest to largest

    return eigenvalues;
};

// electron effective mass from LK eigenvalue
inline double computeEffMassfromEigenval(double eigenval){
    return 1/(2*eigenval);
};

inline double computeEffMassSingleBand(const double kx, const double ky, const double kz, const double m_x, const double m_y, const double m_z){
    if(isEqual(kx,0) && isEqual(ky,0) && isEqual(kz,0)){
        return 1.0;
    }
    else{
        double k_norm = std::sqrt(kx*kx + ky*ky + kz*kz); 
        double kx_norm = kx/k_norm;
        double ky_norm = ky/k_norm;
        double kz_norm = kz/k_norm;

        double inv_mass = kx_norm*kx_norm/m_x + ky_norm*ky_norm/m_y + kz_norm*kz_norm/m_z;
        return (1./inv_mass);
    }
};

// electron dispersion
inline double electronEnergy(const double kx, const double ky, const double kz, const double eff_mass){
    double k_squared = kx*kx + ky*ky + kz*kz;
    return (k_squared/(2*eff_mass));
};

// calc total energy of external phonon, useful at diagram's ends
inline double extPhononEnergy(const int * num_ext_phonons, const double * phonon_modes, int num_phonon_modes){
    double ext_phonon_energy = 0;
    for(int i=0; i < num_phonon_modes; i++){
        ext_phonon_energy += num_ext_phonons[i]*phonon_modes[i];
    }
    return ext_phonon_energy;
}

// single phonon dispersion (constant)
inline double phononEnergy(const double * phonon_modes, int phonon_mode_index){
    return phonon_modes[phonon_mode_index];
};

// vertex evaluation
// strength term (phonon mode-exclusive part)
inline double vertexStrengthTerm(double wx, double wy, double wz, double V_BZ, double V_BvK, double phonon_mode, double born_effective_charge, double dielectric_const){
    return (1./(std::sqrt(wx*wx+wy*wy+wz*wz))*(4.*M_PI/V_BZ)*std::pow(2.*phonon_mode*V_BvK,-1./2)*born_effective_charge/dielectric_const);
};

// overlap term (electron band-dependent part)
// v1
inline double vertexOverlapTerm(Band band_one, Band band_two){
    return (band_one.c1*band_two.c1 + band_one.c2*band_two.c2 + band_one.c3*band_two.c3);
};
// v2
inline double vertexOverlapTerm(Band band_one, double c1_new, double c2_new, double c3_new){
    return (band_one.c1*c1_new + band_one.c2*c2_new + band_one.c3*c3_new);
};
// v3
inline double vertexOverlapTerm(Band band_one, Eigen::Vector3d band_two){
    return (band_one.c1*band_two(0) + band_one.c2*band_two(1) + band_one.c3*band_two(2));
};

// vertex square modulus (|V(q)|^2 real and non-negative, V(q) imaginary (i^2=-1))
inline double calcVertexSquareModulus(double strength, double overlap){
    return std::pow(strength*overlap,2);
};

#endif