#ifndef COMPUTATIONAL_METHODS_HPP
#define COMPUTATIONAL_METHODS_HPP
#include <iostream>
#include <iterator>
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"


// evaluates equality between two double precision values
static inline bool isEqual(long double a, long double b, long double epsilon = 1e-9L) {return std::fabs(a - b) < epsilon;};

// returns eigenvalues and eigenvectors of LK matrix in correct order (band 0, 1 and 2)
void selectionRules(const double A_LK, const double B_LK, const double C_LK,
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
Eigen::Matrix<double, 4, 3> diagonalizeLKHamiltonian(const double kx, const double ky, const double kz,
    const double A_LK, const double B_LK, const double C_LK){
    
    Eigen::Matrix<double, 4, 3> result;        
    
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

double electronEnergy(const double kx, const double ky, const double kz, const double eff_mass){
    double k_squared = kx*kx + ky*ky + kz*kz;
    return (k_squared/(2*eff_mass));
};

double extPhononEnergy(const int * num_ext_phonons, const double * phonon_modes, int num_phonon_modes){
    double ext_phonon_energy = 0;
    for(int i=0; i < num_phonon_modes; i++){
        ext_phonon_energy += num_ext_phonons[i]*phonon_modes[i];
    }
    return ext_phonon_energy;
}

double phononEnergy(const double * phonon_modes, int phonon_mode_index){
    return phonon_modes[phonon_mode_index];
};

#endif