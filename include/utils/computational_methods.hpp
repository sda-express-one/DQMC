#ifndef COMPUTATIONAL_METHODS_HPP
#define COMPUTATIONAL_METHODS_HPP
#include <iostream>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include "MC_data_structures.hpp"


namespace CompMethods{;

    struct choice{
        int position[3] = {0,1,2};
        int sign[3] = {1,1,1};
    };

    // evaluates equality between two double precision values
    inline bool isEqual(long double a, long double b, long double epsilon = 1e-9L) {return std::fabs(a - b) < epsilon;};

    inline void reorderChangeSignEigenvectors(Eigen::Matrix3d& eigenvectors, choice eigenv_gauge){

        Eigen::Matrix3d eigenvectors_new;
        Eigen::RowVector3d eigenvals_new;

        eigenvectors_new = eigenvectors;

        eigenvectors.row(0) = eigenvectors_new.row(eigenv_gauge.position[0]);
        eigenvectors.row(1) = eigenvectors_new.row(eigenv_gauge.position[1]);
        eigenvectors.row(2) = eigenvectors_new.row(eigenv_gauge.position[2]);

        for(int i = 0; i < 3; ++i){
            std::cout << "\n" << eigenvectors(eigenv_gauge.position[i],i) << std::endl;
            std::cout << eigenv_gauge.sign[i] << std::endl;
            std::cout << eigenv_gauge.sign[i] << std::endl;
        }
    
        if(eigenvectors(eigenv_gauge.position[0],0)*eigenv_gauge.sign[0] < 0){eigenvectors.col(0)*=-1;}
        if(eigenvectors(eigenv_gauge.position[1],1)*eigenv_gauge.sign[1] < 0){eigenvectors.col(1)*=-1;}
        if(eigenvectors(eigenv_gauge.position[2],2)*eigenv_gauge.sign[2] < 0){eigenvectors.col(2)*=-1;}
    };


    inline choice selectionRulesLK(const double kx, const double ky, const double kz){
        // decision container
        choice eigenv_gauge;

        if(std::abs(kx) >= std::abs(kz)){
            if(std::abs(kx) >= std::abs(ky)){
                if(std::abs(ky) >= std::abs(kz)){
                    eigenv_gauge.position[0] = 0;
                    eigenv_gauge.position[1] = 1;
                    eigenv_gauge.position[2] = 2;

                    if(kx < 0){eigenv_gauge.sign[0] = -1;}
                    if(ky < 0){eigenv_gauge.sign[1] = -1;}
                    if(kz < 0){eigenv_gauge.sign[2] = -1;}
                }
                else{
                    eigenv_gauge.position[0] = 0;
                    eigenv_gauge.position[1] = 2;
                    eigenv_gauge.position[2] = 1;

                    if(kx < 0){eigenv_gauge.sign[0] = -1;}
                    if(kz < 0){eigenv_gauge.sign[1] = -1;}
                    if(ky < 0){eigenv_gauge.sign[2] = -1;}
                }
            }
            else{
                eigenv_gauge.position[0] = 1;
                eigenv_gauge.position[1] = 0;
                eigenv_gauge.position[2] = 2;

                if(ky < 0){eigenv_gauge.sign[0] = -1;}
                if(kx < 0){eigenv_gauge.sign[1] = -1;}
                if(kz < 0){eigenv_gauge.sign[2] = -1;}
            }
        }
        else{
            if(std::abs(ky) >= std::abs(kz)){
                eigenv_gauge.position[0] = 1;
                eigenv_gauge.position[1] = 2;
                eigenv_gauge.position[2] = 0;

                if(ky < 0){eigenv_gauge.sign[0] = -1;}
                if(kz < 0){eigenv_gauge.sign[1] = -1;}
                if(kx < 0){eigenv_gauge.sign[2] = -1;}
            }
            else{
                if(std::abs(kx) >= std::abs(ky)){
                    eigenv_gauge.position[0] = 2;
                    eigenv_gauge.position[1] = 0;
                    eigenv_gauge.position[2] = 1;

                    if(kz < 0){eigenv_gauge.sign[0] = -1;}
                    if(kx < 0){eigenv_gauge.sign[1] = -1;}
                    if(ky < 0){eigenv_gauge.sign[2] = -1;}
                }
                else{
                    eigenv_gauge.position[0] = 2;
                    eigenv_gauge.position[1] = 1;
                    eigenv_gauge.position[2] = 0;

                    if(kz < 0){eigenv_gauge.sign[0] = -1;}
                    if(ky < 0){eigenv_gauge.sign[1] = -1;}
                    if(kx < 0){eigenv_gauge.sign[2] = -1;}
                }
            }
        }

        return eigenv_gauge;
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
        double k_vector[3] = {kx/k_modulus, ky/k_modulus, kz/k_modulus};
        double k_vector_temp[3] = {kx/k_modulus, ky/k_modulus, kz/k_modulus};

        // check if k vector has kx==ky or kx==kz or ky==kz since it may cause issues in eigenvector ordering
        /*if(isEqual(k_vector[0],ky_n) || isEqual(kx_n,kz_n) || isEqual(ky_n, kz_n)){
            result << -1, -1, -1,
                       0,  0,  0,
                       0,  0,  0,
                        0,  0,  0;
            return result; // return error values
        }*/

        choice eigenv_gauge = selectionRulesLK(kx, ky, kz); // cast k-values into IBZ

        k_vector[0] = std::abs(k_vector_temp[eigenv_gauge.position[0]]);
        k_vector[1] = std::abs(k_vector_temp[eigenv_gauge.position[1]]);
        k_vector[2] = std::abs(k_vector_temp[eigenv_gauge.position[2]]);

        // build LK Hamiltonian matrix (3x3)
        Eigen::Matrix3d LK_matrix;
        LK_matrix << A_LK*k_vector[0]*k_vector[0] + B_LK*(k_vector[1]*k_vector[1] + k_vector[2]*k_vector[2]), C_LK*k_vector[0]*k_vector[1], C_LK*k_vector[0]*k_vector[2],
                     C_LK*k_vector[0]*k_vector[1], A_LK*k_vector[1]*k_vector[1] + B_LK*(k_vector[0]*k_vector[0] + k_vector[2]*k_vector[2]), C_LK*k_vector[1]*k_vector[2],
                     C_LK*k_vector[0]*k_vector[2], C_LK*k_vector[1]*k_vector[2], A_LK*k_vector[2]*k_vector[2] + B_LK*(k_vector[0]*k_vector[0] + k_vector[1]*k_vector[1]);

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
        Eigen::RowVector3d eigenvalues = eigensolver.eigenvalues().transpose(); // eigenvalues row vector (from smallest to largest)

        // reorder eigenvectors and eigenvalues (greatest eigenvalue first band for A_LK > B_LK)
        if(A_LK > B_LK){
            double eigenval_temp;
            Eigen::Vector3d column;

            column = eigenvectors.col(2);
            eigenvectors.col(2) = eigenvectors.col(0);
            eigenvectors.col(0) = column;

            eigenval_temp = eigenvalues(0);
            eigenvalues(0) = eigenvalues(2);
            eigenvalues(2) = eigenval_temp;
        }

        // check all diagonal component > 0
        if(eigenvectors(0,0) < 0){eigenvectors.col(0)*=-1;}
        if(eigenvectors(1,1) < 0){eigenvectors.col(1)*=-1;}
        if(eigenvectors(2,2) < 0){eigenvectors.col(2)*=-1;}

        // transformation matrix from IBZ to generic sector
        Eigen::Matrix3d transformation_matrix;
        transformation_matrix << 0, 0, 0,
                                 0, 0, 0,
                                 0, 0, 0;
    
        // compute transformation matrix
        transformation_matrix(eigenv_gauge.position[0],0) = eigenv_gauge.sign[0];
        transformation_matrix(eigenv_gauge.position[1],1) = eigenv_gauge.sign[1];
        transformation_matrix(eigenv_gauge.position[2],2) = eigenv_gauge.sign[2];

        Eigen::Matrix3d new_eigenvectors;
        Eigen::Vector3d col_zero, col_one, col_two;

        // new eigenvectord
        col_zero = transformation_matrix*eigenvectors.col(0);
        col_one = transformation_matrix*eigenvectors.col(1);
        col_two = transformation_matrix*eigenvectors.col(2);

        // new eigenvector matrix
        new_eigenvectors.col(0) = col_zero;
        new_eigenvectors.col(1) = col_one;
        new_eigenvectors.col(2) = col_two;

        eigenvectors = new_eigenvectors;

        result << eigenvalues, eigenvectors;

        // return packed result
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
        eigensolver.compute(LK_matrix, Eigen::EigenvaluesOnly); // compute eigenvalues

        if(eigensolver.info() != 0){
            eigenvalues << -1, 1, 1;
            return eigenvalues; // return error values
        }

        eigenvalues = eigensolver.eigenvalues().transpose(); // eigenvalues row vector (from smallest to largest)

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
    /*inline double vertexStrengthTerm(double wx, double wy, double wz, double V_BZ, double V_BvK, double phonon_mode, double born_effective_charge, double dielectric_const){
        return (1./(std::sqrt(wx*wx+wy*wy+wz*wz))*(4.*M_PI/V_BZ)*std::pow(2.*phonon_mode*V_BvK,-1./2)*born_effective_charge/dielectric_const);
    };*/

    inline double coupling_strength(double phonon_mode, double dielectric_response, double effective_mass){
        return ((1./dielectric_response)*std::sqrt(effective_mass/(2*phonon_mode)));
    };

    inline double vertexStrengthTerm(double wx, double wy, double wz, double V_BZ, double V_BvK, double phonon_mode, double dielectric_response, double dielectric_const, double effective_mass){
        dielectric_const = 1;
        return(1./std::sqrt(wx*wx+wy*wy+wz*wz)*std::sqrt(2.*std::sqrt(2.)*M_PI*std::pow(phonon_mode,1.5)*coupling_strength(phonon_mode, dielectric_response, effective_mass)/(V_BZ*V_BvK*std::sqrt(effective_mass))))*dielectric_const;
    };

    // overlap term (electron band-dependent part)
    // two bands
    inline double vertexOverlapTerm(Band band_one, Band band_two){
        return (band_one.c1*band_two.c1 + band_one.c2*band_two.c2 + band_one.c3*band_two.c3);
    };
    // band and single values
    inline double vertexOverlapTerm(Band band_one, double c1_new, double c2_new, double c3_new){
        return (band_one.c1*c1_new + band_one.c2*c2_new + band_one.c3*c3_new);
    };
    // band + Eigen::Vector3d
    inline double vertexOverlapTerm(Band band_one, Eigen::Vector3d band_two){
        return (band_one.c1*band_two(0) + band_one.c2*band_two(1) + band_one.c3*band_two(2));
    };

    // vertex square modulus (|V(q)|^2 real and non-negative, V(q) imaginary (i^2=-1))
    inline double calcVertexSquareModulus(double strength, double overlap){
        return std::pow(strength*overlap,2);
    };

    inline unsigned long long int factorial(int n){
        if(n == 0 || n == 1){
            return 1;
        }
        else{
            unsigned long long int result = 1;
            for(int i=2; i != n+1; ++i){
                result *= i;
            }
            return result;
        }
    };

    inline long double computeMean(long double * data, int length){
        long double sum = 0.0L;
        for(int i=0; i<length; i++){
            sum += data[i];
        }
        return sum/static_cast<long double>(length);
    }

    inline double computeMean(double * data, int length){
        double sum = 0.0L;
        for(int i=0; i<length; i++){
            sum += data[i];
        }
        return sum/static_cast<double>(length);
    }

    inline long double computeStdDev(long double * data, long double mean, int length){
        long double sum = 0.0L;
        for(int i=0; i<length; i++){
            sum += (data[i]-mean)*(data[i]-mean);
        }
        return std::sqrt(sum/(static_cast<long double>(length - 1)));
    }

    inline double computeStdDev(double * data, double mean, int length){
        double sum = 0.0L;
        for(int i=0; i<length; i++){
            sum += (data[i]-mean)*(data[i]-mean);
        }
        return std::sqrt(sum/(static_cast<double>(length - 1)));
    }

    inline long double LaguerrePolynomial(int n, long double alpha, long double x){
        if(n == 0){return 1;}
        else if(n == 1){return (1 - alpha*x);}
        else{
            long double L_n_minus_two = 1;
            long double L_n_minus_one = (1 - alpha*x);
            long double L_n;
            for(int i=2; i != n+1; ++i){
                L_n = ((2*i - 1 - alpha * x)*L_n_minus_one - (i - 1)*L_n_minus_two)/static_cast<long double>(i);
                L_n_minus_two = L_n_minus_one;
                L_n_minus_one = L_n;
            }
            return L_n;
        }
    };

    inline long double LaguerrePolynomial(int n, long double alpha, long double L_n_minus_two, long double L_n_minus_one, long double x){
        if(n == 0){return 1;}
        else if(n == 1){return (1 - alpha*x);}
        else{
            long double L_n = ((2*n - 1 - alpha * x)*L_n_minus_one - (n - 1)*L_n_minus_two)/static_cast<long double>(n);
            return L_n;
        }
    };

}
#endif