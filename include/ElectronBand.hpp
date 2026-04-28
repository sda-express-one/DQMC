#ifndef ELECTRONBAND_HPP
#define ELECTRONBAND_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
//#include "ProcessHandler.hpp"

template<typename ElectronBandType>
class ElectronBandEnergy {
    public:
    ElectronBandEnergy() = default;
    ElectronBandEnergy(const std::string& identifier) : _band_type(createElectronBandProcess(identifier)) {}
    inline double computeElectronEnergy(double kx, double ky, double kz) const {return _band_type.compute(kx, ky, kz);}
    inline double computeElectronEnergy(double kx, double ky, double kz, double effective_mass) const {return (kx*kx + ky*ky + kz*kz)/(2.0*effective_mass);}

    private:
    ElectronBandType _band_type;
};

class ParabolicBandSimple {
    public:
    ParabolicBandSimple() = default;
    ParabolicBandSimple(double effective_mass);
    inline double compute(double kx, double ky, double kz) const {
        return (kx*kx + ky*ky + kz*kz)/(2.0*_effective_mass);
    }

    private:
    const double _effective_mass = 1.0;
};

class ParabolicBandAnisotropic {
    public:
    ParabolicBandAnisotropic() = default;
    ParabolicBandAnisotropic(double mx_el, double my_el, double mz_el);

    inline double compute(double kx, double ky, double kz) const {
        return (kx*kx/(2.0*_mx_el) + ky*ky/(2.0*_my_el) + kz*kz/(2.0*_mz_el));
    }

    private:
    double _mx_el = 1.0;
    double _my_el = 1.0;
    double _mz_el = 1.0;
};

class ParabolicBandLK {
    public:
    ParabolicBandLK() = default;
    ParabolicBandLK(double A_LK, double B_LK, double C_LK);

    inline double compute(double kx, double ky, double kz) const {
        return 0;
    }

    private:
    double _A_LK = 0.0;
    double _B_LK = 0.0;
    double _C_LK = 0.0;
};

class FullBand {
    public:
    FullBand() = default;
    FullBand(const std::string& filename);
    FullBand(const std::string& filename, int num_kpoints, int num_bands);

    // destructor
    ~FullBand(){
        delete[] _bands_type;
        delete[] _kpoints_table;
        delete[] _energy_table;
    };
    
    inline double compute(double kx, double ky, double kz) const {
        int index = findClosestKPoint(kx, ky, kz);
        return 0;
    }

    private:
    std::string _cell_type;
    int * _bands_type = nullptr;
    double _lattice_param[3] = {1.0, 1.0, 1.0};
    Eigen::Matrix3d _unit_cell;
    Eigen::Matrix3d _unit_cell_rec;
    Eigen::Vector3d _k_lat_basis;
    int _num_kpoints = 1;
    int _num_bands = 1;
    double * _kpoints_table = nullptr;
    double * _energy_table = nullptr;

    void buildUnitCellRec(); // writes the unit cell matrix based on the cell type and lattice parameters

    inline void changeMomentumBasis(double kx, double ky, double kz) {
        // rewrite in momentum basis of reciprocal lattice vectors
        _k_lat_basis(0) = kx*_unit_cell_rec(0,0) + ky*_unit_cell_rec(1,0) + kz*_unit_cell_rec(2,0);
        _k_lat_basis(1) = kx*_unit_cell_rec(0,1) + ky*_unit_cell_rec(1,1) + kz*_unit_cell_rec(2,1);
        _k_lat_basis(2) = kx*_unit_cell_rec(0,2) + ky*_unit_cell_rec(1,2) + kz*_unit_cell_rec(2,2);

        // transform to first Brillouin zone
        _k_lat_basis(0) = _k_lat_basis(0) - std::round(_k_lat_basis(0));
        _k_lat_basis(1) = _k_lat_basis(1) - std::round(_k_lat_basis(1));
        _k_lat_basis(2) = _k_lat_basis(2) - std::round(_k_lat_basis(2));
    };

    int findClosestKPoint(double kx, double ky, double kz) const;
};

#endif