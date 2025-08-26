#ifndef GREENFUNCNPHBANDS_HPP
#define GREENFUNCNPHBANDS_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "MC_data_structures.hpp"
#include "Diagram.hpp"
#include "progressbar.hpp"
#include "MC_Benchmarking.hpp"
#include <Eigen/Core>        // built with Eigen 3.4.0, download it from https://gitlab.com/libeigen/eigen/-/releases
#include <Eigen/Eigenvalues> // add Eigen directory inside project directory to compile

class GreenFuncNphBands : public Diagram {
    public:

    // constructor
    GreenFuncNphBands() = default;
    GreenFuncNphBands(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int num_bands, int phonon_modes);

    // getters

    // setters
    // electron bands
    void setLongitudinalEffectiveMass(double mass_long_el);
    void setTransversalEffectiveMass(double mass_transv_el);
    void setLuttingerKohnParameters(double A_LK_el, double B_LK_el, double C_LK_el);
    // phonons modes
    void setPhononDispersions(double* phonon_dispersions);
    void setBornEffectiveCharges(double* born_effective_charges);
    // other input quantities
    void set1BZVolume(double V_BZ);
    void setBvKVolume(double V_BvK);
    void setDielectricConstant(double dielectric_const);


    // main simulation method
    void MarkovChainMC();

    private:

    // diagram features
    Band * _bands;

    // simulations features
    long long unsigned int _N0 = 0; // number of diagrams of order 0
    int _D = 3; // dimensions
    int _num_bands = 1; // number of bands
    double _V_BZ = 1.0; // volume of 1BZ
    double _V_BvK = 1.0; // total volume taken into consideration (depends on sampling)
    double _dielectric_const = 1.0; // simple scalar for cubic materials (SC, FCC, BCC, diamond, zincblende, perovskite, etc.)

    // electron band effective mass
    // single band model
    double _mass_long_el = 1.0;
    double _mass_transv_el = 1.0;
    // 3-band model, Luttinger-Kohn parameters
    double _A_LK_el = 1.0;
    double _B_LK_el = 1.0;
    double _C_LK_el = 0.0;

    // phonon longitudinal optical modes
    int _num_phonon_modes = 1;
    double* _phonon_dispersions; // we assume constant phonon dispersion in momentum spaces 
    double* _born_effective_charges; // Born effective charges for each phonon mode

    // transition probabilities
    double _p_length = 1./8.;
    double _p_add_int = 1./8.;
    double _p_rem_int = 1./8.;
    double _p_add_ext = 1./8.;
    double _p_rem_ext = 1./8.;
    double _p_swap = 1./8.;
    double _p_shift = 1./8.;
    double _p_stretch = 1./8.;

    // important variables to keep
    long double _last_vertex = 0.; // last current phonon vertex (0 if no vertices are present)
    int _current_order_int = 0; // internal order of diagram (2*N_{ph_{int}}) 
    int _current_ph_ext = 0; // number of current external phonon lines in diagram (N_{ph_{ext}})

    // vertices computation
    static inline double vertexStrengthTerm(double kx, double ky, double kz, double V_BZ, double V_BvK, double phonon_mode, 
        double born_effective_charge, double dielectric_const){
        return 1/(std::sqrt(kx*kx+ky*ky+kz*kz))*(4*M_PI/V_BZ)*std::pow(2*phonon_mode*V_BvK,-1/2)*born_effective_charge/dielectric_const;
    };
    static inline double vertexOverlapTerm(Band band_one, Band band_two, int num_bands){
        double overlap = band_one.c1*band_two.c1 + band_one.c2*band_two.c2 + band_one.c3*band_two.c3;
        return overlap;
    };
    static inline double calcVertexSquareModulus(double strength, double overlap){return std::pow(strength*overlap,2);}; // -|V(q)|^2, V(q) real (i^2=-1) 

    Flags _flags; // flags for different calculations

    // direct estimator variables
    // renormalized gs energy
    long double _tau_cutoff_energy = _tau_max/10; // cutoff for energy estimator, if tau < tau_cutoff energy estimator is not calculated
    long double _gs_energy = 0.0; // ground state energy estimator of the system
    unsigned long long int _gs_energy_count = 0; // number of times the ground state energy estimator is calculated
    // renormalized effective mass
    long double _tau_cutoff_mass = _tau_max/10; // cutoff for effective mass estimator, if tau < tau_cutoff mass estimator is not calculated
    long double _effective_mass = 0; // effective mass of polaron (in electron mass units)
    unsigned long long int _effective_mass_count = 0; // number of times the effective mass estimator is calculated

    // collect MC statistics
    MC_Benchmarking * _benchmark_sim; // time benchmarking object
    MC_Benchmarking * _benchmark_th; // time benchmarking object for thermalization
    MC_Statistics _mc_statistics; // statistics of the simulation
    long double _tau_cutoff_statistics = 0.; // cutoff for statistics, if tau < tau_cutoff statistics is not calculated

    // manage diagram
    int chooseInternalPhononPropagator();
    int chooseExternalPhononPropagator();
    int findVertexPosition(long double tau);
    int * findVerticesPosition(long double tau_one, long double tau_two);
    void phVertexMakeRoom(int index_one, int index_two);
    void phVertexRemoveRoom(int index_one, int index_two);
    void propagatorArrayMakeRoom(int index_one, int index_two);
    void propagatorArrayRemoveRoom(int index_one, int index_two);

    // MCMC updates
    long double diagramLengthUpdate(long double tau_init);
    void addInternalPhononPropagator();
    void removeInternalPhononPropagator();
    void addExternalPhononPropagator();
    void removeExternalPhononPropagator();
    void swapPhononPropagator();
    void shiftPhononPropagator();
    long double stretchDiagramLength(long double tau_init);

    // estimator methods
    double calcGroundStateEnergy(long double tau_length);
    double calcEffectiveMass(long double tau_length);
    
    // debug methods
    void writeDiagram(std::string filename, int i, double r) const;
    void writeChosenUpdate(std::string filename, int i, double r) const;
};

#endif