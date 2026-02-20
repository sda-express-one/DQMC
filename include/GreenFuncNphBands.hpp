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
#include <omp.h>
#include "utils/MC_data_structures.hpp"
#include "utils/progressbar.hpp"
#include "utils/MC_Benchmarking.hpp"
#include "utils/computational_methods.hpp"
#include <omp.h>
#include "../thirdparty/Eigen/Core"        // built with Eigen 3.4.0, download it from https://gitlab.com/libeigen/eigen/-/releases
#include "../thirdparty/Eigen/Eigenvalues" // add Eigen directory inside project directory to compile
#include "Diagram.hpp"

class GreenFuncNphBands : public Diagram {
    public:

    // constructor
    GreenFuncNphBands() = default;
    GreenFuncNphBands(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
        double chem_potential, int order_int_max, int ph_ext_max, int num_bands, int phonon_modes);
    
    GreenFuncNphBands(Propagator * propagators, Vertex * vertices, Band * bands,
        unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
        double chem_potential, int order_int_max, int ph_ext_max, int num_bands, int phonon_modes);

    // destructor
    ~GreenFuncNphBands(){
        delete[] _phonon_modes;
        delete[] _ext_phonon_type_num;
        delete[] _dielectric_responses;
        delete[] _bands;
        delete[] _new_taus;
        delete[] _bands_init;
        delete[] _bands_fin;

        if(_flags.histo){
            delete[] _histogram;
            delete[] _bin_count;
            delete[] _green_func;
        }
        if(_flags.gf_exact){;
            delete[] _points;
            delete[] _points_gf_exact;
        }
        if(_flags.blocking_analysis){
            if(_flags.gs_energy){
                delete[] _gs_energy_block_array;
            }
            if(_flags.effective_mass){
                delete[] _effective_mass_block_array;
                delete[] _effective_masses_block_array;
            }
        }
    };

    // getters
    // diagram
    Propagator getPropagator(int index) const {return _propagators[index];};
    Vertex getVertex(int index) const {return _vertices[index];};
    Band getBand(int index) const {return _bands[index];};
    int getCurrentOrderInt() const {return _current_order_int;};
    int getCurrentPhExt() const {return _current_ph_ext;};

    // electronic bands
    int getNumBands() const {return _num_bands;};
    double get_m_x_el() const {return _m_x_el;};
    double get_m_y_el() const {return _m_y_el;};
    double get_m_z_el() const {return _m_z_el;};
    double get_A_LK_el() const {return _A_LK_el;};
    double get_B_LK_el() const {return _B_LK_el;};
    double get_C_LK_el() const {return _C_LK_el;};

    // phonon modes
    int getNumPhononModes() const {return _num_phonon_modes;};
    double getPhononMode(int index) const {return _phonon_modes[index];};
    double getDielectricResponse(int index) const {return _dielectric_responses[index];};

    // other input quantities
    double get1BZVolume() const {return _V_BZ;};
    double getBvKVolume() const {return _V_BvK;}
    double getDielectricConst() const {return _dielectric_const;};
    
    // GS energy
    long double getGSEnergy() const {return _gs_energy;};
    long double getGSEnergyVar() const {return _gs_energy_var;};
    long double getTauCutoffEnergy() const {return _tau_cutoff_energy;};

    // effective mass
    long double getEffectiveMass() const {return _effective_mass;};
    long double getEffectiveMassVar() const {return _effective_mass_var;};
    void getEffectiveMasses(long double * effective_masses) const;
    void getEffectiveMassesVar(long double * effective_masses_var) const;
    long double getTauCutoffMass() const {return _tau_cutoff_mass;};

    // exact GF
    int getNumPoints() const {return _num_points;};
    int getGFSelectedOrder() const {return _selected_order;};
    void getGFExactPoints(long double * points, long double * gf_values) const;

    // histogram GF
    int getNumBins() const {return _N_bins;};
    void getHistogram(long double * histogram, long double * green_func) const;

    // blocking method
    int getNumBlocks() const {return _N_blocks;}

    // setters
    // electron bands
    void setEffectiveMasses(double m_x, double m_y, double m_z);
    void setLuttingerKohnParameters(double A_LK_el, double B_LK_el, double C_LK_el);

    // phonons modes
    void setPhononModes(double* phonon_modes);
    void setDielectricResponses(double* dielectric_responses);

    // other input quantities
    void set1BZVolume(double V_BZ);
    void setBvKVolume(double V_BvK);
    void setDielectricConstant(double dielectric_const);

    // set diagram backbone quantities
    void setCurrentOrderInt(int order_int);
    void setCurrentPhExt(int ph_ext);

    // MC updates probability
    void setProbabilities(double * probs);

    // parallelization settings
    void setMaster(bool master_mode = false);
    void setNumNodes(int num_nodes = 1);
    void setNumProcs(int num_procs = 1);

    // calculations performed
    void setCalculations(bool gf_exact, bool histo, bool gs_energy, bool effective_mass, bool Z_factor, bool blocking_analysis, bool fix_tau_value);

    // exact estimator
    // histogram method
    void setN_bins(int N_bins);

    // exact GF estimator method
    void setNumPoints(int num_points);
    void setSelectedOrder(int selected_order);

    // Energy estimator
    void setTauCutoffEnergy(long double tau_cutoff_energy);

    // mass estimator
    void setTauCutoffMass(long double tau_cutoff_mass);

    // blocking method
    void setNumBlocks(int N_blocks);

    // write to file
    // write diagrams
    inline void writeDiagrams(bool write_diagrams = false){_flags.write_diagrams = write_diagrams;};

    // time benchmarking
    inline void setBenchmarking(bool time_benchmark = false){_flags.time_benchmark = time_benchmark;};

    // MC statistics
    inline void setMCStatistics(bool mc_statistics = false){_flags.mc_statistics = mc_statistics;};
    void setTauCutoffStatistics(long double tau_cutoff_statistics);

    // main simulation method
    void markovChainMC();
    void markovChainMCOnlyRelax(); // only relax the system to equilibrium
    void markovChainMCOnlySample(); // only sample the system without relaxation

    // write to file
    void writeHistogram(const std::string& filename) const;
    void writeExactGF(const std::string& filename) const; // write GF with exact method in .txt file
    void writeMCStatistics(std::string filename) const;
    //void writeBoldStatistics(std::string filename, std::string type = "simulation") const;

    protected:

    // support arrays
    long double * _new_taus = nullptr;
    Band * _bands_init = nullptr;
    Band * _bands_fin = nullptr;

    private:

    // diagram features
    Band * _bands = nullptr;

    // cpu features
    int _num_procs = 1;
    int _num_nodes = 1;
    int _autocorr_steps = 0;
    bool _master = false;

    // simulations features
    long long unsigned int _N0 = 0; // number of diagrams of order 0
    int _D = 3; // dimensions
    int _num_bands = 1; // number of bands
    double _V_BZ = 1.0; // volume of 1BZ
    double _V_BvK = 1.0; // total volume taken into consideration (depends on sampling)
    double _dielectric_const = 1.0; // simple scalar for cubic materials (SC, FCC, BCC, diamond, zincblende, perovskite, etc.)

    // electron band effective mass
    // single band model
    double _m_x_el = 1.0;
    double _m_y_el = 1.0;
    double _m_z_el = 1.0;
    // 3-band model, Luttinger-Kohn parameters
    double _A_LK_el = 1.0;
    double _B_LK_el = 1.0;
    double _C_LK_el = 0.0;

    // phonon longitudinal optical modes
    int _num_phonon_modes = 1;
    double* _phonon_modes; // assumption: constant phonon dispersion in momentum spaces
    int* _ext_phonon_type_num;
    double* _dielectric_responses; // Born effective charges for each phonon mode

    // bold diagrammatic Monte Carlo
    int _current_sign = 1; // current sign of the diagram (bold method)
    long long int _num_negative_diagrams[8] = {0,0,0,0,0,0,0,0}; // number of updates that lead to negative diagram weight (bold method)
    long double _ratio_negative_updates = 0.0L; // ratio of negative updates over total updates (bold method)
    long double _average_sign = 1.0L; // average sign of the diagrams (bold method)

    // transition probabilities
    double _p_length = (1./8.);
    double _p_add_int = (1./8.);
    double _p_rem_int = (1./8.);
    double _p_add_ext = (1./8.);
    double _p_rem_ext = (1./8.);
    double _p_swap = (1./8.);
    double _p_shift = (1./8.);
    double _p_stretch = (1./8.);

    // important variables to keep
    long double _last_vertex = 0.; // last current phonon vertex (0 if no vertices are present)
    int _current_order_int = 0; // internal order of diagram (2*N_{ph_{int}}) 
    int _current_ph_ext = 0; // number of current external phonon lines in diagram (N_{ph_{ext}}) 

    Flags _flags; // flags for different calculations

    // block analysis variables
    int _N_blocks = 100;
    unsigned long long int _block_size = 1000000;

    // histogram method
    int _N_bins = 100; // number of bins for histogram
    double _bin_width = _tau_max/_N_bins; // width of each bin
    double _bin_center = _bin_width/2; // center of each bin
    double _bin_width_inv = 1./_bin_width;
    double* _histogram = nullptr; // histogram time points
    unsigned long long int* _bin_count = nullptr; // number of diagrams in each bin
    double* _green_func = nullptr;

    // direct estimator variables
    // renormalized gs energy
    long double _tau_cutoff_energy = _tau_max/10; // cutoff for energy estimator, if tau < tau_cutoff energy estimator is not calculated
    long double _gs_energy = 0.0; // ground state energy estimator of the system
    long double * _gs_energy_block_array = nullptr; // array of mean values of each block
    long double _gs_energy_var = 0.0; // variance of ground state energy estimator (block method)
    unsigned long long int _gs_energy_count = 0; // number of times the ground state energy estimator is calculated

    // renormalized effective mass (polaron)
    long double _tau_cutoff_mass = _tau_max/10; // cutoff for effective mass estimator, if tau < tau_cutoff mass estimator is not calculated
    long double _effective_mass = 0; // isotropic or (1,1,1) direction
    long double * _effective_mass_block_array = nullptr; // array of mean values of each block
    long double _effective_mass_var = 0.0; // variance of effective mass energy estimator
    long double _effective_masses[3] = {0, 0, 0}; // effective masses of polaron (in electron mass units, 1 band model)
    long double * _effective_masses_block_array = nullptr; // array of mean values of each block (mx, my and mz)
    long double _effective_masses_var[3] = {0, 0, 0}; // variance of effective masses of polaron
    Eigen::Matrix3d _effective_masses_bands;
    unsigned long long int _effective_mass_count = 0; // number of times the effective mass estimator is calculated

    // Green function exact estimator
    int _num_points = 100;
    int _selected_order = 0;
    long double _points_step = _tau_max/_num_points;
    long double _points_center = _points_step/2;
    long double* _points = nullptr; // array of evaluated points
    long double* _points_gf_exact = nullptr; // gf values for evaluated points
    unsigned long long int _gf_exact_count = 0;

    // collect MC statistics
    MC_Benchmarking * _benchmark_sim = nullptr; // time benchmarking object
    MC_Benchmarking * _benchmark_th = nullptr; // time benchmarking object for thermalization
    MC_Statistics _mc_statistics; // statistics of the simulation
    long double _tau_cutoff_statistics = 0.; // cutoff for statistics, if tau < tau_cutoff statistics is not calculated


    // manage diagram
    inline void findLastPhVertex(){_last_vertex = _vertices[_current_order_int + 2*_current_ph_ext].tau;};
    int chooseInternalPhononPropagator();
    int chooseExternalPhononPropagator();
    int findVertexPosition(long double tau);
    int * findVerticesPosition(long double tau_one, long double tau_two);
    void phVertexMakeRoom(int index_one, int index_two);
    void phVertexRemoveRoom(int index_one, int index_two);
    void propagatorArrayMakeRoom(int index_one, int index_two);
    void propagatorArrayRemoveRoom(int index_one, int index_two);
    void bandArrayMakeRoom(int index_one, int index_two);
    void bandArrayRemoveRoom(int index_one, int index_two);
    void updateExternalPhononTypes(int index);

    //int choosePhonon();

    // bold methods
    inline void updateSign(){_current_sign = _current_sign*-1;};
    inline void updateNegativeDiagrams(int update_index){if(_current_sign == -1){_num_negative_diagrams[update_index]++;}};
    inline void resetNegativeDiagrams(){for(int i = 0; i < 8; i++){_num_negative_diagrams[i] = 0;}};
    inline void computeRatioNegativeUpdates(long long int num_updates);

    // MCMC updates
    long double diagramLengthUpdate(long double tau_init);
    void addInternalPhononPropagator();
    void removeInternalPhononPropagator();
    void addExternalPhononPropagator();
    void removeExternalPhononPropagator();
    void swapPhononPropagator();
    void shiftPhononPropagator();
    long double stretchDiagramLength(long double tau_init);

    // main MC simulation methods
    long double configSimulation(long double tau_length);
    void configSimulationSilent();
    long double chooseUpdate(long double tau_length, double r, MC_Benchmarking * benchmark);
    void computeQuantities(long double tau_length, double r, int i);
    void computeFinalQuantities();
    void printGFExactEstimator();
    void printhistogramEstimator();
    void printGroundStateEnergyEstimator();
    void printEffectiveMassEstimator();
    void printMCStatistics();
    void printBoldStatistics(std::string type = "simulation");
    
    // histogram methods
    double calcNormConst();
    void normalizeHistogram(double norm_const);

    // exact estimator methods
    void exactEstimatorGF(long double tau_length, int ext_phonon_order);
    double groundStateEnergyExactEstimator(long double tau_length);
    void groundStateEnergyBlockEstimator(long double gs_energy);
    //void calcGroundStateEnergy(std::string filename);
    double effectiveMassExactEstimator(long double tau_length);
    void effectiveMassBlockEstimator(long double avg, long double xP, long double yP, long double zP);
    //void calcEffectiveMasses(std::string filename);

    // debug methods
    void writeDiagram(std::string filename, int i, double r) const;
    void writeChosenUpdate(std::string filename, int i, double r) const;
};

#endif