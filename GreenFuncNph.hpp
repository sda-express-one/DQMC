#ifndef GREENFUNCNPH_HPP
#define GREENFUNCNPH_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>

class GreenFuncNph{
    public:

    struct Propagator{
        double el_propagator_kx;
        double el_propagator_ky;
        double el_propagator_kz;
    };

    struct Vertex{
        long double tau = 0.;
        int type = 0; // +1 outgoing and -1 incoming (internal), +2 outgoing and -2 incoming (external), 0 unassigned (extrema)
        int linked = -1; // describes connection to other vertex of phonon propagator (-1 if not linked for extrema)
        double wx = 0.;
        double wy = 0.;
        double wz = 0.;
    };

    // constructor 
    GreenFuncNph() = default;
    GreenFuncNph(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, double chem_potential, int order_int_max, int ph_ext_max);

    // destructor
    ~GreenFuncNph(){
        if(_gf_exact_written){;
            delete[] _points;
            delete[] _points_gf_exact;
        }
        if(_histogram_calculated){
            delete[] _histogram;
            delete[] _bin_count;
        }
        if(_histogram_normalized){delete[] _green_func;}
        if(_Z_factor_calculated){delete[] _Z_factor;}
        delete[] _vertices;
        delete[] _propagators;
    };

    // getters
    inline long long unsigned int getN() const {return _N_diags;};
    inline long long unsigned int getRelaxSteps() const {return _relax_steps;};
    //inline double * getData() const {return _tau_data;};
    inline double getNormConst() const {return _norm_const;};
    inline double getTauMax() const {return _tau_max;}
    inline double getkx() const {return _kx;};
    inline double getky() const {return _ky;};
    inline double getkz() const {return _kz;};
    inline double getChemPotential() const {return _chem_potential;}
    inline double getAlpha() const {return _alpha;}
    inline double getVolume() const {return _volume;}
    inline int getDimension() const {return _D;};
    inline int getN_bins() const {return _N_bins;};
    inline int getCurrentOrderInt() const {return _current_order_int;};
    inline int getCurrentPhExt() const {return _current_ph_ext;};
    inline long double getTauCutoffEnergy() const {return _tau_cutoff_energy;};
    inline long double getTauCutoffMass() const {return _tau_cutoff_mass;};
    inline double getGSEnergy() const {return _gs_energy;};
    inline double getEffectiveMass() const {return _effective_mass;};
    inline double getnumPoints() const {return _num_points;};
    inline int getSelectedOrder() const {return _selected_order;};
    inline long double getWidthEval()const {return _width_eval;};

    // setters
    void setRelaxSteps(int relax_steps);
    void setAlpha(double alpha);
    void setVolume(double volume);
    void setDimension(int D);
    void setN_bins(int N_bins);
    void setNormConst(double norm_const);
    void setTauCutoffEnergy(long double tau_cutoff_energy);
    void setTauCutoffMass(long double tau_cutoff_mass);
    void setNumPoints(int num_points);
    void setSelectedOrder(int selected_order);
    void setWidthEval(long double width_eval);
    void setProbabilities(double p_length, double p_add_int, double p_rem_int, double p_add_ext, double p_rem_ext, 
        double p_swap, double p_shift, double p_stretch);
    

    // main simulation method
    void markovChainMC(unsigned long long int N_diags, bool gf_exact, bool histo, bool gs_energy, bool effective_mass, bool Z_factor);

    // write to file
    inline void writeDiagrams(bool write_diagrams = false){_write_diagrams = write_diagrams;};
    void writeExactGF(std::string filename) const; // write GF with exact method in .txt file
    void writeHistogram(std::string) const; // write GF with histogram method in .txt file
    void writeZFactor(std::string filename) const; // write Z factor in .txt file

    private:

    std::mt19937 gen; // Mersenne Twister Algorithm, 32-bit

    // initialize seed for random number generator
    static std::mt19937::result_type setSeed(){
        std::mt19937::result_type seed = std::random_device()() ^ std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::system_clock::now().time_since_epoch()).count() 
        ^ std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        return seed;
    };

    // simulations features
    const unsigned long long int _N_diags = 100000000; // number of different generated diagrams
    long long unsigned int _N0 = 0; // number of diagrams of order 0
    const double _tau_max =  50.; // max value for imaginary time
    double _kx = 0.; // x momentum
    double _ky = 0.; // y momentum
    double _kz = 0.; // z momentum
    const double _chem_potential = -2.0; // chemical potential, normalization factor
    const int _order_int_max = 50; // max diagram order
    const int _ph_ext_max = 10; // max number of external phonon lines
    int _D = 3; // dimensions
    double _alpha = 2.0; // coupling strength
    double _volume = 1.0; // volume of 1BZ
    long long unsigned int _relax_steps = 10000000; // steps of DQMC that are not taken into account into the full simulation, useful to relax to equilibrium distrib

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

    // diagram backbone
    Vertex* _vertices; // array  of all possible vertices (also 0 and tau_max)
    Propagator* _propagators; // array of all possible bare electron propagators

    // Green function exact estimator
    int _num_points = 100;
    int _selected_order = 0;
    long double _points_step = _tau_max/_num_points;
    long double _points_center = _points_step/2;
    long double _width_eval = _tau_max/10;
    long double* _points; // array of evaluated points
    long double* _points_gf_exact; // gf values for evaluated points
    int _gf_exact_count = 0; // number of times the exact estimator is calculated

    // histogram method
    int _N_bins = 100; // number of bins for histogram
    double _bin_width = _tau_max/_N_bins; // width of each bin
    double _bin_center = _bin_width/2; // center of each bin
    double* _histogram; // histogram time lengths
    int* _bin_count; // number of diagrams in each bin
    double _norm_const = 1.0;
    double* _green_func;

    // direct estimator variables
    long double _tau_cutoff_energy = _tau_max/10; // cutoff for energy estimator, if tau < tau_cutoff energy estimator is not calculated
    double _gs_energy = 0.0; // ground state energy estimator of the system
    int _gs_energy_count = 0; // number of times the ground state energy estimator is calculated

    long double _tau_cutoff_mass = _tau_max/10; // cutoff for effective mass estimator, if tau < tau_cutoff mass estimator is not calculated
    double _effective_mass = 1.0; // effective mass of polaron (in electron mass units)
    double _effective_mass_count = 0; // number of times the effective mass estimator is calculated

    int* _Z_factor; // Z factor of polaron (overlap between free electron state and polaron state)
    
    // flags for destructor
    bool _gf_exact_written = false;
    bool _histogram_calculated = false;
    bool _histogram_normalized = false;
    bool _Z_factor_calculated = false;

    // fixes errors in input
    static inline int returnEven(int value){if(value%2==0){return value;}
        else{std::cout << "The order of a Diagram must be even, order is " << value << " + 1." << std::endl; return value + 1;}};

    // evaluates equality between two double precision values
    static inline bool isEqual(long double a, long double b, long double epsilon = 1e-9) {return std::fabs(a - b) < epsilon;};

    // computation methods
    // free electron energy
    static inline double calcEnergy(double px, double py, double pz){return (std::pow(px,2) + std::pow(py,2) + std::pow(pz,2))/2;};
     // 3D vertex strength (modulus squared)
    inline double calcVertexStrength(double w_x, double w_y, double w_z){return (2.*std::sqrt(2.)*M_PI*_alpha/_volume)/(std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2));};
     // 2D vertex strength (modulus squared)
    inline double calcVertexStrength(double w_x, double w_y){return (std::sqrt(2)*M_PI*_alpha/_volume)/std::sqrt(pow(w_x,2)+pow(w_y,2));};
    // finds last vertex before diagram end
    inline void findLastPhVertex(){_last_vertex = _vertices[_current_order_int + 2*_current_ph_ext].tau;};
    // general purpose method for generating random numbers between 0 and 1
    inline double drawUniformR(){std::uniform_real_distribution<long double> distrib(0,1); double r = distrib(gen); return r;};

    // manage diagram
    int chooseInternalPhononPropagator();
    int chooseExternalPhononPropagator();
    int findVertexPosition(long double tau);
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

    // histogram methods
    void calcNormConst();
    void normalizeHistogram();

    // estimator methods
    void exactEstimatorGF(double tau_length, int ext_phonon_order);
    double calcGroundStateEnergy(long double tau_length);
    double calcEffectiveMass(long double tau_length);
    void initializeZFactorArray();
    void updateZFactor();

    // debug methods
    void writeDiagram(std::string filename, int i, double r) const;
    void writeChosenUpdate(std::string filename, int i, double r) const;
    bool _write_diagrams = false;
};

#endif