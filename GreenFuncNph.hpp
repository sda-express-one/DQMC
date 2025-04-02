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
        double tau = 0.;
        int type = 0; // +1 outgoing, -1 incoming, -1 unassigned (extrema)
        int linked = -1; // describes connection to other vertex of phonon propagator (-1 if not linked for external ph lines or extrema)
        double wx = 0.;
        double wy = 0.;
        double wz = 0.;
    };

    // constructor 
    GreenFuncNph() = default;
    GreenFuncNph(long long int N_diags, double tau_max, double kx, double ky, double kz, double chem_potential, int order_int_max, int ph_ext_max);

    // getters

    //setters
    void setRelaxSteps(int relax_steps);
    void setAlpha(double alpha);
    void setVolume(double volume);
    void setDimension(int D);
    void setN_bins(int N_bins);
    void setNormConst(double norm_const);

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
    const int _N_diags = 100000000; // number of differently generated diagrams
    int _N0 = 0; // number of diagrams of order 0
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
    int _relax_steps = 10000000; // steps of DQMC that are not taken into account into the full simulation, useful to relax to equilibrium distrib

    // transition probabilities
    double _p_length = 0.2;
    double _p_add_int = 0.2;
    double _p_rem_int = 0.2;
    double _p_add_ext = 0.2;
    double _p_rem_ext = 0.2;

    // important variables to keep
    double _last_vertex = 0.; // last current phonon vertex
    int _current_order_int = 0; // internal order of diagram (2*N_{ph_{int}}) 
    int _current_ph_ext = 0; // number of current external phonon lines in diagram (N_{ph_{ext}})

    // diagram backbone
    Vertex* _vertices; // array  of all possible vertices (also 0 and tau_max)
    Propagator* _propagators; // array of all possible bare electron propagators

    // diagrams data
    double* _tau_data;
    int* _order_data;

    // histogram
    int _N_bins = 100; // number of bins for histogram
    double _bin_width = _tau_max/_N_bins; // width of each bin
    double _bin_center = _bin_width/2; // center of each bin
    double* _histogram; // histogram time lengths
    int* _bin_count; // number of diagrams in each bin
    
    // Green's function
    double _norm_const = 1.0;
    double* _green_func;

    
    // flags for destructor
    bool _data_written = false;
    bool _histogram_calculated = false;
    bool _histogram_normalized = false;

    // fixes errors in input
    static inline int returnEven(int value){if(value%2==0){return value;}
        else{std::cout << "The order of a Diagram must be even, order is " << value << " +1." << std::endl; return value + 1;}};

    // evaluates equality between two double precision values
    static inline bool isEqual(double a, double b, double epsilon = 1e-9) {return std::fabs(a - b) < epsilon;};

    // computation methods
    // free electron energy
    static inline double calcEnergy(double px, double py, double pz){return (std::pow(px,2) + std::pow(py,2) + std::pow(pz,2))/2;}; 
     // 3D vertex strength (modulus squared)
    inline double calcVertexStrength(double w_x, double w_y, double w_z){return (2.*std::sqrt(2.)*M_PI*_alpha/_volume)/(pow(w_x,2)+pow(w_y,2)+pow(w_z,2));};
     // 2D vertex strength (modulus squared)
    inline double calcVertexStrength(double w_x, double w_y){return (std::sqrt(2)*M_PI*_alpha/_volume)/std::sqrt(pow(w_x,2)+pow(w_y,2));};
    // finds last vertex before diagram end
    inline void findLastPhVertex(){_last_vertex = _vertices[_current_order_int + 2*_current_ph_ext].tau;};
    // general purpose method for generating random number between 0 and 1
    inline double drawUniformR(){std::uniform_real_distribution<> distrib(0,1); double r = distrib(gen); return r;};

    // manage diagram
    int chooseInternalPhononPropagator();
    int findVertexPosition(double);
    void phVertexMakeRoom(int, int);
    void phVertexRemoveRoom(int, int);
    void propagatorArrayMakeRoom(int, int);
    void propagatorArrayRemoveRoom(int, int);

    // MCMC updates
    double diagramLengthUpdate(double);
    void addPhononPropagator();
    void removePhononPropagator();

    // debug methods
    void Diagnostic(std::string, int) const;
};

#endif