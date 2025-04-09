#ifndef GREENFUNC_HPP
#define GREENFUNC_HPP 

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>

class GreenFunc{
    public:

    struct Propagator{
        double el_propagator_kx;
        double el_propagator_ky;
        double el_propagator_kz;
    };

    struct Vertex{
        double tau = 0.;
        int type = 0; // +1 outgoing, -1 incoming, 0 unassigned/extrema
        int linked = -1; // describes connection to other vertex of phonon propagator (-1 if not linked)
        double wx = 0.;
        double wy = 0.;
        double wz = 0.;
    };

    // constructor
    GreenFunc() = default;
    GreenFunc(long long int N_diags, double tau_max, double kx, double ky, double kz, double chem_potential, int order_max);

    // destructor
    ~GreenFunc(){
        delete[] _vertices;
        delete[] _propagators;
        if(_data_written){
            delete[] _tau_data;
        }
        if(_histogram_calculated){
            delete[] _histogram;
            delete[] _bin_count;
        }
        if(_histogram_normalized){delete[] _green_func;}
    };

    // getters
    long long int getN() const {return _N_diags;};
    int getRelaxSteps() const {return _relax_steps;};
    double * getData() const {return _tau_data;};
    double getNormConst() const {return _norm_const;};
    double getTauMax() const {return _tau_max;}
    double getkx() const {return _kx;};
    double getky() const {return _ky;};
    double getkz() const {return _kz;};
    double getChemPotential() const {return _chem_potential;}
    double getAlpha() const {return _alpha;}
    double getVolume() const {return _volume;}
    int getDimension() const {return _D;};
    int getOrderMax() const {return _order_max;}
    int getCurrentOrder() const {return _current_order;}
    int getN_bins() const {return _N_bins;}
    double * getHistogram() const {return _histogram;}
    int * getBinCount() const {return _bin_count;}

    // setters
    void setRelaxSteps(int relax_steps);
    void setAlpha(double alpha);
    void setVolume(double volume);
    void setDimension(int D);
    void setN_bins(int N_bins);
    void setNormConst(double norm_const);
    void setProbabilities(double p_length, double p_add);
    void setProbabilities();

    // main simulation
    void markovChain();
    void markovChain(int);

    // write time data to file
    void writeTauData(std::string) const;
    // write histogram to file
    void writeHistogram(std::string) const;

    // test variables
    int test = 0;
    int test_length = 0;
    int test_add = 0;
    int test_remove = 0;
    int test_control = 0;


    protected:

    void calcHistogram(); 
    void calcNormConst(); // calc normalization
    void normalizeHistogram();  // calc GF in histogram method

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
    const double _tau_max =  50.; // max value for imaginary time
    double _kx = 0.; // x momentum
    double _ky = 0.; // y momentum
    double _kz = 0.; // z momentum
    const int _order_max = 50; // max diagram order
    const double _chem_potential = -2.0; // chemical potential, normalization factor
    int _D = 3; // dimensions
    double _alpha = 2.0; // coupling strength
    double _volume = 1.0; // volume of 1BZ
    int _relax_steps = 10000000;

    // transition probabilities
    double _p_length = 1./3;
    double _p_add = 1./3;
    double _p_rem = 1./3;

    // important variables to keep
    double _last_vertex = 0.; // last current phonon vertex
    int _current_order = 0;

    // diagram backbone
    Vertex* _vertices; // array  of all possible vertices (also 0 and tau_max)
    Propagator* _propagators; // array of all possible bare electron propagators

    // diagrams data
    double* _tau_data;

    // histogram
    int _N_bins = 100; // number of bins for histogram
    double _bin_width = _tau_max/_N_bins; // width of each bin
    double _bin_center = _bin_width/2; // center of each bin
    double* _histogram; // histogram time lengths
    int* _bin_count; // number of diagrams in each bin
    int _N0 = 0; // number of diagrams of order 0

    // Green's function
    double _norm_const = 1.0;
    double* _green_func; 

    // flags for destructor
    bool _data_written = false;
    bool _histogram_calculated = false;
    bool _histogram_normalized = false;

    // fixes errors in input
    static inline int returnEven(int value){if(value%2==0){return value;}
    else{std::cout << "The order of a Diagram must be even, order is " << value << " +1." << std::endl; return value + 1;}
    };

    // evaluates equality between two double precision values
    static inline bool isEqual(double a, double b, double epsilon = 1e-9) {return std::fabs(a - b) < epsilon;}; 

    // computation methods
    static inline double calcEnergy(double px, double py, double pz){return (std::pow(px,2) + std::pow(py,2) + std::pow(pz,2))/2;}; // free electron energy
    inline double calcVertexStrength(double w_x, double w_y, double w_z){return (2.*std::sqrt(2.)*M_PI*_alpha/_volume)/(std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2));}; // 3D vertex strength (modulus squared)
    inline double calcVertexStrength(double w_x, double w_y){return (std::sqrt(2)*M_PI*_alpha/_volume)/std::sqrt(pow(w_x,2)+pow(w_y,2));}; // 2D vertex strength (modulus squared)
    inline void findLastPhVertex(){_last_vertex = _vertices[_current_order].tau;};
    inline double drawUniformR(){std::uniform_real_distribution<> distrib(0,1); double r = distrib(gen); return r;};

    // manage diagram
    int choosePhononPropagator();
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