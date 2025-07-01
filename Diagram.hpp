#ifndef DIAGRAM_HPP
#define DIAGRAM_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include "MC_data_structures.hpp"

class Diagram {
    public:
        // constructor
        Diagram() = default;
        Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, 
            double chem_potential, int order_int_max, int ph_ext_max);
        
        // destructor
        ~Diagram(){
            delete[] _vertices;
            delete[] _propagators;
        }

        // evaluates equality between two double precision values
        static inline bool isEqual(long double a, long double b, long double epsilon = 1e-9L) {return std::fabs(a - b) < epsilon;};

        // returns uniform random double precision value between 0 and 1
        inline long double drawUniformR(){std::uniform_real_distribution<long double> distrib(0,1);long double r = distrib(gen); return r;};

        // Metropolis-Hastings Monte Carlo method
        inline int Metropolis(double R){
            double acceptance_ratio = std::min(1.0, R);
            if(drawUniformR() > acceptance_ratio){return 0;} // reject
            else{return 1;} // accept
        };

        // setters
        void setRelaxSteps(int relax_steps);

        // getters
        inline long long unsigned int getNdiags() const {return _N_diags;};
        inline long long unsigned int getRelaxSteps() const {return _N_relax_steps;};
        inline long double getTauMax() const {return _tau_max;};
        inline double getkx() const {return _kx;};
        inline double getky() const {return _ky;};
        inline double getkz() const {return _kz;};
        inline double getChemPotential() const {return _chem_potential;};
        inline int getOrderIntMax() const {return _order_int_max;};
        inline int getPhExtMax() const {return _ph_ext_max;};
        inline int getOrderMax() const {return _order_int_max + 2*_ph_ext_max;};
        inline std::mt19937 getSeed() const {return gen;};

    private:
        // number of diagrams generated
        const unsigned long long int _N_diags = 100000000; // number of diagrams
        unsigned long long int _N_relax_steps = 100000000; // number of diagrams generated for thermalization

        
        // initialize seed for random number generator
        static inline std::mt19937::result_type setSeed(){
            std::mt19937::result_type seed = std::random_device()() ^ std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::system_clock::now().time_since_epoch()).count() 
            ^ std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
            return seed;
        };

        // fixes errors in input
        static inline int returnEven(int value){
            if(value%2==0){return value;}
            else{std::cout << "The order of a Diagram must be even, order is " << value << " + 1." << std::endl; return value + 1;}
        };

    protected:
        // random number generator
        std::mt19937 gen; // Mersenne Twister Algorithm, 32-bit

        // diagram backbone
        Vertex* _vertices;
        Propagator* _propagators;

        // simulation parameters
        long double _tau_max; // maximum value for imaginary time
        double _kx; // x momentum
        double _ky; // y momentum
        double _kz; // z momentum
        double _chem_potential; // chemical potential, normalization factor
        int _order_int_max; // maximum internal order of diagram
        int _ph_ext_max; // maximum number of external phonon lines

};
#endif