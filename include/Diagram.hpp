#ifndef DIAGRAM_HPP
#define DIAGRAM_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
//#include "../thirdparty/pcg-cpp-0.98/include/pcg_random.hpp"
#include "utils/computational_methods.hpp"
#include "utils/MC_data_structures.hpp"

class Diagram {
    public:
        // constructor
        Diagram() = default;
        Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, 
            double chem_potential, int order_int_max, int ph_ext_max);

        Diagram(Propagator * propagators, Vertex * vertices,
            unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, 
            double chem_potential, int order_int_max, int ph_ext_max);

        // destructor
        ~Diagram(){
            delete[] _vertices;
            delete[] _propagators;
        }

        // returns uniform random long double precision value between 0 and 1
        static inline long double drawUniformR(){std::uniform_real_distribution<long double> distrib(0,1);long double r = distrib(gen); return r;};

        // Metropolis-Hastings Monte Carlo method
        static inline int Metropolis(double R){
            double acceptance_ratio = std::min(1.0, R);
            if(drawUniformR() > acceptance_ratio){return 0;} // reject
            else{return 1;} // accept
        };

        // initialize seed for random number generator
        static inline void setSeed(){
            std::mt19937::result_type seed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() 
            ^ std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
            gen.seed(seed);
            //retur seed;
        };

        /*static inline void setSeed(unsigned long int num_thread = 0){
            uint64_t seed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() ^
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
            gen = pcg32(seed, num_thread);        
        };*/

        //pcg32(uint64_t seed = 0x853c49e6748fea9bULL, uint64_t stream = 0xda3e39cb94b95bdbULL){};

        // setters
        void setRelaxSteps(unsigned long long int relax_steps);
        void setAutcorrSteps(unsigned long long int autocorr_steps);

        // getters
        inline long long unsigned int getNdiags() const {return _N_diags;};
        inline long long unsigned int getRelaxSteps() const {return _N_relax_steps;};
        inline long long unsigned int getAutocorrSteps() const {return _N_autocorr_steps;};
        inline long double getTauMax() const {return _tau_max;};
        inline double getkx() const {return _kx;};
        inline double getky() const {return _ky;};
        inline double getkz() const {return _kz;};
        inline double getChemPotential() const {return _chem_potential;};
        inline int getOrderIntMax() const {return _order_int_max;};
        inline int getPhExtMax() const {return _ph_ext_max;};
        inline int getOrderMax() const {return _order_int_max + 2*_ph_ext_max;};
        static inline std::mt19937 getSeed(){return gen;};


        // fix double precision floating errors
        inline void fixDoublePrecisionErrors(int total_order, double error_threshold){
            for(int i = 0; i < total_order+1; i++){
                if(std::abs(_propagators[i].el_propagator_kx - _kx) < error_threshold){
                    _propagators[i].el_propagator_kx = _kx;
                }
                if(std::abs(_propagators[i].el_propagator_ky - _ky) < error_threshold){
                    _propagators[i].el_propagator_ky = _ky;
                }
                if(std::abs(_propagators[i].el_propagator_kz - _kz) < error_threshold){
                    _propagators[i].el_propagator_kz = _kz;
                }
            }
        };

        inline void checkTimeErrors(int total_order){
            for(int i = 1; i < total_order+1; i++)
                if(_vertices[i].tau < _vertices[i-1].tau){
                    std::cerr << "time not valid" << std::endl;
                }
        };

    private:
        // number of diagrams generated
        const unsigned long long int _N_diags = 100000000; // number of diagrams
        unsigned long long int _N_relax_steps = 100000000; // number of diagrams generated for thermalization
        unsigned long long int _N_autocorr_steps = 0; // number of steps for autocorrelation negation

        // fixes errors in input
        static inline int returnEven(int value){
            if(value%2==0){return value;}
            else{std::cout << "The order of a Diagram must be even, order is " << value << " + 1." << std::endl; return value + 1;}
        };

    protected:
        // random number generator
        static thread_local std::mt19937 gen; // Mersenne Twister Algorithm, 32-bit
        //static thread_local pcg32 gen;

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