#ifndef DIAGRAM_HPP
#define DIAGRAM_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
#include "../thirdparty/pcg-cpp-0.98/include/pcg_random.hpp"
#include "../thirdparty/pcg-cpp-0.98/include/pcg_extras.hpp"
#include "utils/computational_methods.hpp"
#include "utils/MC_data_structures.hpp"

class Diagram {
    public:
        // constructor
        Diagram() = default;
        Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, 
            double chem_potential, int order_int_max, int ph_ext_max, int data_type);

        Diagram(FullVertexNode * nodes, FullVertexNodeIndicator * internal_used, FullVertexNodeIndicator * external_used, 
            int current_order, unsigned long long int N_diags, long double tau_max, 
            double kx, double ky, double kz, double chem_potential, int order_int_max, int ph_ext_max, int data_type);

        // destructor
        ~Diagram(){
            if(_data_type == _data_type_array[0]){
                delete[] _vertices;
                delete[] _propagators;
            }
            else if(_data_type == _data_type_array[1]){
                _helper = _head;
                _pointer_one = nullptr;
                while(_helper != nullptr){
                    _pointer_one = _helper->next;
                    _helper->next = _free_list;
                    _helper->prev = nullptr;
                    _free_list = _helper;
                    _helper = _pointer_one;
                }
            
                // reset all pointer
                _tail = nullptr;
                _head = nullptr;
                _diagram_list = nullptr;
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                _helper = nullptr;

                // free array of nodes and move _free_list to nullptr
                delete[] _nodes;
                _free_list = nullptr;
            
                delete[] _internal_used;
                delete[] _external_used;
            }
        };

        // returns uniform random long double precision value between 0 and 1
        static inline long double drawUniformR(){std::uniform_real_distribution<long double> distrib(0,1);long double r = distrib(gen); return r;};

        // Metropolis-Hastings Monte Carlo method
        static inline int Metropolis(double R){
            double acceptance_ratio = std::min(1.0, R);
            if(drawUniformR() > acceptance_ratio){return 0;} // reject
            else{return 1;} // accept
        };

        // initialize seed for random number generator
        static inline void setSeed(uint64_t seed=0x853c49e6748fea9bULL, uint64_t stream=0xda3e39cb94b95bdbULL){
            gen.seed(seed, stream);
        }

        /*static inline void setSeed(){
            std::mt19937::result_type seed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() 
            ^ std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
            gen.seed(seed);
        };*/

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
        //static inline std::mt19937 getSeed(){return gen;};
        static inline pcg32 getSeed(){return gen;};
        void getNodes(FullVertexNode * nodes, FullVertexNodeIndicator * internal_used, FullVertexNodeIndicator * external_used, int size);
        
        // setters
        void setRelaxSteps(unsigned long long int relax_steps);
        void setAutcorrSteps(unsigned long long int autocorr_steps);

        // fix double precision floating errors
        inline void fixDoublePrecisionErrors(double error_threshold){
            _helper = _head;
            while(_helper->next != nullptr){
                if(std::abs(_helper->k[0] - _kx) < error_threshold){
                    _helper->k[0] = _kx;
                }
                if(std::abs(_helper->k[1] - _ky) < error_threshold){
                    _helper->k[1] = _ky;
                }
                if(std::abs(_helper->k[2] - _kz) < error_threshold){
                    _helper->k[2] = _kz;
                }
                _helper = _helper->next;
            }
            _helper = nullptr;
        };

        inline void checkTimeErrors(){
            _helper = _head;
            while(_helper->next != nullptr){
                if(_helper->tau > _helper->next->tau){
                    std::cerr << "time not valid" << std::endl;
                }
                _helper = _helper->next;
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
        //static thread_local std::mt19937 gen; // Mersenne Twister Algorithm, 32-bit
        static thread_local pcg32 gen; // Permuted Congruential Generator, 32-bit

        Vertex* _vertices = nullptr;
        Propagator* _propagators = nullptr;

        // pool of nodes (allocated as array)
        FullVertexNode * _nodes = nullptr;

        // pointer to array of free nodes
        FullVertexNode* _free_list = nullptr;

        // diagram backbone
        FullVertexNode* _diagram_list = nullptr;

        // support arrays for phonon propagator tracing
        FullVertexNodeIndicator* _internal_used = nullptr;
        FullVertexNodeIndicator* _external_used = nullptr;

        // support variables
        FullVertexNode * _head = nullptr;
        FullVertexNode * _tail = nullptr;
        FullVertexNode * _pointer_one = nullptr;
        FullVertexNode * _pointer_two = nullptr;
        FullVertexNode * _helper = nullptr;

        // simulation parameters
        long double _tau_max = 50; // maximum value for imaginary time
        double _kx = 0; // x momentum
        double _ky = 0; // y momentum
        double _kz = 0; // z momentum
        double _chem_potential = -1; // chemical potential, normalization factor
        int _order_int_max = 0; // maximum internal order of diagram
        int _ph_ext_max = 0; // maximum number of external phonon lines

        // data structure type
        const int _data_type = 0;
        int _data_type_array[2] = {0, 1};

        // diagram list methods
        void insertNode(FullVertexNode * node_pointer);
        void deleteNode(FullVertexNode *& node_pointer);
        inline FullVertexNode * findInternalPhononVertex(int index){return _internal_used[index].linked;};
        inline FullVertexNode * findExternalPhononVertex(int index){return _external_used[index].linked;};
};
#endif