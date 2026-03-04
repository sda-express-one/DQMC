#include "../include/Diagram.hpp"

//thread_local std::mt19937 Diagram::gen;

thread_local pcg32 Diagram::gen;

// constructor definition
Diagram::Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) : _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max) {
    
    // assign data structure type
    _data_type = _data_type_array[data_type];

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;
    
    if(_data_type == _data_type_array[0]){
        // initialize array of all possible phonon vertices
        _vertices = new Vertex[_order_int_max + 2*_ph_ext_max + 2];
        for(int i=0; i<_order_int_max + 2*_ph_ext_max + 2; i++){
            _vertices[i].tau = 0;
            _vertices[i].type = 0;
            _vertices[i].linked = -1;
        }

        // initialize array of all possible propagators
        _propagators = new Propagator[_order_int_max + 2*_ph_ext_max + 1];
        for(int i=0; i<(order_int_max + 2*_ph_ext_max + 1); i++){
            _propagators[i].el_propagator_kx = _kx;
            _propagators[i].el_propagator_ky = _ky;
            _propagators[i].el_propagator_kz = _kz;
        }
    }

    else if(_data_type == _data_type_array[1]){
        // initialize array of all possible full vertices
        _full_vertices = new FullVertex * [_order_int_max + 2*_ph_ext_max + 2];
        _indices = new int[_order_int_max + 2*_ph_ext_max + 2];
        for(int i=0; i<(_order_int_max + 2*_ph_ext_max + 2); ++i){
            _full_vertices[i] = new FullVertex;
            _full_vertices[i]->k[0] = _kx;
            _full_vertices[i]->k[1] = _ky;
            _full_vertices[i]->k[2] = _kz;
            _full_vertices[i]->index = 0;
            _indices[i] = -1;
        }
    }
}

Diagram::Diagram(Propagator * propagators, Vertex * vertices, 
    unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) : _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max) {

    _data_type = _data_type_array[data_type];

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    // initialize array of all possible phonon vertices
    _vertices = new Vertex[_order_int_max + 2*_ph_ext_max + 2];
    for(int i=0; i<_order_int_max + 2*_ph_ext_max + 2; i++){
        _vertices[i] = vertices[i];
    }

    // initialize array of all possible propagators
    _propagators = new Propagator[_order_int_max + 2*_ph_ext_max + 1];
    for(int i=0; i<(order_int_max + 2*_ph_ext_max + 1); i++){
        _propagators[i] = propagators[i];
    }
}

Diagram::Diagram(FullVertex * fullvertices, int * indices, unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz, 
    double chem_potential, int order_int_max, int ph_ext_max, int data_type) :  _N_diags(N_diags), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max) {
    
    _data_type = _data_type_array[data_type];

    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    // initialize array of all possible (full) vertices and indices
    _full_vertices = new FullVertex * [order_int_max + 2*ph_ext_max + 2];
    _indices = new int[_order_int_max + 2*_ph_ext_max + 2];

    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; ++i){
        _full_vertices[i] = new FullVertex;
        *_full_vertices[i] = fullvertices[i];
        _indices[i] = indices[i];
    }
};

// setters
    
void Diagram::setRelaxSteps(unsigned long long int relax_steps){_N_relax_steps = relax_steps;}

void Diagram::setAutcorrSteps(unsigned long long int autocorr_steps){_N_autocorr_steps = autocorr_steps;}