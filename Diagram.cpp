#include "Diagram.hpp"

// constructor definition
Diagram::Diagram(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max) : _N_diags(N_diags), gen(setSeed()), _tau_max(tau_max),
    _chem_potential(chem_potential), _order_int_max(returnEven(order_int_max)), _ph_ext_max(ph_ext_max) {

        // assign momentum values
        _kx = kx;
        _ky = ky;
        _kz = kz;

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

    // setters
    
    void Diagram::setRelaxSteps(int relax_steps){
        if(relax_steps < 0){
            std::cerr << "Invalid number of relaxation steps! Number of steps must be >= 0." << std::endl;
            std::cerr << "Reverting to default value of " << _N_relax_steps << "." << std::endl;
            return;
        }
        _N_relax_steps = relax_steps;
    }