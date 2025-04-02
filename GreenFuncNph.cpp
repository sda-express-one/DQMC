#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFuncNph.hpp"

GreenFuncNph::GreenFuncNph(long long int N_diags, double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max) : _N_diags(N_diags), _tau_max(tau_max),
    _order_int_max(_order_int_max), _ph_ext_max(ph_ext_max) {
    // assign momentum values
    _kx = kx;
    _ky = ky;
    _kz = kz;

    // initialize array of all possible phonon vertices
    _vertices = new Vertex[_order_int_max + 2*_ph_ext_max +2];
    for(int i=0; i<(_order_int_max + _ph_ext_max + 2);i++){
        _vertices[i].tau = 0.;
        _vertices[i].type = 0;
        _vertices[i].linked = -1;
    }
};