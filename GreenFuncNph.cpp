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
    for(int i=0; i<(_order_int_max + 2*_ph_ext_max + 2); i++){
        _vertices[i].tau = 0.;
        _vertices[i].type = 0;
        _vertices[i].linked = -1;
    }
    // initialize array of all possible propagators
    _propagators = new Propagator[_order_int_max + 2*_ph_ext_max + 1];
    for(int i=0;i<(order_int_max + 2*_ph_ext_max + 1);i++){
        _propagators[i].el_propagator_kx = _kx;
        _propagators[i].el_propagator_ky = _ky;
        _propagators[i].el_propagator_kz = _kz;
    }
};

void GreenFuncNph::setRelaxSteps(int relax_steps){
    while(relax_steps < 0){
        std::cout << "Invalid number of relaxation steps! Number of steps must be >= 0." << std::endl;
        std::cout << "Enter new number of relaxation steps: ";
        std::cin >> relax_steps;
    }
    _relax_steps = relax_steps;
};

void GreenFuncNph::setAlpha(double alpha){
    while(alpha<=0.){
        std::cout << "Invalid alpha value! Coupling strength must be > 0." << std::endl;
        std::cout << "Enter new alpha value: ";
        std::cin >> alpha;
    }
    _alpha = alpha;
};

void GreenFuncNph::setVolume(double volume){
    while(volume<=0.){
        std::cout << "Invalid volume value! Volume must be > 0." << std::endl;
        std::cout << "Enter new volume value: ";
        std::cin >> volume;
    }
    _volume = volume;
};

void GreenFuncNph::setDimension(int D){
    while(D!=2 && D!=3){
        std::cout << "Invalid dimension! Dimension must be 2 or 3." << std::endl;
        std::cout << "Enter new dimension: ";
        std::cin >> D;
    }
    _D = D;
};

void GreenFuncNph::setN_bins(int N_bins){
    while(N_bins <= 0){
        std::cout << "Invalid number of bins! Number of bins must be > 0." << std::endl;
        std::cout << "Enter new number of bins: ";
        std::cin >> N_bins;
    }
    _N_bins = N_bins;
    _bin_width = _tau_max/_N_bins;
    _bin_center = _bin_width/2;
};

void GreenFuncNph::setNormConst(double norm_const){
    while(norm_const <= 0){
        std::cout << "Invalid normalization constant! Constant must be > 0." << std::endl;
        std::cout << "Enter new normalization constant: ";
        std::cin >> norm_const;
    }
    _norm_const = norm_const;
};

int GreenFuncNph::findVertexPosition(double tau){
    int position = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext +1; i++){
        if(_vertices[i].tau < tau && _vertices[i+1].tau >= tau){
            position = i;
            return position;
        }
    }
    return 0;
};

int GreenFuncNph::chooseInternalPhononPropagator(){
    std::uniform_int_distribution<> distrib_unif(1,int(_current_order_int/2)); // chooses one of the internal phonon propagators at random
    int ph_propagator = distrib_unif(gen);
    int i = 0;
    int counter = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext; i++){
        if(_vertices[i].type == +1){
            counter++;
        }
        if(counter == ph_propagator){ // to be fixed
            return i;
        }
    }
    //test = test + 1;
    return 0;
};

void GreenFuncNph::phVertexMakeRoom(int index_one, int index_two){
    int i = 0;
    while(i < _current_order_int + 2*_current_ph_ext + 1){
        if(_vertices[i].linked > index_two){_vertices[i].linked += 2;}
        else if(_vertices[i].linked > index_one){_vertices[i].linked += 1;}
        i++;
    }
    
    for(int i = _current_order_int + 2*_current_ph_ext + 1; i > index_one; i--){ 
        if(i > index_two){_vertices[i+2] = _vertices[i];}
        else{_vertices[i+1] = _vertices[i];}
  }
};

void GreenFuncNph::phVertexRemoveRoom(int index_one, int index_two){
    int i = 0;
    while(i < _current_order_int + 2*_current_ph_ext + 1){
        if(_vertices[i].linked >= index_two){_vertices[i].linked -= 2;}
        else if(_vertices[i].linked >= index_one){_vertices[i].linked -= 1;}
        i++;
    }

    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; i++){
        if(i < index_two-1){_vertices[i] = _vertices[i+1];}
        else{_vertices[i] = _vertices[i+2];}
    }
};