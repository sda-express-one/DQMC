#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFuncNph.hpp"

GreenFuncNph::GreenFuncNph(long long unsigned int N_diags, double tau_max, double kx, double ky, double kz,
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
        if(_vertices[i].type == +1 && _vertices[i].linked != -1){
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

void GreenFuncNph::propagatorArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order_int + 2*_current_ph_ext; i > index_one-1; i--){
        if(i > index_two - 1){_propagators[i+2] = _propagators[i];}
        else{_propagators[i+1] = _propagators[i];}
    }
    _propagators[index_one+1] = _propagators[index_one];
    _propagators[index_two+1] = _propagators[index_two+2];
};

void GreenFuncNph::propagatorArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; i++){
        if(i < index_two - 2){_propagators[i] = _propagators[i+1];}
        else{_propagators[i] = _propagators[i+2];}
    }
};

double GreenFuncNph::diagramLengthUpdate(double tau_init){
    double tau_fin = _last_vertex - std::log(1-drawUniformR())/(calcEnergy(_kx,_ky,_kz)-_chem_potential);
    if(tau_fin <= _tau_max){
        _vertices[_current_order_int + 2*_current_ph_ext + 1].tau = tau_fin;
        //_vertices[_current_order_int + 2*_current_ph_ext + 1].type = 0;
        return tau_fin;
    }
    else{return tau_init;}
};

void GreenFuncNph::addInternalPhononPropagator(){
    if(_current_order_int+1 >= _order_int_max){return;} // reject if already at max order
    else{
        std::uniform_int_distribution<> distrib_prop(0, _current_order_int + 2*_current_ph_ext);
        int propagator = distrib_prop(gen);
        double tau_init = _vertices[propagator].tau;
        double tau_end = _vertices[propagator+1].tau;
        std::uniform_real_distribution<> distrib_unif(tau_init, tau_end);
        double tau_one = distrib_unif(gen);
        double tau_two = tau_one - std::log(1-drawUniformR())/1;
        if(tau_two >= _vertices[_current_order_int + 2*_current_ph_ext +1].tau){return;} // reject if phonon vertex goes out of bound
        else{
            // sampling momentum values for phonon propagators
            std::normal_distribution<> distrib_norm(0, std::sqrt(1/(tau_two-tau_one))); // may need variance, inserted std dev
            double w_x = distrib_norm(gen);
            double w_y = distrib_norm(gen);
            double w_z = distrib_norm(gen);

            // find position of new tau values
            int index_one = findVertexPosition(tau_one);
            int index_two = findVertexPosition(tau_two);
            
            // create arrays of momentum values
            double* px_init = new double[index_two + 1 - index_one];
            double* py_init = new double[index_two + 1 - index_one];
            double* pz_init = new double[index_two + 1 - index_one];
            double* px_fin = new double[index_two + 1 - index_one];
            double* py_fin = new double[index_two + 1 - index_one];
            double* pz_fin = new double[index_two + 1 - index_one];

            for(int i = index_one; i < index_two +1; i++){
                px_init[i - index_one] = _propagators[i].el_propagator_kx;
                px_fin[i - index_one] = _propagators[i].el_propagator_kx - w_x;
                py_init[i - index_one] = _propagators[i].el_propagator_ky;
                py_fin[i - index_one] = _propagators[i].el_propagator_ky - w_y;
                pz_init[i - index_one] = _propagators[i].el_propagator_kz;
                pz_fin[i - index_one] = _propagators[i].el_propagator_kz - w_z;
            }

            // create arrays of energy values
            double* energy_init = new double[index_two + 1 - index_one];
            double* energy_fin = new double[index_two + 1 - index_one];

            for(int i = 0; i < index_two - index_one + 1; i++){
                energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
                energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
            }

            delete[] px_init;
            delete[] py_init;
            delete[] pz_init;
            delete[] px_fin;
            delete[] py_fin;
            delete[] pz_fin;

            double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext + 1); // may be wrong
            double p_A = _p_add_int*(_current_order_int/2 + 1);

            // initial and final action
            double exponent_init = 0.;
            double exponent_fin = 0.;

            if(index_one == index_two){
                exponent_init = energy_init[0]*(tau_two-tau_one);
                exponent_fin = energy_fin[0]*(tau_two-tau_one);
            }
            else{
                exponent_init = energy_init[0]*(_vertices[index_one+1].tau-tau_one);
                exponent_fin = energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
                for(int i = 1; i < index_two - index_one; i++){
                    exponent_init += energy_init[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                    exponent_fin += energy_fin[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                }
                exponent_init += energy_init[index_two - index_one]*(tau_two-_vertices[index_two].tau);
                exponent_fin += energy_fin[index_two - index_one]*(tau_two-_vertices[index_two].tau);
            }

            // delete array of energies
            delete[] energy_init;
            delete[] energy_fin;

            double numerator = p_B*std::exp(-(exponent_fin - exponent_init + (1.)*(tau_two - tau_one)))*calcVertexStrength(w_x,w_y,w_z)*(tau_end-tau_init);
            double denominator = p_A*std::pow(2*M_PI,_D)*std::exp(-(tau_two-tau_one))*std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_two-tau_one));

            double R_add = numerator/denominator;
            double acceptance_ratio = std::min(1.,R_add);

            if(drawUniformR()>acceptance_ratio){return;}
            else{
                phVertexMakeRoom(index_one, index_two); // make room in vertices array
                propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array

                // assign vertex one values
                _vertices[index_one+1].tau = tau_one;
                _vertices[index_one+1].type = +1;
                _vertices[index_one+1].linked = index_two + 2;
                _vertices[index_one+1].wx = w_x;
                _vertices[index_one+1].wy = w_y;
                _vertices[index_one+1].wz = w_z;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_two;
                _vertices[index_two+2].type = -1;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;

                // update electron propagator energies
                for(int i=index_one+1; i<index_two+2;i++){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                }

                _current_order_int += 2; // update current diagram internal order
                findLastPhVertex();
            }
        } 
    }
};

void GreenFuncNph::removeInternalPhononPropagator(){
    if(_current_order_int < 2){return;} // reject if already at order 0
    else{
        // indexes of initial and final vertices of a random internal phonon propagator
        int index_one = chooseInternalPhononPropagator();
        int index_two = _vertices[index_one].linked;

        double tau_one = _vertices[index_one].tau;
        double tau_two = _vertices[index_two].tau;

        double tau_init = _vertices[index_one-1].tau;
        double tau_end = _tau_max;
        if(index_two != index_one+1){tau_end = _vertices[index_one+1].tau;}
        else{tau_end = _vertices[index_one+2].tau;}

        double w_x = _vertices[index_one].wx;
        double w_y = _vertices[index_one].wy;
        double w_z = _vertices[index_one].wz;

        // create arrays of momentum values
        double* px_init = new double[index_two - index_one];
        double* py_init = new double[index_two - index_one];
        double* pz_init = new double[index_two - index_one];
        double* px_fin = new double[index_two - index_one];
        double* py_fin = new double[index_two - index_one];
        double* pz_fin = new double[index_two - index_one];

        for(int i = index_one; i < index_two; i++){
            px_init[i - index_one] = _propagators[i].el_propagator_kx + w_x;
            px_fin[i - index_one] = _propagators[i].el_propagator_kx;
            py_init[i - index_one] = _propagators[i].el_propagator_ky + w_y;
            py_fin[i - index_one] = _propagators[i].el_propagator_ky;
            pz_init[i - index_one] = _propagators[i].el_propagator_kz + w_z;
            pz_fin[i - index_one] = _propagators[i].el_propagator_kz;
        }

        // create arrays of energy values
        double* energy_init = new double[index_two - index_one];
        double* energy_fin = new double[index_two - index_one];

        for(int i = 0; i < index_two - index_one; i++){
            energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
            energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
        }

        delete[] px_init;
        delete[] px_fin;
        delete[] py_init;
        delete[] py_fin;
        delete[] pz_init;
        delete[] pz_fin;
        
        // initial and final action
        double exponent_fin = 0.;
        double exponent_init = 0.;
        
        if(index_one + 1 == index_two){
            exponent_fin = energy_fin[0]*(tau_two-tau_one);
            exponent_init = energy_init[0]*(tau_two-tau_one);
        }
        else{
            exponent_fin = energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
            exponent_init = energy_init[0]*(_vertices[index_one+1].tau-tau_one);
            for(int i = 1; i < index_two - index_one - 1; i++){
                exponent_fin += energy_fin[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);      
                exponent_init += energy_init[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
            }
            exponent_fin += energy_fin[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
            exponent_init += energy_init[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
        }

        delete[] energy_init;
        delete[] energy_fin;

        double p_A = _p_add_int*((_current_order_int - 2)/2 + 1);
        double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext - 1);
        
        double numerator = p_A*std::pow(2*M_PI,_D)*std::exp(-(tau_two-tau_one))*std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*
        std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_two-tau_one));
        double denominator = p_B*std::exp(-(exponent_fin-exponent_init+(1.)*(tau_two - tau_one)))*calcVertexStrength(w_x,w_y,w_z)*(tau_end-tau_init);

        double acceptance_ratio = std::min(1., numerator/denominator);
        if(drawUniformR() > acceptance_ratio){return;}
        else{
            for(int i=index_one; i<index_two;i++){
                _propagators[i].el_propagator_kx += w_x;
                _propagators[i].el_propagator_ky += w_y;
                _propagators[i].el_propagator_kz += w_z;
            }
            phVertexRemoveRoom(index_one, index_two);
            propagatorArrayRemoveRoom(index_one, index_two);
            _current_order_int -= 2;
            findLastPhVertex();
        }
    }
};

void GreenFuncNph::addExternalPhononPropagator(){
    if(_current_ph_ext >= _ph_ext_max){return;} // return if already at max number of ext phonon propagators
    else{

        int total_order = _current_order_int +2*_current_ph_ext;
        double current_tau = _vertices[total_order + 2].tau; // length of current diagram

        double tau_one = 0. - std::log(1-drawUniformR())/1; // time of left vertex of ext phonon propagator
        if(tau_one >= current_tau){return;} // reject if it goes out of bound

        double tau_two = current_tau + std::log(1-drawUniformR())/1; // right vertex
        if(tau_two <= 0.){return;} // reject if it goes out of bound

        std::normal_distribution<> distrib_norm(0, std::sqrt(1/(current_tau-tau_two+tau_one))); // may need variance, inserted std dev
        double w_x = distrib_norm(gen);
        double w_y = distrib_norm(gen);
        double w_z = distrib_norm(gen);

        if(tau_one <= tau_two){
            int index_one = findVertexPosition(tau_one);
            int index_two = findVertexPosition(tau_two);

            double* px_one_init = new double[index_one + 1];
            double* px_two_init = new double[total_order + 1 - index_two];

            double* py_one_init = new double[index_one + 1];
            double* py_two_init = new double[total_order + 1 - index_two];

            double* pz_one_init = new double[index_one + 1];
            double* pz_two_init = new double[total_order + 1 - index_two];

            double* px_one_fin = new double[index_one + 1];
            double* px_two_fin = new double[total_order + 1 - index_two];

            double* py_one_fin = new double[index_one + 1];
            double* py_two_fin = new double[total_order +1 - index_two];

            double* pz_one_fin = new double[index_one + 1];
            double* pz_two_fin = new double[total_order + 1 - index_two];

            // retrieve momentum values for propagators below first ph vertex
            for(int i = 0; i < index_one + 1; i++){
                px_one_init[i] = _propagators[i].el_propagator_kx;
                px_one_fin[i] = _propagators[i].el_propagator_kx - w_x;
                py_one_init[i] = _propagators[i].el_propagator_ky;
                py_one_fin[i] = _propagators[i].el_propagator_ky - w_y;
                pz_one_init[i] = _propagators[i].el_propagator_kz;
                pz_one_fin[i] = _propagators[i].el_propagator_kz - w_z;
            }
            
            // retrieve momentum values for propagators above second ph vertex
            for(int i = index_two; i < total_order + 1; i++){
                px_one_init[i-index_two] = _propagators[i].el_propagator_kx;
                px_one_fin[i-index_two] = _propagators[i].el_propagator_kx - w_x;
                py_one_init[i-index_two] = _propagators[i].el_propagator_ky;
                py_one_fin[i-index_two] = _propagators[i].el_propagator_ky - w_y;
                pz_one_init[i-index_two] = _propagators[i].el_propagator_kz;
                pz_one_fin[i-index_two] = _propagators[i].el_propagator_kz - w_z;
            }

            double* energy_one_init = new double[index_one+1];
            double* energy_two_init = new double[total_order - index_two];

            double* energy_one_fin = new double[index_one+1];
            double* energy_two_fin = new double[total_order + 1 - index_two];

            // calc energy values for propagators below first ph vertex
            for(int i = 0; i < index_one + 1; i++){
                energy_one_init[i] = calcEnergy(px_one_init[i], py_one_init[i], pz_one_init[i]);
                energy_one_fin[i] = calcEnergy(px_one_fin[i], py_one_fin[i], pz_one_fin[i]);
            }

            // calc energy values for propagators above second ph vertex
            for(int i = 0; i < total_order + 1 - index_two; i++){
                energy_two_init[i] = calcEnergy(px_two_init[i], py_two_init[i], pz_two_init[i]);
                energy_two_fin[i] = calcEnergy(px_two_fin[i], py_two_fin[i], pz_two_fin[i]);
            }

            delete[] px_one_init, px_one_fin, py_one_init, py_one_fin, pz_one_init, pz_one_fin;
            delete[] px_two_init, px_two_fin, py_two_init, py_two_fin, pz_two_init, pz_two_fin;

            // initial and final action
            double exponent_one_init = 0.; 
            double exponent_two_init = 0.;
            double exponent_one_fin = 0.;
            double exponent_two_fin = 0.;

            for(int i = 0; i < index_one; i++){
               exponent_one_init += energy_one_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
               exponent_one_fin += energy_one_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }
            exponent_one_init += energy_one_init[index_one]*(tau_one - _vertices[index_one].tau);
            exponent_one_fin += energy_one_fin[index_one]*(tau_one - _vertices[index_one].tau);

            exponent_two_init = energy_two_init[0]*(_vertices[index_two+1].tau - tau_two);
            exponent_two_fin = energy_two_fin[0]*(_vertices[index_two+1].tau - tau_two);
            for(int i = index_two+1; i < total_order + 1; i++){
                exponent_two_init += energy_two_init[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
                exponent_two_fin += energy_two_fin[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            // delete array of energies
            delete[] energy_one_init, energy_one_fin, energy_two_init, energy_two_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double numerator = p_B*std::exp(-(exponent_two_fin + exponent_one_fin - exponent_two_init - exponent_one_init + 
                1.*(current_tau-tau_two+tau_one)))*calcVertexStrength(w_x, w_y, w_z);

            double denominator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(current_tau-tau_two))*
            std::pow(((current_tau-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(current_tau-tau_two+tau_one));

            double R_add = numerator/denominator;
            double acceptance_ratio = std::min(1.,R_add);
            if(drawUniformR() > acceptance_ratio){return;}
            else{
                phVertexMakeRoom(index_one, index_two); // make room in vertices array
                propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array

                // assign vertex one values
                _vertices[index_one+1].tau = tau_one;
                _vertices[index_one+1].type = -2;
                _vertices[index_one+1].linked = index_two + 2;
                _vertices[index_one+1].wx = w_x;
                _vertices[index_one+1].wy = w_y;
                _vertices[index_one+1].wz = w_z;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_two;
                _vertices[index_two+2].type = +2;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;

                // update electron propagator energies
                for(int i = 0; i < index_one + 2; i++){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                }
                for(int i = index_two+2; i < total_order + 3; i++){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                }

                _current_ph_ext += 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
        else{
            int index_one = findVertexPosition(tau_two);
            int index_two = findVertexPosition(tau_one);

            int total_order = _current_order_int + _current_ph_ext;

            double* px_init = new double[total_order+2];
            double* py_init = new double[total_order+2];
            double* pz_init = new double[total_order+2];
            
            double* px_fin = new double[total_order+4];
            double* py_fin = new double[total_order+4];
            double* pz_fin = new double[total_order+4];

            for(int i = 0; i < total_order+2; i++){
                px_init[i] = _propagators[i].el_propagator_kx;
                py_init[i] = _propagators[i].el_propagator_ky;
                pz_init[i] = _propagators[i].el_propagator_kz;
            }

            for(int i = 0; i < index_one + 1; i++){
                px_fin[i] = _propagators[i].el_propagator_kx - w_x;
                py_fin[i] = _propagators[i].el_propagator_ky - w_y;
                pz_fin[i] = _propagators[i].el_propagator_kz - w_z;
            }
            for(int i = index_one + 1; i < index_two + 2; i++){
                px_fin[i] = _propagators[i-1].el_propagator_kx - 2*w_x;
                py_fin[i] = _propagators[i-1].el_propagator_ky - 2*w_y;
                pz_fin[i] = _propagators[i-1].el_propagator_kz - 2*w_z;
            }
            for(int i = index_two + 2; i < total_order + 4; i++){
                px_fin[i] = _propagators[i-2].el_propagator_kx - w_x;
                py_fin[i] = _propagators[i-2].el_propagator_ky - w_y;
                pz_fin[i] = _propagators[i-2].el_propagator_kz - w_z;
            }

            double* energy_init = new double[total_order+2];
            double* energy_fin = new double[total_order+4];

            for(int i = 0; i < total_order+2; i++){
                energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
            }
            for(int i = 0; i < total_order+4; i++){
                energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
            }

            delete[] px_init, px_fin, py_init;

            double exponent_init = 0.;
            double exponent_fin = 0.;

            for(int i = 0; i < total_order+2; i++){
                exponent_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            for(int i = 0; i < index_one; i++){
                exponent_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }
            exponent_fin += energy_fin[index_one]*(tau_two - _vertices[index_one].tau);
            exponent_fin += energy_fin[index_one+1]*(_vertices[index_one+1].tau - tau_two);
            for(int i = index_one+2; i < index_two + 2; i++){
                exponent_fin += energy_fin[i]*(_vertices[i].tau - _vertices[i-1].tau);
            }
            exponent_fin += energy_fin[index_two+2]*(tau_one - _vertices[index_two].tau);
            exponent_fin += energy_fin[index_two+3]*(_vertices[index_two+1].tau - tau_one);
            for(int i = index_two+3; i < total_order + 4; i++){
                exponent_fin += energy_fin[i]*(_vertices[i-1].tau - _vertices[i-2].tau);
            }
            
            // delete array of energies
            delete[] energy_init, energy_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double numerator = p_B*std::exp(-(exponent_fin - exponent_init + 1.*(current_tau-tau_two+tau_one)))
            *calcVertexStrength(w_x, w_y, w_z);

            double denominator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(current_tau-tau_two))*
            std::pow(((current_tau-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(current_tau-tau_two+tau_one));

            double R_add = numerator/denominator;
            double acceptance_ratio = std::min(1.,R_add);

            if(drawUniformR() < acceptance_ratio){delete px_fin, py_fin, pz_fin; return;}
            else{
                phVertexMakeRoom(index_one, index_two); // make room in vertices array
                propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array

                // assign vertex one values
                _vertices[index_one+1].tau = tau_two;
                _vertices[index_one+1].type = +2;
                _vertices[index_one+1].linked = index_two + 2;
                _vertices[index_one+1].wx = w_x;
                _vertices[index_one+1].wy = w_y;
                _vertices[index_one+1].wz = w_z;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_one;
                _vertices[index_two+2].type = -2;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;

                // update electron propagator energies
                for(int i = 0; i < total_order + 4; i++){
                    _propagators[i].el_propagator_kx = px_fin[i];
                    _propagators[i].el_propagator_ky = py_fin[i];
                    _propagators[i].el_propagator_kz = pz_fin[i];
                }

                delete px_fin, py_fin, pz_fin;

                _current_ph_ext += 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
    }
};

void GreenFuncNph::Diagnostic(std::string filename, int i) const {
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Iteration: " << i << "\n";

    for(int j=0; j<_current_order_int+2*_current_ph_ext+1; j++){
        file << "index: " << j << " time: " << _vertices[j].tau << " wx: " << _vertices[j].wx << " wy: " 
        << _vertices[j].wy << " wz: " << _vertices[j].wz <<" type: " << _vertices[j].type << " linked: " << _vertices[j].linked << "\n";
        file << "propagator: " << j << "          kx: " << _propagators[j].el_propagator_kx << " ky: " << _propagators[j].el_propagator_ky <<
        " kz: " << _propagators[j].el_propagator_kz << "\n";
    }

    file << "index: " << _current_order_int+2*_current_ph_ext+1 << " time: " << _vertices[_current_order_int+2*_current_ph_ext+1].tau << " wx: "
    << _vertices[_current_order_int+2*_current_ph_ext+1].wx << " wy: " << _vertices[_current_order_int+2*_current_ph_ext+1].wy << " wz: " 
    << _vertices[_current_order_int+2*_current_ph_ext+1].wz <<" type: " << _vertices[_current_order_int+2*_current_ph_ext+1].type 
    << " linked: " << _vertices[_current_order_int+2*_current_ph_ext+1].linked << "\n";

    file << std::endl;
    file.close();
};