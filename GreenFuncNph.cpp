#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFuncNph.hpp"

GreenFuncNph::GreenFuncNph(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_int_max, int ph_ext_max) : gen(setSeed()), _N_diags(N_diags), _tau_max(tau_max), 
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
};

void GreenFuncNph::setRelaxSteps(int relax_steps){
    while(relax_steps < 0){
        std::cout << "Invalid number of relaxation steps! Number of steps must be >= 0." << std::endl;
        std::cout << "Enter new number of relaxation steps: ";
        std::cin >> relax_steps;
        std::cout << "\n";
    }
    long long unsigned int temp = relax_steps;
    _relax_steps = temp;
};

void GreenFuncNph::setAlpha(double alpha){
    while(alpha<0){
        std::cout << "Invalid coupling strength value! Coupling strength must be > 0." << std::endl;
        std::cout << "Enter new alpha value: ";
        std::cin >> alpha;
        std::cout << "\n";
    }
    _alpha = alpha;
};

void GreenFuncNph::setVolume(double volume){
    while(volume<=0.){
        std::cout << "Invalid volume value! Volume must be > 0." << std::endl;
        std::cout << "Enter new volume value: ";
        std::cin >> volume;
        std::cout << "\n";
    }
    _volume = volume;
};

void GreenFuncNph::setDimension(int D){
    while(D!=2 && D!=3){
        std::cout << "Invalid dimension! Dimension must be 2 or 3." << std::endl;
        std::cout << "Enter new dimension: ";
        std::cin >> D;
        std::cout << "\n";
    }
    _D = D;
};

void GreenFuncNph::setN_bins(int N_bins){
    while(N_bins <= 0){
        std::cout << "Invalid number of bins! Number of bins must be > 0." << std::endl;
        std::cout << "Enter new number of bins: ";
        std::cin >> N_bins;
        std::cout << "\n";
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
        std::cout << "\n";
    }
    _norm_const = norm_const;
};

void GreenFuncNph::setTauCutoffEnergy(long double tau_cutoff_energy){
    while(tau_cutoff_energy <= 0 || tau_cutoff_energy >= _tau_max){
        std::cout << "Invalid energy cutoff! Cutoff must be > 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "Enter new energy cutoff: ";
        std::cin >> tau_cutoff_energy;
        std::cout << "\n";
    }
    _tau_cutoff_energy = tau_cutoff_energy;
}

void GreenFuncNph::setTauCutoffMass(long double tau_cutoff_mass){
    while(tau_cutoff_mass <= 0 || tau_cutoff_mass >= _tau_max){
        std::cout << "Invalid mass cutoff! Cutoff must be > 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "Enter new mass cutoff: ";
        std::cin >> tau_cutoff_mass;
        std::cout << "\n";
    }
    _tau_cutoff_mass = tau_cutoff_mass;
}

void GreenFuncNph::setNumPoints(int num_points){
    while(num_points <= 0){
        std::cout << "Invalid number of points! Number of points must be > 0." << std::endl;
        std::cout << "Enter new number of points: ";
        std::cin >> num_points;
        std::cout << "\n";
    }
    _num_points = num_points;
    _points_step = _tau_max/(double)_num_points;
    _points_center = _points_step/2;
};

void GreenFuncNph::setSelectedOrder(int selected_order){
    while(selected_order < 0 || selected_order > _ph_ext_max){
        std::cout << "Invalid order! Order must be >= 0 and <= " << _ph_ext_max << "." << std::endl;
        std::cout << "Enter new order: ";
        std::cin >> selected_order;
        std::cout << "\n";
    }
    _selected_order = selected_order;
};

void GreenFuncNph::setProbabilities(double p_length, double p_add_int, double p_rem_int, double p_add_ext, double p_rem_ext, 
    double p_swap, double p_shift, double p_stretch){
    if(!isEqual(p_length + p_add_int + p_rem_int + p_add_ext + p_rem_ext + p_swap + p_shift + p_stretch, 1)){
        std::cerr << "Invalid probabilities, total probability must add to 1.\n";
        double normalization = 1/(p_length + p_add_int + p_rem_int + p_add_ext + p_rem_ext + p_swap + p_shift + p_stretch);
        std::cout << "Probabilities are being riscaled using the value " << normalization <<".\n";
        p_length = p_length*normalization;
        p_add_int = p_add_int*normalization;
        p_rem_int = p_rem_int*normalization;
        p_add_ext = p_add_ext*normalization;
        p_rem_ext = p_rem_ext*normalization;
        p_swap = p_swap*normalization;
        p_shift = p_shift*normalization;
        p_stretch = p_stretch*normalization;
        std::cout << "New probabilities are: " << p_length << " " << p_add_int << " " << p_rem_int << " " 
        << p_add_ext << " " << p_rem_ext << " " << p_swap << " " << p_shift << " " << p_stretch << ".\n";
    }
    _p_length = p_length;
    _p_add_int = p_add_int;
    _p_rem_int = p_rem_int;
    _p_add_ext = p_add_ext;
    _p_rem_ext = p_rem_ext;
    _p_swap = p_swap;
    _p_shift = p_shift;
    _p_stretch = p_stretch;
};

void GreenFuncNph::setProbabilities(double* probs){
    if(!isEqual(probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7], 1)){
        std::cerr << "Invalid probabilities, total probability must add to 1.\n";
        double normalization = 1/(probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7]);
        std::cout << "Probabilities are being riscaled using the value " << normalization <<".\n";
        for(int i = 0; i < 8; i++){
            probs[i] = probs[i]*normalization;
        }
        std::cout << "New probabilities are: " << probs[0] << " " << probs[1] << " " << probs[2] << " " 
        << probs[3] << " " << probs[4] << " " << probs[5] << " " << probs[6] << " " << probs[7] << ".\n";
    }
    _p_length = probs[0];
    _p_add_int = probs[1];
    _p_rem_int = probs[2];
    _p_add_ext = probs[3];
    _p_rem_ext = probs[4];
    _p_swap = probs[5];
    _p_shift = probs[6];
    _p_stretch = probs[7];
};

int GreenFuncNph::findVertexPosition(long double tau){
    int position = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; i++){
        if(_vertices[i].tau < tau && _vertices[i+1].tau >= tau){
            position = i;
            return position;
        }
    }
    return -1; // return -1 if tau is not found in the vertices array
};

int GreenFuncNph::chooseInternalPhononPropagator(){
    std::uniform_int_distribution<int> distrib_unif(1,int(_current_order_int/2)); // chooses one of the internal phonon propagators at random
    int ph_propagator = distrib_unif(gen);
    int counter = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext; i++){
        if(_vertices[i].type == +1){
            counter++;
        }
        if(counter == ph_propagator){
            return i;
        }
    }
    return 0;
};

int GreenFuncNph::chooseExternalPhononPropagator(){
    std::uniform_int_distribution<int> distrib_unif(1, _current_ph_ext); // chooses one of the external phonon propagators at random
    int ph_propagator = distrib_unif(gen);
    int counter = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; i++){
        if(_vertices[i].type == -2){
            counter++;
        }
        if(counter == ph_propagator){ // to be fixed
            return i;
        }
    }
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

long double GreenFuncNph::diagramLengthUpdate(long double tau_init){
    // initialize momentum values for last propagator
    double kx = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_kx;
    double ky = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_ky;
    double kz = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_kz;

    // generate new time value for last vertex, reject if it goes out of bounds
    long double tau_fin = _last_vertex - std::log(1-drawUniformR())/(calcEnergy(kx,ky,kz)-_chem_potential + _current_ph_ext);
    if(tau_fin <= _tau_max){
        _vertices[_current_order_int + 2*_current_ph_ext + 1].tau = tau_fin;
        return tau_fin;
    }
    else{return tau_init;}
};

void GreenFuncNph::addInternalPhononPropagator(){
    if(_current_order_int+1 >= _order_int_max){return;} // reject if already at max order
    else{
        std::uniform_int_distribution<int> distrib_prop(0, _current_order_int + 2*_current_ph_ext);
        int propagator = distrib_prop(gen);
        long double tau_init = _vertices[propagator].tau;
        long double tau_end = _vertices[propagator+1].tau;
        std::uniform_real_distribution<long double> distrib_unif(tau_init, tau_end);
        long double tau_one = distrib_unif(gen);
        long double tau_two = tau_one - std::log(1-drawUniformR())/1;
        if(tau_two >= _vertices[_current_order_int + 2*_current_ph_ext + 1].tau){return;} // reject if phonon vertex goes out of bound
        else{
            // sampling momentum values for phonon propagators
            std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_two-tau_one))); // may need variance, inserted std dev
            double w_x = distrib_norm(gen);
            double w_y = distrib_norm(gen);
            double w_z = distrib_norm(gen);

            // find position of new tau values
            int index_one = findVertexPosition(tau_one);
            int index_two = findVertexPosition(tau_two);

            // control statements to check for floating point errors
            if(index_one == -1 || index_two == -1){return;} // reject if tau values are not found in the vertices array
            if(_current_ph_ext > 0){
                if( tau_one < _vertices[index_one].tau || isEqual(tau_one, _vertices[index_one].tau) 
                    || isEqual(tau_one, _vertices[index_one+1].tau) || tau_one > _vertices[index_one+1].tau){return;}
                if(tau_two < _vertices[index_two].tau || isEqual(tau_two, _vertices[index_two].tau) 
                    || isEqual(tau_two, _vertices[index_two+1].tau) || tau_two > _vertices[index_two+1].tau){return;}
            }
            // create arrays of momentum values
            double* px_init = new double[index_two + 1 - index_one];
            double* py_init = new double[index_two + 1 - index_one];
            double* pz_init = new double[index_two + 1 - index_one];
            double* px_fin = new double[index_two + 1 - index_one];
            double* py_fin = new double[index_two + 1 - index_one];
            double* pz_fin = new double[index_two + 1 - index_one];

            for(int i = index_one; i < index_two + 1; i++){
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

            double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext + 1);
            double p_A = _p_add_int*(_current_order_int/2 + 1);

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

        long double tau_one = _vertices[index_one].tau;
        long double tau_two = _vertices[index_two].tau;

        long double tau_init = _vertices[index_one-1].tau;
        long double tau_end = _tau_max;
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
            exponent_init = energy_init[0]*(tau_two-tau_one);
            exponent_fin = energy_fin[0]*(tau_two-tau_one);
        }
        else{
            exponent_init = energy_init[0]*(_vertices[index_one+1].tau-tau_one);
            exponent_fin = energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
            for(int i = 1; i < index_two - index_one - 1; i++){
                exponent_init += energy_init[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                exponent_fin += energy_fin[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
            }
            exponent_init += energy_init[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
            exponent_fin += energy_fin[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
        }

        delete[] energy_init;
        delete[] energy_fin;

        double p_A = _p_add_int*((_current_order_int - 2)/2 + 1);
        double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext - 1);
        
        double numerator = p_A*std::pow(2*M_PI,_D)*std::exp(-(tau_two-tau_one))*std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*
        std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_two-tau_one));
        double denominator = p_B*std::exp(-(exponent_fin-exponent_init+(1.)*(tau_two - tau_one)))*calcVertexStrength(w_x,w_y,w_z)*(tau_end-tau_init);

        double R_rem = numerator/denominator;
        double acceptance_ratio = std::min(1., R_rem);

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
        int total_order = _current_order_int + 2*_current_ph_ext;
        long double tau_current = _vertices[total_order + 1].tau; // length of current diagram

        long double tau_one = 0. - std::log(1-drawUniformR())/1; // time of ingoing vertex of ext phonon propagator
        if(tau_one >= tau_current){return;} // reject if it goes out of bound

        long double tau_two = tau_current + std::log(1-drawUniformR())/1; // time of outgoing vertex
        if(tau_two <= 0){return;} // reject if it goes out of bound
        if(isEqual(tau_one, tau_two)){return;} // reject if both vertices are equal (should not happen)

        std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_current-tau_two+tau_one))); // may need variance, inserted std dev
        double w_x = distrib_norm(gen);
        double w_y = distrib_norm(gen);
        double w_z = distrib_norm(gen);

        if(tau_one <= tau_two){
            int index_one = findVertexPosition(tau_one);
            int index_two = findVertexPosition(tau_two);

            // control statements to check for floating point errors
            if(index_one == -1 || index_two == -1){return;} // reject if tau values are not found in the vertices array
            if( tau_one < _vertices[index_one].tau || isEqual(tau_one, _vertices[index_one].tau) 
                || isEqual(tau_one, _vertices[index_one+1].tau) || tau_one > _vertices[index_one+1].tau){return;}
            if(tau_two < _vertices[index_two].tau || isEqual(tau_two, _vertices[index_two].tau) 
                || isEqual(tau_two, _vertices[index_two+1].tau) || tau_two > _vertices[index_two+1].tau){return;} 

            double* px_one_init = new double[index_one + 1];
            double* px_two_init = new double[total_order + 1 - index_two];

            double* py_one_init = new double[index_one + 1];
            double* py_two_init = new double[total_order + 1 - index_two];

            double* pz_one_init = new double[index_one + 1];
            double* pz_two_init = new double[total_order + 1 - index_two];

            double* px_one_fin = new double[index_one + 1];
            double* px_two_fin = new double[total_order + 1 - index_two];

            double* py_one_fin = new double[index_one + 1];
            double* py_two_fin = new double[total_order + 1 - index_two];

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
                px_two_init[i-index_two] = _propagators[i].el_propagator_kx;
                px_two_fin[i-index_two] = _propagators[i].el_propagator_kx - w_x;
                py_two_init[i-index_two] = _propagators[i].el_propagator_ky;
                py_two_fin[i-index_two] = _propagators[i].el_propagator_ky - w_y;
                pz_two_init[i-index_two] = _propagators[i].el_propagator_kz;
                pz_two_fin[i-index_two] = _propagators[i].el_propagator_kz - w_z;
            }

            double* energy_one_init = new double[index_one + 1];
            double* energy_two_init = new double[total_order + 1 - index_two];

            double* energy_one_fin = new double[index_one + 1];
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

            delete[] px_one_init; 
            delete[] px_one_fin;
            delete[] py_one_init;
            delete[] py_one_fin;
            delete[] pz_one_init;
            delete[] pz_one_fin;
            delete[] px_two_init;
            delete[] px_two_fin;
            delete[] py_two_init;
            delete[] py_two_fin;
            delete[] pz_two_init;
            delete[] pz_two_fin;

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
            delete[] energy_one_init;
            delete[] energy_one_fin;
            delete[] energy_two_init;
            delete[] energy_two_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double numerator = p_B*std::exp(-(exponent_two_fin + exponent_one_fin - exponent_two_init - exponent_one_init + 
                1.*(tau_current-tau_two+tau_one)))*calcVertexStrength(w_x, w_y, w_z);

            double denominator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(tau_current-tau_two))*
            std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one));

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
                for(int i = 0; i < index_one + 1; i++){
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
            // control statements to check for floating point errors
            if(index_one == -1 || index_two == -1){return;} // reject if tau values are not found in the vertices array
            if( tau_two < _vertices[index_one].tau || isEqual(tau_two, _vertices[index_one].tau) 
                || isEqual(tau_two, _vertices[index_one+1].tau) || tau_two > _vertices[index_one+1].tau){return;}
            if(tau_one < _vertices[index_two].tau || isEqual(tau_one, _vertices[index_two].tau) 
                || isEqual(tau_one, _vertices[index_two+1].tau) || tau_one > _vertices[index_two+1].tau){return;}

            int total_order = _current_order_int + 2*_current_ph_ext;

            double* px_init = new double[total_order+1];
            double* py_init = new double[total_order+1];
            double* pz_init = new double[total_order+1];
            
            double* px_fin = new double[total_order+3];
            double* py_fin = new double[total_order+3];
            double* pz_fin = new double[total_order+3];

            for(int i = 0; i < total_order + 1; i++){
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
            for(int i = index_two + 2; i < total_order + 3; i++){
                px_fin[i] = _propagators[i-2].el_propagator_kx - w_x;
                py_fin[i] = _propagators[i-2].el_propagator_ky - w_y;
                pz_fin[i] = _propagators[i-2].el_propagator_kz - w_z;
            }

            double* energy_init = new double[total_order+1];
            double* energy_fin = new double[total_order+3];

            for(int i = 0; i < total_order + 1; i++){
                energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
            }
            for(int i = 0; i < total_order + 3; i++){
                energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
            }

            delete[] px_init;
            delete[] py_init;
            delete[] pz_init;

            double exponent_init = 0.;
            double exponent_fin = 0.;

            for(int i = 0; i < total_order + 1; i++){
                exponent_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            for(int i = 0; i < index_one; i++){
                exponent_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }
            exponent_fin += energy_fin[index_one]*(tau_two - _vertices[index_one].tau);
            exponent_fin += energy_fin[index_one+1]*(_vertices[index_one+1].tau - tau_two);
            for(int i = index_one + 2; i < index_two + 1; i++){
                exponent_fin += energy_fin[i]*(_vertices[i].tau - _vertices[i-1].tau);
            }
            exponent_fin += energy_fin[index_two+1]*(tau_one - _vertices[index_two].tau);
            exponent_fin += energy_fin[index_two+2]*(_vertices[index_two+1].tau - tau_one);
            for(int i = index_two + 3; i < total_order + 3; i++){
                exponent_fin += energy_fin[i]*(_vertices[i-1].tau - _vertices[i-2].tau);
            }
            
            // delete array of energies
            delete[] energy_init;
            delete[] energy_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double numerator = p_B*std::exp(-(exponent_fin - exponent_init + 1.*(tau_current-tau_two+tau_one)))
            *calcVertexStrength(w_x, w_y, w_z);

            double denominator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(tau_current-tau_two))*
            std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one));

            double R_add = numerator/denominator;
            double acceptance_ratio = std::min(1.,R_add);

            if(drawUniformR() < acceptance_ratio){delete[] px_fin; delete[] py_fin; delete[] pz_fin; return;}
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
                for(int i = 0; i < total_order + 3; i++){
                    _propagators[i].el_propagator_kx = px_fin[i];
                    _propagators[i].el_propagator_ky = py_fin[i];
                    _propagators[i].el_propagator_kz = pz_fin[i];
                }

                delete[] px_fin;
                delete[] py_fin;
                delete[] pz_fin;

                _current_ph_ext += 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
    }
};

void GreenFuncNph::removeExternalPhononPropagator(){
    if(_current_ph_ext <= 0){return;} // reject if already at order 0
    else{
        // indexes of initial and final vertices of a random internal phonon propagator
        int index_one = chooseExternalPhononPropagator();
        int index_two = _vertices[index_one].linked;

        if(index_one >= index_two){
            int temp = index_one;
            index_one = index_two;
            index_two = temp;
        }
        // retrieve time values of the two vertices
        long double tau_one = _vertices[index_one].tau;
        long double tau_two = _vertices[index_two].tau;

        // retrieve phonon momentum
        double w_x = _vertices[index_one].wx;
        double w_y = _vertices[index_one].wy;
        double w_z = _vertices[index_one].wz;

        int total_order = _current_order_int + 2*_current_ph_ext;

        if(_vertices[index_one].type == -2){

            double* px_one_init = new double[index_one + 1];
            double* py_one_init = new double[index_one + 1];
            double* pz_one_init = new double[index_one + 1];

            double* px_one_fin = new double[index_one + 1];
            double* py_one_fin = new double[index_one + 1];
            double* pz_one_fin = new double[index_one + 1];

            double* px_two_init = new double[total_order - index_two + 1];
            double* py_two_init = new double[total_order - index_two + 1];
            double* pz_two_init = new double[total_order - index_two + 1];

            double* px_two_fin = new double[total_order - index_two + 1];
            double* py_two_fin = new double[total_order - index_two + 1];
            double* pz_two_fin = new double[total_order - index_two + 1];

            for(int i = 0; i < index_one + 1; i++){
                px_one_init[i] = _propagators[i].el_propagator_kx + w_x;
                py_one_init[i] = _propagators[i].el_propagator_ky + w_y;
                pz_one_init[i] = _propagators[i].el_propagator_kz + w_z;
                px_one_fin[i] = _propagators[i].el_propagator_kx;
                py_one_fin[i] = _propagators[i].el_propagator_ky;
                pz_one_fin[i] = _propagators[i].el_propagator_kz;
            }

            for(int i = index_two; i < total_order + 1; i++){
                px_two_init[i-index_two] = _propagators[i].el_propagator_kx + w_x;
                py_two_init[i-index_two] = _propagators[i].el_propagator_ky + w_y;
                pz_two_init[i-index_two] = _propagators[i].el_propagator_kz + w_z;
                px_two_fin[i-index_two] = _propagators[i].el_propagator_kx;
                py_two_fin[i-index_two] = _propagators[i].el_propagator_ky;
                pz_two_fin[i-index_two] = _propagators[i].el_propagator_kz;
            }

            double* energy_one_init = new double[index_one + 1];
            double* energy_one_fin = new double[index_one + 1];
            double* energy_two_init = new double[total_order - index_two + 1];
            double* energy_two_fin = new double[total_order - index_two + 1];

            for(int i = 0; i < index_one + 1; i++){
                energy_one_init[i] = calcEnergy(px_one_init[i], py_one_init[i], pz_one_init[i]);
                energy_one_fin[i] = calcEnergy(px_one_fin[i], py_one_fin[i], pz_one_fin[i]);
            }
            
            for(int i = 0; i < total_order - index_two + 1; i++){   
                energy_two_init[i] = calcEnergy(px_two_init[i], py_two_init[i], pz_two_init[i]);
                energy_two_fin[i] = calcEnergy(px_two_fin[i], py_two_fin[i], pz_two_fin[i]);
            }

            delete[] px_one_init;
            delete[] py_one_init;
            delete[] pz_one_init;
            delete[] px_one_fin;
            delete[] py_one_fin;
            delete[] pz_one_fin;
            delete[] px_two_init;
            delete[] py_two_init;
            delete[] pz_two_init;
            delete[] px_two_fin;
            delete[] py_two_fin;
            delete[] pz_two_fin;


            double exponent_one_init = 0.;
            double exponent_one_fin = 0.;
            double exponent_two_init = 0.;
            double exponent_two_fin = 0.;

            for(int i = 0; i < index_one + 1; i++){
                exponent_one_init += energy_one_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                exponent_one_fin += energy_one_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            for(int i = index_two; i < total_order + 1; i++){
                exponent_two_init += energy_two_init[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
                exponent_two_fin += energy_two_fin[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            delete[] energy_one_init;
            delete[] energy_one_fin;
            delete[] energy_two_init;
            delete[] energy_two_fin;

            long double tau_current = _vertices[total_order+1].tau; // length of current diagram

            double p_A = _p_add_ext*_current_ph_ext;
            double p_B = _p_rem_ext;

            double numerator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(tau_current-tau_two))*
            std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one));

            double denominator = p_B*std::exp(-(exponent_two_fin + exponent_one_fin - exponent_two_init - exponent_one_init + 
            1.*(tau_current-tau_two+tau_one)))*calcVertexStrength(w_x, w_y, w_z);
            
            double R_rem = numerator/denominator;
            double acceptance_ratio = std::min(1., R_rem);

            if(drawUniformR() > acceptance_ratio){return;}
            else{
                for(int i=0; i<index_one;i++){
                    _propagators[i].el_propagator_kx += w_x;
                    _propagators[i].el_propagator_ky += w_y;
                    _propagators[i].el_propagator_kz += w_z;
                }
                for(int i=index_two; i<total_order+1;i++){
                    _propagators[i].el_propagator_kx += w_x;
                    _propagators[i].el_propagator_ky += w_y;
                    _propagators[i].el_propagator_kz += w_z;
                }

                phVertexRemoveRoom(index_one, index_two); // remove room in vertices array
                propagatorArrayRemoveRoom(index_one, index_two); // remove room in electron propagators array

                _propagators[total_order-1].el_propagator_kx = 0;
                _propagators[total_order-1].el_propagator_ky = 0;
                _propagators[total_order-1].el_propagator_kz = 0;

                _propagators[total_order].el_propagator_kx = 0;
                _propagators[total_order].el_propagator_ky = 0;
                _propagators[total_order].el_propagator_kz = 0;

                _current_ph_ext -= 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
        else if(_vertices[index_one].type == +2){

                double* px_init = new double[total_order + 1];
                double* py_init = new double[total_order + 1];
                double* pz_init = new double[total_order + 1];
    
                double* px_fin = new double[total_order + 1];
                double* py_fin = new double[total_order + 1];
                double* pz_fin = new double[total_order + 1];

                for(int i = 0; i < index_one; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + w_z;
                }
                for(int i = index_one; i < index_two; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + 2*w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + 2*w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + 2*w_z;
                }
                for(int i = index_two - 1; i < total_order + 1; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + w_z;
                }

                for(int i = 0; i < total_order + 1; i++){
                    px_fin[i] = _propagators[i].el_propagator_kx;
                    py_fin[i] = _propagators[i].el_propagator_ky;
                    pz_fin[i] = _propagators[i].el_propagator_kz;
                }
                
                double* energy_init = new double[total_order+1];
                double* energy_fin = new double[total_order+1];

                for(int i = 0; i < total_order + 1; i++){
                    energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
                    energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
                }

                delete[] px_fin;
                delete[] py_fin;
                delete[] pz_fin;

                double exponent_init = 0.;
                double exponent_fin = 0.;

                for(int i = 0; i < total_order + 1; i++){
                    exponent_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                    exponent_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
                }

                delete[] energy_init;
                delete[] energy_fin;

                long double tau_current = _vertices[total_order+1].tau; // length of current diagram

                double p_A = _p_add_ext*_current_ph_ext;
                double p_B = _p_rem_ext;

                double numerator = p_A*std::pow(2*M_PI,_D)*1*std::exp(-tau_one)*1*std::exp(-(tau_current-tau_two))*
                std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)*
                std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one));

                double denominator = p_B*std::exp(-(exponent_fin - exponent_init + 1.*(tau_current-tau_two+tau_one)))
                *calcVertexStrength(w_x, w_y, w_z);

                double R_rem = numerator/denominator;
                double acceptance_ratio = std::min(1., R_rem);

                if(drawUniformR() > acceptance_ratio){delete[] px_init; delete[] py_init; delete[] pz_init; return;}
                else{
                    for(int i=0; i < total_order+1; i++){
                        _propagators[i].el_propagator_kx = px_init[i];
                        _propagators[i].el_propagator_ky = py_init[i];
                        _propagators[i].el_propagator_kz = pz_init[i];
                    }

                    delete[] px_init;
                    delete[] py_init;
                    delete[] pz_init;

                    phVertexRemoveRoom(index_one, index_two); // remove room in vertices array
                    propagatorArrayRemoveRoom(index_one, index_two); // remove room in electron propagators array

                    _propagators[total_order-1].el_propagator_kx = 0;
                    _propagators[total_order-1].el_propagator_ky = 0;
                    _propagators[total_order-1].el_propagator_kz = 0;

                    _propagators[total_order].el_propagator_kx = 0;
                    _propagators[total_order].el_propagator_ky = 0;
                    _propagators[total_order].el_propagator_kz = 0;

                    _current_ph_ext -= 1; // update current number of external phonons
                    findLastPhVertex();
                    return;
                }            
        }
        else{return;}
    }
};

void GreenFuncNph::swapPhononPropagator(){
    if(_current_order_int < 4){return;} // swap not possible if internal order is less than 4
    else{
        std::uniform_int_distribution<> distrib(1, _current_order_int + 2*_current_ph_ext - 1); // choose random internal propagator
        int index_one = distrib(gen); // choose random internal propagatorS
        
        if(_vertices[index_one].type != +1 && _vertices[index_one].type != -1){return;} // reject if not belonging to internal phonon propagator
        if(_vertices[index_one+1].type != +1 && _vertices[index_one+1].type != -1){return;} // reject if not belonging to internal phonon propagator
        if(_vertices[index_one].linked == index_one + 1 || _vertices[index_one].linked == -1){return;} // reject if the two vertices are linked

        int index_two = index_one +1;

        // get values of first vertex
        int c1 = _vertices[index_one].type;
        long double wx1 = _vertices[index_one].wx;
        long double wy1 = _vertices[index_one].wy;
        long double wz1 = _vertices[index_one].wz;
        long double tau1 = _vertices[index_one].tau;

        // get values of second vertex
        int c2 = _vertices[index_two].type;
        long double wx2 = _vertices[index_two].wx;
        long double wy2 = _vertices[index_two].wy;
        long double wz2 = _vertices[index_two].wz;
        long double tau2 = _vertices[index_two].tau;

        // get momentum of propagator
        long double kx = _propagators[index_one].el_propagator_kx;
        long double ky = _propagators[index_one].el_propagator_ky;
        long double kz = _propagators[index_one].el_propagator_kz;

        double energy_final_el = calcEnergy(kx+c1*wx1-c2*wx2, ky+c1*wy1-c2*wy2, kz+c1*wz1-c2*wz2);
        double energy_initial_el = calcEnergy(kx, ky, kz);

        double R_swap = std::exp(-(energy_final_el - energy_initial_el - 1.*(c1-c2))*(tau2-tau1));
        double acceptance_ratio = std::min(1.,R_swap);

        if(drawUniformR() > acceptance_ratio){return;}
        else{
            _propagators[index_one].el_propagator_kx += c1*wx1-c2*wx2;
            _propagators[index_one].el_propagator_ky += c1*wy1-c2*wy2;
            _propagators[index_one].el_propagator_kz += c1*wz1-c2*wz2;

            int linked1 = _vertices[index_one].linked;
            int linked2 = _vertices[index_two].linked;

            // assing new links to conjugate vertices
            _vertices[linked1].linked = index_two;         
            _vertices[linked2].linked = index_one;

            _vertices[index_one].wx = wx2;
            _vertices[index_one].wy = wy2;
            _vertices[index_one].wz = wz2;
            _vertices[index_one].type = c2;
            _vertices[index_one].linked = linked2;
            _vertices[index_one].tau = tau1;
                
            _vertices[index_two].wx = wx1;
            _vertices[index_two].wy = wy1;
            _vertices[index_two].wz = wz1;
            _vertices[index_two].type = c1;
            _vertices[index_two].linked = linked1;
            _vertices[index_two].tau = tau2;
        }
    }
};

void GreenFuncNph::shiftPhononPropagator(){
    if(_current_order_int + _current_ph_ext == 0){return;} // reject if no vertices are present
    else{
        int total_order = _current_order_int + 2*_current_ph_ext;
        std::uniform_int_distribution<> distrib(1, total_order);
        int vertex_index = distrib(gen); // choose random vertex

        int c = _vertices[vertex_index].type;
        if(c == 2){c = 1;}
        if(c == -2){c = -1;}

        long double tau_init = _vertices[vertex_index - 1].tau;
        long double tau_fin = _vertices[vertex_index + 1].tau;

        long double kx_incoming = _propagators[vertex_index - 1].el_propagator_kx;
        long double ky_incoming = _propagators[vertex_index - 1].el_propagator_ky;
        long double kz_incoming = _propagators[vertex_index - 1].el_propagator_kz;

        long double kx_outgoing = _propagators[vertex_index].el_propagator_kx;
        long double ky_outgoing = _propagators[vertex_index].el_propagator_ky;
        long double kz_outgoing = _propagators[vertex_index].el_propagator_kz;

        double energy_delta = calcEnergy(kx_incoming, ky_incoming, kz_incoming) - calcEnergy(kx_outgoing, ky_outgoing, kz_outgoing) - 1.*c;
        long double tau_new = tau_init - std::log(1 - drawUniformR()*(1 - std::exp(-energy_delta*(tau_fin - tau_init))))/energy_delta;
        if(isEqual(tau_new, tau_init) || isEqual(tau_new, tau_fin) || tau_new > tau_fin){return;}
        _vertices[vertex_index].tau = tau_new; // assign new time value to vertex
    }
};

long double GreenFuncNph::stretchDiagramLength(long double tau_init){
    // initialize momentum values for first electron propagator
    long double kx = _propagators[0].el_propagator_kx;
    long double ky = _propagators[0].el_propagator_ky;
    long double kz = _propagators[0].el_propagator_kz;

    int total_order = _current_order_int + 2*_current_ph_ext;
    long double* new_taus = new long double[total_order+2];
    new_taus[0] = 0; // first vertex time value is always 0

    for(int i = 1; i < total_order+1; i++){
        // assign new time values to every vertex
        int c = _vertices[i].type;
        if(c == 2){c = 1;}
        if(c == -2){c = -1;}

        // value of left vertex is retrieved from new proposed values, right value from old ones (vertex array)
        long double tau_one = new_taus[i-1];
        long double tau_two = _vertices[i+1].tau;

        long double kx_incoming = _propagators[i-1].el_propagator_kx;
        long double ky_incoming = _propagators[i-1].el_propagator_ky;
        long double kz_incoming = _propagators[i-1].el_propagator_kz;

        long double kx_outgoing = _propagators[i].el_propagator_kx;
        long double ky_outgoing = _propagators[i].el_propagator_ky;
        long double kz_outgoing = _propagators[i].el_propagator_kz;

        double energy_delta = calcEnergy(kx_incoming, ky_incoming, kz_incoming) - calcEnergy(kx_outgoing, ky_outgoing, kz_outgoing) - c;
            
        new_taus[i] = tau_one - std::log(1 - drawUniformR()*(1 - std::exp(-energy_delta*(tau_two - tau_one))))/energy_delta;
        if(isEqual(new_taus[i], tau_one) || isEqual(new_taus[i], tau_two) || new_taus[i] > tau_two){delete[] new_taus; return tau_init;}
    }

    new_taus[total_order+1] = new_taus[total_order] - std::log(1-drawUniformR())/(calcEnergy(kx,ky,kz) - _chem_potential + _current_ph_ext);
    
    if(isEqual(new_taus[total_order], new_taus[total_order+1]) || isEqual(new_taus[total_order+1], _tau_max) 
       || new_taus[total_order+1] >= _tau_max){delete[] new_taus; return tau_init;}
    else{
        for(int i = 0; i < total_order + 2; i++){
            _vertices[i].tau = new_taus[i]; // assign new time values to vertices
        }
    }
    for(int i = 0; i < total_order + 2; i++){
        _vertices[i].tau = new_taus[i];} // assign new time values to vertices
    long double tau_fin = new_taus[total_order+1]; // assign new time value to last vertex
    delete[] new_taus;
    return tau_fin; // return new length of diagram
};

void GreenFuncNph::markovChainMC(
    unsigned long long int N_diags = 0,
    bool gf_exact = false,
    bool histo = false,
    bool gs_energy = false,
    bool effective_mass = false,
    bool Z_factor = false
    ){
    if(N_diags == 0){N_diags = _N_diags;}

    if((!(isEqual(_kx,0)) || !(isEqual(_ky,0)) || !(isEqual(_kz,0))) && effective_mass){
        std::cerr << "Warning: kx, ky and kz should be equal to 0 to calculate effective mass." << std::endl;
        std::cerr << "Effective mass calculation is not possible." << std::endl;
        effective_mass = false;
    }

    if(Z_factor && _ph_ext_max == 0){
        std::cerr << "Warning: number of maximum external phonon must be greater than 0 to calculate Z factor." << std::endl;
        std::cerr << "Z factor calculation is not possible." << std::endl;
        Z_factor = false;
    }

    // print simulation parameters
    std::cout <<"Starting simulation..." << std::endl;
    std::cout << "Number of thermalization steps: " << _relax_steps << std::endl;
    std::cout << "Number of diagrams to be generated: " << N_diags << std::endl;
    std::cout << "Maximum length of diagram: " << _tau_max << std::endl;
    std::cout << "Maximum number of internal phonons: " << _order_int_max/2 << std::endl;
    std::cout << "Maximum number of external phonons: " << _ph_ext_max << std::endl;
    std::cout << "Maximum diagram order: " << _order_int_max + 2*_ph_ext_max << std::endl;
    std::cout << "Coupling strength: " << _alpha << ", chemical potential: " << _chem_potential << std::endl;
    std::cout << "Number of dimensions: " << _D << std::endl;
    if(_D == 3){
        std::cout << "Momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
    }
    if(_D == 2){
        std::cout << "Momentum: kx = " << _kx << ", ky = " << _ky << std::endl;
    }
    // print MC update probabilities
    std::cout << "Update probabilities:" << std::endl;
    std::cout << "length update: " << _p_length << ", add internal update: " << _p_add_int <<
    ", remove internal update: " << _p_rem_int << "," << std::endl;
    std::cout << "add external update: " << _p_add_ext << ", remove external update: " << _p_rem_ext <<
    ", swap update: " << _p_swap << "," << std::endl;
    std::cout << "shift update: " << _p_shift << ", stretch update: " << _p_stretch << std::endl; 

    double bin_width_inv = 1./_bin_width;

    if(gf_exact){
        std::cout << "Green Function will be calculated exactly." << std::endl;
        std::cout << "Number of computed points: " << _num_points << std::endl;
        std::cout << "The selected external number of phonons is: " << _selected_order << std::endl;
        _points = new long double[_num_points];
        _points_gf_exact = new long double[_num_points];
        _gf_exact_written = true; // for destructor

        // initialize GF
        for(int i=0; i<_num_points; i++){
            _points[i] = _points_center + i*_points_step;
            _points_gf_exact[i] = 0;
        }
    }

    if(histo){
        std::cout << "Green Function will be computed using the histogram method" << std::endl;
        std::cout << "Number of bins: " << _N_bins << std::endl;
        _histogram = new double[_N_bins];
        _bin_count = new unsigned long long int[_N_bins];
        for(int i=0; i<_N_bins; i++){
            _histogram[i] = _bin_center + i*_bin_width;
            _bin_count[i] = 0;
        }
        _histogram_calculated = true; // for destructor
    }

    if(gs_energy){
        std::cout << "Ground state energy will be calculated using the exact estimator" << std::endl;
        std::cout << "Coupling strength: " << _alpha << ", chemical potential: " << _chem_potential << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
    }

    if(effective_mass){
        std::cout << "Effective mass will be calculated using the exact estimator" << std::endl;
        std::cout << "Coupling strength: " << _alpha << ", chemical potential: " << _chem_potential << ", tau cutoff: " << _tau_cutoff_mass << std::endl;
    }

    if(Z_factor){
        std::cout << "Z factor will be calculated using the exact estimator" << std::endl;
        std::cout << "Coupling strength: " << _alpha << ", chemical potential: " << _chem_potential << ", max number of phonons: " << 
        _ph_ext_max << std::endl;
        initializeZFactorArray();
    }

    if(_write_diagrams){
        if(N_diags > 25000){
            _write_diagrams = false; // if too many diagrams are generated they are not printed to txt file
            std::cerr << "Warning: too many diagrams generated (> 25000), diagrams will not be printed to .txt file." << std::endl;
        }
        else{
            std::cout << "The diagrams generated in the simulation process will be printed in the Diagrams.txt file" << std::endl;
        }    
    }
    // input variables
    std::uniform_real_distribution<long double> distrib(0,_tau_max);
    long double tau_length = distrib(gen);
    double r = 0.5;
    unsigned long long int i = 0;

    std::cout << "Starting thermalization process" << std::endl;
    while(i < _relax_steps){
        r = drawUniformR();
        if(r <= _p_length){
            tau_length = diagramLengthUpdate(tau_length);
        }
        else if(r <= _p_length + _p_add_int){
            addInternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int){
            removeInternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext){
            addExternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext){
            removeExternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap){
            swapPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift){
            shiftPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift + _p_stretch){
            tau_length = stretchDiagramLength(tau_length);
        }

        if(i%(_relax_steps/10) == 0){
            std::cout << "\r" << "Thermalization progress: " << (double)i/(double)_relax_steps*100. << "%";
        }

        i++;
    }
    std::cout << "\r" << "Thermalization progress: 100%";
    std::cout << std::endl;
    std::cout << "Thermalization process finished" << std::endl;

    i = 0;
    std::cout << "Starting simulation process" << std::endl;
    while(i < N_diags){
        r = drawUniformR();
        if(_write_diagrams){
            writeChosenUpdate("Updates.txt", i, r);
        }
        if(r <= _p_length){
            tau_length = diagramLengthUpdate(tau_length);
        }
        else if(r <= _p_length + _p_add_int){
            addInternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int){
            removeInternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext){
            addExternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext){
            removeExternalPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap){
            swapPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift){
            shiftPhononPropagator();
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift + _p_stretch){
            tau_length = stretchDiagramLength(tau_length);
        }

        if(gf_exact && _current_ph_ext == _selected_order){
            exactEstimatorGF(tau_length, _selected_order); // calculate Green function
            _gf_exact_count++; // count number of diagrams for normalization
        }

        if(histo){
            // select correct bin for histogram
            tau_length = (tau_length < 0.) ? 0. : (tau_length >= _tau_max) ? _tau_max - 1e-9 : tau_length;
            int bin = (int)((tau_length - 0.) * bin_width_inv);
            _bin_count[bin]++;
        }

        if(gs_energy){
            _gs_energy += calcGroundStateEnergy(tau_length); // accumulate energy of diagrams
        }

        if(effective_mass){
            _effective_mass += calcEffectiveMass(tau_length); // accumulate effective mass of diagrams
        }

        if(Z_factor){updateZFactor();} // accumulate Z factor data

        if(_write_diagrams){
            writeDiagram("Diagrams.txt", i, r); // debug method to visualize diagram structure
        }

        if(i%(N_diags/100) == 0){
            std::cout << "\r" << "Simulation progress: " << (double)i/(double)N_diags*100. << "%";
        }

        if(_current_order_int == 0 && _current_ph_ext == 0){
                _N0++;
        }

        i++;
    }
    std::cout << "\r" << "Simulation progress: 100%";
    std::cout << std::endl;

    if(gf_exact){
        calcNormConst();
        std::cout << "Exact Green's function computed." << std::endl;
        for(int i=0; i<_num_points; i++){
            //_points_gf_exact[i] = _points_gf_exact[i]/((double)_gf_exact_count); // right normalization
            _points_gf_exact[i] = _points_gf_exact[i]*_norm_const/_N0; // right normalization
        }
        std::string a = "GF^(";
        auto b = std::to_string(_selected_order);
        std::string c = ")";
        std::string d = "exact.txt";
        writeExactGF(a+b+c+d); // write Green function to file
        std::cout << "exaxt GF count = " << _gf_exact_count << std::endl;
    }

    if(histo){
        std::cout << "Histogram computed." << std::endl;
        calcNormConst();
        normalizeHistogram();
        writeHistogram("histo.txt");
    }

    if(gs_energy){
        _gs_energy = _gs_energy/(double)_gs_energy_count; // average energy of diagrams
        std::cout << "Ground state energy of the system is: " << _gs_energy << ". Input parameters are: kx = " << _kx << 
        ", ky = " << _ky << ", kz = " << _kz << " coupling strength = " << _alpha << 
        " chemical potential = " << _chem_potential << std::endl;
        std::cout << " minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;

        std::string filename = "gs_energy.txt";
        std::ofstream file(filename, std::ofstream::app);

        if(!file.is_open()){
            std::cerr << "Could not gs_energy.txt open file " << filename << std::endl;
        }
        else{
            file << "Ground state energy of the system is: " << _gs_energy << " . Input parameters are: kx = " << _kx << 
            ", ky = " << _ky << ", kz = " << _kz << " coupling strength = " << _alpha << " chemical potential = " <<
            _chem_potential << std::endl;
            file << " minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;
            file << std::endl;
            file.close();
        }
    }

    if(effective_mass){
        long double effective_mass_inv = _effective_mass/(double)_effective_mass_count; // average effective mass of diagrams
        _effective_mass = 1./effective_mass_inv; // effective mass is inverse of the value calculated
        std::cout << "Effective mass of system is: " << _effective_mass << ". Input parameters are: coupling strength = " << _alpha << 
        " chemical potential = " << _chem_potential << " minimum length of diagrams for which effective mass is computed = " 
        << _tau_cutoff_mass << ". " << std::endl;

        std::cout << "Inverse effective mass of system is: " << effective_mass_inv << ". Input parameters are: coupling strength = " << _alpha << 
        " chemical potential = " << _chem_potential << " minimum length of diagrams for which effective mass is computed = " 
        << _tau_cutoff_mass << ". " << std::endl;

        std::string filename = "effective_mass.txt";
        std::ofstream file(filename, std::ofstream::app);

        if(!file.is_open()){
            std::cerr << "Could not effective_mass.txt open file " << filename << std::endl;
        }
        else{
            file << "Effective mass of the system is: " << _effective_mass << ". Input parameters are: coupling strength = " << _alpha 
            << " chemical potential = " << _chem_potential << " minimum length of diagrams for which effective mass is computed = " 
            << _tau_cutoff_mass << "." << std::endl;
            file << std::endl;

            file << "Inverse effective mass of the system is: " << effective_mass_inv << ". Input parameters are: coupling strength = " << _alpha 
            << " chemical potential = " << _chem_potential << " minimum length of diagrams for which effective mass is computed = " 
            << _tau_cutoff_mass << "." << std::endl;
            file << std::endl;

            file.close();
        }
    }

    if(Z_factor){
        std::string a = "Z_factor_alpha";
        auto b = std::to_string(_alpha);
        std::string c = ".txt";
        std::string filename = a + b + c;
        writeZFactor(filename);
    }
    std::cout << "Simulation finished!" << std::endl;
};

void GreenFuncNph::calcNormConst(){
    double numerator = 1 - std::exp(-(calcEnergy(_kx,_ky,_kz)-_chem_potential)*_tau_max);
    double denominator = calcEnergy(_kx,_ky,_kz)-_chem_potential;
    std::cout << "Normalization constant calculated." << std::endl;
    _norm_const = numerator/denominator;
};

void GreenFuncNph::normalizeHistogram(){
    _green_func = new double[_N_bins];
    for(int i=0; i<_N_bins; i++){
        _green_func[i] = _bin_count[i]/(_N0*_bin_width);
        _green_func[i] = _green_func[i]*_norm_const;
    }
    _histogram_normalized = true;
    std::cout << "Imaginary time Green's function computed (histogram method)." << std::endl;
};

double GreenFuncNph::calcGroundStateEnergy(long double tau_length){
    if(tau_length <= _tau_cutoff_energy){return 0;} // reject if below cutoff
    else{
        int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
        double electron_action = 0., phonon_action = 0.;

        // compute electron bare propagators action
        for(int i=0; i<current_order+1; i++){
            double electron_energy = calcEnergy(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, _propagators[i].el_propagator_kz);
            electron_action += electron_energy*(_vertices[i+1].tau - _vertices[i].tau);
        }

        // compute phonon bare propagators action
        int i = 0;
        int int_count = 0;
        int ext_count = 0;
        bool int_flag = false;
        bool ext_flag = false;

        while(i < current_order + 1 && !(int_flag && ext_flag)){
            if(_vertices[i].type == +1){
                int index_two = _vertices[i].linked;
                long double tau_one = _vertices[i].tau;  
                long double tau_two = _vertices[index_two].tau;
                phonon_action += tau_two - tau_one;
                int_count++;
                if(int_count == _current_order_int/2){int_flag = true;}
            }
            else if(_vertices[i].type == -2){
                int index_two = _vertices[i].linked;
                long double tau_one = _vertices[i].tau;  
                long double tau_two = _vertices[index_two].tau;
                phonon_action += tau_length + tau_one - tau_two;
                ext_count++;
                if(ext_count == _current_ph_ext){ext_flag = true;}
            }
            i++;
        }
        double diagram_energy = (electron_action + phonon_action - (double)current_order)/tau_length; // energy of current diagram
        _gs_energy_count++;
        return diagram_energy; // return energy of current diagram
    }
};

double GreenFuncNph::calcEffectiveMass(long double tau_length){
    if(tau_length <= _tau_cutoff_mass){return 0;}
    else{
        int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
        double electron_average_kx = 0, electron_average_ky = 0, electron_average_kz = 0;
        for(int i=0; i<current_order+1; i++){
            electron_average_kx += _propagators[i].el_propagator_kx*(_vertices[i+1].tau - _vertices[i].tau);
            electron_average_ky += _propagators[i].el_propagator_ky*(_vertices[i+1].tau - _vertices[i].tau);
            electron_average_kz += _propagators[i].el_propagator_kz*(_vertices[i+1].tau - _vertices[i].tau);
        }
        electron_average_kx = electron_average_kx/tau_length;
        electron_average_ky = electron_average_ky/tau_length;
        electron_average_kz = electron_average_kz/tau_length;

        _effective_mass_count++;
        // return inverse of effective mass, _D dimensionality of the system
        return (1-tau_length*(std::pow(electron_average_kx,2) + std::pow(electron_average_ky,2) + std::pow(electron_average_kz,2))/_D);
    }
};

void GreenFuncNph::initializeZFactorArray(){
    _Z_factor = new unsigned long long int[_ph_ext_max + 1];
    
    for(int i=0; i<_ph_ext_max + 1; i++){
        _Z_factor[i] = 0;
    }
    _Z_factor_calculated = true;
};

void GreenFuncNph::updateZFactor(){
    _Z_factor[_current_ph_ext] += 1;
};

void GreenFuncNph::exactEstimatorGF(long double tau_length, int ext_phonon_order){
    int current_order = _current_order_int + 2*ext_phonon_order; // total order of diagrams (number of phonon vertices)
    double electron_action = 0, phonon_action = 0;

    // compute electron bare propagators action
    for(int i=0; i<current_order+1; i++){
        double electron_energy = calcEnergy(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, _propagators[i].el_propagator_kz);
        electron_action += electron_energy*(_vertices[i+1].tau - _vertices[i].tau);
    }

    // compute phonon bare propagators action
    int i = 0;
    int int_count = 0;
    int ext_count = 0;
    bool int_flag = false;
    bool ext_flag = false;

    while(i < current_order + 1 && !(int_flag && ext_flag)){
        if(_vertices[i].type == +1){
            int index_two = _vertices[i].linked;
            long double tau_one = _vertices[i].tau;  
            long double tau_two = _vertices[index_two].tau;
            phonon_action += tau_two - tau_one;
            int_count++;
            if(int_count == _current_order_int/2){int_flag = true;}
        }
        else if(_vertices[i].type == -2){
            int index_two = _vertices[i].linked;
            long double tau_one = _vertices[i].tau;  
            long double tau_two = _vertices[index_two].tau;
            phonon_action += tau_length + tau_one - tau_two;
            ext_count++;
            if(ext_count == ext_phonon_order){ext_flag = true;}
        }
        i++;
    }

    // compute Green function with exact estimator
    tau_length = (tau_length < 0.) ? 0. : (tau_length >= _tau_max) ? _tau_max - 1e-9 : tau_length;
    int bin = (int)((tau_length - 0.) * 1./_points_step);

    long double prefactor = std::pow(1 + (_points[bin] - tau_length)/tau_length, current_order);
    long double exponential = std::exp(-((_points[bin] - tau_length)/tau_length)*(electron_action + phonon_action));
    long double diagrams_ratio = prefactor*exponential;

    _points_gf_exact[bin] += diagrams_ratio/(_points_step); // accumulate Green function value
};

void GreenFuncNph::writeExactGF(std::string filename) const {
    std::ofstream file;
    if(filename.empty()){
        std::cout << "Enter filename: ";
        std::cin >> filename;
    }

    file.open(filename);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }

    file << "# Exact Green's function calculated for number of external phonons " << _selected_order << " and coupling strength " << _alpha <<".\n";
    file << "# kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << ", chemical potential = " << _chem_potential << "\n";

    for(int i=0; i<_num_points; i++){
        file << _points[i] << " " << _points_gf_exact[i] << "\n";
    }
    file.close();
    std::cout << "Exact Green's function written to file " << filename << "." << std::endl;
}

void GreenFuncNph::writeHistogram(std::string filename) const {
    std::ofstream file;
    if(filename.empty()){
        std::cout << "Enter filename: ";
        std::cin >> filename;
    }

    file.open(filename);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }

    for(int i=0; i<_N_bins; i++){
        file << _histogram[i] << " " << _green_func[i] << "\n";
    }
    file.close();
    std::cout << "Histogram written to file " << filename << "." << std::endl;
};

void GreenFuncNph::writeZFactor(std::string filename) const {
    std::ofstream file;
    file.open(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }
    file << "# Z factor calculated for max number of external phonons " << _ph_ext_max << ", coupling strength " << _alpha << 
    " and chemical potential: " << _chem_potential << ".\n";
    file << "# N_ext Z_factor(N_ext)\n";
    for(int i = 0; i < _ph_ext_max; i++){
        file << i << " " << (double)_Z_factor[i]/(double)_N_diags << "\n";
    }

    file.close();
};

void GreenFuncNph::writeDiagram(std::string filename, int i, double r) const {
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

    file << "ext phonons: " << _current_ph_ext << " int order: " << _current_order_int << " chosen update: " << r <<"\n";

    file << std::endl;
    file.close();
};

void GreenFuncNph::writeChosenUpdate(std::string filename, int i, double r) const {
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Iteration: " << i << " chosen update: " << r << "\n";
    file << "\n";

    file.close();
};