#include <iostream>
#include <fstream> 
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFunc.hpp"

GreenFunc::GreenFunc(long long int N_diags, double tau_max, double kx, double ky, double kz,
    double chem_potential, int order_max) : _N_diags(N_diags), _tau_max(tau_max), _order_max(returnEven(order_max)), 
    _chem_potential(chem_potential), gen(setSeed()) {        
    // initialize array of possible energies
    _kx = kx;
    _ky = ky;
    _kz = kz;

    // initialize array of all possible phonon vertices
    _vertices = new Vertex[_order_max+2];
    for(int i=0; i<(_order_max+2);i++){
        _vertices[i].tau = 0.;
        _vertices[i].type = 0;
        _vertices[i].linked = -1;
    }
    // initialize array of all possible propagators
    _propagators = new Propagator[_order_max+1];
    for(int i=0;i<(order_max+1);i++){
        _propagators[i].el_propagator_kx = _kx;
        _propagators[i].el_propagator_ky = _ky;
        _propagators[i].el_propagator_kz = _kz;
    }
};

void GreenFunc::setRelaxSteps(int relax_steps){
    while(relax_steps < 0){
        std::cout << "Invalid number of relaxation steps! Number of steps must be >= 0." << std::endl;
        std::cout << "Enter new number of relaxation steps: ";
        std::cin >> relax_steps;
    }
    _relax_steps = relax_steps;
};

void GreenFunc::setAlpha(double alpha){
    while(alpha<=0.){
        std::cout << "Invalid alpha value! Coupling strength must be > 0." << std::endl;
        std::cout << "Enter new alpha value: ";
        std::cin >> alpha;
    }
    _alpha = alpha;
};

void GreenFunc::setVolume(double volume){
    while(volume<=0.){
        std::cout << "Invalid volume value! Volume must be > 0." << std::endl;
        std::cout << "Enter new volume value: ";
        std::cin >> volume;
    }
    _volume = volume;
};

void GreenFunc::setDimension(int D){
    while(D!=2 && D!=3){
        std::cout << "Invalid dimension! Dimension must be 2 or 3." << std::endl;
        std::cout << "Enter new dimension: ";
        std::cin >> D;
    }
    _D = D;
};


void GreenFunc::setN_bins(int N_bins){
    while(N_bins <= 0){
        std::cout << "Invalid number of bins! Number of bins must be > 0." << std::endl;
        std::cout << "Enter new number of bins: ";
        std::cin >> N_bins;
    }
    _N_bins = N_bins;
    _bin_width = _tau_max/_N_bins;
    _bin_center = _bin_width/2;
};

void GreenFunc::setNormConst(double norm_const){
    while(norm_const <= 0){
        std::cout << "Invalid normalization constant! Constant must be > 0." << std::endl;
        std::cout << "Enter new normalization constant: ";
        std::cin >> norm_const;
    }
    _norm_const = norm_const;
};

void GreenFunc::setProbabilities(double p_length, double p_add){
    while(p_length + p_add >= 1){
        std::cout << "Invalid probabilities! Probabilities must be < 1 for p_swap to exist." << std::endl;
        std::cout << "Enter new probabilities: ";
        std::cin >> p_length >> p_add;
    }
    _p_length = p_length;
    _p_add = p_add;
    _p_rem = 1 - p_length - p_add;
};

void GreenFunc::setProbabilities(){
    double p_length, p_add, p_rem;
    std::cout << "Enter probabilities for class I update (length change),class II update (add phonon propagator)" << std::endl;
    std::cout << "probabilities for class II update (remove phonon propagator) will be calculated automatically." << std::endl;
    std::cin >> p_length >> p_add;
    while(p_length + p_add >= 1){
        std::cout << "Invalid probabilities! Probabilities must be < 1 for p_swap to exist." << std::endl;
        std::cout << "Enter new probabilities: ";
        std::cin >> p_length >> p_add >> p_rem;
    }
    _p_length = p_length;
    _p_add = p_add;
    _p_rem = 1 - p_length - p_add;
};

int GreenFunc::findVertexPosition(double tau){
    int position = 0;
    for(int i = 0; i < _current_order+1; i++){
        if(_vertices[i].tau < tau && _vertices[i+1].tau >= tau){
            position = i;
            return position;
        }
    }
    //if(isEqual(tau, _vertices[position].tau)){return -1;}
    return position;
};

int GreenFunc::choosePhononPropagator(){
    std::uniform_int_distribution<> distrib_unif(1,int(_current_order/2)); // chooses one of the phonon propagators at random
    int ph_propagator = distrib_unif(gen);
    int i = 0;
    int counter = 0;
    for(int i = 0; i < _current_order; i++){
        if(_vertices[i].type == +1){
            counter++;
        }
        if(counter == ph_propagator){ // to be fixed
            return i;
        }
    }
    test = test + 1;
    return 0;
};


void GreenFunc::phVertexMakeRoom(int index_one, int index_two){
    int i = 0;
    while(i < _current_order + 1){
        if(_vertices[i].linked > index_two){_vertices[i].linked += 2;}
        else if(_vertices[i].linked > index_one){_vertices[i].linked += 1;}
        i++;
    }
    
    for(int i = _current_order+1; i > index_one; i--){ 
        if(i > index_two){_vertices[i+2] = _vertices[i];}
        else{_vertices[i+1] = _vertices[i];}
  }
};

void GreenFunc::phVertexRemoveRoom(int index_one, int index_two){
    int i = 0;
    while(i < _current_order + 1){
        if(_vertices[i].linked >= index_two){_vertices[i].linked -= 2;}
        else if(_vertices[i].linked >= index_one){_vertices[i].linked -= 1;}
        i++;
    }

    for(int i = index_one; i < _current_order; i++){
        if(i < index_two-1){_vertices[i] = _vertices[i+1];}
        else{_vertices[i] = _vertices[i+2];}
    }
};

void GreenFunc::propagatorArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order; i > index_one-1; i--){
        if(i > index_two - 1){_propagators[i+2] = _propagators[i];}
        else{_propagators[i+1] = _propagators[i];}
    }
    //if(index_one == index_two){_propagators[index_one+1]=_propagators[index_one];} // may be needed
    _propagators[index_one+1] = _propagators[index_one];
    _propagators[index_two+1] = _propagators[index_two+2];
};

void GreenFunc::propagatorArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order; i++){
        if(i < index_two - 2){_propagators[i] = _propagators[i+1];}
        else{_propagators[i] = _propagators[i+2];}
    }
};

double GreenFunc::diagramLengthUpdate(double tau_init){
    double tau_fin = _last_vertex - std::log(1-drawUniformR())/(calcEnergy(_kx,_ky,_kz)-_chem_potential);
    if(tau_fin <= _tau_max){
        _vertices[_current_order+1].tau = tau_fin;
        _vertices[_current_order+1].type = 0;
        return tau_fin;
    }
    else{return tau_init;}
};

void GreenFunc::addPhononPropagator(){
    if(_current_order+1 >= _order_max){return;} // reject if already at max order
    else{
        std::uniform_int_distribution<> distrib_prop(0, _current_order);
        int propagator = distrib_prop(gen);
        double tau_init = _vertices[propagator].tau;
        double tau_end = _vertices[propagator+1].tau;
        std::uniform_real_distribution<> distrib_unif(tau_init, tau_end);
        double tau_one = distrib_unif(gen);
        double tau_two = tau_one - std::log(1-drawUniformR())/1;
        if(tau_two >= _vertices[_current_order+1].tau){return;} // reject if phonon vertex goes out of bound
        else{
            // sampling momentum values for phonon propagators
            std::normal_distribution<> distrib_norm(0,std::sqrt(1/(tau_two-tau_one))); // may need variance, inserted std dev
            double w_x = distrib_norm(gen);
            double w_y = distrib_norm(gen);
            double w_z = distrib_norm(gen);

            // find position of new tau values
            int index_1 = findVertexPosition(tau_one);
            int index_2 = findVertexPosition(tau_two);

            double* px_init = new double[index_2 + 1 - index_1];
            double* py_init = new double[index_2 + 1 - index_1];
            double* pz_init = new double[index_2 + 1 - index_1];
            double* px_fin = new double[index_2 + 1 - index_1];
            double* py_fin = new double[index_2 + 1 - index_1];
            double* pz_fin = new double[index_2 + 1 - index_1];

            for(int i = index_1; i < index_2 +1; i++){
                px_init[i - index_1] = _propagators[i].el_propagator_kx;
                px_fin[i - index_1] = _propagators[i].el_propagator_kx - w_x;
                py_init[i - index_1] = _propagators[i].el_propagator_ky;
                py_fin[i - index_1] = _propagators[i].el_propagator_ky - w_y;
                pz_init[i - index_1] = _propagators[i].el_propagator_kz;
                pz_fin[i - index_1] = _propagators[i].el_propagator_kz - w_z;
            }

            double* energy_init = new double[index_2 + 1 - index_1];
            double* energy_fin = new double[index_2 + 1 - index_1];

            for(int i = 0; i < index_2 - index_1 + 1; i++){
                energy_init[i] = calcEnergy(px_init[i], py_init[i], pz_init[i]);
                energy_fin[i] = calcEnergy(px_fin[i], py_fin[i], pz_fin[i]);
            }

            delete[] px_init;
            delete[] py_init;
            delete[] pz_init;
            delete[] px_fin;
            delete[] py_fin;
            delete[] pz_fin;

            double p_B = _p_rem*(_current_order + 1);
            double p_A = _p_add*(_current_order/2 + 1);

            double exponent_fin = 0.;
            double exponent_init = 0.;

            if(index_1 == index_2){
                exponent_fin = energy_fin[0]*(tau_two-tau_one);
                exponent_init = energy_init[0]*(tau_two-tau_one);
            }
            else{
                exponent_fin = energy_fin[0]*(_vertices[index_1+1].tau-tau_one);
                exponent_init = energy_init[0]*(_vertices[index_1+1].tau-tau_one);
                for(int i = 1; i < index_2 - index_1; i++){
                    exponent_fin += energy_fin[i]*(_vertices[index_1+1+i].tau - _vertices[index_1+i].tau);
                    exponent_init += energy_init[i]*(_vertices[index_1+1+i].tau - _vertices[index_1+i].tau);
                }
                exponent_fin += energy_fin[index_2 - index_1]*(tau_two-_vertices[index_2].tau);
                exponent_init += energy_init[index_2 - index_1]*(tau_two-_vertices[index_2].tau);
            }

            delete[] energy_init;
            delete[] energy_fin;

            double numerator = p_B*std::exp(-(exponent_fin - exponent_init + (1.)*(tau_two - tau_one)))*calcVertexStrength(w_x,w_y,w_z)*(tau_end-tau_init);
            double denominator = p_A*std::pow(2*M_PI,_D)*std::exp(-(tau_two-tau_one))*std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*
            std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_two-tau_one));

            double R_add = numerator/denominator;
            double acceptance_ratio = std::min(1.,R_add);

            if(drawUniformR()>acceptance_ratio){return;}
            else{
                phVertexMakeRoom(index_1, index_2); // make room in vertices array
                propagatorArrayMakeRoom(index_1, index_2); // make room in electron propagators array

                // assign vertex one values
                _vertices[index_1+1].tau = tau_one;
                _vertices[index_1+1].type = +1;
                _vertices[index_1+1].linked = index_2 + 2;
                _vertices[index_1+1].wx = w_x;
                _vertices[index_1+1].wy = w_y;
                _vertices[index_1+1].wz = w_z;
                // assign vertex two values
                _vertices[index_2+2].tau = tau_two;
                _vertices[index_2+2].type = -1;
                _vertices[index_2+2].linked = index_1 + 1;
                _vertices[index_2+2].wx = w_x;
                _vertices[index_2+2].wy = w_y;
                _vertices[index_2+2].wz = w_z;

                // update electron propagator energies
                for(int i=index_1+1; i<index_2+2;i++){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                }
                
                _current_order += 2; // update current diagram order
                findLastPhVertex();
            }
        } 
    }
};

void GreenFunc::removePhononPropagator(){
    if(_current_order < 2){return;} // reject if already at order 0
    else{
        //double tau_end = _vertices[_current_order+1].tau;
        int index_one = choosePhononPropagator();
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

        double p_A = _p_add*((_current_order - 2)/2 + 1);
        double p_B = _p_rem*(_current_order - 1);
        
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
            _current_order -= 2;
            findLastPhVertex();
        }
    }
};

// version with no input uses a lot of RAM due to initialization of huge array, use the other version to only obtain the final histogram
void GreenFunc::markovChain(){
    std::cout << "Starting simulation." << std::endl;

    // initialize data array
    _tau_data = new double[getN()];
    _data_written = true; // for destructor

    std::uniform_real_distribution<> distrib(0, _tau_max);
    double tau_init = distrib(gen);
    tau_init = diagramLengthUpdate(tau_init);
    double r = 0.5;
    long long int i = 0;
    while(i < _relax_steps){
        r = drawUniformR();
        if(r < _p_length){
            tau_init = diagramLengthUpdate(tau_init);
        }
        else if(r < _p_length + _p_add){
            addPhononPropagator();
        }
        else{
            removePhononPropagator();
        }
        i++;
    }
    std::cout << "Relaxation completed." << std::endl;
    i = 0;
    while(i < _N_diags){
        r = drawUniformR();
        if(r < _p_length){
            tau_init = diagramLengthUpdate(tau_init);
            //std::cout << "length update" << std::endl;
        }
        else if(r < _p_length + _p_add){
            addPhononPropagator();
            //std::cout << "add propagator" << std::endl;
        }
        else{
            removePhononPropagator();
            //std::cout << "remove propagator" << std::endl;
        }
        if(_current_order == 0){
            _N0++;
        }
        _tau_data[i] = tau_init;
        i++;
        //std::cout << _current_order << std::endl;
        //std::cout << _last_vertex << std::endl;
    }
    calcHistogram();
    calcNormConst();
    normalizeHistogram();
    std::cout << "Simulation completed." << std::endl;
};

// takes as input the number of generated diagrams, -1 to use internal values (can be set), doesn't use a lot of memory
void GreenFunc::markovChain(int N_diags){
    if(N_diags <= 0){
        N_diags = _N_diags;
    }

    // initialize histogram values and bin count
    _histogram = new double[_N_bins];
    _bin_count = new int[_N_bins];
    for(int i=0; i<_N_bins; i++){
        _histogram[i] = _bin_center + i*_bin_width;
        _bin_count[i] = 0;
    }
    _histogram_calculated = true;

    std::cout << "Starting simulation..." << std::endl;
    std::uniform_real_distribution<> distrib(0, getTauMax());
    double tau_init = distrib(gen);
    tau_init = diagramLengthUpdate(tau_init);
    double r = 0.5;
    long long int i = 0;
    while(i < _relax_steps){
        r = drawUniformR();
        if(r < _p_length){
            tau_init = diagramLengthUpdate(tau_init);
        }
        else if(r < _p_length + _p_add){
            addPhononPropagator();
        }
        else{
            removePhononPropagator();
        }
        i++;
    }
    std::cout << "Relaxation completed." << std::endl;
    double bin_width_inv = _N_bins / (_tau_max);
    i = 0;

    while(i < N_diags){
        r = drawUniformR();
        if(r < _p_length){
            tau_init = diagramLengthUpdate(tau_init);
            //std::cout << "length update" << std::endl;
            test_length = test_length + 1;
        }
        else if(r < _p_length + _p_add){
            addPhononPropagator();
            //std::cout << "add propagator" << std::endl;
            test_add = test_add + 1;
        }
        else{
            removePhononPropagator();
            //std::cout << "remove propagator" << std::endl;
            test_remove = test_remove + 1;
        }
        if(_current_order == 0){
            _N0++;
        }
        //bool flag = false;
        tau_init = (tau_init < 0.) ? 0. : (tau_init >= _tau_max) ? _tau_max - 1e-9 : tau_init;
        int bin = (int)((tau_init - 0.) * bin_width_inv);
        _bin_count[bin]++;
        /*int j = 0;
        while(!flag && j<_N_bins){
            if(tau_init < (_tau_max/2)){
                if(tau_init >= _histogram[j] && tau_init < _histogram[j]+_bin_width){
                    _bin_count[j]++;
                    flag = true;
                }
                j++;
            }
            else{
                if(tau_init <= _histogram[_N_bins-1-j] && tau_init > _histogram[_N_bins-1-j]-_bin_width){
                    _bin_count[_N_bins-1-j]++;
                    flag = true;
                }
                j++;
            }
        }*/
        int control = 0;
        for(int m=0; m<_current_order+2; m++){
            control += _vertices[m].type;
        }
        if(control != 0){
            test_control = test_control + 1;
        }
        //Diagnostic("test.txt", i); // debug method to see diagram, comment it for long simulations
        i++;
    }
    std::cout << "Histogram computed." << std::endl;
    /*std::cout << i << std::endl;
    std::cout << test << std::endl;
    std::cout << "length update: " << test_length << std::endl;
    std::cout << "add update: " << test_add << std::endl;
    std::cout << "remove update: " << test_remove << std::endl;
    std::cout << "Controlling nodes: " << test_control << std::endl;*/
    calcNormConst();
    normalizeHistogram();
    std::cout << "Simulation completed." << std::endl;
};

void GreenFunc::calcHistogram(){
    _histogram = new double[_N_bins];
    _bin_count = new int[_N_bins];

    // initialize histogram values and bin count
    for(int i=0; i<_N_bins; i++){
        _histogram[i] = _bin_center + i*_bin_width;
        _bin_count[i] = 0;
    }

    double bin_width_inv = 1/_bin_width;

    for(int i=0; i<_N_diags; i++){
        _tau_data[i] = (_tau_data[i] < 0.) ? 0. : (_tau_data[i] >= _tau_max) ? _tau_max - 1e-9 : _tau_data[i];
        int bin = (int)((_tau_data[i] - 0.) * bin_width_inv);
        _bin_count[bin]++;
    }
    // fill histogram
    /*for(int i=0; i<_N_diags; i++){
        bool flag = false;
        int j = 0;
        while(!flag && j<_N_bins){
            if(_tau_data[i] >= _histogram[j] && _tau_data[i] < _histogram[j]+_bin_width){
                _bin_count[j]++;
                flag = true;
            }
            j++;
        }
    }*/
    _histogram_calculated = true;
    std::cout << "Histogram computed." << std::endl;
};

void GreenFunc::calcNormConst(){
    double numerator = 1 - std::exp(-(calcEnergy(_kx,_ky,_kz)-_chem_potential)*_tau_max);
    double denominator = calcEnergy(_kx,_ky,_kz)-_chem_potential;
    std::cout << "Normalization constant calculated." << std::endl;
    _norm_const = numerator/denominator;
};

void GreenFunc::normalizeHistogram(){
    _green_func = new double[_N_bins];
    for(int i=0; i<_N_bins; i++){
        _green_func[i] = _bin_count[i]/(_N0*_bin_width);
        _green_func[i] = _green_func[i]*_norm_const;
    }
    _histogram_normalized = true;
    std::cout << "Imaginary time Green's function computed." << std::endl;
};

void GreenFunc::writeTauData(std::string filename) const {
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

    
    for(int i=0; i<_N_diags; i++){
        file << _tau_data[i] << " " <<"\n";
    }
    file.close();
    std::cout << "Data written to file " << filename << "." << std::endl;
};

void GreenFunc::writeHistogram(std::string filename) const {
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

void GreenFunc::Diagnostic(std::string filename, int i) const {
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Iteration: " << i << "\n";

    for(int j=0; j<_current_order+1; j++){
        file << "index: " << j << " time: " << _vertices[j].tau << " wx: " << _vertices[j].wx << " wy: " 
        << _vertices[j].wy << " wz: " << _vertices[j].wz <<" type: " << _vertices[j].type << " linked: " << _vertices[j].linked << "\n";
        file << "propagator: " << j << "          kx: " << _propagators[j].el_propagator_kx << " ky: " << _propagators[j].el_propagator_ky <<
        " kz: " << _propagators[j].el_propagator_kz << "\n";
    }

    file << "index: " << _current_order+1 << " time: " << _vertices[_current_order+1].tau << " wx: " << _vertices[_current_order+1].wx 
    << " wy: " << _vertices[_current_order+1].wy << " wz: " << _vertices[_current_order+1].wz <<" type: " << _vertices[_current_order+1].type 
    << " linked: " << _vertices[_current_order+1].linked << "\n";

    file << std::endl;
    file.close();
};