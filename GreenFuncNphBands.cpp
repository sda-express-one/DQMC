#include "GreenFuncNphBands.hpp"

// constructor
GreenFuncNphBands::GreenFuncNphBands(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
double chem_potential, int order_int_max, int ph_ext_max, int num_bands, int phonon_modes) 
: Diagram(N_diags, tau_max, kx, ky, kz, chem_potential, order_int_max, ph_ext_max) {
    
    // initialize flags
    _flags.gs_energy = false;
    _flags.effective_mass = false;
    _flags.write_diagrams = false;
    _flags.time_benchmark = false;

    // set number of bands and phonon modes
    num_bands = num_bands > 0 ? num_bands : 1; // ensure at least one band
    phonon_modes = phonon_modes > 0 ? phonon_modes : 1; // ensure at least one phonon mode
    _num_bands = num_bands;
    _num_phonon_modes = phonon_modes;

    // initialize arrays
    _phonon_modes = new double[_num_phonon_modes];
    _ext_phonon_type_num = new int[_num_phonon_modes];
    _born_effective_charges = new double[_num_phonon_modes];

    for(int i=0; i < _num_phonon_modes; i++){
        _ext_phonon_type_num[i] = 0;
    }

    // initialize bands
    _bands = new Band[_order_int_max + 2*_ph_ext_max + 1];  // needs further development

    if(_num_bands == 1){
        for(int i = 0; i <_order_int_max + 2*_ph_ext_max + 1; i++){
            _bands[i].band_number = 0;
            _bands[i].c1 = 1;
            _bands[i].c2 = 0;
            _bands[i].c3 = 0;
        }
    }
    else if(_num_bands == 3){
        for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 1; i++){
            _bands[i].band_number = -1; // unassigned
            _bands[i].c1 = 1/3;
            _bands[i].c2 = 1/3;
            _bands[i].c3 = 1/3;
        }
    }
};


// electron bands setters
void GreenFuncNphBands::setLongitudinalEffectiveMass(double mass_long_el){
    _mass_long_el = mass_long_el;
};

void GreenFuncNphBands::setTransversalEffectiveMass(double mass_transv_el){
    _mass_transv_el = mass_transv_el;
};

void GreenFuncNphBands::setLuttingerKohnParameters(double A_LK_el, double B_LK_el, double C_LK_el){
    _A_LK_el = A_LK_el;
    _B_LK_el = B_LK_el;
    _C_LK_el = C_LK_el;
};


// phonon modes setters
void GreenFuncNphBands::setPhononDispersions(double * phonon_modes){
    _phonon_modes = phonon_modes;
};

void GreenFuncNphBands::setBornEffectiveCharges(double * born_effective_charges){
    _born_effective_charges = born_effective_charges;
}

// other constant quantities setters
void GreenFuncNphBands::set1BZVolume(double V_BZ){
    _V_BZ = V_BZ;
};

void GreenFuncNphBands::setBvKVolume(double V_BvK){
    _V_BvK = V_BvK;
};

void GreenFuncNphBands::setDielectricConstant(double dielectric_const){
    _dielectric_const = dielectric_const;
}

int GreenFuncNphBands::findVertexPosition(long double tau){
    int position = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; i++){
        if(_vertices[i].tau < tau && _vertices[i+1].tau >= tau){
            position = i;
            return position;
        }
    }
    return -1; // return -1 if tau is not found in the vertices array
};

int * GreenFuncNphBands::findVerticesPosition(long double tau_one, long double tau_two){
    int* positions = new int[2];
    int counts = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; i++){
        if(_vertices[i].tau < tau_one && _vertices[i+1].tau >= tau_one){
            positions[0] = i;
            counts+=1;
        }
        if(counts > 0 && _vertices[i].tau < tau_two && _vertices[i+1].tau >= tau_two){
            positions[1] = i;
            counts+=1;
        }
    }

    // return [-1,-1] if tau_one or tau_two is not found in the vertices array (or multiple positions are somehow found)
    if(counts != 2){
        positions[0] = -1;
        positions[1] = -1;
        return positions;
    }
    else{
        return positions;
    }
};

int GreenFuncNphBands::chooseInternalPhononPropagator(){
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

int GreenFuncNphBands::chooseExternalPhononPropagator(){
    std::uniform_int_distribution<int> distrib_unif(1, _current_ph_ext); // chooses one of the external phonon propagators at random
    int ph_propagator = distrib_unif(gen);
    int counter = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; i++){
        if(_vertices[i].type == -2){
            counter++;
        }
        if(counter == ph_propagator){
            return i;
        }
    }
    return 0;
};

void GreenFuncNphBands::phVertexMakeRoom(int index_one, int index_two){
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

void GreenFuncNphBands::phVertexRemoveRoom(int index_one, int index_two){
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

void GreenFuncNphBands::propagatorArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order_int + 2*_current_ph_ext; i > index_one-1; i--){
        if(i > index_two - 1){_propagators[i+2] = _propagators[i];}
        else{_propagators[i+1] = _propagators[i];}
    }
    _propagators[index_one+1] = _propagators[index_one];
    _propagators[index_two+1] = _propagators[index_two+2];
};

void GreenFuncNphBands::propagatorArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; i++){
        if(i < index_two - 2){_propagators[i] = _propagators[i+1];}
        else{_propagators[i] = _propagators[i+2];}
    }
};

void GreenFuncNphBands::bandArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order_int + 2*_current_ph_ext; i > index_one-1; i--){
        if(i > index_two - 1){_bands[i+2] = _bands[i];}
        else{_bands[i+1] = _bands[i];}
    }
};

void GreenFuncNphBands::bandArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; i++){
        if(i < index_two - 2){_bands[i] = _bands[i+1];}
        else{_bands[i] = _bands[i+2];}
    }
};

void GreenFuncNphBands::updateExternalPhononTypes(int index){
    _ext_phonon_type_num[index] += 1;
};

long double GreenFuncNphBands::diagramLengthUpdate(long double tau_init){
    // initialize momentum values for last propagator
    int total_order = _current_order_int + 2*_current_ph_ext;
    double kx = _propagators[total_order].el_propagator_kx;
    double ky = _propagators[total_order].el_propagator_ky;
    double kz = _propagators[total_order].el_propagator_kz;

    // generate new time value for last vertex, reject if it goes out of bounds
    long double tau_fin = _last_vertex - std::log(1-drawUniformR())/(electronEnergy(kx,ky,kz,_bands[total_order].effective_mass)
        -_chem_potential + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes));
    if(tau_fin <= _tau_max){
        _vertices[_current_order_int + 2*_current_ph_ext + 1].tau = tau_fin;
        return tau_fin;
    }
    else{return tau_init;}
};

void GreenFuncNphBands::addInternalPhononPropagator(){
    if(_current_order_int+1 >= _order_int_max){return;} // reject if already at max order
    else{
        // choose random electron propagator for new vertex
        std::uniform_int_distribution<int> distrib_prop(0, _current_order_int + 2*_current_ph_ext);
        int propagator = distrib_prop(gen);
        long double tau_init = _vertices[propagator].tau;
        long double tau_end = _vertices[propagator+1].tau;

        // choose phonon index
        std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes);
        int phonon_index = distrib_phon(gen);

        // choose time value of new vertex (between ends of chosen propagator)
        std::uniform_real_distribution<long double> distrib_unif(tau_init, tau_end);
        long double tau_one = distrib_unif(gen);

        // choose time value of second vertex, different energy for different phonon modes
        long double tau_two = tau_one - std::log(1-drawUniformR())/_phonon_modes[phonon_index]; // may be on a different propagator

        if(tau_two >= _vertices[_current_order_int + 2*_current_ph_ext + 1].tau){return;} // reject if phonon vertex goes out of bound
        else{
            // sampling momentum values for phonon propagators
            std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_two-tau_one))); // may need to specify phonon mode
            double w_x = distrib_norm(gen); 
            double w_y = distrib_norm(gen);
            double w_z = distrib_norm(gen);
            
            // find position of new tau values
            int index_one = propagator;
            int index_two = findVertexPosition(tau_two);

            // control statements to check for double precision point errors
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

            // create array of band values
            Band* bands_init = new Band[index_two + 1 - index_one];
            Band* bands_fin = new Band[index_two + 1 - index_one];

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            double c1_new = 1;
            double c2_new = 0;
            double c3_new = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;

            for(int i = index_one; i < index_two + 1; i++){
                px_init[i - index_one] = _propagators[i].el_propagator_kx;
                px_fin[i - index_one] = _propagators[i].el_propagator_kx - w_x;
                py_init[i - index_one] = _propagators[i].el_propagator_ky;
                py_fin[i - index_one] = _propagators[i].el_propagator_ky - w_y;
                pz_init[i - index_one] = _propagators[i].el_propagator_kz;
                pz_fin[i - index_one] = _propagators[i].el_propagator_kz - w_z;

                if(_num_bands == 3){
                    if(i == index_one || i == index_two){
                        bands_init[i - index_one] = _bands[i];

                        std::uniform_int_distribution<int> unif(0, _num_bands-1);
                        chosen_band = unif(gen);
                        new_values_matrix = diagonalizeLKHamiltonian(_propagators[i].el_propagator_kx-w_x, 
                                                                _propagators[i].el_propagator_ky-w_y, 
                                                                _propagators[i].el_propagator_kz-w_z, 
                                                                _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);

                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_fin[i - index_one].effective_mass = computeEffMassfromEigenval(eigenval);
                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_fin[i - index_one].c1 = new_overlap(0); 
                        bands_fin[i - index_one].c2 = new_overlap(1);
                        bands_fin[i - index_one].c3 = new_overlap(2);
                        
                        // compute vertex terms
                        if(i == index_one){
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_init[0], new_overlap);
                            if(index_one == index_two){
                                prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_init[0], new_overlap);
                            }
                        }
                        else{
                            prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[index_two - index_one - 1], bands_init[index_two - index_one]);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_init[index_two - index_one], new_overlap);
                        }
                    }
                    else{
                        bands_init[i - index_one] = _bands[i];

                        chosen_band = _bands[i].band_number;
                        new_values_matrix = diagonalizeLKHamiltonian(_propagators[i].el_propagator_kx-w_x, 
                                                                    _propagators[i].el_propagator_ky-w_y, 
                                                                    _propagators[i].el_propagator_kz-w_z, 
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);
                        
                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_fin[i - index_one].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_fin[i - index_one].c1 = new_overlap(0); 
                        bands_fin[i - index_one].c2 = new_overlap(1); 
                        bands_fin[i - index_one].c3 = new_overlap(2);

                        // compute vertex terms
                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[i-1], bands_init[i]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], bands_fin[i]);
                    }
                }
            }

            // create arrays of energy values
            double* energy_init = new double[index_two + 1 - index_one];
            double* energy_fin = new double[index_two + 1 - index_one];

            for(int i = 0; i < index_two - index_one + 1; i++){
                energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);
                energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);
            }

            delete[] px_init;
            delete[] py_init;
            delete[] pz_init;
            delete[] px_fin;
            delete[] py_fin;
            delete[] pz_fin;

            // initial and final action
            double action_init = 0.;
            double action_fin = 0.;

            if(index_one == index_two){
                action_init = energy_init[0]*(tau_two-tau_one);
                action_fin = energy_fin[0]*(tau_two-tau_one);
            }
            else{
                action_init = energy_init[0]*(_vertices[index_one+1].tau-tau_one);
                action_fin = energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
                for(int i = 1; i < index_two - index_one; i++){
                    action_init += energy_init[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                    action_fin += energy_fin[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                }
                action_init += energy_init[index_two - index_one]*(tau_two-_vertices[index_two].tau);
                action_fin += energy_fin[index_two - index_one]*(tau_two-_vertices[index_two].tau);
            }

            // delete array of energies
            delete[] energy_init;
            delete[] energy_fin;

            double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext + 1);
            double p_A = _p_add_int*(_current_order_int/2 + 1);

            double numerator = p_B*std::exp(-(action_fin - action_init + (phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one))))*
                prefactor_fin*(tau_end-tau_init)*_num_bands*_num_bands*_num_phonon_modes;
            double denominator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)
                *std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_two-tau_one))
                *std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)
                *(tau_two-tau_one));

            double R_add = numerator/denominator;
            
            if(!(Metropolis(R_add))){delete[] bands_init; delete[] bands_fin; return;}
            else{
                phVertexMakeRoom(index_one, index_two); // make room in vertices array
                propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array
                bandArrayMakeRoom(index_one, index_two); // make room in bands array

                // assign vertex one values
                _vertices[index_one+1].tau = tau_one;
                _vertices[index_one+1].type = +1;
                _vertices[index_one+1].linked = index_two + 2;
                _vertices[index_one+1].wx = w_x;
                _vertices[index_one+1].wy = w_y;
                _vertices[index_one+1].wz = w_z;
                _vertices[index_one+1].index = phonon_index;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_two;
                _vertices[index_two+2].type = -1;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;
                _vertices[index_two+2].index = phonon_index;

                // update electron propagator
                for(int i=index_one+1; i<index_two+2;i++){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;

                    _bands[i] = bands_fin[i];
                }

                delete[] bands_init; delete[] bands_fin;
                _current_order_int += 2; // update current diagram internal order
                findLastPhVertex();
            }
        } 
    }
};

void GreenFuncNphBands::removeInternalPhononPropagator(){
    if(_current_order_int < 2){return;} // reject if already at order 0
    else{
        // indexes of initial and final vertices of a random internal phonon propagator
        int index_one = chooseInternalPhononPropagator();
        if(index_one == 0){return;}
        int index_two = _vertices[index_one].linked;

        // vertices' time values
        long double tau_one = _vertices[index_one].tau;
        long double tau_two = _vertices[index_two].tau;

        long double tau_init = _vertices[index_one-1].tau;
        long double tau_end = _tau_max;
        if(index_two != index_one+1){tau_end = _vertices[index_one+1].tau;}
        else{tau_end = _vertices[index_one+2].tau;}

        double w_x = _vertices[index_one].wx;
        double w_y = _vertices[index_one].wy;
        double w_z = _vertices[index_one].wz;

        // phonon mode
        int phonon_index = _vertices[index_one].index;

        // create arrays of momentum values
        double* px_init = new double[index_two - index_one];
        double* py_init = new double[index_two - index_one];
        double* pz_init = new double[index_two - index_one];
        double* px_fin = new double[index_two - index_one];
        double* py_fin = new double[index_two - index_one];
        double* pz_fin = new double[index_two - index_one];

        Band* bands_init = new Band[index_two - index_one];
        Band* bands_fin = new Band[index_two - index_one];

        // vertices weights of two diagrams
        double prefactor_init = 1;
        double prefactor_fin = 1;

        // temporary variables for new proposed diagram
        int chosen_band = 0;
        double c1_new = 1;
        double c2_new = 0;
        double c3_new = 0;
        double eigenval = 1.0;
        Eigen::Matrix<double,4,3> new_values_matrix;
        Eigen::Vector3d new_overlap;

        for(int i = index_one; i < index_two; i++){
            px_init[i - index_one] = _propagators[i].el_propagator_kx + w_x;
            px_fin[i - index_one] = _propagators[i].el_propagator_kx;
            py_init[i - index_one] = _propagators[i].el_propagator_ky + w_y;
            py_fin[i - index_one] = _propagators[i].el_propagator_ky;
            pz_init[i - index_one] = _propagators[i].el_propagator_kz + w_z;
            pz_fin[i - index_one] = _propagators[i].el_propagator_kz;

            if(_num_bands == 3){
                if(i == index_one || i == index_two){
                    bands_fin[i - index_one] = _bands[i];

                    new_values_matrix = diagonalizeLKHamiltonian(_propagators[i].el_propagator_kx + w_x, 
                                                            _propagators[i].el_propagator_ky + w_y, 
                                                            _propagators[i].el_propagator_kz + w_z, 
                                                            _A_LK_el, _B_LK_el, _C_LK_el);

                    // check if it is a free propagator
                    if(new_values_matrix(0,0) == -2){
                        bands_init[i - index_one].band_number = -1;
                        bands_init[i - index_one].effective_mass = 1;
                        bands_init[i - index_one].c1 = 1/3;
                        bands_init[i - index_one].c2 = 1/3;
                        bands_init[i - index_one].c3 = 1/3;

                        new_overlap << (1/3), (1/3), (1/3);
                    }
                    else{
                        eigenval = new_values_matrix(0,chosen_band);

                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_init[i - index_one].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_init[i - index_one].c1 = new_overlap(0); 
                        bands_init[i - index_one].c2 = new_overlap(1);
                        bands_init[i - index_one].c3 = new_overlap(2);
                    }

                    // compute vertex terms
                    if(index_one + 1 == index_two){
                        prefactor_init = prefactor_init*vertexOverlapTerm(_bands[index_one-1], new_overlap)
                            *vertexOverlapTerm(_bands[index_one + 1], new_overlap);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands[index_one-1], bands_fin[0])
                            *vertexOverlapTerm(bands_fin[0], _bands[index_one+1]);
                    }
                    else if(i == index_one){
                        prefactor_init = prefactor_init*vertexOverlapTerm(_bands[index_one-1], new_overlap);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands[index_one-1], bands_fin[0]);
                    }
                    else{
                        prefactor_init = prefactor_init*vertexOverlapTerm(_bands[index_two+1], new_overlap);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[index_two], _bands[index_two+1]);
                    }
                }
                else{
                    bands_fin[i - index_one] = _bands[i];

                    chosen_band = _bands[i].band_number;
                    new_values_matrix = diagonalizeLKHamiltonian(_propagators[i].el_propagator_kx-w_x, 
                                                                _propagators[i].el_propagator_ky-w_y, 
                                                                _propagators[i].el_propagator_kz-w_z, 
                                                                _A_LK_el, _B_LK_el, _C_LK_el);
                    
                    // check if it is a free propagator
                    if(new_values_matrix(0,0) == -2){
                        bands_init[i - index_one].band_number = -1;
                        bands_init[i - index_one].effective_mass = 1;
                        bands_init[i - index_one].c1 = 1/3;
                        bands_init[i - index_one].c2 = 1/3;
                        bands_init[i - index_one].c3 = 1/3;

                        new_overlap << (1/3), (1/3), (1/3);
                    }
                    else{                        
                        eigenval = new_values_matrix(0,chosen_band);

                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_init[i - index_one].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_fin[i - index_one].c1 = new_overlap(0); 
                        bands_fin[i - index_one].c2 = new_overlap(1); 
                        bands_fin[i - index_one].c3 = new_overlap(2);
                    }
                        
                    // compute vertex terms
                    prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[i-1], bands_init[i]);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], bands_fin[i]);
                }
            }
        }

        // create arrays of energy values
        double* energy_init = new double[index_two - index_one];
        double* energy_fin = new double[index_two - index_one];

        for(int i = 0; i < index_two - index_one; i++){
            energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);
            energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);
        }

        delete[] px_init;
        delete[] px_fin;
        delete[] py_init;
        delete[] py_fin;
        delete[] pz_init;
        delete[] pz_fin;
        
        // initial and final action
        double action_fin = 0;
        double action_init = 0;
        
        if(index_one + 1 == index_two){
            action_init = energy_init[0]*(tau_two-tau_one);
            action_fin = energy_fin[0]*(tau_two-tau_one);
        }
        else{
            action_init = energy_init[0]*(_vertices[index_one+1].tau-tau_one);
            action_fin = energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
            for(int i = 1; i < index_two - index_one - 1; i++){
                action_init += energy_init[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
                action_fin += energy_fin[i]*(_vertices[index_one+1+i].tau - _vertices[index_one+i].tau);
            }
            action_init += energy_init[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
            action_fin += energy_fin[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
        }

        delete[] energy_init;
        delete[] energy_fin;

        double p_A = _p_add_int*((_current_order_int - 2)/2 + 1);
        double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext - 1);
        
        double numerator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)
                *std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_two-tau_one))
                *std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)
                *(tau_two-tau_one));
        double denominator = p_B*std::exp(-(action_fin - action_init + (phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one))))*
                prefactor_fin*(tau_end-tau_init)*_num_bands*_num_bands*_num_phonon_modes;

        double R_rem = numerator/denominator;

        if(!(Metropolis(R_rem))){delete[] bands_init; delete[] bands_fin; return;}
        else{
            for(int i=index_one; i<index_two;i++){
                _propagators[i].el_propagator_kx += w_x;
                _propagators[i].el_propagator_ky += w_y;
                _propagators[i].el_propagator_kz += w_z;
                _bands[i] = bands_init[i];
            }

            phVertexRemoveRoom(index_one, index_two);
            propagatorArrayRemoveRoom(index_one, index_two);
            bandArrayRemoveRoom(index_one, index_two);

            delete[] bands_init; delete[] bands_fin;

            _current_order_int -= 2;
            findLastPhVertex();
        }
    }
};


void GreenFuncNphBands::swapPhononPropagator(){
    if(_current_order_int < 4){return;} // swap not possible if internal order is less than 4
    else{
        std::uniform_int_distribution<int> distrib(1, _current_order_int + 2*_current_ph_ext - 1); // choose random internal propagator
        int index_one = distrib(gen); // choose random internal propagatorS
        
        if(_vertices[index_one].type != +1 && _vertices[index_one].type != -1){return;} // reject if vertex does not belong to internal phonon propagator
        if(_vertices[index_one+1].type != +1 && _vertices[index_one+1].type != -1){return;} // reject if vertex does not belong to internal phonon propagator
        if(_vertices[index_one].linked == index_one + 1 || _vertices[index_one].linked == -1){return;} // reject if the two vertices are linked

        int index_two = index_one +1;

        // get values of first vertex
        int v1 = _vertices[index_one].type;
        long double wx1 = _vertices[index_one].wx;
        long double wy1 = _vertices[index_one].wy;
        long double wz1 = _vertices[index_one].wz;
        long double tau1 = _vertices[index_one].tau;
        int phonon_index1 = _vertices[index_one].index;

        // get values of second vertex
        int v2 = _vertices[index_two].type;
        long double wx2 = _vertices[index_two].wx;
        long double wy2 = _vertices[index_two].wy;
        long double wz2 = _vertices[index_two].wz;
        long double tau2 = _vertices[index_two].tau;
        int phonon_index2 = _vertices[index_two].index;

        // get momentum of propagator
        long double kx = _propagators[index_one].el_propagator_kx;
        long double ky = _propagators[index_one].el_propagator_ky;
        long double kz = _propagators[index_one].el_propagator_kz;

        // retrieve value of initial electron effective mass and propose new value of final effective mass
        double eff_mass_el_initial, eff_mass_el_final;
        eff_mass_el_initial = _bands[index_one].effective_mass;

        // new proposed band for propagator
        int chosen_band = 0;
        
        // new proposed eigenfunction
        Eigen::Vector3d new_overlap;
        new_overlap << 1, 0 ,0;

        // vertex terms that go into R_swap evaluation
        double prefactor_num = 1;
        double prefactor_den = 1;
        
        if(_num_bands == 3){
            Eigen::Matrix<double,4,3> new_values_matrix = diagonalizeLKHamiltonian(kx+v1*wx1-v2*wx2, ky+v1*wy1-v2*wy2, kz+v1*wz1-v2*wz2, _A_LK_el, _B_LK_el, _C_LK_el);

            // choose at random one of the bands
            std::uniform_int_distribution<int> distrib_unif(0,2);
            chosen_band = distrib_unif(gen);

            double eigenval = new_values_matrix(0,chosen_band);
            eff_mass_el_final = computeEffMassfromEigenval(eigenval); // computing new proposed electron effective mass from chosen eigenvalue
            new_overlap = new_values_matrix.block<3,1>(1,chosen_band); // new proposed band eigenstate

            // compute only overlap term for vertices that go into R_swap evaluation,
            // the strength term does not change
            prefactor_num = vertexOverlapTerm(_bands[index_one-1], new_overlap)*vertexOverlapTerm(_bands[index_two], new_overlap);
            prefactor_den = vertexOverlapTerm(_bands[index_one-1], _bands[index_one])*vertexOverlapTerm(_bands[index_two], _bands[index_one]);
        }

        // compute energies
        double energy_final_el = electronEnergy(kx+v1*wx1-v2*wx2, ky+v1*wy1-v2*wy2, kz+v1*wz1-v2*wz2, eff_mass_el_final);
        double energy_initial_el = electronEnergy(kx, ky, kz, eff_mass_el_initial);
        double energy_phonons = phononEnergy(_phonon_modes, phonon_index1)*v1-phononEnergy(_phonon_modes, phonon_index2)*v2;
        
        // need to address negative values issue
        // compute transition probability
        double R_swap = (prefactor_num/prefactor_den)*std::exp(-(energy_final_el-energy_initial_el-energy_phonons)*(tau2-tau1))*_num_bands; // to be fixed
        double acceptance_ratio = std::min(1.,R_swap);

        if(drawUniformR() > acceptance_ratio){return;}
        else{
            // assign new momentum values to propagator
            _propagators[index_one].el_propagator_kx += v1*wx1-v2*wx2;
            _propagators[index_one].el_propagator_ky += v1*wy1-v2*wy2;
            _propagators[index_one].el_propagator_kz += v1*wz1-v2*wz2;
            
            // assign new band values
            _bands[index_one].band_number = chosen_band;
            _bands[index_one].effective_mass = eff_mass_el_final;
            _bands[index_one].c1 = new_overlap[0];
            _bands[index_one].c2 = new_overlap[1];
            _bands[index_one].c3 = new_overlap[2];

            int linked1 = _vertices[index_one].linked;
            int linked2 = _vertices[index_two].linked;

            // assing new links to conjugate vertices
            _vertices[linked1].linked = index_two;         
            _vertices[linked2].linked = index_one;

            _vertices[index_one].wx = wx2;
            _vertices[index_one].wy = wy2;
            _vertices[index_one].wz = wz2;
            _vertices[index_one].type = v2;
            _vertices[index_one].linked = linked2;
            _vertices[index_one].tau = tau1;
            _vertices[index_one].index = phonon_index2;
                
            _vertices[index_two].wx = wx1;
            _vertices[index_two].wy = wy1;
            _vertices[index_two].wz = wz1;
            _vertices[index_two].type = v1;
            _vertices[index_two].linked = linked1;
            _vertices[index_two].tau = tau2;
            _vertices[index_two].index = phonon_index1;
        }
    }
};

void GreenFuncNphBands::shiftPhononPropagator(){
    if(_current_order_int + _current_ph_ext == 0){return;} // reject if no vertices are present
    else{
        int total_order = _current_order_int + 2*_current_ph_ext;
        std::uniform_int_distribution<int> distrib(1, total_order);
        int vertex_index = distrib(gen); // choose random vertex

        // necessary step to address phonon type into evaluation
        int c = _vertices[vertex_index].type;
        if(c == 2){c = 1;}
        else if(c == -2){c = -1;}
        int phonon_index = _vertices[vertex_index].index;  // index of phonon mode in el-phonon vertex

        long double tau_init = _vertices[vertex_index - 1].tau;
        long double tau_fin = _vertices[vertex_index + 1].tau;

        // incoming electron momentum
        long double kx_incoming = _propagators[vertex_index - 1].el_propagator_kx;
        long double ky_incoming = _propagators[vertex_index - 1].el_propagator_ky;
        long double kz_incoming = _propagators[vertex_index - 1].el_propagator_kz;
        double el_eff_mass_incoming = _bands[vertex_index - 1].effective_mass;

        // outgoing electron momentum
        long double kx_outgoing = _propagators[vertex_index].el_propagator_kx;
        long double ky_outgoing = _propagators[vertex_index].el_propagator_ky;
        long double kz_outgoing = _propagators[vertex_index].el_propagator_kz;
        double el_eff_mass_outgoing = _bands[vertex_index].effective_mass;

        double energy_delta = electronEnergy(kx_incoming, ky_incoming, kz_incoming, el_eff_mass_incoming) 
            - electronEnergy(kx_outgoing, ky_outgoing, kz_outgoing, el_eff_mass_outgoing) - phononEnergy(_phonon_modes, phonon_index)*c;
        
        long double tau_new = tau_init - std::log(1 - drawUniformR()*(1 - std::exp(-energy_delta*(tau_fin - tau_init))))/energy_delta;

        if(isEqual(tau_new, tau_init) || isEqual(tau_new, tau_fin) || tau_new > tau_fin){return;} // check for possible double precision errors
        
        _vertices[vertex_index].tau = tau_new; // assign new time value to vertex
    }
};

long double GreenFuncNphBands::stretchDiagramLength(long double tau_init){
    // initialize momentum values for first electron propagator
    long double kx = _propagators[0].el_propagator_kx;
    long double ky = _propagators[0].el_propagator_ky;
    long double kz = _propagators[0].el_propagator_kz;

    // create complete array of new vertex time values
    int total_order = _current_order_int + 2*_current_ph_ext;
    long double* new_taus = new long double[total_order+2];
    new_taus[0] = 0; // first vertex time value is always 0

    // declare variables used in for loop
    int c = 0;
    int phonon_index = 0;
    long double tau_one, tau_two;
    long double kx_incoming, ky_incoming, kz_incoming;
    long double kx_outgoing, ky_outgoing, kz_outgoing;
    double el_eff_mass_incoming, el_eff_mass_outgoing;
    double energy_delta = 0;

    // assign new time values to every vertex
    for(int i = 1; i < total_order+1; i++){
        // necessary step to address phonon type into evaluation
        c = _vertices[i].type;
        if(c == 2){c = 1;}
        else if(c == -2){c = -1;}
        phonon_index = _vertices[i].index;  // index of phonon mode in el-phonon vertex

        // left vertex value is retrieved from new proposed values, right vertex value from old ones (vertex array)
        tau_one = new_taus[i-1];
        tau_two = _vertices[i+1].tau;

        kx_incoming = _propagators[i-1].el_propagator_kx;
        ky_incoming = _propagators[i-1].el_propagator_ky;
        kz_incoming = _propagators[i-1].el_propagator_kz;
        el_eff_mass_incoming = _bands[i-1].effective_mass;

        kx_outgoing = _propagators[i].el_propagator_kx;
        ky_outgoing = _propagators[i].el_propagator_ky;
        kz_outgoing = _propagators[i].el_propagator_kz;
        el_eff_mass_outgoing = _bands[i].effective_mass;

        energy_delta = electronEnergy(kx_incoming, ky_incoming, kz_incoming, el_eff_mass_incoming) - 
        electronEnergy(kx_outgoing, ky_outgoing, kz_outgoing, el_eff_mass_outgoing) - phononEnergy(_phonon_modes, phonon_index)*c;
            
        new_taus[i] = tau_one - std::log(1 - drawUniformR()*(1 - std::exp(-energy_delta*(tau_two - tau_one))))/energy_delta;

        // check for possible double precision errors
        if(isEqual(new_taus[i], tau_one) || isEqual(new_taus[i], tau_two) || new_taus[i] > tau_two){delete[] new_taus; return tau_init;}
    }

    new_taus[total_order+1] = new_taus[total_order] - std::log(1-drawUniformR())/(electronEnergy(kx,ky,kz,_bands[total_order].effective_mass) 
        - _chem_potential + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes));
    
    // check for possible double precision errors
    if(isEqual(new_taus[total_order], new_taus[total_order+1]) || isEqual(new_taus[total_order+1], _tau_max) 
       || new_taus[total_order+1] >= _tau_max){delete[] new_taus; return tau_init;}
    else{
        for(int i = 0; i < total_order + 2; i++){
            _vertices[i].tau = new_taus[i]; // assign new time values to vertices
        }
    }
    long double tau_fin = new_taus[total_order+1]; // assign new time value to last vertex
    delete[] new_taus;
    return tau_fin; // return new length of diagram
};