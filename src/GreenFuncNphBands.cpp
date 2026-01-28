#include "../include/GreenFuncNphBands.hpp"

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
    _dielectric_responses = new double[_num_phonon_modes];


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
            _bands[i].c1 = (1./3);
            _bands[i].c2 = (1./3);
            _bands[i].c3 = (1./3);
        }
    }
};

GreenFuncNphBands::GreenFuncNphBands(Propagator * propagators, Vertex * vertices, Band * bands,
unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
double chem_potential, int order_int_max, int ph_ext_max, int num_bands, int phonon_modes) 
: Diagram(propagators, vertices, N_diags, tau_max, kx, ky, kz, chem_potential, order_int_max, ph_ext_max) {
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
    _dielectric_responses = new double[_num_phonon_modes];


    for(int i=0; i < _num_phonon_modes; i++){
        _ext_phonon_type_num[i] = 0;
    }

    // initialize bands
    _bands = new Band[_order_int_max + 2*_ph_ext_max + 1];

    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 1; i++){
        _bands[i] = bands[i];
    }
};


void GreenFuncNphBands::getEffectiveMasses(long double * effective_masses) const {
    effective_masses[0] = _effective_masses[0];
    effective_masses[1] = _effective_masses[1];
    effective_masses[2] = _effective_masses[2];
};

void GreenFuncNphBands::getGFExactPoints(long double * points, long double * gf_values) const {
    for(int i=0; i<_num_points; i++){
        points[i] = _points[i];
        gf_values[i] = _points_gf_exact[i];
    }
};

void GreenFuncNphBands::getHistogram(long double * histogram, long double * green_func) const {
    for(int i=0; i<_N_bins; i++){
        histogram[i] = _histogram[i];
        green_func[i] = _green_func[i];
    }
};

// electron bands setters
void GreenFuncNphBands::setEffectiveMasses(double m_x, double m_y, double m_z){
    _m_x_el = m_x;
    _m_y_el = m_y;
    _m_z_el = m_z;
};

void GreenFuncNphBands::setLuttingerKohnParameters(double A_LK_el, double B_LK_el, double C_LK_el){
    _A_LK_el = A_LK_el;
    _B_LK_el = B_LK_el;
    _C_LK_el = C_LK_el;
};

// phonon modes setters
void GreenFuncNphBands::setPhononModes(double * phonon_modes){
    for(int i = 0; i < _num_phonon_modes; i++){
        _phonon_modes[i] = phonon_modes[i];
    }
};

void GreenFuncNphBands::setDielectricResponses(double * dielectric_responses){
    for(int i = 0; i < _num_phonon_modes; i++){
        _dielectric_responses[i] = dielectric_responses[i];
    }
}

// other constant quantities setters
void GreenFuncNphBands::set1BZVolume(double V_BZ){_V_BZ = V_BZ;};

void GreenFuncNphBands::setBvKVolume(double V_BvK){_V_BvK = V_BvK;};

void GreenFuncNphBands::setDielectricConstant(double dielectric_const){_dielectric_const = dielectric_const;};

void GreenFuncNphBands::setCurrentOrderInt(int order_int){_current_order_int = (order_int >= 0 ? order_int : 0);};

void GreenFuncNphBands::setCurrentPhExt(int ph_ext){_current_ph_ext = (ph_ext >= 0 ? ph_ext : 0);};

// MC updates probability setter
void GreenFuncNphBands::setProbabilities(double* probs){
    if(!isEqual(probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7], 1)){
        std::cerr << "Invalid probabilities, total probability must add to 1.\n";
        double normalization = 1/(probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7]);
        std::cerr << "Probabilities are being riscaled using the value " << normalization <<".\n";
        for(int i = 0; i < 8; i++){
            probs[i] = probs[i]*normalization;
        }
        std::cerr << "New probabilities are: " << probs[0] << " " << probs[1] << " " << probs[2] << " " 
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

// parallelization settings
void GreenFuncNphBands::setMaster(bool master_mode){_master = master_mode;};

// calculations setters
void GreenFuncNphBands::setCalculations(bool gf_exact, bool histo, bool gs_energy, bool effective_mass, bool Z_factor, bool fix_tau_value){
    _flags.gf_exact = gf_exact;
    _flags.histo = histo;
    _flags.gs_energy = gs_energy;
    _flags.effective_mass = effective_mass;
    _flags.Z_factor = Z_factor;
    _flags.fix_tau_value = fix_tau_value;
};

// histogram setters
void GreenFuncNphBands::setN_bins(int N_bins){
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

// exact GF setters
void GreenFuncNphBands::setNumPoints(int num_points){
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

void GreenFuncNphBands::setSelectedOrder(int selected_order){
    if(selected_order > _ph_ext_max){
        std::cerr << "Invalid order! Order must be <= " << _ph_ext_max << "." << std::endl;
        std::cerr << "Selected order for exact GF calculation is set to 0." << std::endl;
        selected_order = 0;
    }
    _selected_order = selected_order;
};

// exact energy estimator setters
void GreenFuncNphBands::setTauCutoffEnergy(long double tau_cutoff_energy){
    while(tau_cutoff_energy <= 0 || tau_cutoff_energy >= _tau_max){
        std::cout << "Invalid energy cutoff! Cutoff must be > 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "Enter new energy cutoff: ";
        std::cin >> tau_cutoff_energy;
        std::cout << "\n";
    }
    _tau_cutoff_energy = tau_cutoff_energy;
}

// exact mass estimator setters
void GreenFuncNphBands::setTauCutoffMass(long double tau_cutoff_mass){
    while(tau_cutoff_mass <= 0 || tau_cutoff_mass >= _tau_max){
        std::cout << "Invalid mass cutoff! Cutoff must be > 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "Enter new mass cutoff: ";
        std::cin >> tau_cutoff_mass;
        std::cout << "\n";
    }
    _tau_cutoff_mass = tau_cutoff_mass;
}

// MC statistics setter
void GreenFuncNphBands::setTauCutoffStatistics(long double tau_cutoff_statistics){
    while(tau_cutoff_statistics < 0 || tau_cutoff_statistics >= _tau_max){
        std::cout << "Invalid statistics cutoff! Cutoff must be >= 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "The default value is set to 0.\n";
        tau_cutoff_statistics = 0;
    }
    _tau_cutoff_statistics = tau_cutoff_statistics;
};

int GreenFuncNphBands::findVertexPosition(long double tau){
    int position = 0;
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; ++i){
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
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; ++i){
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
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext; ++i){
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
    for(int i = 0; i < _current_order_int + 2*_current_ph_ext + 1; ++i){
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
        ++i;
    }
    
    for(int i = _current_order_int + 2*_current_ph_ext + 1; i > index_one; --i){ 
        if(i > index_two){_vertices[i+2] = _vertices[i];}
        else{_vertices[i+1] = _vertices[i];}
  }
};

void GreenFuncNphBands::phVertexRemoveRoom(int index_one, int index_two){
    int i = 0;
    while(i < _current_order_int + 2*_current_ph_ext + 1){
        if(_vertices[i].linked >= index_two){_vertices[i].linked -= 2;}
        else if(_vertices[i].linked >= index_one){_vertices[i].linked -= 1;}
        ++i;
    }

    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; ++i){
        if(i < index_two-1){_vertices[i] = _vertices[i+1];}
        else{_vertices[i] = _vertices[i+2];}
    }
};

void GreenFuncNphBands::propagatorArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order_int + 2*_current_ph_ext; i > index_one-1; --i){
        if(i > index_two - 1){_propagators[i+2] = _propagators[i];}
        else{_propagators[i+1] = _propagators[i];}
    }
    _propagators[index_one+1] = _propagators[index_one];
    _propagators[index_two+1] = _propagators[index_two+2];
};

void GreenFuncNphBands::propagatorArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; ++i){
        if(i < index_two - 2){_propagators[i] = _propagators[i+1];}
        else{_propagators[i] = _propagators[i+2];}
    }
};

void GreenFuncNphBands::bandArrayMakeRoom(int index_one, int index_two){
    for(int i = _current_order_int + 2*_current_ph_ext; i > index_one-1; --i){
        if(i > index_two - 1){_bands[i+2] = _bands[i];}
        else{_bands[i+1] = _bands[i];}
    }
};

void GreenFuncNphBands::bandArrayRemoveRoom(int index_one, int index_two){
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext; ++i){
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
        std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes-1);
        int phonon_index = distrib_phon(gen);

        // choose time value of new vertex (between ends of chosen propagator)
        std::uniform_real_distribution<long double> distrib_unif(tau_init, tau_end);
        long double tau_one = distrib_unif(gen);

        // choose time value of second vertex, different energy for different phonon modes
        long double tau_two = tau_one - std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // may be on a different propagator

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

            // create arrays of energy values
            double* energy_init = new double[index_two + 1 - index_one];
            double* energy_fin = new double[index_two + 1 - index_one];

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            //double c1_new = 1;
            //double c2_new = 0;
            //double c3_new = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> unif(0, _num_bands-1);

            // initial and final action
            double action_init = 0.;
            double action_fin = 0.;

            int j = 0;
            
            for(int i = index_one; i < index_two + 1; ++i){
                j = i - index_one;

                // initial diagram momenta
                px_init[j] = _propagators[i].el_propagator_kx;
                py_init[j] = _propagators[i].el_propagator_ky;
                pz_init[j] = _propagators[i].el_propagator_kz;

                // proposed diagram momenta
                px_fin[j] = _propagators[i].el_propagator_kx - w_x;
                py_fin[j] = _propagators[i].el_propagator_ky - w_y;
                pz_fin[j] = _propagators[i].el_propagator_kz - w_z;

                bands_init[j] = _bands[i];

                if(_num_bands == 3){
                    chosen_band = unif(gen);
                    bands_fin[j].band_number = chosen_band;

                    new_values_matrix = diagonalizeLKHamiltonian(px_fin[j], py_fin[j], pz_fin[j], _A_LK_el, _B_LK_el, _C_LK_el);
                    eigenval = new_values_matrix(0,chosen_band);

                    // computing new proposed electron effective mass from chosen eigenvalue
                    bands_fin[j].effective_mass = computeEffMassfromEigenval(eigenval);
                    // new proposed band eigenstate
                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                    bands_fin[j].c1 = new_overlap(0); 
                    bands_fin[j].c2 = new_overlap(1);
                    bands_fin[j].c3 = new_overlap(2);
                        
                    // compute vertex terms
                    if(i == index_one){
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_init[0], new_overlap);
                    }
                    else{
                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[j-1], bands_init[j]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[j-1], bands_fin[j]);                        
                    }
                    if (i == index_two){
                        // prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[index_two - index_one - 1], bands_init[index_two - index_one]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_init[index_two - index_one], new_overlap);
                    }
                }
                else if(_num_bands == 1){
                    bands_fin[j].effective_mass = computeEffMassSingleBand(px_fin[j], py_fin[j], pz_fin[j],
                                                                        _m_x_el, _m_y_el, _m_z_el);
                }

                energy_init[j] = electronEnergy(px_init[j], py_init[j], pz_init[j], bands_init[j].effective_mass);
                energy_fin[j] = electronEnergy(px_fin[j], py_fin[j], pz_fin[j], bands_fin[j].effective_mass);

                if(index_one == index_two){
                    action_init += energy_init[0]*(tau_two-tau_one);
                    action_fin += energy_fin[0]*(tau_two-tau_one);
                }
                else if(i == index_one){
                    action_init += energy_init[0]*(_vertices[index_one+1].tau-tau_one);
                    action_fin += energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
                }
                else if(i == index_two){
                    action_init += energy_init[index_two - index_one]*(tau_two-_vertices[index_two].tau);
                    action_fin += energy_fin[index_two - index_one]*(tau_two-_vertices[index_two].tau);
                }
                else{
                    action_init += energy_init[j]*(_vertices[i+1].tau - _vertices[i].tau);
                    action_fin += energy_fin[j]*(_vertices[i+1].tau - _vertices[i].tau);
                }
            }

            /*for(int i = 0; i < index_two - index_one + 1; i++){
                energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);
                energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);
            }*/

            /*if(index_one == index_two){
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
            }*/
            
            delete[] px_init;
            delete[] py_init;
            delete[] pz_init;
            delete[] px_fin;
            delete[] py_fin;
            delete[] pz_fin;

            // delete array of energies
            delete[] energy_init;
            delete[] energy_fin;

            double mass_q = 1;
            if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}

            // multiply prefactor final by strength term of 2 extrema
            prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                         _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                         *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                         _dielectric_responses[phonon_index], _dielectric_const, mass_q);

            double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext + 1);
            double p_A = _p_add_int*(_current_order_int/2 + 1);

            double numerator = p_B*std::exp(-(action_fin - action_init + (phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one))))*
                (tau_end-tau_init)*prefactor_fin*_V_BZ; //*_num_bands*_num_bands
            double denominator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)
                *std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_two-tau_one))
                *std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)
                *(tau_two-tau_one))*prefactor_init;

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

                    _bands[i] = bands_fin[i-(index_one+1)];
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

        // create arrays of energy values
        double* energy_init = new double[index_two - index_one];
        double* energy_fin = new double[index_two - index_one];

        // initial and final action
        double action_fin = 0;
        double action_init = 0;

        // vertices weights of two diagrams
        double prefactor_init = 1;
        double prefactor_fin = 1;

        // temporary variables for new proposed diagram
        int chosen_band = 0;
        double eigenval = 1.0;
        Eigen::Matrix<double,4,3> new_values_matrix;
        Eigen::Vector3d new_overlap;
        std::uniform_int_distribution<int> band_number(0, _num_bands - 1);

        int j = 0;
        for(int i = index_one; i < index_two; ++i){
            j = i - index_one;

            // proposed new momentum values
            px_init[j] = _propagators[i].el_propagator_kx + w_x;
            py_init[j] = _propagators[i].el_propagator_ky + w_y;
            pz_init[j] = _propagators[i].el_propagator_kz + w_z;

            // current momentum values
            px_fin[j] = _propagators[i].el_propagator_kx;
            py_fin[j] = _propagators[i].el_propagator_ky;
            pz_fin[j] = _propagators[i].el_propagator_kz;

            bands_fin[j] = _bands[i];

            if(_num_bands == 3){
                /*if(i == index_one || i == index_two){
                    chosen_band = band_number(gen);
                    new_values_matrix = diagonalizeLKHamiltonian(_propagators[i].el_propagator_kx + w_x, 
                                                            _propagators[i].el_propagator_ky + w_y, 
                                                            _propagators[i].el_propagator_kz + w_z, 
                                                            _A_LK_el, _B_LK_el, _C_LK_el);

                    // check if it is a free propagator
                    if(new_values_matrix(0,0) == -2){
                        bands_init[j].band_number = -1;
                        bands_init[j].effective_mass = 1;
                        bands_init[j].c1 = (1./3);
                        bands_init[j].c2 = (1./3);
                        bands_init[j].c3 = (1./3);

                        new_overlap << (1./3), (1./3), (1./3);
                    }
                    else{
                        bands_init[j].band_number = chosen_band;
                        eigenval = new_values_matrix(0,chosen_band);

                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_init[j].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_init[j].c1 = new_overlap(0); 
                        bands_init[j].c2 = new_overlap(1);
                        bands_init[j].c3 = new_overlap(2);
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
                }*/
                //else{
                chosen_band = band_number(gen);
                new_values_matrix = diagonalizeLKHamiltonian(px_init[j], py_init[j], pz_init[j], _A_LK_el, _B_LK_el, _C_LK_el);
                
                // check if it is a free propagator
                if(new_values_matrix(0,0) == -2){
                    bands_init[j].band_number = -1;
                    bands_init[j].effective_mass = 1;
                    bands_init[j].c1 = (1./3);
                    bands_init[j].c2 = (1./3);
                    bands_init[j].c3 = (1./3);

                    new_overlap << (1./3), (1./3), (1./3);
                }
                else{                        
                    eigenval = new_values_matrix(0,chosen_band);

                    // computing new proposed electron effective mass from chosen eigenvalue
                    bands_init[j].effective_mass = computeEffMassfromEigenval(eigenval);

                    // new proposed band eigenstate
                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                    bands_init[j].c1 = new_overlap(0); 
                    bands_init[j].c2 = new_overlap(1); 
                    bands_init[j].c3 = new_overlap(2);
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
                else if(i == index_two - 1){
                    prefactor_init = prefactor_init*vertexOverlapTerm(_bands[index_two], new_overlap);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[j-1], _bands[index_two]);
                }
                else{
                    prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[j-1], bands_init[j]);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[j-1], bands_fin[j]);
                }
            }
            //}
            else if(_num_bands == 1){
                bands_init[j].effective_mass = computeEffMassSingleBand(px_init[j], py_init[j], pz_init[j], 
                                                                    _m_x_el, _m_y_el, _m_z_el);
            }
            energy_init[j] = electronEnergy(px_init[j], py_init[j], pz_init[j], bands_init[j].effective_mass);
            energy_fin[j] = electronEnergy(px_fin[j], py_fin[j], pz_fin[j], bands_fin[j].effective_mass);

            if(index_one + 1 == index_two){
                action_init += energy_init[0]*(tau_two-tau_one);
                action_fin += energy_fin[0]*(tau_two-tau_one);
            }
            else if(i == index_one){
                action_init += energy_init[0]*(_vertices[index_one+1].tau-tau_one);
                action_fin += energy_fin[0]*(_vertices[index_one+1].tau-tau_one);
            }
            else if(i == index_two - 1){
                action_init += energy_init[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
                action_fin += energy_fin[index_two - index_one - 1]*(tau_two-_vertices[index_two-1].tau);
            }
            else{
                action_init += energy_init[j]*(_vertices[i+1].tau - _vertices[i].tau);
                action_fin += energy_fin[j]*(_vertices[i+1].tau - _vertices[i].tau);
            }
        }

        /*for(int i = 0; i < index_two - index_one; i++){
            energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);
            energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);
        }*/


        /*if(index_one + 1 == index_two){
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
        }*/

        delete[] px_init; delete[] py_init; delete[] pz_init;
        delete[] px_fin; delete[] py_fin; delete[] pz_fin;

        delete[] energy_init; delete[] energy_fin;

        double mass_q = 1;
        if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}

        // multiply prefactor final by strength term of 2 extrema
        prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                    _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                    *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                    _dielectric_responses[phonon_index], _dielectric_const, mass_q);

        double p_A = _p_add_int*((_current_order_int - 2)/2 + 1);
        double p_B = _p_rem_int*(_current_order_int + 2*_current_ph_ext - 1);
        
        double numerator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)
                *std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_two-tau_one))
                *std::pow(((tau_two-tau_one)/(2*M_PI)),double(_D)/2.)*std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)
                *(tau_two-tau_one))*prefactor_init;
        double denominator = p_B*std::exp(-(action_fin - action_init + (phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one))))*
                (tau_end-tau_init)*prefactor_fin*_V_BZ; //*_num_bands*_num_bands

        double R_rem = numerator/denominator;

        if(!(Metropolis(R_rem))){delete[] bands_init; delete[] bands_fin; return;}
        else{
            for(int i=index_one; i<index_two; ++i){
                _propagators[i].el_propagator_kx += w_x;
                _propagators[i].el_propagator_ky += w_y;
                _propagators[i].el_propagator_kz += w_z;
                _bands[i] = bands_init[i-index_one];
            }

            phVertexRemoveRoom(index_one, index_two);
            propagatorArrayRemoveRoom(index_one, index_two);
            bandArrayRemoveRoom(index_one, index_two);

            int total_order = _current_order_int + 2*_current_ph_ext;
            
            _propagators[total_order-1].el_propagator_kx = 0;
            _propagators[total_order-1].el_propagator_ky = 0;
            _propagators[total_order-1].el_propagator_kz = 0;
            _bands[total_order-1].effective_mass = 1.0;
                    
            _propagators[total_order].el_propagator_kx = 0;
            _propagators[total_order].el_propagator_ky = 0;
            _propagators[total_order].el_propagator_kz = 0;
            _bands[total_order].effective_mass = 1.0;
            
            if(_num_bands == 3){
                _bands[total_order-1].band_number = -1;
                _bands[total_order-1].c1 = (1./3);
                _bands[total_order-1].c2 = (1./3);
                _bands[total_order-1].c3 = (1./3);

                _bands[total_order].band_number = -1;
                _bands[total_order].c1 = (1./3);
                _bands[total_order].c2 = (1./3);
                _bands[total_order].c3 = (1./3);
            }

            delete[] bands_init; delete[] bands_fin;

            _current_order_int -= 2;
            findLastPhVertex();
        }
    }
};

void GreenFuncNphBands::addExternalPhononPropagator(){
    if(_current_ph_ext >= _ph_ext_max){return;} // return if already at max number of ext phonon propagators
    else{
        int total_order = _current_order_int + 2*_current_ph_ext;
        long double tau_current = _vertices[total_order + 1].tau; // length of current diagram

        // choose phonon index
        std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes-1);
        int phonon_index = distrib_phon(gen);
        //int phonon_index = choosePhonon();

        // time of ingoing vertex of ext phonon propagator
        long double tau_one = 0 - std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // time of ingoing vertex of ext phonon propagator
        if(isEqual(tau_one, tau_current) || tau_one >= tau_current){return;} // reject if it goes out of bound

        long double tau_two = tau_current + std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // time of outgoing vertex

        if(isEqual(tau_two,0) || tau_two <= 0){return;} // reject if it goes out of bound

        if(isEqual(tau_one, tau_two)){return;} // reject if both vertices are equal (should not happen)
        
        // sampling momentum values for phonon propagators
        std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_current-tau_two+tau_one)));
        double w_x = distrib_norm(gen);
        double w_y = distrib_norm(gen);
        double w_z = distrib_norm(gen);

        if(tau_one <= tau_two){
            int * indexes = findVerticesPosition(tau_one, tau_two);
            int index_one = indexes[0];
            int index_two = indexes[1];
            delete[] indexes;
            //int index_one = findVertexPosition(tau_one);
            //int index_two = findVertexPosition(tau_two);

            // control statements to check for floating point errors
            if(index_one == -1 || index_two == -1){return;} // reject if tau values are not found in the vertices array
            if(tau_one < _vertices[index_one].tau || isEqual(tau_one, _vertices[index_one].tau) 
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

            Band* bands_one_init = new Band[index_one + 1];
            Band* bands_two_init = new Band[total_order + 1 - index_two];

            Band* bands_one_fin = new Band[index_one + 1];
            Band* bands_two_fin = new Band[total_order + 1 - index_two];

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            //double c1_new = 1;
            //double c2_new = 0;
            //double c3_new = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> band_number(0, _num_bands-1);

            // energy arrays
            double* energy_one_init = new double[index_one + 1];
            double* energy_two_init = new double[total_order + 1 - index_two];
            double* energy_one_fin = new double[index_one + 1];
            double* energy_two_fin = new double[total_order + 1 - index_two];

            // initial and final action
            double action_one_init = 0.; 
            double action_two_init = 0.;
            double action_one_fin = 0.;
            double action_two_fin = 0.;

            
            for(int i = 0; i < index_one + 1; ++i){
                // retrieve momentum values for propagators below first ph vertex
                px_one_init[i] = _propagators[i].el_propagator_kx;
                px_one_fin[i] = _propagators[i].el_propagator_kx - w_x;

                py_one_init[i] = _propagators[i].el_propagator_ky;
                py_one_fin[i] = _propagators[i].el_propagator_ky - w_y;

                pz_one_init[i] = _propagators[i].el_propagator_kz;
                pz_one_fin[i] = _propagators[i].el_propagator_kz - w_z;

                bands_one_init[i] = _bands[i];

                if(_num_bands == 3){
                    chosen_band = band_number(gen);
                    bands_one_fin[i].band_number = chosen_band;

                    new_values_matrix = diagonalizeLKHamiltonian(px_one_fin[i], py_one_fin[i], pz_one_fin[i], _A_LK_el, _B_LK_el, _C_LK_el);    
                    eigenval = new_values_matrix(0,chosen_band);
                
                    // computing new proposed electron effective mass from chosen eigenvalue
                    bands_one_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                    // new proposed band eigenstate
                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                    bands_one_fin[i].c1 = new_overlap(0); 
                    bands_one_fin[i].c2 = new_overlap(1); 
                    bands_one_fin[i].c3 = new_overlap(2);

                    if(i != 0){
                    // compute vertex terms
                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_one_init[i-1], bands_one_init[i]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_one_fin[i-1], new_overlap);
                    }
                }
                else if(_num_bands == 1){
                    bands_one_fin[i].effective_mass = computeEffMassSingleBand(px_one_fin[i], py_one_fin[i], pz_one_fin[i],
                                                                                _m_x_el, _m_y_el, _m_z_el);   
                }

                // compute energies
                energy_one_init[i] = electronEnergy(px_one_init[i], py_one_init[i], pz_one_init[i], bands_one_init[i].effective_mass);
                energy_one_fin[i] = electronEnergy(px_one_fin[i], py_one_fin[i], pz_one_fin[i], bands_one_fin[i].effective_mass);

                // compute action
                if(i != index_one){
                    action_one_init += energy_one_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                    action_one_fin += energy_one_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
                }
                else{
                    action_one_init += energy_one_init[index_one]*(tau_one - _vertices[index_one].tau);
                    action_one_fin += energy_one_fin[index_one]*(tau_one - _vertices[index_one].tau);
                }
            }
            
            // new vertex (left)
            if(_num_bands == 3){
                prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_one_fin[index_one], _bands[index_one]);
            }
            
            int j = 0;
            
            for(int i = index_two; i < total_order + 1; ++i){
                j = i-index_two;

                // retrieve momentum values for propagators above second ph vertex
                px_two_init[j] = _propagators[i].el_propagator_kx;
                px_two_fin[j] = _propagators[i].el_propagator_kx - w_x;

                py_two_init[j] = _propagators[i].el_propagator_ky;
                py_two_fin[j] = _propagators[i].el_propagator_ky - w_y;

                pz_two_init[j] = _propagators[i].el_propagator_kz;
                pz_two_fin[j] = _propagators[i].el_propagator_kz - w_z;

                bands_two_init[j] = _bands[i];

                if(_num_bands == 3){
                    if(i != total_order){
                        chosen_band = band_number(gen);
                        bands_two_fin[j].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_two_fin[j], py_two_fin[j], pz_two_fin[j], 
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);

                        // computing new proposed electron effective mass from chosen eigenvalue
                        bands_two_fin[j].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_two_fin[j].c1 = new_overlap(0); 
                        bands_two_fin[j].c2 = new_overlap(1); 
                        bands_two_fin[j].c3 = new_overlap(2);
                    }
                    else{
                        // last propagator must be the same as the first one (conservation of four-momentum)
                        bands_two_fin[total_order - index_two] = bands_one_fin[0];
                        new_overlap(0) = bands_one_fin[0].c1;
                        new_overlap(1) = bands_one_fin[0].c2;
                        new_overlap(2) = bands_one_fin[0].c3;
                    }

                    // compute vertex terms
                    if(i != index_two){
                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_two_init[j-1], bands_two_init[j]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_two_fin[j-1], new_overlap);
                    }
                    else{
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_two_init[0], new_overlap);
                    }
                }
                else if(_num_bands == 1){
                    bands_two_fin[j].effective_mass = computeEffMassSingleBand(px_two_fin[j], py_two_fin[j], pz_two_fin[j], 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                }

                // calc energy values for propagators above second ph vertex
                energy_two_init[j] = electronEnergy(px_two_init[j], py_two_init[j], pz_two_init[j], bands_two_init[j].effective_mass);
                energy_two_fin[j] = electronEnergy(px_two_fin[j], py_two_fin[j], pz_two_fin[j], bands_two_fin[j].effective_mass);
                
                if(i == index_two){
                    action_two_init += energy_two_init[0]*(_vertices[index_two+1].tau - tau_two);
                    action_two_fin += energy_two_fin[0]*(_vertices[index_two+1].tau - tau_two);
                }
                else{
                    action_two_init += energy_two_init[j]*(_vertices[i+1].tau - _vertices[i].tau);
                    action_two_fin += energy_two_fin[j]*(_vertices[i+1].tau - _vertices[i].tau);
                }
            }

            // delete momentum arrays
            delete[] px_one_init; delete[] px_one_fin; 
            delete[] py_one_init; delete[] py_one_fin; 
            delete[] pz_one_init; delete[] pz_one_fin;
            delete[] px_two_init; delete[] px_two_fin; 
            delete[] py_two_init; delete[] py_two_fin; 
            delete[] pz_two_init; delete[] pz_two_fin;

            // delete array of energies
            delete[] energy_one_init;
            delete[] energy_one_fin;
            delete[] energy_two_init;
            delete[] energy_two_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double mass_q = 1;
            if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}

            // multiply prefactor final by strength term of 2 extrema
            prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                        *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q);

            double numerator = p_B*std::exp(-(action_two_fin + action_one_fin - action_two_init - action_one_init + 
                phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two+tau_one)))*prefactor_fin*_V_BZ; // *_num_bands*_num_bands

            double denominator = p_A*std::pow(2*M_PI,_D)
                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*tau_one)
                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two))
                *std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)
                *std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one))*prefactor_init;

            double R_add = numerator/denominator;

            if(!(Metropolis(R_add))){
                delete[] bands_one_init; delete[] bands_one_fin;
                delete[] bands_two_init; delete[] bands_two_fin;
                return;
            }
            else{
                phVertexMakeRoom(index_one, index_two); // make room in vertices array
                propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array
                bandArrayMakeRoom(index_one, index_two); // make room in bands array

                // assign vertex one values
                _vertices[index_one+1].tau = tau_one;
                _vertices[index_one+1].type = -2;
                _vertices[index_one+1].linked = index_two + 2;
                _vertices[index_one+1].index = phonon_index;
                _vertices[index_one+1].wx = w_x;
                _vertices[index_one+1].wy = w_y;
                _vertices[index_one+1].wz = w_z;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_two;
                _vertices[index_two+2].type = +2;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].index = phonon_index;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;

                // update electron propagator energies
                for(int i = 0; i < index_one + 1; ++i){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                    _bands[i] = bands_one_fin[i];
                }
                for(int i = index_two+2; i < total_order + 3; ++i){
                    _propagators[i].el_propagator_kx -= w_x;
                    _propagators[i].el_propagator_ky -= w_y;
                    _propagators[i].el_propagator_kz -= w_z;
                    _bands[i] = bands_two_fin[i - index_two - 2];
                }

                delete[] bands_one_init; delete[] bands_one_fin;
                delete[] bands_two_init; delete[] bands_two_fin;

                _current_ph_ext += 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
        else{
            int * indexes = findVerticesPosition(tau_two, tau_one);
            int index_one = indexes[0];
            int index_two = indexes[1];
            delete[] indexes;

            //int index_one = findVertexPosition(tau_two);
            //int index_two = findVertexPosition(tau_one);
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

            Band* bands_init = new Band[total_order+1];
            Band* bands_fin = new Band[total_order+3];

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            //double c1_new = 1;
            //double c2_new = 0;
            //double c3_new = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> band_number(0, _num_bands-1);

            double* energy_init = new double[total_order+1];
            double* energy_fin = new double[total_order+3];

            double action_init = 0.;
            double action_fin = 0.;

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            /*for(int i = 0; i < total_order + 1; i++){
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
            }*/

            for(int i = 0; i < total_order + 3; ++i){

                // initial diagram
                if(i < total_order + 1){
                    px_init[i] = _propagators[i].el_propagator_kx;
                    py_init[i] = _propagators[i].el_propagator_ky;
                    pz_init[i] = _propagators[i].el_propagator_kz;

                    bands_init[i] = _bands[i];

                    if(_num_bands == 3 && i > 0){    
                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[i-1], bands_init[i]);
                    }

                    energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);

                    action_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                }

                // final diagram
                if(i < index_one + 1){
                    px_fin[i] = _propagators[i].el_propagator_kx - w_x;
                    py_fin[i] = _propagators[i].el_propagator_ky - w_y;
                    pz_fin[i] = _propagators[i].el_propagator_kz - w_z;

                    if(_num_bands == 3){
                        chosen_band = band_number(gen);
                        bands_fin[i].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_fin[i], py_fin[i], pz_fin[i],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);
                        bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_fin[i].c1 = new_overlap(0); 
                        bands_fin[i].c2 = new_overlap(1); 
                        bands_fin[i].c3 = new_overlap(2);

                        if(i > 0){
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], new_overlap);
                        }
                    }
                    else if(_num_bands == 1){
                        bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin[i], py_fin[i], pz_fin[i],
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);

                    if(i < index_one){
                        action_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
                    }
                    else{
                        action_fin += energy_fin[index_one]*(tau_two - _vertices[index_one].tau);
                    }
                }
                else if(i < index_two + 2){
                    px_fin[i] = _propagators[i-1].el_propagator_kx - 2*w_x;
                    py_fin[i] = _propagators[i-1].el_propagator_ky - 2*w_y;
                    pz_fin[i] = _propagators[i-1].el_propagator_kz - 2*w_z;

                    if(_num_bands == 3){
                        chosen_band = band_number(gen);
                        bands_fin[i].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_fin[i], py_fin[i], pz_fin[i],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);
                        bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        bands_fin[i].c1 = new_overlap(0); 
                        bands_fin[i].c2 = new_overlap(1); 
                        bands_fin[i].c3 = new_overlap(2);

                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], new_overlap);
                    }
                    else if(_num_bands == 1){
                        bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin[i], py_fin[i], pz_fin[i], 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);

                    if(i == index_one + 1){
                        action_fin += energy_fin[index_one+1]*(_vertices[index_one+1].tau - tau_two);
                    }
                    else{
                        action_fin += energy_fin[i]*(_vertices[i].tau - _vertices[i-1].tau);
                    }
                }
                else{
                    px_fin[i] = _propagators[i-2].el_propagator_kx - w_x;
                    py_fin[i] = _propagators[i-2].el_propagator_ky - w_y;
                    pz_fin[i] = _propagators[i-2].el_propagator_kz - w_z;

                    if(_num_bands == 3){
                        if(i != total_order + 2){
                            chosen_band = band_number(gen);
                            bands_fin[i].band_number = chosen_band;

                            new_values_matrix = diagonalizeLKHamiltonian(px_fin[i], py_fin[i], pz_fin[i],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            eigenval = new_values_matrix(0,chosen_band);
                            bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                            new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                            bands_fin[i].c1 = new_overlap(0); 
                            bands_fin[i].c2 = new_overlap(1); 
                            bands_fin[i].c3 = new_overlap(2);
                        }
                        else{
                            bands_fin[total_order + 2] = bands_fin[0];
                        }
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], new_overlap);
                    }
                    else if(_num_bands == 1){
                        bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin[i], py_fin[i], pz_fin[i], 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);

                    if(i == index_two + 2){
                        action_fin += energy_fin[index_two+2]*(_vertices[index_two+1].tau - tau_one);
                    }
                    else{
                        action_fin += energy_fin[i]*(_vertices[i-1].tau - _vertices[i-2].tau);
                    }
                }
            }

            /*for(int i = 0; i < total_order + 1; i++){
                energy_init[i] = electronDispersion(px_init[i], py_init[i], pz_init[i], _el_eff_mass);
            }*/
            /*for(int i = 0; i < total_order + 3; i++){
                energy_fin[i] = electronDispersion(px_fin[i], py_fin[i], pz_fin[i], _el_eff_mass);
            }*/
            
            /*for(int i = 0; i < total_order + 1; i++){
                action_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }*/

            /*for(int i = 0; i < index_one; i++){
                action_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }
            action_fin += energy_fin[index_one]*(tau_two - _vertices[index_one].tau);
            action_fin += energy_fin[index_one+1]*(_vertices[index_one+1].tau - tau_two);
            for(int i = index_one + 2; i < index_two + 1; i++){
                action_fin += energy_fin[i]*(_vertices[i].tau - _vertices[i-1].tau);
            }
            action_fin += energy_fin[index_two+1]*(tau_one - _vertices[index_two].tau);
            action_fin += energy_fin[index_two+2]*(_vertices[index_two+1].tau - tau_one);
            for(int i = index_two + 3; i < total_order + 3; i++){
                action_fin += energy_fin[i]*(_vertices[i-1].tau - _vertices[i-2].tau);
            }*/

            delete[] px_init; delete[] py_init; delete[] pz_init;
            
            // delete array of energies
            delete[] energy_init;
            delete[] energy_fin;
            
            double p_B = _p_rem_ext;
            double p_A = _p_add_ext*(_current_ph_ext+1);

            double mass_q = 1;
            if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}

            // multiply prefactor final by strength term of 2 extrema
            prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                        *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q);

            double numerator = p_B*std::exp(-(action_fin - action_init + phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two+tau_one)))
                                *prefactor_fin*_V_BZ; // *_num_bands*_num_bands

            double denominator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*tau_one)
                                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two))
                                *std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)
                                *std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one))*prefactor_init;

            double R_add = numerator/denominator;

            if(!(Metropolis(R_add))){
                delete[] px_fin; delete[] py_fin; delete[] pz_fin; 
                delete[] bands_init; delete[] bands_fin;
                return;
            }
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
                _vertices[index_one+1].index = phonon_index;

                // assign vertex two values
                _vertices[index_two+2].tau = tau_one;
                _vertices[index_two+2].type = -2;
                _vertices[index_two+2].linked = index_one + 1;
                _vertices[index_two+2].wx = w_x;
                _vertices[index_two+2].wy = w_y;
                _vertices[index_two+2].wz = w_z;
                _vertices[index_two+2].index = phonon_index;

                // update electron propagator energies
                for(int i = 0; i < total_order + 3; ++i){
                    _propagators[i].el_propagator_kx = px_fin[i];
                    _propagators[i].el_propagator_ky = py_fin[i];
                    _propagators[i].el_propagator_kz = pz_fin[i];
                    _bands[i] = bands_fin[i];
                }

                delete[] px_fin; delete[] py_fin; delete[] pz_fin;
                delete[] bands_init; delete[] bands_fin;

                _current_ph_ext += 1; // update current number of external phonons
                findLastPhVertex();
                return;
            }
        }
    }
};

void GreenFuncNphBands::removeExternalPhononPropagator(){
    if(_current_ph_ext <= 0){return;} // reject if already at order 0
    else{
        // indexes of initial and final vertices of a random internal phonon propagator
        int index_one = chooseExternalPhononPropagator();
        if(index_one == 0){return;}
        int index_two = _vertices[index_one].linked;

        if(index_one >= index_two){
            int temp = index_one;
            index_one = index_two;
            index_two = temp;
        }   
        int total_order = _current_order_int + 2*_current_ph_ext;

        // retrieve time values of the two vertices
        long double tau_one = _vertices[index_one].tau;
        long double tau_two = _vertices[index_two].tau;

        // retrieve phonon momentum
        double w_x = _vertices[index_one].wx;
        double w_y = _vertices[index_one].wy;
        double w_z = _vertices[index_one].wz;

        // retrieve phonon mode
        int phonon_index = _vertices[index_one].index;


        // vertices weights of two diagrams
        double prefactor_fin = 1;
        double prefactor_init = 1;

        // temporary variables for new proposed diagram
        int chosen_band = 0;
        //double c1_new = 1;
        //double c2_new = 0;
        //double c3_new = 0;
        double eigenval = 1.0;
        Eigen::Matrix<double,4,3> new_values_matrix;
        Eigen::Vector3d new_overlap;
        new_overlap << 1,0,0;
        std::uniform_int_distribution<int> band_number(0, _num_bands-1);

        if(_vertices[index_one].type == -2){

            // momentum arrays
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

            // band arrays
            Band* bands_one_init = new Band[index_one + 1];
            Band* bands_one_fin = new Band[index_one + 1];

            Band* bands_two_init = new Band[total_order - index_two + 1];
            Band* bands_two_fin = new Band[total_order - index_two + 1];

            // energy arrays
            double* energy_one_init = new double[index_one + 1];
            double* energy_one_fin = new double[index_one + 1];
            double* energy_two_init = new double[total_order - index_two + 1];
            double* energy_two_fin = new double[total_order - index_two + 1];

            // initial and final action
            double action_one_init = 0.;
            double action_one_fin = 0.;
            double action_two_init = 0.;
            double action_two_fin = 0.;

            for(int i = 0; i < index_one + 1; ++i){
                px_one_init[i] = _propagators[i].el_propagator_kx + w_x;
                py_one_init[i] = _propagators[i].el_propagator_ky + w_y;
                pz_one_init[i] = _propagators[i].el_propagator_kz + w_z;

                px_one_fin[i] = _propagators[i].el_propagator_kx;
                py_one_fin[i] = _propagators[i].el_propagator_ky;
                pz_one_fin[i] = _propagators[i].el_propagator_kz;

                bands_one_fin[i] = _bands[i];

                if(_num_bands == 3){
                    if(i != index_one){
                        chosen_band = band_number(gen);

                        new_values_matrix = diagonalizeLKHamiltonian(px_one_init[i], py_one_init[i], pz_one_init[i],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        // check if it is a free propagator
                        if(new_values_matrix(0,0) == -2){
                            bands_one_init[i].band_number = -1;
                            bands_one_init[i].effective_mass = 1;
                            bands_one_init[i].c1 = (1./3);
                            bands_one_init[i].c2 = (1./3);
                            bands_one_init[i].c3 = (1./3);

                            new_overlap << (1./3), (1./3), (1./3);
                        }
                        else{
                            bands_one_init[i].band_number = chosen_band;

                            eigenval = new_values_matrix(0,chosen_band);
                            bands_one_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                            new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                            bands_one_init[i].c1 = new_overlap(0); 
                            bands_one_init[i].c2 = new_overlap(1); 
                            bands_one_init[i].c3 = new_overlap(2);
                        }
                        if(i != 0){
                            prefactor_init = prefactor_init*vertexOverlapTerm(bands_one_init[i-1], new_overlap);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_one_fin[i-1], bands_one_fin[i]);
                        }
                    }
                    else{
                        bands_one_init[index_one] = _bands[index_one + 1];
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_one_fin[index_one-1], bands_one_fin[index_one]); 
                    }
                }
                else if(_num_bands == 1){
                    bands_one_init[i].effective_mass = computeEffMassSingleBand(px_one_init[i], py_one_init[i], pz_one_init[i], 
                                                                                _m_x_el, _m_y_el, _m_z_el);
                }

                // compute energies
                energy_one_init[i] = electronEnergy(px_one_init[i], py_one_init[i], pz_one_init[i], bands_one_init[i].effective_mass);
                energy_one_fin[i] = electronEnergy(px_one_fin[i], py_one_fin[i], pz_one_fin[i], bands_one_fin[i].effective_mass);

                // compute actions
                action_one_init += energy_one_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                action_one_fin += energy_one_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            int j = 0;
            for(int i = index_two; i < total_order + 1; ++i){
                j = i - index_two;
                px_two_init[j] = _propagators[i].el_propagator_kx + w_x;
                py_two_init[j] = _propagators[i].el_propagator_ky + w_y;
                pz_two_init[j] = _propagators[i].el_propagator_kz + w_z;
                px_two_fin[j] = _propagators[i].el_propagator_kx;
                py_two_fin[j] = _propagators[i].el_propagator_ky;
                pz_two_fin[j] = _propagators[i].el_propagator_kz;

                bands_two_fin[j] = _bands[i];

                if(_num_bands == 3){
                    if(i != index_two && i != total_order){
                        chosen_band = band_number(gen);

                        new_values_matrix = diagonalizeLKHamiltonian(px_two_init[j], py_two_init[j], pz_two_init[j],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        if(new_values_matrix(0,0) == -2){
                            bands_two_init[j].band_number = -1;
                            bands_two_init[j].effective_mass = 1;
                            bands_two_init[j].c1 = (1./3);
                            bands_two_init[j].c2 = (1./3);
                            bands_two_init[j].c3 = (1./3);

                            new_overlap << (1./3), (1./3), (1./3);
                        }
                        else{
                            bands_two_init[j].band_number = chosen_band;
                            
                            eigenval = new_values_matrix(0,chosen_band);
                            bands_two_init[j].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                            new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                            bands_two_init[j].c1 = new_overlap(0); 
                            bands_two_init[j].c2 = new_overlap(1); 
                            bands_two_init[j].c3 = new_overlap(2);
                        }

                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_two_init[j-1], new_overlap);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_two_fin[j-1], bands_two_fin[j]);
                    }
                    else if( i == index_two){
                        bands_two_init[0] = _bands[index_two-1];
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands[index_two-1], bands_two_fin[0]);
                    }
                    else{
                        bands_two_init[j] = bands_one_init[0];

                        prefactor_init = prefactor_init*vertexOverlapTerm(bands_two_init[j], bands_two_init[j-1]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_two_fin[j], bands_two_fin[j-1]);
                    }
                }
                else if(_num_bands == 1){
                    bands_two_init[j].effective_mass = computeEffMassSingleBand(px_two_init[j], py_two_init[j], pz_two_init[j],
                                                                                _m_x_el, _m_y_el, _m_z_el);
                }

                energy_two_init[j] = electronEnergy(px_two_init[j], py_two_init[j], pz_two_init[j], bands_two_init[j].effective_mass);
                energy_two_fin[j] = electronEnergy(px_two_fin[j], py_two_fin[j], pz_two_fin[j], bands_two_fin[j].effective_mass);

                action_two_init += energy_two_init[j]*(_vertices[i+1].tau - _vertices[i].tau);
                action_two_fin += energy_two_fin[j]*(_vertices[i+1].tau - _vertices[i].tau);
            }


            /*for(int i = 0; i < index_one + 1; i++){
                energy_one_init[i] = electronDispersion(px_one_init[i], py_one_init[i], pz_one_init[i], _el_eff_mass);
                energy_one_fin[i] = electronDispersion(px_one_fin[i], py_one_fin[i], pz_one_fin[i], _el_eff_mass);
            }*/
            
            /*for(int i = 0; i < total_order - index_two + 1; i++){   
                energy_two_init[i] = electronDispersion(px_two_init[i], py_two_init[i], pz_two_init[i], _el_eff_mass);
                energy_two_fin[i] = electronDispersion(px_two_fin[i], py_two_fin[i], pz_two_fin[i], _el_eff_mass);
            }*/

            /*for(int i = 0; i < index_one + 1; i++){
                action_one_init += energy_one_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                action_one_fin += energy_one_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
            }

            for(int i = index_two; i < total_order + 1; i++){
                action_two_init += energy_two_init[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
                action_two_fin += energy_two_fin[i-index_two]*(_vertices[i+1].tau - _vertices[i].tau);
            }*/

            delete[] px_one_init; delete[] py_one_init; delete[] pz_one_init;
            delete[] px_one_fin; delete[] py_one_fin; delete[] pz_one_fin;
            delete[] px_two_init; delete[] py_two_init; delete[] pz_two_init;
            delete[] px_two_fin; delete[] py_two_fin; delete[] pz_two_fin;

            delete[] energy_one_init; delete[] energy_one_fin;
            delete[] energy_two_init; delete[] energy_two_fin;

            long double tau_current = _vertices[total_order+1].tau; // length of current diagram

            double p_A = _p_add_ext*_current_ph_ext;
            double p_B = _p_rem_ext;

            double mass_q = 1;
            if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}

            // multiply prefactor final by strength term of 2 extrema
            prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                        *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q);

            double numerator = p_A*std::pow(2*M_PI,_D)
                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*tau_one)
                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two))
                *std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)
                *std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one))*prefactor_init;

            double denominator = p_B*std::exp(-(action_two_fin + action_one_fin - action_two_init - action_one_init + 
                phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two+tau_one)))*prefactor_fin*_V_BZ; // *_num_bands*_num_bands
            
            double R_rem = numerator/denominator;

            if(!Metropolis(R_rem)){
                delete[] bands_one_init; delete[] bands_one_fin;
                delete[] bands_two_init; delete[] bands_two_fin;
                return;
            }
            else{
                for(int i=0; i<index_one;++i){
                    _propagators[i].el_propagator_kx += w_x;
                    _propagators[i].el_propagator_ky += w_y;
                    _propagators[i].el_propagator_kz += w_z;
                    _bands[i] = bands_one_init[i];
                }
                for(int i=index_two; i<total_order+1;++i){
                    _propagators[i].el_propagator_kx += w_x;
                    _propagators[i].el_propagator_ky += w_y;
                    _propagators[i].el_propagator_kz += w_z;
                    _bands[i] = bands_two_init[i-index_two];
                }

                phVertexRemoveRoom(index_one, index_two); // remove room in vertices array
                propagatorArrayRemoveRoom(index_one, index_two); // remove room in electron propagators array
                bandArrayRemoveRoom(index_one, index_two); // remove room in bands array

                _propagators[total_order-1].el_propagator_kx = 0;
                _propagators[total_order-1].el_propagator_ky = 0;
                _propagators[total_order-1].el_propagator_kz = 0;
                _bands[total_order-1].effective_mass = 1.0;         

                _propagators[total_order].el_propagator_kx = 0;
                _propagators[total_order].el_propagator_ky = 0;
                _propagators[total_order].el_propagator_kz = 0;
                _bands[total_order].effective_mass = 1.0;
                
                if(_num_bands == 3){
                    _bands[total_order-1].band_number = -1;
                    _bands[total_order-1].c1 = (1./3);
                    _bands[total_order-1].c2 = (1./3);
                    _bands[total_order-1].c3 = (1./3);
                    _bands[total_order].band_number = -1;
                    _bands[total_order].c1 = (1./3);
                    _bands[total_order].c2 = (1./3);
                    _bands[total_order].c3 = (1./3);
                }

                delete[] bands_one_init; delete[] bands_one_fin;
                delete[] bands_two_init; delete[] bands_two_fin;

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

                Band* bands_init = new Band[total_order + 1];
                Band* bands_fin = new Band[total_order + 1];
                
                double* energy_init = new double[total_order + 1];
                double* energy_fin = new double[total_order + 1];

                double action_init = 0.;
                double action_fin = 0.;

                for(int i = 0; i < total_order + 1; ++i){
                    px_fin[i] = _propagators[i].el_propagator_kx;
                    py_fin[i] = _propagators[i].el_propagator_ky;
                    pz_fin[i] = _propagators[i].el_propagator_kz;

                    bands_fin[i] = _bands[i];

                    if(i < index_one || i >= index_two){
                        px_init[i] = _propagators[i].el_propagator_kx + w_x;
                        py_init[i] = _propagators[i].el_propagator_ky + w_y;
                        pz_init[i] = _propagators[i].el_propagator_kz + w_z;

                        if(_num_bands == 3){
                            chosen_band = band_number(gen);

                            new_values_matrix = diagonalizeLKHamiltonian(px_init[i], py_init[i], pz_init[i],
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            if(new_values_matrix(0,0) == -2){
                                bands_init[i].band_number = -1;
                                bands_init[i].effective_mass = 1;
                                bands_init[i].c1 = (1./3);
                                bands_init[i].c2 = (1./3);
                                bands_init[i].c3 = (1./3);

                                new_overlap << (1./3), (1./3), (1./3);
                            }
                            else{
                                bands_init[i].band_number = chosen_band;
                            
                                eigenval = new_values_matrix(0,chosen_band);
                                bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                bands_init[i].c1 = new_overlap(0); 
                                bands_init[i].c2 = new_overlap(1); 
                                bands_init[i].c3 = new_overlap(2);
                            }

                            if(i != 0 && i != index_two){
                                prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[i-1], new_overlap);
                                prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], bands_fin[i]);
                            }
                            else if(i == index_two){
                                prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[index_two-1], bands_fin[index_two]);
                            }

                        }
                        else if(_num_bands == 1){
                            bands_init[i].effective_mass = computeEffMassSingleBand(px_init[i], py_init[i], pz_init[i], 
                                                                                    _m_x_el, _m_y_el, _m_z_el);
                        }
                    }
                    else{
                        px_init[i] = _propagators[i].el_propagator_kx + 2*w_x;
                        py_init[i] = _propagators[i].el_propagator_ky + 2*w_y;
                        pz_init[i] = _propagators[i].el_propagator_kz + 2*w_z;

                        if(_num_bands == 3){
                            if(i != total_order){
                                chosen_band = band_number(gen);
                                new_values_matrix = diagonalizeLKHamiltonian(px_init[i], py_init[i], pz_init[i],
                                                                            _A_LK_el, _B_LK_el, _C_LK_el);
                                if(new_values_matrix(0,0) == -2){
                                    bands_init[i].band_number = -1;
                                    bands_init[i].effective_mass = 1;
                                    bands_init[i].c1 = (1./3);
                                    bands_init[i].c2 = (1./3);
                                    bands_init[i].c3 = (1./3);

                                    new_overlap << (1./3), (1./3), (1./3);
                                }
                                else{
                                    bands_init[i].band_number = chosen_band;
                            
                                    eigenval = new_values_matrix(0,chosen_band);
                                    bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                    bands_init[i].c1 = new_overlap(0); 
                                    bands_init[i].c2 = new_overlap(1); 
                                    bands_init[i].c3 = new_overlap(2);
                                }
                                
                                if(i != index_one){
                                    prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[i-1], new_overlap);
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], bands_fin[i]);
                                }
                                else{
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[i-1], bands_fin[i]);
                                }
                            }
                            else{
                                bands_init[total_order] = bands_init[0]; // first and last propagators, 4-momentum conservation

                                prefactor_init = prefactor_init*vertexOverlapTerm(bands_init[total_order-1], new_overlap);
                                prefactor_fin = prefactor_fin*vertexOverlapTerm(bands_fin[total_order-1], bands_fin[total_order]);
                            }
                        }
                        else if(_num_bands == 1){
                            bands_init[i].effective_mass = computeEffMassSingleBand(px_init[i], py_init[i], pz_init[i], 
                                                                                    _m_x_el, _m_y_el, _m_z_el);
                        }
                    }
                    energy_init[i] = electronEnergy(px_init[i], py_init[i], pz_init[i], bands_init[i].effective_mass);
                    energy_fin[i] = electronEnergy(px_fin[i], py_fin[i], pz_fin[i], bands_fin[i].effective_mass);

                    action_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                    action_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
                }

                /*for(int i = 0; i < index_one; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + w_z;
                }
                for(int i = index_one; i < index_two; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + 2*w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + 2*w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + 2*w_z;
                }
                for(int i = index_two; i < total_order + 1; i++){
                    px_init[i] = _propagators[i].el_propagator_kx + w_x;
                    py_init[i] = _propagators[i].el_propagator_ky + w_y;
                    pz_init[i] = _propagators[i].el_propagator_kz + w_z;
                }

                for(int i = 0; i < total_order + 1; i++){
                    px_fin[i] = _propagators[i].el_propagator_kx;
                    py_fin[i] = _propagators[i].el_propagator_ky;
                    pz_fin[i] = _propagators[i].el_propagator_kz;
                }*/
                
                

                /*for(int i = 0; i < total_order + 1; i++){
                    energy_init[i] = electronDispersion(px_init[i], py_init[i], pz_init[i], _el_eff_mass);
                    energy_fin[i] = electronDispersion(px_fin[i], py_fin[i], pz_fin[i], _el_eff_mass);
                }*/

                /*for(int i = 0; i < total_order + 1; i++){
                    action_init += energy_init[i]*(_vertices[i+1].tau - _vertices[i].tau);
                    action_fin += energy_fin[i]*(_vertices[i+1].tau - _vertices[i].tau);
                }*/

                delete[] px_fin;
                delete[] py_fin;
                delete[] pz_fin;

                delete[] energy_init;
                delete[] energy_fin;

                long double tau_current = _vertices[total_order+1].tau; // length of current diagram

                double p_A = _p_add_ext*_current_ph_ext;
                double p_B = _p_rem_ext;

                double mass_q = 1;
                if(_num_bands == 1){mass_q = computeEffMassSingleBand(w_x, w_y, w_z, _m_x_el, _m_y_el, _m_z_el);}
                
                prefactor_fin = prefactor_fin*vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q)
                                        *vertexStrengthTerm(w_x, w_y, w_z, _V_BZ, _V_BvK, _phonon_modes[phonon_index],
                                        _dielectric_responses[phonon_index], _dielectric_const, mass_q);

                double numerator = p_A*std::pow(2*M_PI,_D)*phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*tau_one)
                                *phononEnergy(_phonon_modes, phonon_index)*std::exp(-phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two))
                                *std::pow(((tau_current-tau_two+tau_one)/(2*M_PI)), (double)_D/2.)
                                *std::exp(-((std::pow(w_x,2)+std::pow(w_y,2)+std::pow(w_z,2))/2)*(tau_current-tau_two+tau_one))*prefactor_init;
                
                double denominator = p_B*std::exp(-(action_fin - action_init + phononEnergy(_phonon_modes, phonon_index)*(tau_current-tau_two+tau_one)))
                                *prefactor_fin*_V_BZ; // *_num_bands*_num_bands

                double R_rem = numerator/denominator;

                if(!(Metropolis(R_rem))){
                    delete[] px_init; delete[] py_init; delete[] pz_init; 
                    delete[] bands_init; delete[] bands_fin;
                    return;
                }
                else{
                    for(int i=0; i < total_order+1; ++i){
                        _propagators[i].el_propagator_kx = px_init[i];
                        _propagators[i].el_propagator_ky = py_init[i];
                        _propagators[i].el_propagator_kz = pz_init[i];
                        _bands[i] = bands_init[i];
                    }

                    delete[] px_init; delete[] py_init; delete[] pz_init;
                    delete[] bands_init; delete[] bands_fin;

                    phVertexRemoveRoom(index_one, index_two); // remove room in vertices array
                    propagatorArrayRemoveRoom(index_one, index_two); // remove room in electron propagators array
                    bandArrayRemoveRoom(index_one, index_two); // remove room in bands array

                    _propagators[total_order-1].el_propagator_kx = 0;
                    _propagators[total_order-1].el_propagator_ky = 0;
                    _propagators[total_order-1].el_propagator_kz = 0;
                    _bands[total_order-1].effective_mass = 1.0;

                    _propagators[total_order].el_propagator_kx = 0;
                    _propagators[total_order].el_propagator_ky = 0;
                    _propagators[total_order].el_propagator_kz = 0;
                    _bands[total_order].effective_mass = 1.0;

                    if(_num_bands == 3){
                        _bands[total_order-1].band_number = -1;
                        _bands[total_order-1].c1 = (1./3);
                        _bands[total_order-1].c2 = (1./3);
                        _bands[total_order-1].c3 = (1./3);

                        _bands[total_order].band_number = -1;
                        _bands[total_order].c1 = (1./3);
                        _bands[total_order].c2 = (1./3);
                        _bands[total_order].c3 = (1./3);
                    }

                    _current_ph_ext -= 1; // update current number of external phonons
                    findLastPhVertex();
                    return;
                }            
        }
        else{return;}
    }
};

void GreenFuncNphBands::swapPhononPropagator(){
    if(_current_order_int < 4){return;} // swap not possible if internal order is less than 4
    else{
        std::uniform_int_distribution<int> distrib(1, _current_order_int + 2*_current_ph_ext - 1);
        int index_one = distrib(gen); // choose random internal propagator
        
        if((_vertices[index_one].type != +1) && (_vertices[index_one].type != -1)){return;} // reject if vertex does not belong to internal phonon propagator
        if((_vertices[index_one+1].type != +1) && (_vertices[index_one+1].type != -1)){return;} // reject if vertex does not belong to internal phonon propagator
        if((_vertices[index_one].linked == index_one + 1) || (_vertices[index_one].linked == -1)){return;} // reject if the two vertices are linked

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
        eff_mass_el_final = 1.0;

        // new proposed band for propagator
        int chosen_band = 0;
        
        // new proposed eigenfunction
        Eigen::Vector3d new_overlap;
        new_overlap << 1, 0 ,0;

        // vertex terms that go into R_swap evaluation
        double prefactor_fin = 1;
        double prefactor_init = 1;
        
        if(_num_bands == 3){
            Eigen::Matrix<double,4,3> new_values_matrix = diagonalizeLKHamiltonian(kx+v1*wx1-v2*wx2, ky+v1*wy1-v2*wy2, kz+v1*wz1-v2*wz2, _A_LK_el, _B_LK_el, _C_LK_el);

            // check if free propagator
            if(new_values_matrix(0,0) == -2){
                chosen_band = -1;
                eff_mass_el_final = 1.0;
                new_overlap << (1./3), (1./3), (1./3);
            }
            else{
                // choose at random one of the bands
                std::uniform_int_distribution<int> distrib_unif(0,2);
                chosen_band = distrib_unif(gen);

                double eigenval = new_values_matrix(0,chosen_band);
                eff_mass_el_final = computeEffMassfromEigenval(eigenval); // computing new proposed electron effective mass from chosen eigenvalue
                new_overlap = new_values_matrix.block<3,1>(1,chosen_band); // new proposed band eigenstate
            }

            // compute only overlap term for vertices that go into R_swap evaluation,
            // the strength term does not change
            prefactor_fin = vertexOverlapTerm(_bands[index_one-1], new_overlap)*vertexOverlapTerm(_bands[index_two], new_overlap);
            prefactor_init = vertexOverlapTerm(_bands[index_one-1], _bands[index_one])*vertexOverlapTerm(_bands[index_two], _bands[index_one]);
        }
        else if(_num_bands == 1){
            eff_mass_el_final = computeEffMassSingleBand(kx+v1*wx1-v2*wx2, ky+v1*wy1-v2*wy2, kz+v1*wz1-v2*wz2,
                                                        _m_x_el, _m_y_el, _m_z_el);
        }

        // compute energies
        double energy_final_el = electronEnergy(kx+v1*wx1-v2*wx2, ky+v1*wy1-v2*wy2, kz+v1*wz1-v2*wz2, eff_mass_el_final);
        double energy_initial_el = electronEnergy(kx, ky, kz, eff_mass_el_initial);
        double energy_phonons = ((phononEnergy(_phonon_modes, phonon_index1)*v1)-(phononEnergy(_phonon_modes, phonon_index2)*v2));
        
        // need to address negative values issue
        // compute transition probability
        double R_swap = (prefactor_fin/prefactor_init)*std::exp(-(energy_final_el-energy_initial_el-energy_phonons)*(tau2-tau1)); //*_num_bands // to be fixed
        double acceptance_ratio = std::min(1.,R_swap);

        if(drawUniformR() > acceptance_ratio){return;}
        else{
            // assign new momentum values to propagator
            _propagators[index_one].el_propagator_kx += v1*wx1-v2*wx2;
            _propagators[index_one].el_propagator_ky += v1*wy1-v2*wy2;
            _propagators[index_one].el_propagator_kz += v1*wz1-v2*wz2;
            
            // assign new band values
            _bands[index_one].effective_mass = eff_mass_el_final;
            if(_num_bands == 3){
                _bands[index_one].band_number = chosen_band;
                _bands[index_one].c1 = new_overlap[0];
                _bands[index_one].c2 = new_overlap[1];
                _bands[index_one].c3 = new_overlap[2];
            }

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
    int total_order = _current_order_int + 2*_current_ph_ext;
    if(total_order <= 0){return;} // reject if no vertices are present
    else{
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

        if(isEqual(tau_new, tau_init) || isEqual(tau_new, tau_fin) || tau_new < tau_init || tau_new > tau_fin){return;} // check for possible double precision errors
        
        _vertices[vertex_index].tau = tau_new; // assign new time value to vertex

        findLastPhVertex();
        return;
    }
};

long double GreenFuncNphBands::stretchDiagramLength(long double tau_init){
    int total_order = _current_order_int + 2*_current_ph_ext;
    double kx = 0, ky = 0, kz = 0;

    long double * new_taus = new long double[total_order+2];
    new_taus[0] = 0.L; // first vertex time value is always 0

    int c = 0;
    int phonon_index = -1;
    double phonon_lines_energies = 0;

    for(int i=1; i < total_order+2; ++i){

        kx = _propagators[i-1].el_propagator_kx;
        ky = _propagators[i-1].el_propagator_ky;
        kz = _propagators[i-1].el_propagator_kz;

        c = _vertices[i-1].type;
        phonon_index = _vertices[i-1].index;
        
        if(c == +1 || c == +2){
            phonon_lines_energies += phononEnergy(_phonon_modes, phonon_index);
        }
        else if( c == -1 || c == -2){
            phonon_lines_energies -= phononEnergy(_phonon_modes, phonon_index);
        }

        new_taus[i] = new_taus[i-1] - std::log(1-drawUniformR())/(electronEnergy(kx,ky,kz,_bands[i-1].effective_mass)
        -_chem_potential + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes) + phonon_lines_energies);

        if(new_taus[i] < new_taus[i-1] || isEqual(new_taus[i], new_taus[i-1])){
            delete[] new_taus; 
            return tau_init;
        }
    }

    if(isEqual(new_taus[total_order+1], _tau_max) || new_taus[total_order+1] > _tau_max){
        delete[] new_taus; 
        return tau_init;
    }

    for(int i = 0; i < total_order+2; ++i){
        _vertices[i].tau = new_taus[i];
    }
    findLastPhVertex();
    delete[] new_taus;
    return _vertices[total_order+1].tau;

    // initialize momentum values for first electron propagator
    /*double kx = _propagators[0].el_propagator_kx;
    double ky = _propagators[0].el_propagator_ky;
    double kz = _propagators[0].el_propagator_kz;

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

    */
};

long double GreenFuncNphBands::configSimulation(long double tau_length = 1.0L){

    if((!(isEqual(_kx,0)) || !(isEqual(_ky,0)) || !(isEqual(_kz,0))) && _flags.effective_mass){
        std::cerr << "Warning: kx, ky and kz should be equal to 0 to calculate effective mass." << std::endl;
        std::cerr << "Effective mass calculation is not possible." << std::endl;
        _flags.effective_mass = false;
    }

    if(_flags.Z_factor && _ph_ext_max == 0){
        std::cerr << "Warning: number of maximum external phonon must be greater than 0 to calculate Z factor." << std::endl;
        std::cerr << "Z factor calculation is not possible." << std::endl;
        _flags.Z_factor = false;
    }

    // print simulation parameters
    std::cout <<"Starting simulation..." << std::endl;
    std::cout << std::endl;
    std::cout << "Number of thermalization steps: " << getRelaxSteps() << std::endl;
    std::cout << "Number of diagrams to be generated: " << getNdiags() << std::endl;
    if(_master){
        std::cout <<  "Number of autocorrelation steps to be generated: " << getAutocorrSteps() << std::endl;
    }
    std::cout << "Maximum length of diagram: " << _tau_max << std::endl;
    std::cout << "Maximum number of internal phonons: " << _order_int_max/2 << std::endl;
    std::cout << "Maximum number of external phonons: " << _ph_ext_max << std::endl;
    std::cout << "Maximum diagram order: " << _order_int_max + 2*_ph_ext_max << std::endl;
    std::cout << std::endl;

    std::cout << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands: " << _num_bands << std::endl;
    if(_num_bands == 1){
        std::cout << "electronic effective masses: mx_el = " << _m_x_el << ", my_el = " << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }
    else if(_num_bands == 3){
        std::cout << "electronic Luttinger-Kohn parameters: A_LK_el = " 
        << _A_LK_el << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }
    if(_D == 3){
        std::cout << "free electron momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
    }
    std::cout << std::endl;

    std::cout << "1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " << _dielectric_const << std::endl;
    std::cout << "number of phonon modes: " << _num_phonon_modes << std::endl;
    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }
    std::cout << "Number of dimensions: " << _D << std::endl;
    std::cout << std::endl;

    // print MC update probabilities
    std::cout << "Update probabilities:" << std::endl;
    std::cout << "length update: " << _p_length << ", add internal update: " << _p_add_int <<
    ", remove internal update: " << _p_rem_int << "," << std::endl;
    std::cout << "add external update: " << _p_add_ext << ", remove external update: " << _p_rem_ext <<
    ", swap update: " << _p_swap << "," << std::endl;
    std::cout << "shift update: " << _p_shift << ", stretch update: " << _p_stretch << std::endl; 
    std::cout << std::endl;

    if(_flags.gf_exact){
        std::cout << "Green Function will be calculated exactly." << std::endl;
        std::cout << "Number of computed points: " << _num_points << std::endl;

        if(_selected_order >= 0){
            std::cout << "The selected number of external phonons is: " << _selected_order << std::endl;
        }
        else {
            std::cout << "GF will be calculated for every number of external phonons." << std::endl;
        }

        _points = new long double[_num_points];
        _points_gf_exact = new long double[_num_points];

        // initialize GF
        for(int i=0; i<_num_points; i++){
            _points[i] = _points_center + i*_points_step;
            _points_gf_exact[i] = 0;
        }
        std::cout << std::endl;
    }

    if(_flags.histo){
        std::cout << "Green Function will be computed using the histogram method" << std::endl;
        std::cout << "Number of bins: " << _N_bins << std::endl;
        _bin_width_inv = 1./_bin_width;
        _histogram = new double[_N_bins];
        _bin_count = new unsigned long long int[_N_bins];
        _green_func = new double[_N_bins];
        for(int i=0; i<_N_bins; i++){
            _histogram[i] = _bin_center + i*_bin_width;
            _bin_count[i] = 0;
            _green_func[i] = 0;
        }
        std::cout << std::endl;
    }

    if(_flags.gs_energy){
        std::cout << "Ground state energy will be calculated using the exact estimator" << std::endl;
        std::cout << "Free electron momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
        std::cout << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        if(_num_bands == 1){
            std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                    << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                    << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;

        for(int i=0; i<_num_phonon_modes; i++){
            std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
            << _dielectric_responses[i] << std::endl;
        }
        std::cout << std::endl;
    }

    if(_flags.effective_mass){
        std::cout << "Effective mass will be calculated using the exact estimator" << std::endl;
        std::cout << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        if(_num_bands == 1){
            std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                    << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                    << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
            _effective_masses_bands << 0, 0, 0,
                                       0, 0, 0,
                                       0, 0, 0;
        }
        std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;

        for(int i=0; i<_num_phonon_modes; i++){
            std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
            << _dielectric_responses[i] << std::endl;
        }
        std::cout << std::endl;
    }

    /*if(_flags.Z_factor){
        std::cout << "Z factor will be calculated using the exact estimator" << std::endl;
        std::cout << "Coupling strength: " << _alpha << ", chemical potential: " << _chem_potential << ", max number of phonons: " << 
        _ph_ext_max << std::endl;
        initializeZFactorArray();
        std::cout << std::endl;
    }*/

    if(_flags.write_diagrams){
        if(getNdiags() > 25000){
            _flags.write_diagrams = false; // if too many diagrams are generated they are not printed to txt file
            std::cerr << "Warning: too many diagrams generated (> 25000), diagrams will not be printed to .txt file." << std::endl;
            std::cerr << std::endl;
        }
        else{
            std::cout << "The diagrams generated in the simulation process will be printed in the Diagrams.txt file" << std::endl;
            std::cout << std::endl;
        }    
    }

    if(_flags.time_benchmark){
        _benchmark_sim = new MC_Benchmarking(getNdiags(), 8);
        _benchmark_th = new MC_Benchmarking(getRelaxSteps(), 8);
        std::cout << "Time benchmark will be performed." << std::endl;
        std::cout << std::endl;
    }

    if(_flags.mc_statistics){
        std::cout << "Monte Carlo statistics will be calculated (simulation)." << std::endl;
        std::cout << "Average length of diagram, average order, average number of internal" << 
        " and external phonons and number of order 0 diagrams will be calculated." << std::endl;
        std::cout << "Cutoff for diagram length is set to: " << _tau_cutoff_statistics << std::endl;
        std::cout << std::endl;
    }

    if(_flags.fix_tau_value){
        std::cout << "Length of diagrams is fixed to: " << _tau_max << std::endl;
        tau_length = _tau_max - 1e-7L; // fix length of diagrams to tau_max
        _vertices[1].tau = tau_length; // assign fixed time value to first vertex
        std::cout << "All diagrams will have the same length (tau_max)." << std::endl;
        std::cout << std::endl;
        if(_p_length > 0 || _p_stretch > 0){
            std::cerr << "Warning: probabilities for length and stretch updates are set to non-zero values, but they will not be used." << std::endl;
            std::cerr << "Length and stretch updates will not be performed." << std::endl;

            double probs[8] = {0., _p_add_int, _p_rem_int, _p_add_ext, _p_rem_ext, _p_swap, _p_shift, 0.};
            setProbabilities(probs);

            std::cout << "Update probabilities have been adjusted to have a fixed diagram length." << std::endl;
            std::cout << "p_length = 0, p_stretch = 0" << std::endl;
            std::cout << "New update probabilities: add internal update = " << _p_add_int << ", remove internal update = " << _p_rem_int <<
            ", add external update = " << _p_add_ext << ", remove external update = " << _p_rem_ext <<
            ", swap update = " << _p_swap << ", shift update = " << _p_shift << std::endl;
            std::cout << std::endl;
        }
    }
    return tau_length;
};

void GreenFuncNphBands::configSimulationSilent(){
    if((!(isEqual(_kx,0)) || !(isEqual(_ky,0)) || !(isEqual(_kz,0))) && _flags.effective_mass){
        _flags.effective_mass = false;
    }

    if(_flags.Z_factor && _ph_ext_max == 0){
        _flags.Z_factor = false;
    }

    if(_flags.gf_exact){

        _points = new long double[_num_points];
        _points_gf_exact = new long double[_num_points];

        // initialize GF
        for(int i=0; i<_num_points; i++){
            _points[i] = _points_center + i*_points_step;
            _points_gf_exact[i] = 0;
        }
    }

    if(_flags.histo){
        _bin_width_inv = 1./_bin_width;
        _histogram = new double[_N_bins];
        _bin_count = new unsigned long long int[_N_bins];
        _green_func = new double[_N_bins];
        for(int i=0; i<_N_bins; i++){
            _histogram[i] = _bin_center + i*_bin_width;
            _bin_count[i] = 0;
            _green_func[i] = 0;
        }
    }

    if(_flags.effective_mass){
        if(_num_bands == 3){
            _effective_masses_bands << 0, 0, 0,
                                       0, 0, 0,
                                       0, 0, 0;
        }
    }

    /*if(_flags.Z_factor){
        initializeZFactorArray();
    }*/

    if(_flags.write_diagrams){
        if(getNdiags() > 25000){
            _flags.write_diagrams = false; // if too many diagrams are generated they are not printed to txt file
        }
    }

    if(_flags.time_benchmark){
        _benchmark_sim = new MC_Benchmarking(getNdiags(), 8);
        _benchmark_th = new MC_Benchmarking(getAutocorrSteps(), 8);
    }

    if(_flags.fix_tau_value){
        if(_p_length > 0 || _p_stretch > 0){

            double probs[8] = {0., _p_add_int, _p_rem_int, _p_add_ext, _p_rem_ext, _p_swap, _p_shift, 0.};
            setProbabilities(probs);
            std::cout << std::endl;
        }
    }
};

long double GreenFuncNphBands::chooseUpdate(long double tau_length, double r , MC_Benchmarking * benchmark){

    if(r <= _p_length){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        tau_length = diagramLengthUpdate(tau_length);
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(0);}
    }
    else if(r <= _p_length + _p_add_int){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        addInternalPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(1);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
            removeInternalPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(2);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        addExternalPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(3);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        removeExternalPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(4);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        swapPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(5);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        shiftPhononPropagator();
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(6);}
    }
    else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift + _p_stretch){
        if(_flags.time_benchmark){benchmark->startUpdateTimer();}
        tau_length = stretchDiagramLength(tau_length);
        if(_flags.time_benchmark){benchmark->stopUpdateTimer(7);}
    }
    return tau_length;
};

void GreenFuncNphBands::computeQuantities(long double tau_length, double r, int i){
    if(_flags.gf_exact && (_current_ph_ext == _selected_order || _selected_order < 0)){
        exactEstimatorGF(tau_length, _selected_order); // calculate Green function
        _gf_exact_count++; // count number of diagrams for normalization
    }

    if(_flags.histo){
        // select correct bin for histogram
        tau_length = (tau_length < 0.) ? 0. : (tau_length >= _tau_max) ? _tau_max - 1e-9 : tau_length;
        int bin = (int)((tau_length - 0.) * _bin_width_inv);
        _bin_count[bin]++;
    }

    if(_flags.gs_energy){
        _gs_energy += groundStateEnergyExactEstimator(tau_length); // accumulate energy of diagrams
    }

    if(_flags.effective_mass){
        _effective_mass += effectiveMassExactEstimator(tau_length); // accumulate effective mass of diagrams
    }

    // if(_flags.Z_factor){updateZFactor();} // accumulate Z factor data

    if(_flags.write_diagrams){
        writeDiagram("Diagrams.txt", i, r); // method to visualize diagram structure
    }

    if(_current_order_int == 0 && _current_ph_ext == 0){
            _N0++;
    }

    if(_flags.mc_statistics && (tau_length >= _tau_cutoff_statistics)){
        _mc_statistics.num_diagrams++; // accumulate number of diagrams
        _mc_statistics.avg_tau += tau_length; // accumulate average length of diagrams
        _mc_statistics.avg_tau_squared += tau_length*tau_length; // accumulate average squared length of diagrams
        _mc_statistics.avg_order += _current_order_int + 2*_current_ph_ext; // accumulate average order of diagrams
        _mc_statistics.avg_order_squared += (_current_order_int + 2*_current_ph_ext)*(_current_order_int + 2*_current_ph_ext); // accumulate average squared order of diagrams
        _mc_statistics.avg_ph_ext += _current_ph_ext; // accumulate average number of external phonons
        _mc_statistics.avg_ph_ext_squared += _current_ph_ext*_current_ph_ext; // accumulate average squared number of external phonons
        _mc_statistics.avg_ph_int += _current_order_int/2; // accumulate average number of internal phonons
        _mc_statistics.avg_ph_int_squared += (_current_order_int/2)*(_current_order_int/2); // accumulate average squared number of internal phonons
        if(_current_order_int + 2*_current_ph_ext == 0){
            _mc_statistics.zero_order_diagrams++; // accumulate number of zero order diagrams
        }
    }
};

void GreenFuncNphBands::computeFinalQuantities(){
    if(_flags.gf_exact){
        double norm_const = calcNormConst();
        for(int i=0; i<_num_points; i++){
            _points_gf_exact[i] = _points_gf_exact[i]*norm_const/_N0; // right normalization
        }
    }

    if(_flags.histo){
        double norm_const = calcNormConst();
        normalizeHistogram(norm_const);
    }

    if(_flags.gs_energy){
        _gs_energy = _gs_energy/static_cast<long double>(_gs_energy_count); // average energy of diagrams
    }

    if(_flags.effective_mass){
        long double effective_mass_inv = ((3.L/(static_cast<long double>(_m_x_el)+static_cast<long double>(_m_y_el))+static_cast<long double>(_m_z_el)) - _effective_mass/static_cast<long double>(_effective_mass_count)); // average effective mass of diagrams
        _effective_mass = 1.L/effective_mass_inv; // effective mass is inverse of the value calculated

        if(_num_bands == 1){
            long double effective_masses_inv[3] =  {0., 0., 0.};
            effective_masses_inv[0] = _effective_masses[0]/static_cast<long double>(_effective_mass_count);
            effective_masses_inv[1] = _effective_masses[1]/static_cast<long double>(_effective_mass_count);
            effective_masses_inv[2] = _effective_masses[2]/static_cast<long double>(_effective_mass_count);
            _effective_masses[0] = (1.L)/effective_masses_inv[0];
            _effective_masses[1] = (1.L)/effective_masses_inv[1];
            _effective_masses[2] = (1.L)/effective_masses_inv[2];
        }

        else if (_num_bands == 3){
            // stuff
        }
    }

    /*if(_flags.Z_factor){
        std::string a = "Z_factor_alpha";
        auto b = std::to_string(_alpha);
        std::string c = ".txt";
        std::string filename = a + b + c;
        writeZFactor(filename);
        std::cout << "Z factor calculated." << std::endl;
        std::cout << std::endl;
    }*/

    if(_flags.mc_statistics){
        _mc_statistics.avg_tau /= static_cast<long double>(_mc_statistics.num_diagrams); // average length of diagrams
        _mc_statistics.avg_tau_squared /= static_cast<long double>(_mc_statistics.num_diagrams); // average squared length of diagrams
    }
};

void GreenFuncNphBands::printGFExactEstimator(){
    double norm_const = calcNormConst();
    for(int i=0; i<_num_points; i++){
        //_points_gf_exact[i] = _points_gf_exact[i]/((double)_gf_exact_count);
        _points_gf_exact[i] = _points_gf_exact[i]*norm_const/_N0; // right normalization
    }
    std::string a = "GF_";
    auto b = std::to_string(_selected_order);
    if(_selected_order < 0){
        b = "total";
    }
    std::cout << "Exact Green's function computed." << std::endl;
    std::string c = "_exact.txt";
    writeExactGF(a+b+c); // write Green function to file
    std::cout << std::endl;
};

void GreenFuncNphBands::printhistogramEstimator(){
    double norm_const = calcNormConst();
    normalizeHistogram(norm_const);
    std::cout << "Histogram computed." << std::endl;
    writeHistogram("histo.txt");
    std::cout << std::endl;
};

void GreenFuncNphBands::printGroundStateEnergyEstimator(){
    _gs_energy = _gs_energy/static_cast<long double>(_gs_energy_count); // average energy of diagrams
    std::cout << "Ground state energy of the system is: " << _gs_energy << ". Input parameters are: kx = " << _kx << 
    ", ky = " << _ky << ", kz = " << _kz << std::endl;
    std::cout << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
    std::cout << "minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;
    std::cout << "Number of diagrams used for ground state energy calculation: " << _gs_energy_count << std::endl;
    if(_num_bands == 1){
        std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
            << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
    else if(_num_bands == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
            << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }
    std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
    std::cout << std::endl;

    std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;
    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }

    std::string filename = "gs_energy.txt";
    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not gs_energy.txt open file " << filename << std::endl;
    }
    else{
        file << "Ground state energy of the system is: " << _gs_energy << " . Input parameters are: kx = " << _kx << 
            ", ky = " << _ky << ", kz = " << _kz << std::endl;
        file << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        file << " minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;
        file << "Number of diagrams used for ground state energy calculation: " << _gs_energy_count << std::endl;

        if(_num_bands == 1){
            file << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        file << "1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        file << std::endl;

        file << "Number of phonon modes: " << _num_phonon_modes << std::endl;
        for(int i=0; i<_num_phonon_modes; i++){
            file << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
            << _dielectric_responses[i] << std::endl;
        }
        file << std::endl;
        file.close();
    }
    std::cout << std::endl;
};

void GreenFuncNphBands::printEffectiveMassEstimator(){
    long double effective_mass_inv = ((3.L/(static_cast<long double>(_m_x_el)+static_cast<long double>(_m_y_el))+static_cast<long double>(_m_z_el)) - _effective_mass/static_cast<long double>(_effective_mass_count)); // average effective mass of diagrams
    _effective_mass = 1.L/effective_mass_inv; // effective mass is inverse of the value calculated
    std::cout << "Input parameters are: chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
    std::cout << "Number of diagrams used for effective mass calculation: " << _effective_mass_count << std::endl;
    if(_num_bands == 1){
        long double effective_masses_inv[3] =  {0., 0., 0.};
        effective_masses_inv[0] = _effective_masses[0]/static_cast<long double>(_effective_mass_count);
        effective_masses_inv[1] = _effective_masses[1]/static_cast<long double>(_effective_mass_count);
        effective_masses_inv[2] = _effective_masses[2]/static_cast<long double>(_effective_mass_count);
        _effective_masses[0] = (1.L)/effective_masses_inv[0];
        _effective_masses[1] = (1.L)/effective_masses_inv[1];
        _effective_masses[2] = (1.L)/effective_masses_inv[2];
        std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
            << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }

    else if(_num_bands == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }

    std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_mass << std::endl;
    std::cout << std::endl;

    std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;
    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }
    std::cout << std::endl;
    if(_num_bands == 1){
        std::cout << "Polaronic effective masses are: mx_pol = " << (_effective_masses[0]) << ", my_pol = " 
            << _effective_masses[1] << ", mz_pol = " << _effective_masses[2] << std::endl; 
    }
    else if(_num_bands == 3){
        // stuff
    }
    std::cout << std::endl;
    std::cout << "Average effective mass of diagrams is: " << /*static_cast<long double>(_m_x_el/3.+_m_y_el/3.+_m_z_el/3.)**/_effective_mass << "." << std::endl;
    std::cout << "Average inverse effective mass of system is: " << effective_mass_inv << "." << std::endl;
    std::cout << std::endl;

    std::string filename = "effective_mass.txt";
    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not effective_mass.txt open file " << filename << std::endl;
    }
    else{
        file << "Effective mass of the system is: " << _effective_mass << "." << std::endl;
        file << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        file << "Number of diagrams used for effective mass calculation: " << _effective_mass_count << std::endl;
        if(_num_bands == 1){
            file << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                    << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        file <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_mass << std::endl;
        file << std::endl;

        file << "Number of phonon modes: " << _num_phonon_modes << std::endl;
        for(int i=0; i<_num_phonon_modes; i++){
            file << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
                << _dielectric_responses[i] << std::endl;
        }
        file << std::endl;

        if(_num_bands == 1){
            file << "Polaronic effective masses are: mx_pol = " << _effective_masses[0] << ", my_pol = " 
                << _effective_masses[1] << ", mz_pol = " << _effective_masses[2] << std::endl; 
        }
       else if(_num_bands == 3){
            // stuff
        }
        file << std::endl;
        file << "Average effective mass of diagrams is: " << _effective_mass << "." << std::endl;
        file << "Average inverse effective mass of the system is: " << effective_mass_inv << "." << std::endl;
        file << std::endl;

        file.close();
    }
};

void GreenFuncNphBands::printMCStatistics(){
    _mc_statistics.avg_tau /= static_cast<long double>(_mc_statistics.num_diagrams); // average length of diagrams
    _mc_statistics.avg_tau_squared /= static_cast<long double>(_mc_statistics.num_diagrams); // average squared length of diagrams

    std::cout << "Monte Carlo statistics:" << std::endl;
       
    std::cout << "chemical potential: " << _chem_potential << ", number of degenerate electronic bands: " << _num_bands << std::endl;
    if(_num_bands == 1){
        std::cout << "electronic effective masses: mx_el = " << _m_x_el << ", my_el = " << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }
    else if(_num_bands == 3){
        std::cout << "electronic Luttinger-Kohn parameters: A_LK_el = " << _A_LK_el << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }

    std::cout << "total momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
    std::cout << std::endl;
    std::cout << "1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " << _dielectric_const << std::endl;
    std::cout << std::endl;
    std::cout << "number of phonon modes: " << _num_phonon_modes << std::endl;
    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Cutoff for statistics: " << _tau_cutoff_statistics << std::endl;
    std::cout << "Number of diagrams (taken into account): " << _mc_statistics.num_diagrams << std::endl;
    std::cout << "Average length of diagrams: " << _mc_statistics.avg_tau << std::endl;
    std::cout << "Std dev length of diagrams: " << std::sqrt(_mc_statistics.avg_tau_squared - _mc_statistics.avg_tau*_mc_statistics.avg_tau) << std::endl;
    std::cout << "Average order of diagrams: " << static_cast<long double>(_mc_statistics.avg_order)/static_cast<long double>(_mc_statistics.num_diagrams) << std::endl;
    std::cout << "Std dev order of diagrams: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_order_squared)/static_cast<long double>(_mc_statistics.num_diagrams)
        - static_cast<long double>(_mc_statistics.avg_order*_mc_statistics.avg_order)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << std::endl;
    std::cout << "Average number of internal phonons: " << static_cast<long double>(_mc_statistics.avg_ph_int)/static_cast<long double>(_mc_statistics.num_diagrams) << std::endl;
    std::cout << "Std dev number of internal phonons: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_ph_int_squared)/static_cast<long double>(_mc_statistics.num_diagrams)
        - static_cast<long double>(_mc_statistics.avg_ph_int*_mc_statistics.avg_ph_int)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << std::endl;
    std::cout << "Average number of external phonons: " << static_cast<long double>(_mc_statistics.avg_ph_ext)/static_cast<long double>(_mc_statistics.num_diagrams) << std::endl;
    std::cout << "Std dev number of external phonons: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_ph_ext_squared)/static_cast<long double>(_mc_statistics.num_diagrams) 
        - static_cast<long double>(_mc_statistics.avg_ph_ext*_mc_statistics.avg_ph_ext)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << std::endl;
    std::cout << "Number of zero order diagrams: " << _mc_statistics.zero_order_diagrams << std::endl;
    std::cout << std::endl;
    writeMCStatistics("MC_Statistics.txt");
};

void GreenFuncNphBands::markovChainMC(){

    // input variables
    long double tau_length = diagramLengthUpdate(_tau_max/100);
    double r = 0.5;
    unsigned long long int i = 0;

    unsigned long long int N_diags = getNdiags(); // number of diagrams to be generated
    unsigned long long int N_relax = getRelaxSteps(); // number of thermalization steps

    tau_length = configSimulation(tau_length); // print simulation parameters

    //double bin_width_inv = 1./_bin_width;

    std::cout << "Starting thermalization process" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        std::cout << "Benchmarking thermalization time..." << std::endl;
        std::cout << std::endl;
        _benchmark_th->startTimer();
    }

    ProgressBar bar(N_relax, 70);
    for(i = 0; i < N_relax; ++i){
        r = drawUniformR();

        tau_length = chooseUpdate(tau_length, r, _benchmark_th);

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(_current_order_int + 2*_current_ph_ext, 1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors(_current_order_int+2*_current_ph_ext);}

        if(static_cast<int>(i%100000) == 0 && _flags.write_diagrams){
            writeDiagram("Diagrams.txt", i, r); // method to visualize diagram structure
        }

        if(static_cast<int>(i%(N_relax/100)) == 0){bar.update(i);}
    }
    bar.finish();

    if(_flags.time_benchmark){_benchmark_th->stopTimer();}

    std::cout << std::endl;
    std::cout << "Thermalization process finished" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        _benchmark_th->printResults();
        _benchmark_th->writeResultsToFile("thermalization_benchmark.txt");
        delete _benchmark_th;
        std::cout << std::endl;
    }

    i = 0;
    std::cout << "Starting simulation process" << std::endl;
    std::cout << std::endl;

    bar.setTotal(N_diags);

    if(_flags.time_benchmark){
        std::cout << "Benchmarking simulation time..." << std::endl;
        std::cout << std::endl;
        _benchmark_sim->startTimer();
    }

    for(i = 0; i < N_diags; ++i){
        r = drawUniformR();
        if(_flags.write_diagrams){writeChosenUpdate("Updates.txt", i, r);}

        tau_length = chooseUpdate(tau_length, r, _benchmark_sim);

        computeQuantities(tau_length, r, i); // compute desired quantities from MC computation

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(_current_order_int + 2*_current_ph_ext, 1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors(_current_order_int+2*_current_ph_ext);}

        if(static_cast<int>(i%(N_diags/100)) == 0){bar.update(i);}
    }

    bar.finish();

    std::cout << std::endl;

    if(_flags.time_benchmark){
        _benchmark_sim->stopTimer();
        std::cout << "Simulation time benchmark finished." << std::endl;
        _benchmark_sim->printResults(); 
        _benchmark_sim->writeResultsToFile("simulation_benchmark.txt");
        delete _benchmark_sim;
        std::cout << std::endl;
    }

    if(_flags.gf_exact){
        printGFExactEstimator();
    }

    if(_flags.histo){
        printhistogramEstimator();
    }

    if(_flags.gs_energy){
        printGroundStateEnergyEstimator();
    }

    if(_flags.effective_mass){
        printEffectiveMassEstimator();
    }

    /*if(_flags.Z_factor){
        std::string a = "Z_factor_alpha";
        auto b = std::to_string(_alpha);
        std::string c = ".txt";
        std::string filename = a + b + c;
        writeZFactor(filename);
        std::cout << "Z factor calculated." << std::endl;
        std::cout << std::endl;
    }*/

    if(_flags.mc_statistics){
        printMCStatistics();
    }
};

void GreenFuncNphBands::markovChainMCOnlyRelax(){
    // input variables
    long double tau_length = diagramLengthUpdate(_tau_max/100);
    double r = 0.5;
    unsigned long long int i = 0;

    unsigned long long int N_relax = getRelaxSteps(); // number of thermalization steps

    tau_length = configSimulation(tau_length); // print simulation parameters

    std::cout << "Starting thermalization process" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        std::cout << "Benchmarking thermalization time..." << std::endl;
        std::cout << std::endl;
        _benchmark_th->startTimer();
    }

    ProgressBar bar(N_relax, 70);
    for(i = 0; i < N_relax; ++i){
        r = drawUniformR();

        tau_length = chooseUpdate(tau_length, r, _benchmark_th);

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(_current_order_int + 2*_current_ph_ext, 1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors(_current_order_int+2*_current_ph_ext);}

        if(static_cast<int>(i%100000) == 0 && _flags.write_diagrams){
            writeDiagram("Diagrams.txt", i, r); // method to visualize diagram structure
        }

        if(static_cast<int>(i%(N_relax/100)) == 0){bar.update(i);}
    }
    bar.finish();

    if(_flags.time_benchmark){_benchmark_th->stopTimer();}

    std::cout << std::endl;
    std::cout << "Thermalization process finished" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        _benchmark_th->printResults();
        _benchmark_th->writeResultsToFile("thermalization_benchmark.txt");
        delete _benchmark_th;
        delete _benchmark_sim;
        std::cout << std::endl;
    }
};

void GreenFuncNphBands::markovChainMCOnlySample(){
    long double tau_length = _vertices[_current_order_int + 2*_current_ph_ext + 1].tau; // get current length of diagram
    double r = 0.5;
    unsigned long long int i = 0;

    configSimulationSilent();

    unsigned long long int N_autocorr = getAutocorrSteps(); // number of autocorrelation steps
    unsigned long long int N_diags = getNdiags(); // number of diagrams to be generated

    ProgressBar bar(N_autocorr, 70);

    if(_master){
        std::cout << "Starting autocorrelation process (to obtain independent results)" << std::endl;
        std::cout << std::endl;

        bar.setTotal(N_autocorr);

        if(_flags.time_benchmark){
            std::cout << "Benchmarking autocorrelation process time..." << std::endl;
            std::cout << std::endl;
            _benchmark_th->startTimer();
        }
    }

    for(i = 0; i < N_autocorr; ++i){
        r = drawUniformR();

        tau_length = chooseUpdate(tau_length, r, _benchmark_th);

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(_current_order_int + 2*_current_ph_ext, 1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors(_current_order_int+2*_current_ph_ext);}

        if(static_cast<int>(i%(N_autocorr/100)) == 0 && _master){bar.update(i);}

    }

    if(_master){
        bar.finish();
        std::cout << std::endl;
        std::cout << "Autocorrelation process finished." << std::endl;
        std::cout << std::endl;
        if(_flags.time_benchmark){
            _benchmark_th->stopTimer();
            _benchmark_th->printResults(); 
            _benchmark_th->writeResultsToFile("autocorrelation_benchmark.txt");
        }
    }

    if(_flags.time_benchmark){delete _benchmark_th;}
    
    _flags.write_diagrams = false; 

    if(_master){
        std::cout << "Starting simulation process" << std::endl;
        std::cout << std::endl;

        _flags.write_diagrams = true; // set true for master process

        bar.setTotal(N_diags);
        
        if(_flags.time_benchmark){
            std::cout << "Benchmarking simulation time..." << std::endl;
            std::cout << std::endl;
            _benchmark_sim->startTimer();
        }
    }

    for(i = 0; i < N_diags; ++i){
        r = drawUniformR();
        if(_flags.write_diagrams){writeChosenUpdate("Updates.txt", i, r);}

        tau_length = chooseUpdate(tau_length, r, _benchmark_sim);

        computeQuantities(tau_length, r, i); // compute desired quantities from MC computation

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(_current_order_int + 2*_current_ph_ext, 1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors(_current_order_int+2*_current_ph_ext);}
    
        if(_flags.write_diagrams && _master){writeDiagram("Diagrams.txt", i, r);}
        
        if(static_cast<int>(i%(N_diags/100)) == 0 && _master){bar.update(i);}
    }

    if(_master){
        bar.finish();

        std::cout << std::endl;

        if(_flags.time_benchmark){
            _benchmark_sim->stopTimer();
            std::cout << "Simulation time benchmark finished." << std::endl;
            _benchmark_sim->printResults(); 
            _benchmark_sim->writeResultsToFile("simulation_benchmark.txt");
            std::cout << std::endl;
        }
    }

    if(_flags.time_benchmark){delete _benchmark_sim;}

    computeFinalQuantities(); // compute final quantities for all processes (no write to console or files)
    
};

double GreenFuncNphBands::groundStateEnergyExactEstimator(long double tau_length){
    if(tau_length <= _tau_cutoff_energy){return 0;} // reject if below cutoff
    else{

        // compute electron bare propagators action
        int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
        double electron_energy = 0;
        double electron_action = 0, phonon_action = 0;

        // phonon bare propagators action
        //int i = 0;
        int index_two = 0;
        int phonon_index = -1;
        long double tau_one = 0;
        long double tau_two = 0;
        // int int_count = 0;
        // int ext_count = 0;  
        // bool int_flag = false;
        // bool ext_flag = false;

        // main loop
        for(int i=0; i<current_order+1; ++i){
            electron_energy = electronEnergy(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, 
                _propagators[i].el_propagator_kz, _bands[i].effective_mass);
            electron_action += electron_energy*(_vertices[i+1].tau - _vertices[i].tau);

            if(_vertices[i].type == +1){
                index_two = _vertices[i].linked;
                phonon_index = _vertices[i].index;
                tau_one = _vertices[i].tau;  
                tau_two = _vertices[index_two].tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one);
            }
            else if(_vertices[i].type == -2){
                index_two = _vertices[i].linked;
                phonon_index = _vertices[i].index;
                tau_one = _vertices[i].tau;  
                tau_two = _vertices[index_two].tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_length + tau_one - tau_two);
            }
        }

        /*while(i < current_order + 1 && !(int_flag && ext_flag)){
            if(_vertices[i].type == +1){
                index_two = _vertices[i].linked;
                phonon_index = _vertices[i].index;
                tau_one = _vertices[i].tau;  
                tau_two = _vertices[index_two].tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one);
                int_count++;
                if(int_count == _current_order_int/2){int_flag = true;}
            }
            else if(_vertices[i].type == -2){
                index_two = _vertices[i].linked;
                tau_one = _vertices[i].tau;  
                tau_two = _vertices[index_two].tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_length + tau_one - tau_two);
                ext_count++;
                if(ext_count == _current_ph_ext){ext_flag = true;}
            }
            i++;
        }*/

        double diagram_energy = (electron_action + phonon_action - static_cast<double>(current_order))/tau_length; // energy of current diagram
        _gs_energy_count++; // update number of computed exact energy estimators
        return diagram_energy; // return energy of current diagram
    }
};

void GreenFuncNphBands::exactEstimatorGF(long double tau_length, int ext_phonon_order){

    if(ext_phonon_order < 0){ext_phonon_order = _current_ph_ext;}

    // compute electron bare propagators action
    int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
    double electron_energy = 0;
    double electron_action = 0., phonon_action = 0.;

    // phonon bare propagators action
    //int i = 0;
    int index_two = 0;
    int phonon_index = -1;
    long double tau_one = 0;
    long double tau_two = 0;
    // int int_count = 0;
    // int ext_count = 0;  
    // bool int_flag = false;
    // bool ext_flag = false;

    // main loop
    for(int i=0; i<current_order+1; ++i){
        electron_energy = electronEnergy(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, 
            _propagators[i].el_propagator_kz, _bands[i].effective_mass);
        electron_action += electron_energy*(_vertices[i+1].tau - _vertices[i].tau);

        if(_vertices[i].type == +1){
            index_two = _vertices[i].linked;
            phonon_index = _vertices[i].index;
            tau_one = _vertices[i].tau;  
            tau_two = _vertices[index_two].tau;
            phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one);
        }
        else if(_vertices[i].type == -2){
            index_two = _vertices[i].linked;
            phonon_index = _vertices[i].index;
            tau_one = _vertices[i].tau;  
            tau_two = _vertices[index_two].tau;
            phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_length + tau_one - tau_two);
        }
    }

    // compute Green function with exact estimator
    tau_length = (tau_length < 1e-9) ? 1e-9 : (tau_length >= _tau_max) ? _tau_max - 1e-9 : tau_length;
    int bin = static_cast<int>((tau_length - 0.) * 1./_points_step);

    long double prefactor = std::pow(1 + (_points[bin] - tau_length)/tau_length, current_order);
    long double exponential = std::exp(-((_points[bin] - tau_length)/tau_length)*(electron_action + phonon_action));
    long double diagrams_ratio = prefactor*exponential;

    _points_gf_exact[bin] += diagrams_ratio/(_points_step); // accumulate Green function value
};

double GreenFuncNphBands::effectiveMassExactEstimator(long double tau_length){
    if(tau_length <= _tau_cutoff_mass){return 0;}
    else{
        int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
        double mass_average_inv_x = 0; double mass_average_inv_y = 0; double mass_average_inv_z = 0;
        double electron_average_kx = 0; double electron_average_ky = 0; double electron_average_kz = 0;
        double mx, my, mz;

        if(_num_bands == 1){
            for(int i=0; i< current_order+1; ++i){
                /*if(isEqual(_propagators[i].el_propagator_kx, 0) && isEqual(_propagators[i].el_propagator_ky, 0) && isEqual(_propagators[i].el_propagator_kz, 0)){
                    mx = _m_x_el; my = _m_y_el; mz = _m_z_el;
                }
                else{*/
                mx = _m_x_el; my = _m_y_el; mz = _m_z_el;
                //}
                mass_average_inv_x += (_vertices[i+1].tau - _vertices[i].tau)/mx;
                mass_average_inv_y += (_vertices[i+1].tau - _vertices[i].tau)/my;
                mass_average_inv_z += (_vertices[i+1].tau - _vertices[i].tau)/mz;

                electron_average_kx += _propagators[i].el_propagator_kx*(_vertices[i+1].tau - _vertices[i].tau)/mx;
                electron_average_ky += _propagators[i].el_propagator_ky*(_vertices[i+1].tau - _vertices[i].tau)/my;
                electron_average_kz += _propagators[i].el_propagator_kz*(_vertices[i+1].tau - _vertices[i].tau)/mz;
            }

            mass_average_inv_x = (1./tau_length)*mass_average_inv_x;
            mass_average_inv_y = (1./tau_length)*mass_average_inv_y;
            mass_average_inv_z = (1./tau_length)*mass_average_inv_z;

            electron_average_kx = (1./tau_length)*electron_average_kx*electron_average_kx;
            electron_average_ky = (1./tau_length)*electron_average_ky*electron_average_ky;
            electron_average_kz = (1./tau_length)*electron_average_kz*electron_average_kz;

            _effective_masses[0] += static_cast<long double>(mass_average_inv_x - electron_average_kx); // x component
            _effective_masses[1] += static_cast<long double>(mass_average_inv_y - electron_average_ky); // y component
            _effective_masses[2] += static_cast<long double>(mass_average_inv_z - electron_average_kz); // z component
        }
        else if (_num_bands == 3){

        }

        /*for(int i=0; i<current_order+1; ++i){
            if(isEqual(_propagators[i].el_propagator_kx, 0) && isEqual(_propagators[i].el_propagator_ky, 0) && isEqual(_propagators[i].el_propagator_kz, 0)){
                mx = _m_x_el; my = _m_y_el; mz = _m_z_el;
            }
            else{
                mx = computeEffMassSingleBand(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, _propagators[i].el_propagator_kz, _m_x_el, _m_y_el, _m_z_el);
                my = mx;
                mz = mx;
            }
            electron_average_kx += ((_vertices[i+1].tau - _vertices[i].tau)/mx - std::pow(_propagators[i].el_propagator_kx,2)*std::pow(_vertices[i+1].tau - _vertices[i].tau,2)/(mx*mx));
            electron_average_ky += ((_vertices[i+1].tau - _vertices[i].tau)/my - std::pow(_propagators[i].el_propagator_ky,2)*std::pow(_vertices[i+1].tau - _vertices[i].tau,2)/(my*my));
            electron_average_kz += ((_vertices[i+1].tau - _vertices[i].tau)/mz- std::pow(_propagators[i].el_propagator_kz,2)*std::pow(_vertices[i+1].tau - _vertices[i].tau,2)/(mz*mz));
        }

        for(int i=0; i< current_order+1; ++i){
            electron_average_kx += _propagators[i].el_propagator_kx*(_vertices[i+1].tau - _vertices[i].tau);
            electron_average_ky += _propagators[i].el_propagator_ky*(_vertices[i+1].tau - _vertices[i].tau);
            electron_average_kz += _propagators[i].el_propagator_kz*(_vertices[i+1].tau - _vertices[i].tau);
        }

        electron_average_kx = electron_average_kx/tau_length;  //static_cast<long double>(_m_x_el*_m_x_el);
        electron_average_ky = electron_average_ky/tau_length;  //static_cast<long double>(_m_y_el*_m_y_el);
        electron_average_kz = electron_average_kz/tau_length;  //static_cast<long double>(_m_z_el*_m_z_el);

        if(_num_bands == 3){
            Eigen::RowVector3d unit;
            unit << 1, 1, 1;

            _effective_masses_bands.row(0) += (unit - tau_length*std::pow(electron_average_kx,2)
                            *diagonalizeLKHamiltonianEigenval(1,0,0,_A_LK_el,_B_LK_el,_C_LK_el)); // (100) direction
            _effective_masses_bands.row(1) += (unit - tau_length*(std::pow(electron_average_kx,2)+std::pow(electron_average_ky,2))
                            *diagonalizeLKHamiltonianEigenval(1,1,0,_A_LK_el,_B_LK_el,_C_LK_el)/(2.)); // (110) direction
            _effective_masses_bands.row(2) += (unit - tau_length*(std::pow(electron_average_kx,2)+std::pow(electron_average_ky,2)+std::pow(electron_average_kz,2))
                            *diagonalizeLKHamiltonianEigenval(1,1,1,_A_LK_el,_B_LK_el,_C_LK_el)/(2.)); // (111) direction
        }
        else if(_num_bands == 1){
            _effective_masses[0] += ((1.L/static_cast<long double>(_m_x_el)) - (tau_length*electron_average_kx*electron_average_kx)/static_cast<long double>(_m_x_el*_m_x_el));
            _effective_masses[1] += ((1.L/static_cast<long double>(_m_y_el)) - (tau_length*electron_average_ky*electron_average_ky)/static_cast<long double>(_m_y_el*_m_y_el));
            _effective_masses[2] += ((1.L/static_cast<long double>(_m_z_el)) - (tau_length*electron_average_kz*electron_average_kz)/static_cast<long double>(_m_z_el*_m_z_el));
            //_effective_masses[0] += 1.L/static_cast<long double>(_m_x_el) - electron_average_kx/tau_length; // x comp
            //_effective_masses[1] += 1.L/static_cast<long double>(_m_y_el) - electron_average_ky/tau_length; // y comp
            //_effective_masses[2] += 1.L/static_cast<long double>(_m_z_el) - electron_average_kz/tau_length; // z comp
            //_effective_masses[0] += (1.L - tau_length*std::pow(electron_average_kx,2))/static_cast<long double>(_m_x_el); // x component
            //_effective_masses[1] += (1.L - tau_length*std::pow(electron_average_ky,2))/static_cast<long double>(_m_y_el); // y component
            //_effective_masses[2] += (1.L - tau_length*std::pow(electron_average_kz,2))/static_cast<long double>(_m_z_el); // z component
        }
        
        */
       
        _effective_mass_count++;

        // return inverse of effective mass, _D dimensionality of the system
        //static_cast<long double>(_m_x_el/3.+_m_y_el/3.+_m_z_el/3.)
        return (tau_length*(std::pow(electron_average_kx,2) + std::pow(electron_average_ky,2) + std::pow(electron_average_kz,2))/_D);
        //return (1.L/(static_cast<long double>(_m_x_el/3.+_m_y_el/3.+_m_z_el/3.))-tau_length*(std::pow(electron_average_kx,2) + std::pow(electron_average_ky,2) + std::pow(electron_average_kz,2))/_D);
    }
};

void GreenFuncNphBands::calcGroundStateEnergy(std::string filename){
    _gs_energy = _gs_energy/(double)_gs_energy_count; // average energy of diagrams
    std::cout << "Ground state energy of the system is: " << _gs_energy << ". Input parameters are: kx = " << _kx << 
    ", ky = " << _ky << ", kz = " << _kz << std::endl;
    std::cout << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;

    if(_num_bands == 1){
        std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }
    else if(_num_bands == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }

    std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
    << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
    std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;

    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }

    std::cout << " minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;

    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not gs_energy.txt open file " << filename << std::endl;
    }
    else{
        file << "Ground state energy of the system is: " << _gs_energy << " . Input parameters are: kx = " << _kx << 
        ", ky = " << _ky << ", kz = " << _kz << std::endl;
        file << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        if(_num_bands == 1){
            file << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        file << "1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
            << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        file << "Number of phonon modes: " << _num_phonon_modes << std::endl;

        for(int i=0; i<_num_phonon_modes; i++){
            file << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
            << _dielectric_responses[i] << std::endl;
        }
        file << " minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;
        file << std::endl;
        file.close();
    }
    std::cout << std::endl;
};

/*void GreenFuncNphBands::calcEffectiveMasses(std::string filename){
    long double effective_mass_inv = _effective_mass/(double)_effective_mass_count; // average effective mass of diagrams
    _effective_mass = 1./effective_mass_inv; // effective mass is inverse of the value calculated
    //std::cout << "Effective mass of system is: " << _effective_mass << "." << std::endl;
    std::cout << "Input parameters are: chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;

    if(_num_bands == 1){
        std::cout << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
            << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }
    else if(_num_bands == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
            << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }

    std::cout <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
    << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
    std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;

    for(int i=0; i<_num_phonon_modes; i++){
        std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _born_effective_charges[i] << std::endl;
    }

    if(_num_bands == 1){
        std::cout << "Polaronic effective masses are: mx_pol = " << (1./_effective_masses[0]) << " my_pol = " 
            << (1./_effective_masses[1]) << "mz_pol = " << (1./_effective_masses[2]) << std::endl; 
    }
    else if(_num_bands == 3){
        double A_LK_pol = 0;
        double B_LK_pol = 0;
        double C_LK_pol = 0;
        std::cout << "Computed polaronic inverse masses for (100), (110) and (111) high symmetry directions." << std::endl;
        std::cout << "(100)" << std::endl;

        if(_A_LK_el >= _B_LK_el){
            std::cout << "2A = " << _effective_masses_bands(0,0) << ", 2B = " << _effective_masses_bands(0,1) 
                << ", 2B = " << _effective_masses_bands(0,2) << std::endl;
            A_LK_pol = _effective_masses_bands(0,0)/2;
            B_LK_pol = _effective_masses_bands(0,1)/2;

            std::cout << "(110)" << std::endl;

            if(_C_LK_el >= 0){
                if(_A_LK_el - _C_LK_el >= _B_LK_el){
                    std::cout << "A + B + C = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                    << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,1))/2;
                }
                else{
                    std::cout << "A + B + C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                    << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,2))/2;
                }
                std::cout << "(111)" << std::endl;
                std::cout << "(2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,2) << std::endl;
            }
            else{
                if(_A_LK_el + _C_LK_el > _B_LK_el){
                    std::cout << "A + B - C = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                    << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,0))/2;
                }
            else{
                    std::cout << "A + B - C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                    << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                }
                std::cout << "(111)" << std::endl;
                std::cout << "(2/3)(A + 2B - C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                << ", (2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,2) << std::endl;
            }
        }
        else{
            std::cout << "2B = " << _effective_masses_bands(0,0) << ", 2B = " << _effective_masses_bands(0,1) 
            << ", 2A = " << _effective_masses_bands(0,2) << std::endl;
            A_LK_pol = _effective_masses_bands(0,2)/2;
            B_LK_pol = _effective_masses_bands(0,0)/2;

            std::cout << "(110)" << std::endl;
        
            if(_C_LK_el >= 0){
                if(_A_LK_el - _C_LK_el >= _B_LK_el){
                    std::cout << "A + B + C = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                    << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,1))/2;
                }
                else if(_A_LK_el + _C_LK_el >= _B_LK_el){
                    std::cout << "A + B + C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                    << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,2))/2;
                }
                else{
                    std::cout << "2B = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                    << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,2))/2;
                }

                std::cout << "(111)" << std::endl;
                std::cout << "(2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,2) << std::endl;
            }
            else{
                if(_A_LK_el + _C_LK_el >= _B_LK_el){
                    std::cout << "A + B - C = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                    << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,0))/2;
                }
                else if(_A_LK_el - _C_LK_el >= _B_LK_el){
                    std::cout << "A + B - C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                    << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                }
                else{
                    std::cout << "2B = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                    << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                    C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                }

                std::cout << "(111)" << std::endl;
                std::cout << "(2/3)(A + 2B - C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                << ", (2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,2) << std::endl;
            }
        }
        std::cout << "Polaronic Luttinger-Kohn parameters are:" << std::endl;
        std::cout << "A_LK_pol = " << A_LK_pol << ", B_LK_pol = " << B_LK_pol << ", C_LK_pol = " << C_LK_pol << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Inverse effective mass of system is: " << effective_mass_inv << "." << std::endl;

    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not effective_mass.txt open file " << filename << std::endl;
    }
    else{
        file << "Effective mass of the system is: " << _effective_mass << "." << std::endl;
        file << "Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;

        if(_num_bands == 1){
            file << "Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }

        file <<"1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
            << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        file << "Number of phonon modes: " << _num_phonon_modes << std::endl;

        for(int i=0; i<_num_phonon_modes; i++){
            file << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
                << _born_effective_charges[i] << std::endl;
        }

        if(_num_bands == 1){
            file << "Polaronic effective masses are: mx_pol = " << (1./_effective_masses[0]) << " my_pol = " 
                << (1./_effective_masses[1]) << "mz_pol = " << (1./_effective_masses[2]) << std::endl; 
        }
        else if(_num_bands == 3){
            double A_LK_pol = 0;
            double B_LK_pol = 0;
            double C_LK_pol = 0;
            file << "Computed polaronic inverse masses for (100), (110) and (111) high symmetry directions." << std::endl;
            file << "(100)" << std::endl;

            if(_A_LK_el >= _B_LK_el){
                file << "2A = " << _effective_masses_bands(0,0) << ", 2B = " << _effective_masses_bands(0,1) 
                    << ", 2B = " << _effective_masses_bands(0,2) << std::endl;
                A_LK_pol = _effective_masses_bands(0,0)/2;
                B_LK_pol = _effective_masses_bands(0,1)/2;

                file << "(110)" << std::endl;

                if(_C_LK_el >= 0){
                    if(_A_LK_el - _C_LK_el >= _B_LK_el){
                        file << "A + B + C = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                        << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,1))/2;
                    }
                    else{
                        file << "A + B + C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                        << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,2))/2;
                    }
                file << "(111)" << std::endl;
                file << "(2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,2) << std::endl;
                }
                else{
                    if(_A_LK_el + _C_LK_el > _B_LK_el){
                        file << "A + B - C = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                        << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,0))/2;
                    }
                    else{
                        file << "A + B - C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                        << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                    }
                    file << "(111)" << std::endl;
                    file << "(2/3)(A + 2B - C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                    << ", (2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,2) << std::endl;
                }
            }
            else{
                file << "2B = " << _effective_masses_bands(0,0) << ", 2B = " << _effective_masses_bands(0,1) 
                << ", 2A = " << _effective_masses_bands(0,2) << std::endl;
                A_LK_pol = _effective_masses_bands(0,2)/2;
                B_LK_pol = _effective_masses_bands(0,0)/2;

                file << "(110)" << std::endl;
        
                if(_C_LK_el >= 0){
                    if(_A_LK_el - _C_LK_el >= _B_LK_el){
                        file << "A + B + C = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                        << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,1))/2;
                    }
                    else if(_A_LK_el + _C_LK_el >= _B_LK_el){
                        file << "A + B + C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                        << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,0) - _effective_masses_bands(1,2))/2;
                    }
                    else{
                        file << "2B = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                        << ", A + B - C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,2))/2;
                    }

                    file << "(111)" << std::endl;
                    file << "(2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                    << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,2) << std::endl;
                }
                else{
                    if(_A_LK_el + _C_LK_el >= _B_LK_el){
                        file << "A + B - C = " << _effective_masses_bands(1,0) << ", A + B + C = " << _effective_masses_bands(1,1) 
                        << ", 2B = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,1) - _effective_masses_bands(1,0))/2;
                    }
                    else if(_A_LK_el - _C_LK_el >= _B_LK_el){
                        file << "A + B - C = " << _effective_masses_bands(1,0) << ", 2B = " << _effective_masses_bands(1,1) 
                        << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                    }
                    else{
                        file << "2B = " << _effective_masses_bands(1,0) << ", A + B - C = " << _effective_masses_bands(1,1) 
                        << ", A + B + C = " << _effective_masses_bands(1,2) << std::endl;
                        C_LK_pol = (_effective_masses_bands(1,2) - _effective_masses_bands(1,0))/2;
                    }

                    file << "(111)" << std::endl;
                    file << "(2/3)(A + 2B - C) = " << _effective_masses_bands(2,0) << ", (2/3)(A + 2B - C) = " << _effective_masses_bands(2,1)
                    << ", (2/3)(A + 2B + 2C) = " << _effective_masses_bands(2,2) << std::endl;
                }
            }

            file << "Polaronic Luttinger-Kohn parameters are:" << std::endl;
            file << "A_LK_pol = " << A_LK_pol << ", B_LK_pol = " << B_LK_pol << ", C_LK_pol = " << C_LK_pol << std::endl;
        }
        

        file << std::endl;

        file << "Inverse effective mass of the system is: " << effective_mass_inv << "." << std::endl;
        file << std::endl;

        file.close();
    }

    std::cout << std::endl;
};*/

/*{
    //int current_order = _current_order_int + 2*ext_phonon_order; // total order of diagrams (number of phonon vertices)
    //double electron_action = 0, phonon_action = 0;
*/
    // compute electron bare propagators action
    /*for(int i=0; i<current_order+1; i++){
        double electron_energy = electronEnergy(_propagators[i].el_propagator_kx, _propagators[i].el_propagator_ky, 
            _propagators[i].el_propagator_kz, _bands[i].effective_mass);
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
            phonon_action += phononEnergy(_phonon_modes,_vertices[i].index)*(tau_two - tau_one);
            int_count++;
            if(int_count == _current_order_int/2){int_flag = true;}
        }
        else if(_vertices[i].type == -2){
            int index_two = _vertices[i].linked;
            long double tau_one = _vertices[i].tau;  
            long double tau_two = _vertices[index_two].tau;
            phonon_action += phononEnergy(_phonon_modes, _vertices[i].index)*(tau_length + tau_one - tau_two);
            ext_count++;
            if(ext_count == ext_phonon_order){ext_flag = true;}
        }
        i++;
    }

}
}*/

    /*int total_order = _current_order_int + 2*_current_ph_ext;
    double kx = 0, ky = 0, kz = 0;

    long double * new_taus = new long double[total_order+2];
    new_taus[0] = 0.L; // first vertex time value is always 0

    int c = 0;
    int phonon_index = -1;
    double phonon_lines_energies = 0;

    for(int i=1; i < total_order+2; ++i){

        kx = _propagators[i-1].el_propagator_kx;
        ky = _propagators[i-1].el_propagator_ky;
        kz = _propagators[i-1].el_propagator_kz;

        c = _vertices[i-1].type;
        phonon_index = _vertices[i-1].index;
        
        if(c == +1 || c == +2){
            phonon_lines_energies += phononEnergy(_phonon_modes, phonon_index);
        }
        else if( c == -1 || c == -2){
            phonon_lines_energies -= phononEnergy(_phonon_modes, phonon_index);
        }

        new_taus[i] = new_taus[i-1] - std::log(1-drawUniformR())/(electronEnergy(kx,ky,kz,_bands[i-1].effective_mass)
        -_chem_potential + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes) + phonon_lines_energies);

        if(new_taus[i] < new_taus[i-1] || isEqual(new_taus[i], new_taus[i-1])){
            delete[] new_taus; 
            return tau_init;
        }
    }

    if(isEqual(new_taus[total_order+1], _tau_max) || new_taus[total_order+1] > _tau_max){
        delete[] new_taus; 
        return tau_init;
    }

    for(int i = 0; i < total_order+2; ++i){
        _vertices[i].tau = new_taus[i];
    }
    

    delete[] new_taus;
    return _vertices[total_order+1].tau;*/

double GreenFuncNphBands::calcNormConst(){
    double eff_mass_electron = computeEffMassSingleBand(_kx, _ky, _kz, _m_x_el, _m_y_el, _m_z_el);
    double numerator = 1 - std::exp(-(electronEnergy(_kx,_ky,_kz, eff_mass_electron)-_chem_potential)*_tau_max);
    double denominator = electronEnergy(_kx,_ky,_kz,eff_mass_electron)-_chem_potential;
    //std::cout << "Normalization constant calculated." << std::endl;
    return (numerator/denominator);
};

void GreenFuncNphBands::normalizeHistogram(double norm_const){
    for(int i=0; i<_N_bins; i++){
        _green_func[i] = _bin_count[i]/(_N0*_bin_width);
        _green_func[i] = _green_func[i]*norm_const;
    }
    //std::cout << "Imaginary time Green's function computed (histogram method)." << std::endl;
};

void GreenFuncNphBands::writeHistogram(const std::string& filename) const {
    std::ofstream file;

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

void GreenFuncNphBands::writeExactGF(const std::string& filename) const {
    std::ofstream file;

    file.open(filename);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }

    if(_selected_order < 0){
        file << "# Exact Green's function calculated for all orders of external phonons.\n";
    }
    else{
        file << "# Exact Green's function calculated for number of external phonons " << _selected_order << ".\n";
    }
    file << "# kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << ", chemical potential = " << _chem_potential << "\n";

    for(int i=0; i<_num_points; i++){
        file << _points[i] << " " << _points_gf_exact[i] << "\n";
    }
    file.close();
    std::cout << "Exact Green's function written to file " << filename << "." << std::endl;
}

void GreenFuncNphBands::writeDiagram(std::string filename, int i, double r) const {
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Iteration: " << i << "\n";

    file << "index: " << 0 << " time: " << _vertices[0].tau << " wx: " << _vertices[0].wx << " wy: " 
        << _vertices[0].wy << " wz: " << _vertices[0].wz <<" type: " << _vertices[0].type << " linked: " << _vertices[0].linked << "\n";
        file << "\n";
        file << "propagator: " << 0 << "          kx: " << _propagators[0].el_propagator_kx << " ky: " << _propagators[0].el_propagator_ky 
        << " kz: " << _propagators[0].el_propagator_kz << "\n";
        file << "band number: " << _bands[0].band_number << " el eff mass: " << _bands[0].effective_mass 
        << " (c1, c2, c3): (" << _bands[0].c1 << ", " << _bands[0].c2 << ", " << _bands[0].c3 << ")\n";
        file << "\n";

    for(int j=1; j<_current_order_int+2*_current_ph_ext+1; j++){
        file << "index: " << j << " time: " << _vertices[j].tau << " wx: " << _vertices[j].wx << " wy: " 
        << _vertices[j].wy << " wz: " << _vertices[j].wz <<" type: " << _vertices[j].type << " linked: " << _vertices[j].linked << 
        " phonon mode: " << _vertices[j].index << " phonon energy: " << _phonon_modes[_vertices[j].index] << "\n";
        file << "\n";
        file << "propagator: " << j << "          kx: " << _propagators[j].el_propagator_kx << " ky: " << _propagators[j].el_propagator_ky 
        << " kz: " << _propagators[j].el_propagator_kz << "\n";
        file << "band number: " << _bands[j].band_number << " el eff mass: " << _bands[j].effective_mass 
        << " (c1, c2, c3): (" << _bands[j].c1 << ", " << _bands[j].c2 << ", " << _bands[j].c3 << ")\n";
        file << "\n";
    }
    int final_vertex = _current_order_int+2*_current_ph_ext+1;

    file << "index: " << final_vertex << " time: " << _vertices[final_vertex].tau << " wx: " << _vertices[final_vertex].wx 
    << " wy: " << _vertices[final_vertex].wy << " wz: " << _vertices[final_vertex].wz 
    << " type: " << _vertices[final_vertex].type << " linked: " << _vertices[final_vertex].linked 
    <<  "\n";
    file << "\n";

    file << "ext phonons: " << _current_ph_ext << " int order: " << _current_order_int << " chosen update: " << r <<"\n";
    file << "\n";

    file << "\n";
    file << "\n";
    file << std::endl;
    file.close();
};

void GreenFuncNphBands::writeChosenUpdate(std::string filename, int i, double r) const {
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Iteration: " << i << " chosen update: " << r << "\n";
    file << "\n";

    file.close();
};

void GreenFuncNphBands::writeMCStatistics(std::string filename) const {
    std::ofstream file(filename, std::ofstream::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }
    file << "Monte Carlo statistics:\n";
    file << "chemical potential: " << _chem_potential << ", number of degenerate electronic bands: " << _num_bands << "\n";
    if(_num_bands == 1){
        file << "electronic effective masses: mx_el = " << _m_x_el << ", my_el = " << _m_y_el << ", mz_el = " << _m_z_el << "\n";
    }
    else if(_num_bands == 3){
        file << "electronic Luttinger-Kohn parameters: A_LK_el = " << _A_LK_el << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << "\n";
    }
    file << "total momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << "\n";
    file << "\n";
    file << "1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " << _dielectric_const << "\n";
    file << "\n";
    file << "number of phonon modes: " << _num_phonon_modes << "\n";
    for(int i=0; i<_num_phonon_modes; i++){
        file << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << "\n";
    }
    file << "\n";
    file << "Cutoff for statistics: " << _tau_cutoff_statistics << "\n";
    file << "Number of diagrams (taken into account): " << _mc_statistics.num_diagrams << "\n";
    file << "\n";
    file << "Average length of diagrams: " << _mc_statistics.avg_tau << "\n";
    file << "Std dev length of diagrams: " << std::sqrt(_mc_statistics.avg_tau_squared - _mc_statistics.avg_tau*_mc_statistics.avg_tau) << "\n";
    file << "\n";
    file << "Average order of diagrams: " << static_cast<long double>(_mc_statistics.avg_order)/static_cast<long double>(_mc_statistics.num_diagrams) << "\n";
    file << "Std dev order of diagrams: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_order_squared)/static_cast<long double>(_mc_statistics.num_diagrams)
    - static_cast<long double>(_mc_statistics.avg_order*_mc_statistics.avg_order)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << "\n";
    file << "\n";
    file << "Average number of internal phonons: " << static_cast<long double>(_mc_statistics.avg_ph_int)/static_cast<long double>(_mc_statistics.num_diagrams) << "\n";
    file << "Std dev number of internal phonons: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_ph_int_squared)/static_cast<long double>(_mc_statistics.num_diagrams)
    - static_cast<long double>(_mc_statistics.avg_ph_int*_mc_statistics.avg_ph_int)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << "\n";
    file << "\n";
    file << "Average number of external phonons: " << static_cast<long double>(_mc_statistics.avg_ph_ext)/static_cast<long double>(_mc_statistics.num_diagrams) << "\n";
    file << "Std dev number of external phonons: " << std::sqrt(static_cast<long double>(_mc_statistics.avg_ph_ext_squared)/static_cast<long double>(_mc_statistics.num_diagrams)
    - static_cast<long double>(_mc_statistics.avg_ph_ext*_mc_statistics.avg_ph_ext)/static_cast<long double>(_mc_statistics.num_diagrams*_mc_statistics.num_diagrams)) << "\n";
    file << "\n";
    file << "Number of zero order diagrams: " << _mc_statistics.zero_order_diagrams << "\n";
    file << "\n";

    file.close();
    std::cout << "Values' statistics written to file " << filename << "." << std::endl;
};


/*
int GreenFuncNphBands::choosePhonon(){
    double total_phonons=0;

    for(int i = 0; i < _num_phonon_modes; i++){
        total_phonons += 1./_phonon_modes[i];
    }

    double prob_phonon[_num_phonon_modes];
    for(int i=0; i < _num_phonon_modes; i++){
        prob_phonon[i] = (1./_phonon_modes[i])/(total_phonons);
    }

    // choose phonon index
    //std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes-1);
    int phonon_index;
    double r_phonon = drawUniformR();
    double cumulative_prob = 0;

    for(int i=0; i < _num_phonon_modes; i++){
        cumulative_prob += (5.84/5.98);
        if(r_phonon <= cumulative_prob){
            phonon_index = i;
            return phonon_index;
        }
    }
    return 0; // should never reach this point
};
*/