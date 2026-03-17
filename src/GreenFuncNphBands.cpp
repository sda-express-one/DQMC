#include "../include/GreenFuncNphBands.hpp"

// constructor
GreenFuncNphBands::GreenFuncNphBands(unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
double chem_potential, int order_int_max, int ph_ext_max, int data_type = 1, int num_bands = 1, int phonon_modes = 1) 
: Diagram(N_diags, tau_max, kx, ky, kz, chem_potential, order_int_max, ph_ext_max, data_type) {
    
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


    for(int i=0; i < _num_phonon_modes; i++){_ext_phonon_type_num[i] = 0;}

    // initialize support arrays
    _new_taus = new long double[_order_int_max + 2*_ph_ext_max + 2];
    _bands_init = new Band[_order_int_max + 2*_ph_ext_max + 2];
    _bands_fin = new Band[_order_int_max + 2*_ph_ext_max + 2];
    
    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; i++){
        _new_taus[i] = 0.L;
        if(_num_bands == 1){ 
            _bands_init[i].band_number = 0;
            _bands_init[i].c1 = 1;
            _bands_init[i].c2 = 0;
            _bands_init[i].c3 = 0;

            _bands_fin[i].band_number = 0;
            _bands_fin[i].c1 = 1;
            _bands_fin[i].c2 = 0;
            _bands_fin[i].c3 = 0;

            _nodes[i].electronic_band.band_number = 0;
            _nodes[i].electronic_band.c1 = 1;
            _nodes[i].electronic_band.c2 = 0;
            _nodes[i].electronic_band.c3 = 0;
        }
        else if(_num_bands == 3){
            _bands_init[i].band_number = -1;
            _bands_init[i].c1 = (1./3);
            _bands_init[i].c2 = (1./3);
            _bands_init[i].c3 = (1./3);

            _bands_fin[i].band_number = -1;
            _bands_fin[i].c1 = (1./3);
            _bands_fin[i].c2 = (1./3);
            _bands_fin[i].c3 = (1./3);

            _nodes[i].electronic_band.band_number = -1;
            _nodes[i].electronic_band.c1 = (1./3);
            _nodes[i].electronic_band.c2 = (1./3);
            _nodes[i].electronic_band.c3 = (1./3);
        }
    }

    findLastPhVertex();
};

GreenFuncNphBands::GreenFuncNphBands(FullVertexNode* nodes, int current_order,
unsigned long long int N_diags, long double tau_max, double kx, double ky, double kz,
double chem_potential, int order_int_max, int ph_ext_max, int data_type = 1, int num_bands = 1, int phonon_modes = 1) 
: Diagram(nodes, current_order, N_diags, tau_max, kx, ky, kz, chem_potential, order_int_max, ph_ext_max, data_type) {
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
    _bands_init = new Band[_order_int_max + 2*_ph_ext_max + 1];
    _bands_fin = new Band[_order_int_max + 2*_ph_ext_max + 1];

    // initialize support arrays
    _new_taus = new long double[_order_int_max + 2*_ph_ext_max + 2];
    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 2; i++){_new_taus[i] = 0.L;}

    // initialize bands
    _bands = new Band[_order_int_max + 2*_ph_ext_max + 1];

    for(int i = 0; i < _order_int_max + 2*_ph_ext_max + 1; i++){

        if(_num_bands == 1){
            _bands_init[i].band_number = 0;
            _bands_init[i].c1 = 1;
            _bands_init[i].c2 = 0;
            _bands_init[i].c3 = 0;

            _bands_fin[i].band_number = 0;
            _bands_fin[i].c1 = 1;
            _bands_fin[i].c2 = 0;
            _bands_fin[i].c3 = 0;
        }
        else if(_num_bands == 3){
            _bands_init[i].band_number = -1; // unassigned
            _bands_init[i].c1 = (1./3);
            _bands_init[i].c2 = (1./3);
            _bands_init[i].c3 = (1./3);

            _bands_fin[i].band_number = -1;
            _bands_fin[i].c1 = (1./3);
            _bands_fin[i].c2 = (1./3);
            _bands_fin[i].c3 = (1./3);
        }
    }

    findLastPhVertex();
};

void GreenFuncNphBands::getEffectiveMasses(long double * effective_masses) const {
    effective_masses[0] = _effective_masses[0];
    effective_masses[1] = _effective_masses[1];
    effective_masses[2] = _effective_masses[2];
};

void GreenFuncNphBands::getEffectiveMassesVar(long double * effective_masses_var) const {
    effective_masses_var[0] = _effective_masses_var[0];
    effective_masses_var[1] = _effective_masses_var[1];
    effective_masses_var[2] = _effective_masses_var[2];
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
        if(_master){std::cerr << "Invalid probabilities, total probability must add to 1.\n";}
        double normalization = 1/(probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7]);
        if(_master){std::cerr << "Probabilities are being riscaled using the value " << normalization <<".\n";}
        for(int i = 0; i < 8; i++){
            probs[i] = probs[i]*normalization;
        }
        if(_master){std::cerr << "New probabilities are: " << probs[0] << " " << probs[1] << " " << probs[2] << " " 
        << probs[3] << " " << probs[4] << " " << probs[5] << " " << probs[6] << " " << probs[7] << ".\n\n";}
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

void GreenFuncNphBands::setNumNodes(int num_nodes){_num_nodes = num_nodes;};

void GreenFuncNphBands::setNumProcs(int num_procs){_num_procs = num_procs;};

// calculations setters
void GreenFuncNphBands::setCalculations(bool gf_exact, bool histo, bool gs_energy, bool effective_mass, bool Z_factor, bool blocking_analysis, bool fix_tau_value){
    _flags.gf_exact = gf_exact;
    _flags.histo = histo;
    _flags.gs_energy = gs_energy;
    _flags.effective_mass = effective_mass;
    _flags.Z_factor = Z_factor;
    _flags.blocking_analysis = blocking_analysis;
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

void GreenFuncNphBands::setNumBlocks(int N_blocks){if(N_blocks > 0){_N_blocks = N_blocks;}};

// MC statistics setter
void GreenFuncNphBands::setTauCutoffStatistics(long double tau_cutoff_statistics){
    while(tau_cutoff_statistics < 0 || tau_cutoff_statistics >= _tau_max){
        std::cout << "Invalid statistics cutoff! Cutoff must be >= 0 and < " << _tau_max << " (max tau value).\n";
        std::cout << "The default value is set to 0.\n";
        tau_cutoff_statistics = 0;
    }
    _tau_cutoff_statistics = tau_cutoff_statistics;
};

// NOT OKAY???
FullVertexNode * GreenFuncNphBands::findVertexPosition(long double tau){
    // ?????????????????????????
    if(tau > _tail->tau || isEqual(tau, _tail->tau)){return nullptr;}
    else if(tau < _head->tau_next){
        if(!isEqual(0, tau) && !isEqual(tau, _head->tau_next)){return _head;}
        else{return nullptr;}
    }
    else{
        double tau_inf = 0;
        double tau_sup = _tau_max;
        int iterations = std::max(_current_order_int, 2*_current_ph_ext);
        for(int i=0; i<iterations; ++i){
            if(i < _current_order_int){
                tau_inf = _internal_used[i].linked->tau;
                tau_sup = _internal_used[i].linked->tau_next;
                // return pointer if conditions are satisfied
                if(tau > tau_inf && tau < tau_sup && !isEqual(tau, tau_inf) && !isEqual(tau, tau_sup)){return _internal_used[i].linked;}
            }
            if(i < 2*_current_ph_ext){
                tau_inf = _external_used[i].linked->tau;
                tau_sup = _external_used[i].linked->tau_next;
                // return pointer if condition is satisfied
                if(tau > tau_inf && tau < tau_sup && !isEqual(tau, tau_inf) && !isEqual(tau, tau_sup)){return _external_used[i].linked;}
            }
        }
        return nullptr;
    }
    return nullptr;
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

FullVertexNodeIndicator GreenFuncNphBands::chooseInternalPhononPropagator(){
    std::uniform_int_distribution<int> distrib_unif(0,_current_order_int-1); // chooses one of the internal phonon propagators at random
    int ph_vertex_internal = distrib_unif(gen);
    return _internal_used[ph_vertex_internal];
};

FullVertexNodeIndicator GreenFuncNphBands::chooseExternalPhononPropagator(){
    std::uniform_int_distribution<int> distrib_unif(0, 2*_current_ph_ext-1); // chooses one of the external phonon propagators at random
    int ph_vertex_external = distrib_unif(gen);
    return _external_used[ph_vertex_external];
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
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext - 1; ++i){
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
    for(int i = index_one; i < _current_order_int + 2*_current_ph_ext - 1; ++i){
        if(i < index_two - 2){_bands[i] = _bands[i+1];}
        else{_bands[i] = _bands[i+2];}
    }
};

void GreenFuncNphBands::updateExternalPhononTypes(int index){_ext_phonon_type_num[index] += 1;};

long double GreenFuncNphBands::diagramLengthUpdate(long double tau_init){
    // initialize momentum values from last node
    _helper = _tail->prev;
    double kx = _helper->k[0];
    double ky = _helper->k[1];
    double kz = _helper->k[2];
    double effective_mass = _helper->electronic_band.effective_mass;
    _helper = nullptr;

    // generate new time value for last vertex, reject if it goes out of bounds
    long double tau_fin = _last_vertex - std::log(1-drawUniformR())/(electronEnergy(kx, ky, kz, effective_mass)
        -_chem_potential + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes));
    
    if(_num_bands > 1){updateNegativeDiagrams(0);}

    // reject if it goes out of bounds
    if(tau_fin <= _tau_max){
        _tail->tau = tau_fin;
        _tail->prev->tau_next = tau_fin;
        findLastPhVertex();
        return tau_fin;
    }
    else{return tau_init;}
};

void GreenFuncNphBands::addInternalPhononPropagator(){
    if(_current_order_int+1 >= _order_int_max){
        if(_num_bands > 1){updateNegativeDiagrams(1);}
        return; // reject if already at max order
    } 
    else{
        // choose random electron propagator for new vertex
        std::uniform_int_distribution<int> distrib_prop(0, _current_order_int + 2*_current_ph_ext);
        int propagator = distrib_prop(gen);
        long double tau_init = 0;
        long double tau_end = 0;
        
        // find right propagator to insert new proposed tau_{1} value
        // 0 = first propagator, 1 <= x <= order_int_current propagator starting with internal phonon vertex, otherwise propagator with external phonon vertex
        if(propagator > 0){
            int temp = propagator;
            if(propagator < _current_order_int + 1){
                propagator = propagator - 1;
                _pointer_one = findInternalPhononVertex(propagator);
                tau_init = _pointer_one->tau;
                tau_end = _pointer_one->next->tau;
            }
            else{
                propagator -= (_current_order_int + 1);
                _pointer_one = findExternalPhononVertex(propagator);
                tau_init = _pointer_one->tau;
                tau_end = _pointer_one->next->tau;
            }
            propagator = temp;
        }
        else{
            tau_init = _head->tau;
            tau_end = _head->next->tau;
            _pointer_one = _head;
        }

        // choose phonon index
        std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes-1);
        int phonon_index = distrib_phon(gen);

        // choose time value of new vertex (between ends of chosen propagator)
        std::uniform_real_distribution<long double> distrib_unif(tau_init, tau_end);
        long double tau_one = distrib_unif(gen);

        // choose time value of second vertex, different energy for different phonon modes
        long double tau_two = tau_one - std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // may be on a different propagator

        if(tau_two >= _tail->tau){
            if(_num_bands > 1){updateNegativeDiagrams(1);}
            return; // reject if phonon vertex goes out of bound
        }  
        else{
            // sampling momentum values for phonon propagators
            std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_two-tau_one))); // may need to specify phonon mode
            double w_x = distrib_norm(gen); 
            double w_y = distrib_norm(gen);
            double w_z = distrib_norm(gen);
            
            // find position of new tau values
            _pointer_two = findVertexPosition(tau_two);

            // control statements to check for double precision point errors
            if(_pointer_one == nullptr || _pointer_two == nullptr){
                if(_num_bands > 1){updateNegativeDiagrams(1);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return; // reject if tau values are not found in the vertices array
            } 
            if(_current_ph_ext > 0){
                if( tau_one <  _pointer_one->tau || isEqual(tau_one, _pointer_one->tau) 
                    || isEqual(tau_one, _pointer_one->tau_next) || tau_one > _pointer_one->tau_next){
                    if(_num_bands > 1){updateNegativeDiagrams(1);}
                    _pointer_one = nullptr;
                    _pointer_two = nullptr;
                    return;
                }
                if(tau_two < _pointer_two->tau || isEqual(tau_two, _pointer_two->tau) 
                    || isEqual(tau_two, _pointer_two->tau_next) || tau_two > _pointer_two->tau_next){
                    if(_num_bands > 1){updateNegativeDiagrams(1);}
                    _pointer_one = nullptr;
                    _pointer_two = nullptr;
                    return;
                }
            }

            // momentum values
            double px_init = 0;
            double py_init = 0;
            double pz_init = 0;
            double px_fin = 0;
            double py_fin = 0;
            double pz_fin = 0;

            // energy values
            double energy_init = 0;
            double energy_fin = 0;

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> unif(0, _num_bands-1);

            // initial and final action
            double action_init = 0.;
            double action_fin = 0.;

            int i = 0;
            //int current_order = _current_order_int + 2*_current_ph_ext + 1;

            _helper = _pointer_one;

            while(_helper != _pointer_two->next){
                
                // initial diagram momenta
                px_init = _helper->k[0];
                py_init = _helper->k[1];
                pz_init = _helper->k[2];

                // proposed diagram momenta
                px_fin = px_init - w_x;
                py_fin = py_init - w_y;
                pz_fin = pz_init - w_z;

                _bands_init[i] = _helper->electronic_band;

                if(_num_bands == 3){
                    chosen_band = unif(gen);
                    _bands_fin[i].band_number = chosen_band;

                    new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin, _A_LK_el, _B_LK_el, _C_LK_el);
                    eigenval = new_values_matrix(0,chosen_band);

                    // computing new proposed electron effective mass from chosen eigenvalue
                    _bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    // new proposed band eigenstate
                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                    _bands_fin[i].c1 = new_overlap(0); 
                    _bands_fin[i].c2 = new_overlap(1);
                    _bands_fin[i].c3 = new_overlap(2);
                        
                    // compute vertex terms
                    if(_helper == _pointer_one){
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_init[0], new_overlap);
                    }
                    else{
                        prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], _bands_init[i]);
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);                        
                    }
                    if (_helper == _pointer_two){
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_init[i], new_overlap);
                    }
                }
                else if(_num_bands == 1){
                    _bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin,
                                                                        _m_x_el, _m_y_el, _m_z_el);
                }

                energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                if(_pointer_one == _pointer_two){
                    action_init += energy_init*(tau_two-tau_one);
                    action_fin += energy_fin*(tau_two-tau_one);
                }
                else if(_helper == _pointer_one){
                    action_init += energy_init*(_helper->tau_next-tau_one);
                    action_fin += energy_fin*(_helper->tau_next-tau_one);
                }
                else if(_helper == _pointer_two){
                    action_init += energy_init*(tau_two-_helper->tau);
                    action_fin += energy_fin*(tau_two-_helper->tau);
                }
                else{
                    action_init += energy_init*(_helper->tau_next - _helper->tau);
                    action_fin += energy_fin*(_helper->tau_next - _helper->tau);
                }

                _helper = _helper->next;
                ++i;
            }
            _helper = nullptr;

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
            
            // check for sign problem
            if(R_add < 0){
                R_add = std::abs(R_add);
            } 
            
            if(!(Metropolis(R_add))){
                if(_num_bands > 1){updateNegativeDiagrams(1);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            else{
                // update is accepted, build new diagram
                insertNode(_pointer_one);
                insertNode(_pointer_two);
                _pointer_one = _pointer_one->next;
                if(_pointer_two == _pointer_one->prev){_pointer_two = _pointer_one->next;}
                else{_pointer_two = _pointer_two->next;}

                // assign vertex one values
                _pointer_one->tau = tau_one;
                _pointer_one->type = +1;
                _pointer_one->w[0] = w_x;
                _pointer_one->w[1] = w_y;
                _pointer_one->w[2] = w_z;
                _pointer_one->k[0] = _pointer_one->prev->k[0];
                _pointer_one->k[1] = _pointer_one->prev->k[1];
                _pointer_one->k[2] = _pointer_one->prev->k[2];
                _pointer_one->index = phonon_index;

                // assign vertex two values
                _pointer_two->tau = tau_two;
                _pointer_two->type = -1;
                _pointer_two->w[0] = w_x;
                _pointer_two->w[1] = w_y;
                _pointer_two->w[2] = w_z;
                // ???????????????????
                if(_pointer_two->prev != _pointer_one){
                    _pointer_two->k[0] = _pointer_two->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->k[2];
                    _pointer_two->electronic_band = _pointer_two->prev->electronic_band;
                }
                else{
                    _pointer_two->k[0] = _pointer_two->prev->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->prev->k[2];
                    _pointer_two->electronic_band = _pointer_two->prev->prev->electronic_band;
                }
                _pointer_two->index = phonon_index;

                _pointer_one->tau_next = _pointer_one->next->tau;
                _pointer_one->prev->tau_next = tau_one;
                _pointer_two->tau_next = _pointer_two->next->tau;
                _pointer_two->prev->tau_next = tau_two;

                // update electron propagators
                _helper = _pointer_one;
                i = 0;
                while(_helper != _pointer_two){
                    _helper->k[0] -= w_x;
                    _helper->k[1] -= w_y;
                    _helper->k[2] -= w_z;
                    _helper->electronic_band = _bands_fin[i];
                    ++i;
                    _helper = _helper->next;
                }
                _helper = nullptr;

                int current_order_int = _current_order_int;

                // update support array (links)
                _internal_used[current_order_int].linked = _pointer_one;
                _internal_used[current_order_int].used = true;
                _internal_used[current_order_int+1].linked = _pointer_two;
                _internal_used[current_order_int+1].used = true;

                // link two phonon vertices
                _internal_used[current_order_int].conjugated = &_internal_used[current_order_int+1];
                _internal_used[current_order_int+1].conjugated = &_internal_used[current_order_int];

                // update current order and find new last ph vertex
                _current_order_int += 2;
                findLastPhVertex();

                // update sign and negative diagram count if necessary
                if(_num_bands > 1){
                    double ratio = numerator/denominator;
                    if(ratio < 0 && !isEqual(ratio, 0)){
                        updateSign();
                    }
                    updateNegativeDiagrams(1);
                }
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
        } 
    }
};

void GreenFuncNphBands::removeInternalPhononPropagator(){
    if(_current_order_int < 2){
        if(_num_bands > 1){updateNegativeDiagrams(2);}
        return; // reject if already at order 0
    }
    else{
        // indexes of initial and final vertices of a random internal phonon propagator
        _pointer_one = nullptr;
        _pointer_two = nullptr;

        FullVertexNodeIndicator vertex_data = chooseInternalPhononPropagator();
        _pointer_one = vertex_data.linked;
        _pointer_two = vertex_data.conjugated->linked;
        if(_pointer_one == nullptr || _pointer_two == nullptr){
            if(_num_bands > 1){updateNegativeDiagrams(2);}
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
        if(_pointer_one->type == -1){
            _helper = _pointer_one;
            _pointer_one = _pointer_two;
            _pointer_two = _helper;
            _helper = nullptr;
        }

        // vertices' time values
        long double tau_one = _pointer_one->tau;
        long double tau_two = _pointer_two->tau;

        long double tau_init = _pointer_one->prev->tau;
        long double tau_end = _tau_max;
        
        if(_pointer_two != _pointer_one->next){tau_end = _pointer_one->tau_next;}
        else{tau_end = _pointer_two->tau_next;}

        double w_x = _pointer_one->w[0];
        double w_y = _pointer_one->w[1];
        double w_z = _pointer_one->w[2];

        // phonon mode
        int phonon_index = _pointer_one->index;

        // momentum values
        double px_init = 0;
        double py_init = 0;
        double pz_init = 0;
        double px_fin = 0;
        double py_fin = 0;
        double pz_fin = 0;


        // energy values
        double energy_init = 0;
        double energy_fin = 0;

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

        _helper = _pointer_one; // AAAAAAAAAAAAAH DIOCANE
        int i = 0;

        while(_helper != _pointer_two){
            // current momentum values
            px_fin = _helper->k[0];
            py_fin = _helper->k[1];
            pz_fin = _helper->k[2];

            // proposed new momentum values
            px_init = px_fin + w_x;
            py_init = py_fin + w_y;
            pz_init = pz_fin + w_z;
            
            _bands_fin[i] = _helper->electronic_band;

            if(_num_bands == 3){
                chosen_band = band_number(gen);
                new_values_matrix = diagonalizeLKHamiltonian(px_init, py_init, pz_init, _A_LK_el, _B_LK_el, _C_LK_el);
                
                // check if it is a free propagator
                if(new_values_matrix(0,0) == -2){
                    _bands_init[i].band_number = -1;
                    _bands_init[i].effective_mass = 1;
                    _bands_init[i].c1 = (1./3);
                    _bands_init[i].c2 = (1./3);
                    _bands_init[i].c3 = (1./3);

                    new_overlap << (1./3), (1./3), (1./3);
                }
                else{                        
                    eigenval = new_values_matrix(0,chosen_band);

                    // computing new proposed electron effective mass from chosen eigenvalue
                    _bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);

                    // new proposed band eigenstate
                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                    _bands_init[i].c1 = new_overlap(0); 
                    _bands_init[i].c2 = new_overlap(1); 
                    _bands_init[i].c3 = new_overlap(2);
                }

                // compute vertex terms
                if(_pointer_one->next == _pointer_two){
                    prefactor_init = prefactor_init*vertexOverlapTerm(_pointer_one->prev->electronic_band, new_overlap)
                        *vertexOverlapTerm(_pointer_one->next->electronic_band, new_overlap);
                    //################# TO BE CHECKED, something sketchy here diocane
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_pointer_one->prev->electronic_band, _bands_fin[i])
                                    *vertexOverlapTerm(_pointer_one->next->electronic_band, _bands_fin[i]);
                }
                else if(_helper == _pointer_one){
                    prefactor_init = prefactor_init*vertexOverlapTerm(_pointer_one->prev->electronic_band, new_overlap);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_pointer_one->prev->electronic_band, _bands_fin[0]);
                }
                else if(_helper->next == _pointer_two){
                    prefactor_init = prefactor_init*vertexOverlapTerm(_pointer_two->electronic_band, new_overlap);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i], _pointer_two->electronic_band);
                }
                else{
                    prefactor_init = prefactor_init*vertexOverlapTerm(_helper->prev->electronic_band, _helper->electronic_band);
                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                }
            }
            else if(_num_bands == 1){
                _bands_init[i].effective_mass = computeEffMassSingleBand(px_init, py_init, pz_init, 
                                                                    _m_x_el, _m_y_el, _m_z_el);
            }
            energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
            energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

            if(_pointer_one->next == _pointer_two){
                action_init += energy_init*(tau_two-tau_one);
                action_fin += energy_fin*(tau_two-tau_one);
            }
            else if(_helper == _pointer_one){
                action_init += energy_init*(_helper->tau_next-tau_one);
                action_fin += energy_fin*(_helper->tau_next-tau_one);
            }
            else if(_helper->next == _pointer_two){
                action_init += energy_init*(tau_two-_helper->tau);
                action_fin += energy_fin*(tau_two-_helper->tau);
            }
            else{
                action_init += energy_init*(_helper->tau_next - _helper->tau);
                action_fin += energy_fin*(_helper->tau_next - _helper->tau);
            }
            _helper = _helper->next;
            ++i;
        }
        _helper = nullptr;

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

        // check for sign problem
        if(R_rem < 0){
            R_rem = std::abs(R_rem);
        } 

        if(!(Metropolis(R_rem))){
            if(_num_bands > 1){updateNegativeDiagrams(2);}
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
        else{
            // possibly needs a fix
            int index_one = vertex_data.position;
            int index_two = vertex_data.conjugated->position;

            // move array pointers in suitable positions
            _internal_used[index_one].linked = _internal_used[_current_order_int-2].linked;
            _internal_used[index_two].linked = _internal_used[_current_order_int-1].linked;
            _internal_used[index_one].conjugated = &_internal_used[index_two];
            _internal_used[index_two].conjugated = &_internal_used[index_one];

            // free suitable array positions
            _internal_used[_current_order_int - 2].used = false;
            _internal_used[_current_order_int - 2].linked = nullptr;
            _internal_used[_current_order_int - 2].conjugated = nullptr;
            _internal_used[_current_order_int - 1].used = false;
            _internal_used[_current_order_int - 1].linked = nullptr;
            _internal_used[_current_order_int - 1].conjugated = nullptr;

            // update is accepted, build new diagram
            _helper = _pointer_one;
            i = 0;
            while(_helper != _pointer_two){
                _helper->k[0] += w_x;
                _helper->k[1] += w_y;
                _helper->k[2] += w_z;
                _helper->electronic_band = _bands_init[i];
                _helper = _helper->next;
                ++i;
            }
            _helper = nullptr;

            deleteNode(_pointer_one);
            deleteNode(_pointer_two);
            if(_pointer_one == _pointer_two){_pointer_one->tau_next = _pointer_one->next->tau;}
            else{
                _pointer_one->tau_next = _pointer_one->next->tau;
                _pointer_two->tau_next = _pointer_two->next->tau;
            }
            
            if(_pointer_one != _pointer_two){}

            _helper = _free_list;
            i = 0;
            while(i<2){
                _helper->k[0] = 0;
                _helper->k[1] = 0;
                _helper->k[2] = 0;
                _helper->electronic_band.effective_mass = 1.0;
                if(_num_bands == 3){
                    _helper->electronic_band.band_number = -1;
                    _helper->electronic_band.c1 = (1./3);
                    _helper->electronic_band.c2 = (1./3);
                    _helper->electronic_band.c3 = (1./3);
                }
                _helper = _helper->next;
                ++i;
            }
            _helper = nullptr;

            


            // update current order and find new last ph vertex
            _current_order_int -= 2;
            findLastPhVertex();

            // update sign and negative diagram count if necessary
            if(_num_bands > 1){
                double ratio = numerator/denominator;
                if(ratio < 0 && !isEqual(ratio, 0)){
                    updateSign();
                }
                updateNegativeDiagrams(2);
            }
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
    }
};

void GreenFuncNphBands::addExternalPhononPropagator(){
    if(_current_ph_ext >= _ph_ext_max){
        if(_num_bands > 1){updateNegativeDiagrams(3);}
        return; // return if already at max number of ext phonon propagators
    } 
    else{
        _pointer_one = nullptr;
        _pointer_two = nullptr;
        int total_order = _current_order_int + 2*_current_ph_ext;
        long double tau_current = _tail->tau; // length of current diagram

        // choose phonon index
        std::uniform_int_distribution<int> distrib_phon(0, _num_phonon_modes-1);
        int phonon_index = distrib_phon(gen);
        //int phonon_index = choosePhonon();

        // time of ingoing vertex of ext phonon propagator
        long double tau_one = 0 - std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // time of ingoing vertex of ext phonon propagator
        if(isEqual(tau_one, tau_current) || tau_one >= tau_current){
            if(_num_bands > 1){updateNegativeDiagrams(3);}
            return; // reject if it goes out of bound
        }

        long double tau_two = tau_current + std::log(1-drawUniformR())/phononEnergy(_phonon_modes, phonon_index); // time of outgoing vertex

        if(isEqual(tau_two,0) || tau_two <= 0){
            if(_num_bands > 1){updateNegativeDiagrams(3);}
            return; // reject if it goes out of bound
        } 

        if(isEqual(tau_one, tau_two)){
            if(_num_bands > 1){updateNegativeDiagrams(3);}
            return; // reject if both vertices are equal (should not happen)
        } 
        
        // sampling momentum values for phonon propagators
        std::normal_distribution<double> distrib_norm(0, std::sqrt(1/(tau_current-tau_two+tau_one)));
        double w_x = distrib_norm(gen);
        double w_y = distrib_norm(gen);
        double w_z = distrib_norm(gen);

        if(tau_one <= tau_two){

            _pointer_one = findVertexPosition(tau_one);
            _pointer_two = findVertexPosition(tau_two);

            // control statements to check for floating point errors
            if(_pointer_one == nullptr || _pointer_two == nullptr){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return; // reject if tau values are not found in the vertices array
            } 
            if(tau_one < _pointer_one->tau || isEqual(tau_one, _pointer_one->tau) 
                || isEqual(tau_one, _pointer_one->tau_next) || tau_one > _pointer_one->tau_next){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            if(tau_two < _pointer_two->tau || isEqual(tau_two, _pointer_two->tau) 
                || isEqual(tau_two, _pointer_two->tau_next) || tau_two > _pointer_two->tau_next){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            } 

            // momentum values
            double px_init = 0;
            double px_fin = 0;
            double py_init = 0;
            double py_fin = 0;
            double pz_init = 0;
            double pz_fin = 0;

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> band_number(0, _num_bands-1);

            // energy values
            double energy_init = 0;
            double energy_fin = 0;

            // initial and final action
            double action_one_init = 0.; 
            double action_two_init = 0.;
            double action_one_fin = 0.;
            double action_two_fin = 0.;

            int i = 0;
            bool first_part = false, second_part = true;
            _helper = _head;

            while(_helper->next != nullptr){ 
                while(!first_part){
                    // retrieve momentum values for propagators below first ph vertex
                    px_init = _helper->k[0];
                    px_fin = px_init - w_x;

                    py_init = _helper->k[1];
                    py_fin = py_init - w_y;

                    pz_init = _helper->k[2];
                    pz_fin = pz_init - w_z;

                    _bands_init[i] = _helper->electronic_band;

                    if(_num_bands == 3){
                        chosen_band = band_number(gen);
                        _bands_fin[i].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin, _A_LK_el, _B_LK_el, _C_LK_el);    
                        eigenval = new_values_matrix(0,chosen_band);
                
                        // computing new proposed electron effective mass from chosen eigenvalue
                        _bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                        // new proposed band eigenstate
                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        _bands_fin[i].c1 = new_overlap(0); 
                        _bands_fin[i].c2 = new_overlap(1);
                        _bands_fin[i].c3 = new_overlap(2);

                        if(_helper != _head){
                            // compute vertex terms
                            prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], _bands_init[i]);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], new_overlap);
                        }
                    }
                    else if(_num_bands == 1){
                        _bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin,
                                                                                _m_x_el, _m_y_el, _m_z_el);   
                    }

                    // compute energies
                    energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    // compute action
                    if(_helper != _pointer_one){
                        action_one_init += energy_init*(_helper->tau_next - _helper->tau);
                        action_one_fin += energy_fin*(_helper->tau_next - _helper->tau);
                    }
                    else{
                        action_one_init += energy_init*(tau_one - _helper->tau);
                        action_one_fin += energy_fin*(tau_one - _helper->tau);
                        first_part = true;
                        second_part = false;
                    }
                    _helper = _helper->next;
                    ++i;
                }

                // new vertex (left)
                if(_num_bands == 3){
                     prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _pointer_one->electronic_band);
                }
                _helper = _pointer_two;

                while(!second_part){
                    // retrieve momentum values for propagators above second ph vertex
                    px_init = _helper->k[0];
                    px_fin = px_init - w_x;

                    py_init = _helper->k[1];
                    py_fin = py_init - w_y;

                    pz_init = _helper->k[2];
                    pz_fin = pz_init - w_z;

                    _bands_init[i] = _helper->electronic_band;

                    if(_num_bands == 3){
                        if(_helper->next != _tail){
                            chosen_band = band_number(gen);
                            _bands_fin[i].band_number = chosen_band;

                            new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin, 
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            eigenval = new_values_matrix(0,chosen_band);

                            // computing new proposed electron effective mass from chosen eigenvalue
                            _bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                            // new proposed band eigenstate
                            new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                            _bands_fin[i].c1 = new_overlap(0); 
                            _bands_fin[i].c2 = new_overlap(1); 
                            _bands_fin[i].c3 = new_overlap(2);
                        }
                        else{
                            // last propagator must be the same as the first one (conservation of four-momentum)
                            _bands_fin[i] = _bands_fin[0];
                            new_overlap(0) = _bands_fin[0].c1;
                            new_overlap(1) = _bands_fin[0].c2;
                            new_overlap(2) = _bands_fin[0].c3;
                        }

                        // compute vertex terms
                        if(_helper != _pointer_two){
                            prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], _bands_init[i]);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], new_overlap);
                        }
                        else{
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_init[i], new_overlap);
                        }
                    }
                    else if(_num_bands == 1){
                        _bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin, 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    // calc energy values for propagators above second ph vertex
                    energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);
                
                    if(_helper == _pointer_two){
                        action_two_init += energy_init*(_helper->tau_next - tau_two);
                        action_two_fin += energy_fin*(_helper->tau_next - tau_two);
                    }
                    else{
                        action_two_init += energy_init*(_helper->tau_next - _helper->tau);
                        action_two_fin += energy_fin*(_helper->tau_next - _helper->tau);
                    }
                    if(_helper->next == _tail){second_part = true;}

                    _helper = _helper->next;
                    ++i;
                }
            }
            
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

            if(R_add < 0){
                R_add = std::abs(R_add);
            } 

            if(!(Metropolis(R_add))){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            else{
                // update is accepted, build new diagram
                insertNode(_pointer_one);
                insertNode(_pointer_two);
                _pointer_one = _pointer_one->next;
                if(_pointer_two == _pointer_one->prev){_pointer_two = _pointer_one->next;}
                else{_pointer_two = _pointer_two->next;}

                // assign vertex one values
                _pointer_one->tau = tau_one;
                _pointer_one->type = -2;
                _pointer_one->index = phonon_index;
                _pointer_one->w[0] = w_x;
                _pointer_one->w[1] = w_y;
                _pointer_one->w[2] = w_z;
                _pointer_one->k[0] = _pointer_one->prev->k[0];
                _pointer_one->k[1] = _pointer_one->prev->k[1]; // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                _pointer_one->k[2] = _pointer_one->prev->k[2];
                _pointer_one->electronic_band = _pointer_one->prev->electronic_band;

                // assign vertex two values
                _pointer_two->tau = tau_two;
                _pointer_two->type = +2;
                _pointer_two->index = phonon_index;
                _pointer_two->w[0] = w_x;
                _pointer_two->w[1] = w_y;
                _pointer_two->w[2] = w_z;
                if(_pointer_two->prev != _pointer_one){
                    _pointer_two->k[0] = _pointer_two->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->k[2];
                }
                else{
                    _pointer_two->k[0] = _pointer_two->prev->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->prev->k[2];
                }

                _pointer_one->prev->tau_next = tau_one;
                _pointer_one->tau_next = _pointer_one->next->tau;
                _pointer_two->prev->tau_next = tau_two;
                _pointer_two->tau_next = _pointer_two->next->tau;

                _helper = _head;
                int i = 0;
                while(_helper != _pointer_one){
                    _helper->k[0] -= w_x;
                    _helper->k[1] -= w_y;
                    _helper->k[2] -= w_z;
                    _helper->electronic_band = _bands_fin[i];
                    ++i;
                    _helper = _helper->next;
                }
                _helper = _pointer_two;
                while(_helper->next != nullptr){
                    _helper->k[0] -= w_x;
                    _helper->k[1] -= w_y;
                    _helper->k[2] -= w_z;
                    _helper->electronic_band = _bands_fin[i];
                    ++i;
                    _helper = _helper->next;
                }
                _helper = nullptr;

                _external_used[2*_current_ph_ext].linked = _pointer_one;
                _external_used[2*_current_ph_ext].used = true;
                _external_used[2*_current_ph_ext + 1].linked = _pointer_two;
                _external_used[2*_current_ph_ext + 1].used = true;
                _external_used[2*_current_ph_ext].conjugated = &_external_used[2*_current_ph_ext + 1];
                _external_used[2*_current_ph_ext + 1].conjugated = &_external_used[2*_current_ph_ext];

                _pointer_one = nullptr;
                _pointer_two = nullptr;

                _current_ph_ext += 1; // update current number of external phonons
                _ext_phonon_type_num[phonon_index]++;
                findLastPhVertex();

                // update sign and negative diagram count if necessary
                if(_num_bands > 1){
                    double ratio = numerator/denominator;
                    if(ratio < 0 && !isEqual(ratio, 0)){
                        updateSign();
                    }
                    updateNegativeDiagrams(3);
                }
                return;
            }
        }
        else{
            _pointer_one = findVertexPosition(tau_two);
            _pointer_two = findVertexPosition(tau_one);

            // control statements to check for floating point errors
            if(_pointer_one == nullptr || _pointer_two == nullptr){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return; // reject if tau values are not found in the vertices array
            }
            if( tau_two < _pointer_one->tau || isEqual(tau_two, _pointer_one->tau) 
                || isEqual(tau_two, _pointer_one->tau_next) || tau_two > _pointer_one->tau_next){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            if(tau_one < _pointer_two->tau || isEqual(tau_one, _pointer_two->tau) 
                || isEqual(tau_one, _pointer_two->tau_next) || tau_one > _pointer_two->tau_next){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }

            // momentum values
            double px_init = 0;
            double px_fin = 0;
            double py_init = 0;
            double py_fin = 0;
            double pz_init = 0;
            double pz_fin = 0;

            // temporary variables for new proposed diagram
            int chosen_band = 0;
            double eigenval = 1.0;
            Eigen::Matrix<double,4,3> new_values_matrix;
            Eigen::Vector3d new_overlap;
            std::uniform_int_distribution<int> band_number(0, _num_bands-1);

            // energy values
            double energy_init = 0;
            double energy_fin = 0;

            double action_init = 0.;
            double action_fin = 0.;

            // vertices weights of two diagrams
            double prefactor_fin = 1;
            double prefactor_init = 1;
            bool first_part = false, second_part = true, third_part = true;

            _helper = _head;
            int i = 0;

            while(_helper->next != nullptr){
                // initial diagram
                px_init = _helper->k[0];
                py_init = _helper->k[1];
                pz_init = _helper->k[2];

                _bands_init[i] = _helper->electronic_band;

                if(_num_bands == 3 && _helper != _head){    
                    prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], _bands_init[i]);
                }

                energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                action_init += energy_init*(_helper->tau_next - _helper->tau);

                // final diagram
                if(!first_part){
                    px_fin = px_init - w_x;
                    py_fin = py_init - w_y;
                    pz_fin = pz_init - w_z;

                    if(_num_bands == 3){
                        chosen_band = band_number(gen);
                        _bands_fin[i].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin,
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);
                        _bands_fin[i].effective_mass = computeEffMassfromEigenval(eigenval);

                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        _bands_fin[i].c1 = new_overlap(0); 
                        _bands_fin[i].c2 = new_overlap(1); 
                        _bands_fin[i].c3 = new_overlap(2);

                        if(_helper != _head){
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], new_overlap);
                        }
                    }
                    else if(_num_bands == 1){
                        _bands_fin[i].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin,
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    if(_helper != _pointer_one){
                        action_fin += energy_fin*(_helper->tau_next - _helper->tau);
                    }
                    else{
                        action_fin += energy_fin*(tau_two - _helper->tau);
                        first_part = true;
                        second_part = false;
                    }
                }
                if(!second_part){
                    px_fin = px_init - 2*w_x;
                    py_fin = py_init - 2*w_y;
                    pz_fin = pz_init - 2*w_z;

                    if(_num_bands == 3){
                        chosen_band = band_number(gen);
                        _bands_fin[i+1].band_number = chosen_band;

                        new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin,
                                                                _A_LK_el, _B_LK_el, _C_LK_el);
                        eigenval = new_values_matrix(0,chosen_band);
                        _bands_fin[i+1].effective_mass = computeEffMassfromEigenval(eigenval);

                        new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                        _bands_fin[i+1].c1 = new_overlap(0); 
                        _bands_fin[i+1].c2 = new_overlap(1); 
                        _bands_fin[i+1].c3 = new_overlap(2);

                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i], new_overlap);
                    }
                    else if(_num_bands == 1){
                        _bands_fin[i+1].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin, 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i+1].effective_mass);

                    if(_helper == _pointer_one && _helper == _pointer_two){
                        action_fin += energy_fin*(tau_one - tau_two);
                        second_part = true;
                        third_part = false;
                    }
                    else if(_helper == _pointer_one){
                        action_fin += energy_fin*(_helper->tau_next - tau_two);
                    }
                    else if(_helper == _pointer_two){
                        action_fin += energy_fin*(tau_one - _helper->tau);
                        second_part = true;
                        third_part = false;
                    }
                    else{
                        action_fin += energy_fin*(_helper->tau_next - _helper->tau);
                    }
                }
                if(!third_part){
                    px_fin = px_init - w_x;
                    py_fin = py_init - w_y;
                    pz_fin = pz_init - w_z;

                    if(_num_bands == 3){
                        if(_helper->next != _tail){
                            chosen_band = band_number(gen);
                            _bands_fin[i+2].band_number = chosen_band;

                            new_values_matrix = diagonalizeLKHamiltonian(px_fin, py_fin, pz_fin,
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            eigenval = new_values_matrix(0,chosen_band);
                            _bands_fin[i+2].effective_mass = computeEffMassfromEigenval(eigenval);

                            new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                            _bands_fin[i+2].c1 = new_overlap(0); 
                            _bands_fin[i+2].c2 = new_overlap(1); 
                            _bands_fin[i+2].c3 = new_overlap(2);
                        }
                        else{
                            _bands_fin[total_order+2] = _bands_fin[0];
                        }
                        prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i+1], new_overlap);
                    }
                    else if(_num_bands == 1){
                        _bands_fin[i+2].effective_mass = computeEffMassSingleBand(px_fin, py_fin, pz_fin, 
                                                                            _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    if(_helper == _pointer_two){
                        action_fin += energy_fin*(_helper->tau_next - tau_one);
                    }
                    else{
                        action_fin += energy_fin*(_helper->tau_next - _helper->tau);
                    }

                    if(_helper->next == _tail){third_part = true;}
                }
                _helper = _helper->next;
                ++i;
            }
            _helper = nullptr;

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

            if(R_add < 0){
                R_add = std::abs(R_add);
            } 

            if(!(Metropolis(R_add))){
                if(_num_bands > 1){updateNegativeDiagrams(3);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            else{
                //phVertexMakeRoom(index_one, index_two); // make room in vertices array
                //propagatorArrayMakeRoom(index_one, index_two); // make room in electron propagators array

                // insert new nodes
                insertNode(_pointer_one);
                insertNode(_pointer_two);
                _pointer_one = _pointer_one->next;
                if(_pointer_two == _pointer_one->prev){_pointer_two = _pointer_one->next;}
                else{_pointer_two = _pointer_two->next;}

                // assign vertex one values
                _pointer_one->tau = tau_two;
                _pointer_one->type = +2;
                _pointer_one->w[0] = w_x;
                _pointer_one->w[1] = w_y;
                _pointer_one->w[2] = w_z;
                _pointer_one->k[0] = _pointer_one->prev->k[0];
                _pointer_one->k[1] = _pointer_one->prev->k[1];
                _pointer_one->k[2] = _pointer_one->prev->k[2];
                _pointer_one->index = phonon_index;

                // assign vertex two values
                _pointer_two->tau = tau_one;
                _pointer_two->type = -2;
                _pointer_two->w[0] = w_x;
                _pointer_two->w[1] = w_y;
                _pointer_two->w[2] = w_z;
                if(_pointer_two->prev != _pointer_one){
                    _pointer_two->k[0] = _pointer_two->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->k[2];
                }
                else{
                    _pointer_two->k[0] = _pointer_two->prev->prev->k[0];
                    _pointer_two->k[1] = _pointer_two->prev->prev->k[1];
                    _pointer_two->k[2] = _pointer_two->prev->prev->k[2];
                }
                _pointer_two->index = phonon_index;

                _pointer_one->tau_next =_pointer_one->next->tau;
                _pointer_one->prev->tau_next = tau_two;
                _pointer_two->tau_next = _pointer_two->next->tau;
                _pointer_two->prev->tau_next = tau_one;

                

                // update electron propagator energies
                _helper = _head;
                i = 0;
                //bool external = false;
                while(_helper->next != nullptr){
                    if(_helper == _pointer_one || (_helper->tau > _pointer_one->tau && _helper->tau < _pointer_two->tau && _helper != _pointer_two)){
                        _helper->k[0] -= 2*w_x;
                        _helper->k[1] -= 2*w_y;
                        _helper->k[2] -= 2*w_z;
                    }
                    else{
                        _helper->k[0] -= w_x;
                        _helper->k[1] -= w_y;
                        _helper->k[2] -= w_z;
                    }
                    _helper->electronic_band = _bands_fin[i];

                    _helper = _helper->next;
                    ++i;
                }
                _helper = nullptr;

                // update phonon vertices counts
                _external_used[2*_current_ph_ext].linked = _pointer_one;
                _external_used[2*_current_ph_ext].used = true;
                _external_used[2*_current_ph_ext+1].linked = _pointer_two;
                _external_used[2*_current_ph_ext+1].used = true;
                // link two phonons to each other
                _external_used[2*_current_ph_ext].conjugated = &_external_used[2*_current_ph_ext+1];
                _external_used[2*_current_ph_ext+1].conjugated = &_external_used[2*_current_ph_ext];
                
                _pointer_one = nullptr;
                _pointer_two = nullptr;

                _current_ph_ext += 1; // update current number of external phonons
                _ext_phonon_type_num[phonon_index]++;
                findLastPhVertex();

                if(_num_bands > 1){
                    double ratio = numerator/denominator;
                    if(ratio < 0 && !isEqual(ratio, 0)){
                        updateSign();
                    }
                    updateNegativeDiagrams(3);
                }
                return;
            }
        }
    }
};

void GreenFuncNphBands::removeExternalPhononPropagator(){
    if(_current_ph_ext <= 0){
        if(_num_bands > 1){updateNegativeDiagrams(4);}
        return; // reject if already at order 0
    }
    else{
        _pointer_one = nullptr;
        _pointer_two = nullptr;
        // indexes of initial and final vertices of a random internal phonon propagator
        FullVertexNodeIndicator vertex_data = chooseExternalPhononPropagator();
        
        if(vertex_data.linked == nullptr){
            if(_num_bands > 1){updateNegativeDiagrams(4);}
            return;
        }

        _pointer_one = vertex_data.linked;
        _pointer_two = vertex_data.conjugated->linked;

        if(isEqual(_pointer_one->tau, _pointer_two->tau)){
            _pointer_one = nullptr; _pointer_two = nullptr;
            if(_num_bands > 1){updateNegativeDiagrams(4);}
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }

        if(_pointer_one->tau > _pointer_two->tau){
            FullVertexNode * temp;
            temp = _pointer_one;
            _pointer_one = _pointer_two;
            _pointer_two = temp;
        }

        // retrieve time values of the two vertices
        long double tau_one = _pointer_one->tau;
        long double tau_two = _pointer_two->tau;

        // retrieve phonon momentum
        double w_x = _pointer_one->w[0];
        double w_y = _pointer_one->w[1];
        double w_z = _pointer_one->w[2];

        // retrieve phonon mode
        int phonon_index = _pointer_one->index;


        // vertices weights of two diagrams
        double prefactor_fin = 1;
        double prefactor_init = 1;

        // temporary variables for new proposed diagram
        int chosen_band = 0;
        double eigenval = 1.0;
        Eigen::Matrix<double,4,3> new_values_matrix;
        Eigen::Vector3d new_overlap;
        new_overlap << 1,0,0;
        std::uniform_int_distribution<int> band_number(0, _num_bands-1);

        if(_pointer_one->type == -2){

            // momentum values
            double px_init = 0;
            double px_fin = 0;
            double py_init = 0;
            double py_fin = 0;
            double pz_init = 0;
            double pz_fin = 0;

            // energy values
            double energy_init = 0;
            double energy_fin = 0;

            // initial and final action
            double action_one_init = 0.;
            double action_one_fin = 0.;
            double action_two_init = 0.;
            double action_two_fin = 0.;

            _helper = _head;
            int i = 0;

            while(_helper->next != nullptr){
                if(_helper->tau < _pointer_one->tau && _helper != _pointer_one){
                    px_fin = _helper->k[0];
                    py_fin = _helper->k[1];
                    pz_fin = _helper->k[2];

                    _bands_fin[i] = _helper->electronic_band;

                    px_init = px_fin + w_x;
                    py_init = py_fin + w_y;
                    pz_init = pz_fin + w_z;

                    if(_num_bands == 3){
                        if(_helper->next != _pointer_one){

                            chosen_band = band_number(gen);
                            new_values_matrix = diagonalizeLKHamiltonian(px_init, py_init, pz_init,
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            // check if it is a free propagator
                            if(new_values_matrix(0,0) == -2){
                                _bands_init[i].band_number = -1;
                                _bands_init[i].effective_mass = 1;
                                _bands_init[i].c1 = (1./3);
                                _bands_init[i].c2 = (1./3);
                                _bands_init[i].c3 = (1./3);

                                new_overlap << (1./3), (1./3), (1./3);
                            }
                            else{
                                _bands_init[i].band_number = chosen_band;

                                eigenval = new_values_matrix(0,chosen_band);
                                _bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                _bands_init[i].c1 = new_overlap(0); 
                                _bands_init[i].c2 = new_overlap(1); 
                                _bands_init[i].c3 = new_overlap(2);
                            }
                            if(_helper != _head){
                                prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], new_overlap);
                                prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                            }
                        }
                        else{
                            _bands_init[i] = _pointer_one->electronic_band;
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]); 
                        }
                    }
                    else if(_num_bands == 1){
                        _bands_init[i].effective_mass = computeEffMassSingleBand(px_init, py_init, pz_init, 
                                                                                _m_x_el, _m_y_el, _m_z_el);
                    }

                    // compute energies
                    energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    // compute actions
                    action_one_init += energy_init*(_helper->tau_next - _helper->tau);
                    action_one_fin += energy_fin*(_helper->tau_next - _helper->tau);
                }

                
                if(_helper->tau > _pointer_two->tau || _helper == _pointer_two){

                    px_fin = _helper->k[0];
                    py_fin = _helper->k[1];
                    pz_fin = _helper->k[2];

                    px_init = px_fin + w_x;
                    py_init = py_fin + w_y;
                    pz_init = pz_fin + w_z;
                    
                    _bands_fin[i] = _helper->electronic_band;

                    if(_num_bands == 3){
                        if(_helper != _pointer_two && _helper->next != _tail /*i != index_two && i != total_order*/){
                            chosen_band = band_number(gen);

                            new_values_matrix = diagonalizeLKHamiltonian(px_init, py_init, pz_init,
                                                                    _A_LK_el, _B_LK_el, _C_LK_el);
                            if(new_values_matrix(0,0) == -2){
                                _bands_init[i].band_number = -1;
                                _bands_init[i].effective_mass = 1;
                                _bands_init[i].c1 = (1./3);
                                _bands_init[i].c2 = (1./3);
                                _bands_init[i].c3 = (1./3);

                                new_overlap << (1./3), (1./3), (1./3);
                            }
                            else{
                                _bands_init[i].band_number = chosen_band;
                            
                                eigenval = new_values_matrix(0,chosen_band);
                                _bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                _bands_init[i].c1 = new_overlap(0); 
                                _bands_init[i].c2 = new_overlap(1); 
                                _bands_init[i].c3 = new_overlap(2);
                            }

                            prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], new_overlap);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                        }
                        else if(_helper == _pointer_two){
                            _bands_init[i] = _pointer_two->prev->electronic_band;
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_pointer_two->prev->electronic_band, _bands_fin[i]);
                        }
                        else{
                            _bands_init[i] = _bands_init[0];

                            prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i], _bands_init[i-1]);
                            prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i], _bands_fin[i-1]);
                        }
                    }
                    else if(_num_bands == 1){
                        _bands_init[i].effective_mass = computeEffMassSingleBand(px_init, py_init, pz_init,
                                                                                _m_x_el, _m_y_el, _m_z_el);
                    }

                    energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    action_two_init += energy_init*(_helper->tau_next - _helper->tau);
                    action_two_fin += energy_fin*(_helper->tau_next - _helper->tau);
                }

                _helper = _helper->next;
                ++i;
                if(_helper == _pointer_one){_helper = _pointer_two;}
            }
            _helper = nullptr;

            long double tau_current = _tail->tau; // length of current diagram

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

            // check for sign problems
            if(R_rem < 0){R_rem = std::abs(R_rem);} 

            if(!Metropolis(R_rem)){
                if(_num_bands > 1){updateNegativeDiagrams(4);}
                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
            else{
                _helper = _head;
                i = 0;
                while(_helper->next != nullptr){
                    if(_helper->tau < _pointer_one->tau_next && _helper != _pointer_one->next){
                        _helper->k[0] += w_x;
                        _helper->k[1] += w_y;
                        _helper->k[2] += w_z;
                        _helper->electronic_band = _bands_init[i];
                    }
                    if(_helper->tau > _pointer_two->tau || _helper == _pointer_two){
                        _helper->k[0] += w_x;
                        _helper->k[1] += w_y;
                        _helper->k[2] += w_z;
                        _helper->electronic_band = _bands_init[i];
                    }
                    _helper = _helper->next;
                    ++i;
                    if(_helper == _pointer_one){_helper = _pointer_two;}
                }
                _helper = nullptr;

                deleteNode(_pointer_one);
                deleteNode(_pointer_two);

                _pointer_one->tau_next = _pointer_one->next->tau;
                if(_pointer_one != _pointer_two){_pointer_two->tau_next = _pointer_two->next->tau;}

                _helper = _free_list;
                i = 0;
                while(i<2){
                    _helper->k[0] = 0;
                    _helper->k[1] = 0;
                    _helper->k[2] = 0;
                    _helper->electronic_band.effective_mass = 1.0;
                    if(_num_bands == 3){
                        _helper->electronic_band.band_number = -1;
                        _helper->electronic_band.c1 = (1./3);
                        _helper->electronic_band.c2 = (1./3);
                        _helper->electronic_band.c3 = (1./3);
                    }
                    _helper = _helper->next;
                    ++i;
                }
                _helper = nullptr;

                int index_one = vertex_data.position;
                int index_two = vertex_data.conjugated->position;

                // move array pointers in suitable positions
                _external_used[index_one].linked = _external_used[2*_current_ph_ext-2].linked;
                _external_used[index_two].linked = _external_used[2*_current_ph_ext-1].linked;
                _external_used[index_one].conjugated = &_external_used[index_two];
                _external_used[index_two].conjugated = &_external_used[index_one];

                // free suitable array positions
                _external_used[2*_current_ph_ext - 2].used = false;
                _external_used[2*_current_ph_ext - 2].linked = nullptr;
                _external_used[2*_current_ph_ext - 2].conjugated = nullptr;
                _external_used[2*_current_ph_ext - 1].used = false;
                _external_used[2*_current_ph_ext - 1].linked = nullptr;
                _external_used[2*_current_ph_ext - 1].conjugated = nullptr;

                _current_ph_ext -= 1; // update current number of external phonons
                _ext_phonon_type_num[phonon_index]--;
                findLastPhVertex();

                // update sign and negative diagram count if necessary
                if(_num_bands > 1){
                    double ratio = numerator/denominator;
                    if(ratio < 0 && !isEqual(ratio, 0)){
                        updateSign();
                    }
                    updateNegativeDiagrams(4);
                }

                _pointer_one = nullptr;
                _pointer_two = nullptr;
                return;
            }
        }
        else if(_pointer_one->type == 2){

                // momentum values
                double px_init = 0;
                double px_fin = 0;
                double py_init = 0;
                double py_fin = 0;
                double pz_init = 0;
                double pz_fin = 0;

                // energy values
                double energy_init = 0;
                double energy_fin = 0;

                double action_init = 0.;
                double action_fin = 0.;

                _helper = _head;
                int i = 0;
                bool external = true;

                while(_helper != _tail){
                    px_fin = _helper->k[0];
                    py_fin = _helper->k[1];
                    pz_fin = _helper->k[2];

                    _bands_fin[i] = _helper->electronic_band;

                    if(external){
                        px_init = px_fin + w_x;
                        py_init = py_fin + w_y;
                        pz_init = pz_fin + w_z;

                        if(_num_bands == 3){
                            // first and last electron propagators must be the same
                            if(_helper->next == _tail){
                                _bands_init[i] = _bands_init[0];
                            }
                            else if(_helper == _pointer_two){
                                _bands_init[i] = _bands_init[i-1];
                            }
                            else{
                                chosen_band = band_number(gen);
                                new_values_matrix = diagonalizeLKHamiltonian(px_init, py_init, pz_init,
                                                                        _A_LK_el, _B_LK_el, _C_LK_el);
                                if(new_values_matrix(0,0) == -2){
                                    _bands_init[i].band_number = -1;
                                    _bands_init[i].effective_mass = 1;
                                    _bands_init[i].c1 = (1./3);
                                    _bands_init[i].c2 = (1./3);
                                    _bands_init[i].c3 = (1./3);

                                    new_overlap << (1./3), (1./3), (1./3);
                                }
                                else{
                                    _bands_init[i].band_number = chosen_band;
                            
                                    eigenval = new_values_matrix(0,chosen_band);
                                    _bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                    _bands_init[i].c1 = new_overlap(0); 
                                    _bands_init[i].c2 = new_overlap(1); 
                                    _bands_init[i].c3 = new_overlap(2);
                                }
                                
                                if(_helper->next == _pointer_one || _helper == _pointer_two){
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                                }
                                else if(_helper != _head){
                                    prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], new_overlap);
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                                }
                            }

                        }
                        else if(_num_bands == 1){
                            _bands_init[i].effective_mass = computeEffMassSingleBand(px_init, py_init, pz_init, 
                                                                                    _m_x_el, _m_y_el, _m_z_el);
                        }

                        if(_helper->next == _pointer_one){external = false;}
                    }
                    else{
                        px_init = px_fin + 2*w_x;
                        py_init = py_fin + 2*w_y;
                        pz_init = pz_fin + 2*w_z;

                        if(_num_bands == 3){
                            if(_helper == _pointer_one){
                                _bands_init[i] = _bands_init[i-1];
                            }
                            else{
                                chosen_band = band_number(gen);
                                new_values_matrix = diagonalizeLKHamiltonian(px_init, py_init, pz_init,
                                                                            _A_LK_el, _B_LK_el, _C_LK_el);
                                if(new_values_matrix(0,0) == -2){
                                    _bands_init[i].band_number = -1;
                                    _bands_init[i].effective_mass = 1;
                                    _bands_init[i].c1 = (1./3);
                                    _bands_init[i].c2 = (1./3);
                                    _bands_init[i].c3 = (1./3);

                                    new_overlap << (1./3), (1./3), (1./3);
                                }
                                else{
                                    _bands_init[i].band_number = chosen_band;
                            
                                    eigenval = new_values_matrix(0,chosen_band);
                                    _bands_init[i].effective_mass = computeEffMassfromEigenval(eigenval);
                    
                                    new_overlap = new_values_matrix.block<3,1>(1,chosen_band);
                                    _bands_init[i].c1 = new_overlap(0); 
                                    _bands_init[i].c2 = new_overlap(1); 
                                    _bands_init[i].c3 = new_overlap(2);
                                }
                                
                                if(_helper != _pointer_one){
                                    prefactor_init = prefactor_init*vertexOverlapTerm(_bands_init[i-1], new_overlap);
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                                }
                                else{
                                    prefactor_fin = prefactor_fin*vertexOverlapTerm(_bands_fin[i-1], _bands_fin[i]);
                                }
                            }
                        }  
                        else if(_num_bands == 1){
                            _bands_init[i].effective_mass = computeEffMassSingleBand(px_init, py_init, pz_init, 
                                                                                _m_x_el, _m_y_el, _m_z_el);
                        }
                        if(_helper->next == _pointer_two){external = true;}
                    }
                    energy_init = electronEnergy(px_init, py_init, pz_init, _bands_init[i].effective_mass);
                    energy_fin = electronEnergy(px_fin, py_fin, pz_fin, _bands_fin[i].effective_mass);

                    action_init += energy_init*(_helper->tau_next - _helper->tau);
                    action_fin += energy_fin*(_helper->tau_next - _helper->tau);

                    _helper = _helper->next;
                    ++i;
                }
                _helper = nullptr;

                long double tau_current = _tail->tau; // length of current diagram

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
                
                // check for sign problems
                if(R_rem < 0){R_rem = std::abs(R_rem);}   
                    
                if(!(Metropolis(R_rem))){
                    if(_num_bands > 1){updateNegativeDiagrams(4);}
                    _pointer_one = nullptr;
                    _pointer_two = nullptr;
                    return;
                }
                else{
                    external = true;
                    _helper = _head;
                    i = 0;
                    while(_helper->next != nullptr){
                        if(external){
                            _helper->k[0] += w_x;
                            _helper->k[1] += w_y;
                            _helper->k[2] += w_z;
                        }
                        else{
                            _helper->k[0] += 2*w_x;
                            _helper->k[1] += 2*w_y;
                            _helper->k[2] += 2*w_z;
                        }
                        _helper->electronic_band = _bands_init[i];

                        if(_helper->next == _pointer_one){external = false;}
                        else if(_helper->next == _pointer_two){external = true;}

                        _helper = _helper->next;
                        ++i;
                    }
                    _helper = nullptr;

                    deleteNode(_pointer_one);
                    deleteNode(_pointer_two);

                    _pointer_one->tau_next = _pointer_one->next->tau;
                    if(_pointer_one != _pointer_two){_pointer_two->tau_next = _pointer_two->next->tau;}

                    _helper = _free_list;
                    i = 0;
                    while(i<2){
                        _helper->k[0] = 0;
                        _helper->k[1] = 0;
                        _helper->k[2] = 0;
                        _helper->electronic_band.effective_mass = 1.0;
                        if(_num_bands == 3){
                            _helper->electronic_band.band_number = -1;
                            _helper->electronic_band.c1 = (1./3);
                            _helper->electronic_band.c2 = (1./3);
                            _helper->electronic_band.c3 = (1./3);
                        }
                        _helper = _helper->next;
                        ++i;
                    }
                    _helper = nullptr;

                    int index_one = vertex_data.position;
                    int index_two = vertex_data.conjugated->position;

                    // move array pointers in suitable positions
                    _external_used[index_one].linked = _external_used[2*_current_ph_ext-2].linked;
                    _external_used[index_two].linked = _external_used[2*_current_ph_ext-1].linked;
                    _external_used[index_one].conjugated = &_external_used[index_two];
                    _external_used[index_two].conjugated = &_external_used[index_one];

                    // free suitable array positions
                    _external_used[2*_current_ph_ext - 2].used = false;
                    _external_used[2*_current_ph_ext - 2].linked = nullptr;
                    _external_used[2*_current_ph_ext - 2].conjugated = nullptr;
                    _external_used[2*_current_ph_ext - 1].used = false;
                    _external_used[2*_current_ph_ext - 1].linked = nullptr;
                    _external_used[2*_current_ph_ext - 1].conjugated = nullptr;
                    
                    _current_ph_ext -= 1; // update current number of external phonons
                    _ext_phonon_type_num[phonon_index]--;
                    findLastPhVertex();

                    // update sign and negative diagram count if necessary
                    if(_num_bands > 1){
                        double ratio = numerator/denominator;
                        if(ratio < 0 && !isEqual(ratio, 0)){
                            updateSign();
                        }
                        updateNegativeDiagrams(4);
                    }
                    _pointer_one = nullptr;
                    _pointer_two = nullptr;
                    return;
                }            
        }
        else{
            if(_num_bands > 1){updateNegativeDiagrams(4);}
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
    }
};

void GreenFuncNphBands::swapPhononPropagator(){
    if(_current_order_int < 4){
        if(_num_bands > 1){updateNegativeDiagrams(5);}
        return;  // swap not possible if internal order is less than 4
    }
    else{
        FullVertexNodeIndicator data_vertex = chooseInternalPhononPropagator();

        int index_one = data_vertex.position;
        int index_two = data_vertex.conjugated->position;

        if(data_vertex.linked->next == _internal_used[index_two].linked){
            if(_num_bands > 1){updateNegativeDiagrams(5);}
            return; // reject if the two vertices are linked
        }

        _pointer_one = data_vertex.linked;
        _pointer_two = _pointer_one->next;

        if(_pointer_two->type != 1 || _pointer_two->type != -1){
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            if(_num_bands > 1){updateNegativeDiagrams(5);}
            return; // reject if they are not both internal vertices
        }

        // get values of first vertex
        int v1 = _pointer_one->type;
        long double wx1 = _pointer_one->w[0];
        long double wy1 = _pointer_one->w[1];
        long double wz1 = _pointer_one->w[2];
        long double tau1 = _pointer_one->tau;
        int phonon_index1 = _pointer_one->index;

        // get values of second vertex
        int v2 = _pointer_two->type;
        long double wx2 = _pointer_two->w[0];
        long double wy2 = _pointer_two->w[1];
        long double wz2 = _pointer_two->w[2];
        long double tau2 = _pointer_two->tau;
        int phonon_index2 = _pointer_two->index;

        // get momentum of propagator
        long double kx = _pointer_one->k[0];
        long double ky = _pointer_one->k[1];
        long double kz = _pointer_one->k[2];

        // retrieve value of initial electron effective mass and propose new value of final effective mass
        double eff_mass_el_initial, eff_mass_el_final;
        eff_mass_el_initial = _pointer_one->electronic_band.effective_mass;
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
                std::uniform_int_distribution<int> distrib_unif(0,_num_bands-1);
                chosen_band = distrib_unif(gen);

                double eigenval = new_values_matrix(0,chosen_band);
                eff_mass_el_final = computeEffMassfromEigenval(eigenval); // computing new proposed electron effective mass from chosen eigenvalue
                new_overlap = new_values_matrix.block<3,1>(1,chosen_band); // new proposed band eigenstate
            }

            // compute only overlap term for vertices that go into R_swap evaluation,
            // the strength term does not change
            prefactor_fin = vertexOverlapTerm(_pointer_one->prev->electronic_band, new_overlap)*vertexOverlapTerm(_pointer_two->electronic_band, new_overlap);
            prefactor_init = vertexOverlapTerm(_pointer_one->prev->electronic_band, _pointer_one->electronic_band)*vertexOverlapTerm(_pointer_two->electronic_band, _pointer_one->electronic_band);
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

        // check for sign problems
        double ratio = R_swap;
        if(R_swap < 0){
            R_swap = std::abs(R_swap);
        } 

        if(!(Metropolis(R_swap))){
            if(_num_bands > 1){updateNegativeDiagrams(5);}
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
        else{
            // assign new momentum values to propagator
            _pointer_one->k[0] += v1*wx1-v2*wx2;
            _pointer_one->k[1] += v1*wy1-v2*wy2;
            _pointer_one->k[2] += v1*wz1-v2*wz2;
            
            // assign new band values
            _pointer_one->electronic_band.effective_mass = eff_mass_el_final;
            if(_num_bands == 3){
                _pointer_one->electronic_band.band_number = chosen_band;
                _pointer_one->electronic_band.c1 = new_overlap[0];
                _pointer_one->electronic_band.c2 = new_overlap[1];
                _pointer_one->electronic_band.c3 = new_overlap[2];
            }

            FullVertexNodeIndicator * indicator = _internal_used[index_one].conjugated;
            _internal_used[index_one].conjugated = _internal_used[index_two].conjugated;
            _internal_used[index_two].conjugated = indicator;

            _pointer_one->w[0] = wx2;
            _pointer_one->w[1] = wy2;
            _pointer_one->w[2] = wz2;
            _pointer_one->type = v2;
            _pointer_one->tau = tau1;
            _pointer_one->index = phonon_index2;
            
            _pointer_two->w[0] = wx1;
            _pointer_two->w[1] = wy1;
            _pointer_two->w[2] = wz1;
            _pointer_two->type = v1;
            _pointer_two->tau = tau2;
            _pointer_two->index = phonon_index1;

            if(_num_bands > 1){
                if(ratio < 0 && !isEqual(ratio, 0)){
                    updateSign();
                }
                updateNegativeDiagrams(5);
            }
            _pointer_one = nullptr;
            _pointer_two = nullptr;
            return;
        }
    }
};

void GreenFuncNphBands::shiftPhononPropagator(){
    int total_order = _current_order_int + 2*_current_ph_ext;
    if(total_order <= 0){
        if(_num_bands > 1){updateNegativeDiagrams(6);}
        return; // reject if no vertices are present
    }
    else{
        std::uniform_int_distribution<int> distrib(0, total_order-1);
        int vertex_index = distrib(gen); // choose random vertex

        if(vertex_index < _current_order_int){
            _pointer_one = findInternalPhononVertex(vertex_index);
        }
        else{
            _pointer_one = findExternalPhononVertex(vertex_index - _current_order_int);
        }

        // necessary step to address phonon type into evaluation
        int c = _pointer_one->type;
        if(c == 2){c = 1;}
        else if(c == -2){c = -1;}
        int phonon_index = _pointer_one->index;  // index of phonon mode in el-phonon vertex

        long double tau_init = _pointer_one->prev->tau;
        long double tau_fin = _pointer_one->tau_next;

        // incoming electron momentum
        long double kx_incoming = _pointer_one->prev->k[0];
        long double ky_incoming = _pointer_one->prev->k[1];
        long double kz_incoming = _pointer_one->prev->k[2];
        double el_eff_mass_incoming = _pointer_one->prev->electronic_band.effective_mass;

        // outgoing electron momentum
        long double kx_outgoing = _pointer_one->k[0];
        long double ky_outgoing = _pointer_one->k[1];
        long double kz_outgoing = _pointer_one->k[2];
        double el_eff_mass_outgoing = _pointer_one->electronic_band.effective_mass;

        double energy_delta = electronEnergy(kx_incoming, ky_incoming, kz_incoming, el_eff_mass_incoming) 
            - electronEnergy(kx_outgoing, ky_outgoing, kz_outgoing, el_eff_mass_outgoing) - phononEnergy(_phonon_modes, phonon_index)*c;
        
        long double tau_new = tau_init - std::log(1 - drawUniformR()*(1 - std::exp(-energy_delta*(tau_fin - tau_init))))/energy_delta;

        // check for possible double precision errors
        if(isEqual(tau_new, tau_init) || isEqual(tau_new, tau_fin) || tau_new < tau_init || tau_new > tau_fin){
            updateNegativeDiagrams(6);
            _pointer_one = nullptr;
            return;
        } 
        
        _pointer_one->tau = tau_new; // assign new time value to vertex
        _pointer_one->prev->tau_next = tau_new;

        findLastPhVertex();

        if(_num_bands > 1){updateNegativeDiagrams(6);}

        _pointer_one = nullptr;
        return;
    }
};

long double GreenFuncNphBands::stretchDiagramLength(long double tau_init){
    int total_order = _current_order_int + 2*_current_ph_ext;
    double kx = 0, ky = 0, kz = 0;
    double effective_mass = 1;

    _new_taus[0] = 0; // first vertex time value is always 0

    int c = 0;
    int phonon_index = -1;
    double phonon_lines_energies = 0;

    _helper = _head;
    int i = 1;

    while(_helper != _tail){
        kx = _helper->k[0];
        ky = _helper->k[1];
        kz = _helper->k[2];
        effective_mass = _helper->electronic_band.effective_mass;

        c = _helper->type;
        phonon_index = _helper->index;
        
        if(c == +1 || c == +2){
            phonon_lines_energies += phononEnergy(_phonon_modes, phonon_index);
        }
        else if( c == -1 || c == -2){
            phonon_lines_energies -= phononEnergy(_phonon_modes, phonon_index);
        }

        _new_taus[i] = _new_taus[i-1] - std::log(1-drawUniformR())/(electronEnergy(kx,ky,kz,effective_mass) - _chem_potential 
            + extPhononEnergy(_ext_phonon_type_num, _phonon_modes, _num_phonon_modes)  + phonon_lines_energies);

        if(_new_taus[i] < _new_taus[i-1] || isEqual(_new_taus[i], _new_taus[i-1])){
            if(_num_bands > 1){updateNegativeDiagrams(7);}
            _helper = nullptr;
            return tau_init;
        }
        _helper = _helper->next;
        ++i;
    }
    _helper = nullptr;

    if(isEqual(_new_taus[total_order+1], _tau_max) || _new_taus[total_order+1] > _tau_max){
        if(_num_bands > 1){updateNegativeDiagrams(7);}
        _helper = nullptr;
        return tau_init;
    }

    _helper = _head;
    i = 0;
    while(_helper != nullptr){
        _helper->tau = _new_taus[i];
        if(_helper != _tail){_helper->tau_next = _new_taus[i+1];}
        _helper = _helper->next;
        ++i;
    }

    findLastPhVertex();

    if(_num_bands > 1){updateNegativeDiagrams(7);}
    _helper = nullptr;
    return _tail->tau;
};

long double GreenFuncNphBands::configSimulation(long double tau_length = 1.0L){

    if((!(isEqual(_kx,0)) || !(isEqual(_ky,0)) || !(isEqual(_kz,0))) && _flags.effective_mass){
        std::cerr << "Warning: kx, ky and kz should be equal to 0 to calculate effective mass." << std::endl;
        std::cerr << "Effective mass calculation is not possible." << std::endl;
        _flags.effective_mass = false;
    }

    if(_flags.Z_factor && _ph_ext_max <= 0){
        std::cerr << "Warning: number of maximum external phonon must be greater than 0 to calculate Z factor." << std::endl;
        std::cerr << "Z factor calculation is not possible." << std::endl;
        _flags.Z_factor = false;
    }

    // print simulation parameters
    std::cout <<"Starting simulation..." << std::endl;
    std::cout << std::endl;
    if(_master){
        std::cout << "Employing parallelized version of the program" << std::endl;
        std::cout << "Number of nodes employed: " << _num_nodes << std::endl;
        std::cout << "Number of parallelized processes (cpus) per node: " << _num_procs  << std::endl;
        std::cout << std::endl;
    }
    std::cout << "Number of thermalization steps: " << getRelaxSteps() << std::endl;
    if(_master){std::cout << "Number of steps to perform to avoid correlations between different parallel processes: " << getAutocorrSteps() << std::endl;}
    std::cout << "Number of diagrams to be generated: " << getNdiags() << std::endl;
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

        if(_flags.blocking_analysis){
            _gs_energy_block_array = new long double[_N_blocks];

            for(int i=0; i<_N_blocks; ++i){
                _gs_energy_block_array[i] = 0;
            }
        }
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

        if(_flags.blocking_analysis){
            _effective_mass_block_array = new long double[_N_blocks];
            _effective_masses_block_array = new long double[3*_N_blocks];

            for(int i=0; i<_N_blocks; ++i){
                _effective_mass_block_array[i] = 0;
                _effective_masses_block_array[3*i] = 0;
                _effective_masses_block_array[3*i+1] = 0;
                _effective_masses_block_array[3*i+2] = 0;
            }
        }
    }

    if(_flags.Z_factor){
        std::cout << "Z factor will be calculated using the exact estimator" << std::endl;
        std::cout << "Free electron momentum: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
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
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_Z << std::endl;
        std::cout << "Number of phonon modes: " << _num_phonon_modes << std::endl;

        for(int i=0; i<_num_phonon_modes; i++){
            std::cout << "phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
            << _dielectric_responses[i] << std::endl;
        }
        std::cout << std::endl;
        _Z_factor_array = new int[_ph_ext_max+1];
        for(int i=0; i<_ph_ext_max+1; i++){_Z_factor_array[i] = 0;}
    }

    if(_flags.blocking_analysis){
        std::cout << "Blocking analysis method will be employed to compute variance of quantities." << std::endl;
        if(getNdiags() % _N_blocks != 0){
            _N_blocks = _N_blocks - 1;
            std::cerr << "Warning! Size of each block non compatible with the total number of computed diagrams." << std::endl;
            std::cerr << "Blocking analysis will be performed only on the first " << _N_blocks << "blocks." << std::endl;
        }
        std::cout << "Number of blocks computed: " << _N_blocks << std::endl;
        _block_size = static_cast<long long int>(getNdiags()/_N_blocks);
        std::cout << "Size of each block is: " << _block_size << "." << std::endl;
    }

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
        _tail->tau = tau_length; // assign fixed time value to first (last) vertex
        _tail->prev->tau_next = tau_length;
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

    if(_flags.blocking_analysis){
        _block_size = static_cast<long long int>(getNdiags()/_N_blocks);

        if(_flags.gs_energy){
            _gs_energy_block_array = new long double[_N_blocks];
            for(int i=0; i<_N_blocks; ++i){
                _gs_energy_block_array[i] = 0;
            }
        }
    }

    if(_flags.effective_mass){
        if(_num_bands == 3){
            _effective_masses_bands << 0, 0, 0,
                                       0, 0, 0,
                                       0, 0, 0;
        }
        if(_flags.blocking_analysis){
            _effective_mass_block_array = new long double[_N_blocks];
            _effective_masses_block_array = new long double[3*_N_blocks];

            for(int i=0; i<_N_blocks; ++i){
                _effective_mass_block_array[i] = 0;
                _effective_masses_block_array[3*i] = 0;
                _effective_masses_block_array[3*i+1] = 0;
                _effective_masses_block_array[3*i+2] = 0;
            }
        }
    }

    if(_flags.Z_factor){
        _Z_factor_array = new int[_ph_ext_max+1];
        for(int i=0; i<_ph_ext_max+1; i++){_Z_factor_array[i] = 0;}
    }

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
            //std::cout << std::endl;
        }
    }
};

long double GreenFuncNphBands::chooseUpdate(long double tau_length, double r , MC_Benchmarking * benchmark){

    if(_flags.time_benchmark){
        if(r <= _p_length){
            benchmark->startUpdateTimer();
            tau_length = diagramLengthUpdate(tau_length);
            benchmark->stopUpdateTimer(0);
        }
        else if(r <= _p_length + _p_add_int){
            benchmark->startUpdateTimer();
            addInternalPhononPropagator();
            benchmark->stopUpdateTimer(1);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int){
            benchmark->startUpdateTimer();
            removeInternalPhononPropagator();
            benchmark->stopUpdateTimer(2);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext){
            benchmark->startUpdateTimer();
            addExternalPhononPropagator();
            benchmark->stopUpdateTimer(3);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext){
            benchmark->startUpdateTimer();
            removeExternalPhononPropagator();
            benchmark->stopUpdateTimer(4);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap){
            benchmark->startUpdateTimer();
            swapPhononPropagator();
            benchmark->stopUpdateTimer(5);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift){
            benchmark->startUpdateTimer();
            shiftPhononPropagator();
            benchmark->stopUpdateTimer(6);
        }
        else if(r <= _p_length + _p_add_int + _p_rem_int + _p_add_ext + _p_rem_ext + _p_swap + _p_shift + _p_stretch){
            benchmark->startUpdateTimer();
            tau_length = stretchDiagramLength(tau_length);
            benchmark->stopUpdateTimer(7);
        }
    }
    else{
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

    if(_flags.gs_energy){_gs_energy += groundStateEnergyExactEstimator(tau_length); } // accumulate energy of diagrams

    if(_flags.effective_mass){_effective_mass += effectiveMassExactEstimator(tau_length);} // accumulate effective mass of diagrams

    if(_flags.Z_factor){ZFactorExactEstimator(tau_length);} // accumulate Z factor data

    if(_flags.write_diagrams){writeDiagram("Diagrams.txt", i, r);} // method to visualize diagram structure

    if(_current_order_int == 0 && _current_ph_ext == 0){_N0++;} // count number of order 0 diagrams for normalization

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
    // method used to compute final quantities at the end of the simulation in the pararellized version of the program (no output to console)

    if(_num_bands > 1){computeRatioNegativeUpdates(getNdiags());} // compute ratio of negative updates to total updates for final normalization of quantities

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
        if(_num_bands > 1){_gs_energy = _gs_energy/((1-_ratio_negative_updates)-_ratio_negative_updates);}
        if(_flags.blocking_analysis){
            long double squared_sum = 0;
            int count = 0;

            // compute variance using blocking analysis procedure
            for(int i = 0; i < _N_blocks; ++i){
                if(!isEqual(_gs_energy_block_array[i],0)){
                    squared_sum += (_gs_energy_block_array[i]-_gs_energy)*(_gs_energy_block_array[i]-_gs_energy);
                    ++count;
                }
            }
            _gs_energy_var = static_cast<long double>(1.0)/(count*(count-1))*squared_sum;
        }
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

        if(_flags.blocking_analysis){
            long double squared_sum_avg = 0;
            long double squared_sum_xP = 0;
            long double squared_sum_yP = 0;
            long double squared_sum_zP = 0;
            int count = 0;

            // compute variance using blocking analysis procedure
            for(int i = 0; i < _N_blocks; ++i){
                if(!isEqual(_gs_energy_block_array[i],0)){
                    squared_sum_avg += (_effective_mass_block_array[i]-_effective_mass)*(_effective_mass_block_array[i]-_effective_mass);
                    squared_sum_xP += (_effective_masses_block_array[3*i]-_effective_masses[0])*(_effective_masses_block_array[3*i]-_effective_masses[0]);
                    squared_sum_yP += (_effective_masses_block_array[3*i+1]-_effective_masses[1])*(_effective_masses_block_array[3*i+1]-_effective_masses[1]);
                    squared_sum_zP += (_effective_masses_block_array[3*i+2]-_effective_masses[2])*(_effective_masses_block_array[3*i+2]-_effective_masses[2]);
                    ++count;
                }
            }
            _effective_mass_var = static_cast<long double>(1.L)/(count*(count-1))*squared_sum_avg;
            _effective_masses_var[0] = static_cast<long double>(1.L)/(count*(count-1))*squared_sum_xP;
            _effective_masses_var[1] = static_cast<long double>(1.L)/(count*(count-1))*squared_sum_yP;
            _effective_masses_var[2] = static_cast<long double>(1.L)/(count*(count-1))*squared_sum_zP;
        }
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
    if(_num_bands > 1){_gs_energy = _gs_energy/((1-_ratio_negative_updates)-_ratio_negative_updates);}

    if(_flags.blocking_analysis){
        long double squared_sum = 0;
        int count = 0;

        // compute variance using blocking analysis procedure
        for(int i = 0; i < _N_blocks; ++i){
            if(!isEqual(_gs_energy_block_array[i],0)){
                squared_sum += (_gs_energy_block_array[i]-_gs_energy)*(_gs_energy_block_array[i]-_gs_energy);
                ++count;
            }
        }
        _gs_energy_var = static_cast<long double>(1.0)/(count*(count-1))*squared_sum;
    }

    std::cout << "Ground state energy of the system is: " << _gs_energy;
    if(_flags.blocking_analysis){std::cout << " +\\- "  << _gs_energy_var;}
    std::cout << "." << std::endl;

    if(_flags.blocking_analysis){std::cout << "Number of blocks used for blocking analysis: " << _N_blocks << std::endl;}

    std::cout << "Input parameters are: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
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
        file << "# Ground state energy of the system is: " << _gs_energy;
        if(_flags.blocking_analysis){file << " +\\- " << _gs_energy_var;}
        file << "." << std::endl;

        if(_flags.blocking_analysis){file << "# Number of blocks used for blocking analysis: " << _N_blocks << std::endl;}

        file << "# Input parameters are: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
        file << "# Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        file << "# Minimum length of diagrams for which gs energy is computed = " << _tau_cutoff_energy << "." << std::endl;
        file << "# Number of diagrams used for ground state energy calculation: " << _gs_energy_count << std::endl;

        if(_num_bands == 1){
            file << "# Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        file << "# 1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_energy << std::endl;
        file << std::endl;

        file << "# Number of phonon modes: " << _num_phonon_modes << std::endl;
        for(int i=0; i<_num_phonon_modes; i++){
            file << "# phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
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

    if(_flags.blocking_analysis){std::cout << "Number of blocks used for blocking analysis: " << _N_blocks << std::endl;}
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

    if(_flags.blocking_analysis){
        long double squared_sum_avg = 0;
        long double squared_sum_xP = 0;

        long double squared_sum_yP = 0;
        long double squared_sum_zP = 0;
        int count = 0;

        // compute variance using blocking analysis procedure
        for(int i = 0; i < _N_blocks; ++i){
            if(!isEqual(_gs_energy_block_array[i],0)){
                squared_sum_avg += (_effective_mass_block_array[i]-_effective_mass)*(_effective_mass_block_array[i]-_effective_mass);
                squared_sum_xP += (_effective_masses_block_array[3*i]-_effective_masses[0])*(_effective_masses_block_array[3*i]-_effective_masses[0]);
                squared_sum_yP += (_effective_masses_block_array[3*i+1]-_effective_masses[1])*(_effective_masses_block_array[3*i+1]-_effective_masses[1]);
                squared_sum_zP += (_effective_masses_block_array[3*i+2]-_effective_masses[2])*(_effective_masses_block_array[3*i+2]-_effective_masses[2]);
                ++count;
            }
        }

        
        _effective_mass_var = static_cast<long double>(1.0)/(count*(count-1))*squared_sum_avg;
        _effective_masses_var[0] = static_cast<long double>(1.0)/(count*(count-1))*squared_sum_xP;
        _effective_masses_var[1] = static_cast<long double>(1.0)/(count*(count-1))*squared_sum_yP;
        _effective_masses_var[2] = static_cast<long double>(1.0)/(count*(count-1))*squared_sum_zP;
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
        std::cout << "Polaronic effective masses are: mx_pol = " << (_effective_masses[0]); 
        if(_flags.blocking_analysis){std::cout << " +\\- " << std::sqrt(_effective_masses_var[0]);}
        std::cout << ", my_pol = " << _effective_masses[1];
        if(_flags.blocking_analysis){std::cout << " +\\- " << std::sqrt(_effective_masses_var[1]);}
        std::cout << ", mz_pol = " << _effective_masses[2];
        if(_flags.blocking_analysis){std::cout << " +\\- " << std::sqrt(_effective_masses_var[2]);}
        std::cout << "." << std::endl; 
    }
    else if(_num_bands == 3){
        // stuff
    }
    std::cout << std::endl;
    std::cout << "Average effective mass of diagrams is: " << _effective_mass;
    if(_flags.blocking_analysis){std::cout << " +\\- " << _effective_mass_var;}
    std::cout << "." << std::endl;
    
    std::cout << "Average inverse effective mass of system is: " << effective_mass_inv << "." << std::endl;
    std::cout << std::endl;

    std::string filename = "effective_mass.txt";
    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not effective_mass.txt open file " << filename << std::endl;
    }
    else{
        file << "# Average effective mass of the system is: " << _effective_mass; 
        if(_flags.blocking_analysis){file << " +\\- " << _effective_mass_var;}
        file << "." << std::endl;

        if(_flags.blocking_analysis){file << "# Number of blocks used for blocking analysis: " << _N_blocks << std::endl;}

        file << "# Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
        file << "# Number of diagrams used for effective mass calculation: " << _effective_mass_count << std::endl;
        if(_num_bands == 1){
            file << "# Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
                << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
        }
        else if(_num_bands == 3){
            file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
                    << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
        }
        file <<"# 1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
        << _dielectric_const << ", tau cutoff: " << _tau_cutoff_mass << std::endl;
        file << std::endl;

        file << "# Number of phonon modes: " << _num_phonon_modes << std::endl;
        for(int i=0; i<_num_phonon_modes; i++){
            file << "# phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
                << _dielectric_responses[i] << std::endl;
        }
        file << std::endl;

        if(_num_bands == 1){
            file << "# Polaronic effective masses are: mx_pol = " << _effective_masses[0];
            if(_flags.blocking_analysis){file << " +\\- " << _effective_masses_var[0];}
            file << ", my_pol = " << _effective_masses[1];
            if(_flags.blocking_analysis){file << " +\\- " << _effective_masses_var[1];}
            file << ", mz_pol = " << _effective_masses[2];
            if(_flags.blocking_analysis){file << " +\\- " << _effective_masses_var[2];}
            file << "." << std::endl; 
        }
       else if(_num_bands == 3){
            // stuff
        }
        file << std::endl;

        file.close();
    }
};

void GreenFuncNphBands::printZFactor(){
    std::cout << "Z factor calculated." << std::endl;
    std::string filename = "quasiparticle_weights.txt";
    writeZFactor(filename);
    std::cout << "Z factor written to file " << filename << "." << std::endl;
    std::cout << std::endl;
};

void GreenFuncNphBands::printMCStatistics(){
    _mc_statistics.avg_tau /= static_cast<long double>(_mc_statistics.num_diagrams); // average length of diagrams
    _mc_statistics.avg_tau_squared /= static_cast<long double>(_mc_statistics.num_diagrams); // average squared length of diagrams
    long double avg_order = static_cast<long double>(_mc_statistics.avg_order) / static_cast<long double>(_mc_statistics.num_diagrams); // average order of diagrams
    long double avg_order_sq = static_cast<long double>(_mc_statistics.avg_order_squared) / static_cast<long double>(_mc_statistics.num_diagrams); // average squared order of diagrams
    long double avg_ph_ext = static_cast<long double>(_mc_statistics.avg_ph_ext) / static_cast<long double>(_mc_statistics.num_diagrams); // average number of external phonons
    long double avg_ph_ext_sq = static_cast<long double>(_mc_statistics.avg_ph_ext_squared) / static_cast<long double>(_mc_statistics.num_diagrams); // average squared number of external phonons
    long double avg_ph_int = static_cast<long double>(_mc_statistics.avg_ph_int) / static_cast<long double>(_mc_statistics.num_diagrams); // average number of internal phonons
    long double avg_ph_int_sq = static_cast<long double>(_mc_statistics.avg_ph_int_squared)/static_cast<long double>(_mc_statistics.num_diagrams); // average squared number of internal phonons

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
    std::cout << "Average order of diagrams: " << avg_order << std::endl;
    std::cout << "Std dev order of diagrams: " << std::sqrt(avg_order_sq - avg_order*avg_order) << std::endl;
    std::cout << "Average number of internal phonons: " << avg_ph_int << std::endl;
    std::cout << "Std dev number of internal phonons: " << std::sqrt(avg_ph_int_sq - avg_ph_int*avg_ph_int) << std::endl;
    std::cout << "Average number of external phonons: " << avg_ph_ext << std::endl;
    std::cout << "Std dev number of external phonons: " << std::sqrt(avg_ph_ext_sq - avg_ph_ext*avg_ph_ext) << std::endl;
    std::cout << "Number of zero order diagrams: " << _mc_statistics.zero_order_diagrams << std::endl;
    std::cout << std::endl;

    if(_num_bands > 1){
        printBoldStatistics("simulation");
    }

    writeMCStatistics("MC_Statistics.txt");

    std::cout << std::endl;
};

void GreenFuncNphBands::printBoldStatistics(std::string type){
    long long int N = 1;
    if(type == "thermalization"){
        N = getRelaxSteps();
    }
    else if(type == "simulation"){
        N = getNdiags();
    }
    else if(type == "autocorrelation"){
        N = getAutocorrSteps();
    }

    if(type != "simulation"){computeRatioNegativeUpdates(N);}

    std::cout << "Bold Diagrammatic Monte Carlo analysis" << std::endl;
    for(int i = 0; i < 8; ++i){std::cout << "Update " << i << ", negative diagram updates: " << _num_negative_diagrams[i] << std::endl;}
    std::cout << std::endl;
    std::cout << "Total ratio (" << type << "): " << _ratio_negative_updates << "." << std::endl;
    std::cout << std::endl;
}

void GreenFuncNphBands::markovChainMC(){

    // input variables
    long double tau_length = diagramLengthUpdate(_tail->tau);
    double r = 0.5;
    unsigned long long int i = 0;

    unsigned long long int N_diags = getNdiags(); // number of diagrams to be generated
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

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors();}

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

    if(_num_bands > 1){
        if(_flags.mc_statistics){printBoldStatistics("thermalization");}
        resetNegativeDiagrams(); // reset negative diagram count after thermalization
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

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors();}

        if(static_cast<int>(i%(N_diags/100)) == 0){bar.update(i);}
    }

    bar.finish();

    if(_flags.time_benchmark){
        _benchmark_sim->stopTimer();
    }

    std::cout << std::endl;
    std::cout << "Simulation process finished" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        std::cout << "Simulation time benchmark finished." << std::endl;
        _benchmark_sim->printResults(); 
        _benchmark_sim->writeResultsToFile("simulation_benchmark.txt");
        delete _benchmark_sim;
        std::cout << std::endl;
    }

    if(_num_bands > 1){computeRatioNegativeUpdates(getNdiags());}

    if(_flags.mc_statistics){
        printMCStatistics();
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

    if(_flags.Z_factor){
        printZFactor();
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

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors();}

        /*if(static_cast<int>(i%100000) == 0 && _flags.write_diagrams){
            writeDiagram("Diagrams.txt", i, r); // method to visualize diagram structure
        }*/

        if(static_cast<int>(i%(N_relax/100)) == 0){bar.update(i);}
    }
    bar.finish();

    if(_flags.time_benchmark){_benchmark_th->stopTimer();}

    std::cout << std::endl;
    std::cout << "Thermalization process finished" << std::endl;
    std::cout << std::endl;

    if(_flags.time_benchmark){
        std::cout << "Thermalization time benchmark finished." << std::endl;
        _benchmark_th->printResults();
        _benchmark_th->writeResultsToFile("thermalization_benchmark.txt");
        delete _benchmark_th;
        delete _benchmark_sim;
        std::cout << std::endl;
    }

    if(_num_bands > 1){
        if(_flags.mc_statistics){printBoldStatistics("thermalization");}
        resetNegativeDiagrams(); // reset negative diagram count after thermalization
    }
};

void GreenFuncNphBands::markovChainMCOnlySample(){
    long double tau_length = _tail->tau;
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

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors();}

        if(static_cast<int>(i%(N_autocorr/100)) == 0 && _master){bar.update(i);}

    }

    if(_master){
        bar.finish();

        if(_flags.time_benchmark){
            _benchmark_th->stopTimer();
        }

        std::cout << std::endl;
        std::cout << "Autocorrelation process finished." << std::endl;
        std::cout << std::endl;

        if(_flags.time_benchmark){
            std::cout << "Autocorrelation time benchmark finished." << std::endl;
            _benchmark_th->printResults(); 
            _benchmark_th->writeResultsToFile("autocorrelation_benchmark.txt");
            std::cout << std::endl;
        }

        if(_num_bands > 1 && _flags.mc_statistics){
            printBoldStatistics("autocorrelation");
        }
    }

    if(_num_bands > 1){resetNegativeDiagrams();} // compute ratio of negative updates to total updates for autocorrelation process

    if(_flags.time_benchmark){delete _benchmark_th;}
    
    _flags.write_diagrams = false; 

    if(_master){
        std::cout << "Starting simulation process" << std::endl;
        std::cout << std::endl;

        if(N_diags <= 15000){
            _flags.write_diagrams = true; // set true for master process
        }
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

        if(static_cast<int>(i%100000) == 0){fixDoublePrecisionErrors(1e-9);}

        if(static_cast<int>(i%100000) == 0){checkTimeErrors();}
    
        if(_flags.write_diagrams && _master){writeDiagram("Diagrams.txt", i, r);}
        
        if(static_cast<int>(i%(N_diags/100)) == 0 && _master){bar.update(i);}
    }

    if(_master){
        bar.finish();
        if(_flags.time_benchmark){
            _benchmark_sim->stopTimer();
        }

        std::cout << std::endl;
        std::cout << "Simulation process finished." << std::endl;
        std::cout << std::endl;

        if(_flags.time_benchmark){
            std::cout << "Simulation time benchmark finished." << std::endl;
            _benchmark_sim->printResults(); 
            _benchmark_sim->writeResultsToFile("simulation_benchmark.txt");
            std::cout << std::endl;
        }
    }

    if(_flags.time_benchmark){delete _benchmark_sim;}

    computeFinalQuantities(); // compute final quantities for all processes (no write to console or files) 

    if(_master && _flags.mc_statistics){printMCStatistics();}
};

void GreenFuncNphBands::exactEstimatorGF(long double tau_length, int ext_phonon_order){

    if(ext_phonon_order < 0){ext_phonon_order = _current_ph_ext;}

    // compute electron bare propagators action
    int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
    double electron_energy = 0;
    double electron_action = 0., phonon_action = 0.;

    // phonon bare propagators action
    int phonon_index = -1;
    long double tau_one = 0;
    long double tau_two = 0;
    _helper = nullptr;
    FullVertexNode * second_helper = nullptr;

    electron_energy = electronEnergy(_head->k[0], _head->k[1], _head->k[2], _head->electronic_band.effective_mass);
    electron_action = electron_energy*(_head->next->tau - _head->tau);

    // main loop
    for(int i = 0; i < std::max(_current_order_int, 2*_current_ph_ext); ++i){
        if(i < _current_order_int){
            _helper = _internal_used[i].linked;
            tau_one = _helper->tau;
            electron_energy = electronEnergy(_helper->k[0], _helper->k[1], _helper->k[2], _helper->electronic_band.effective_mass);
            electron_action += electron_energy*(_helper->next->tau - tau_one);

            if(_helper->type == 1){
                phonon_index = _helper->index;
                second_helper = _internal_used[i].conjugated->linked;
                tau_two = second_helper->tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one);
            }
        }
        _helper = nullptr;
        second_helper = nullptr;
        if(i < 2*_current_ph_ext){
            _helper = _external_used[i].linked;
            tau_one = _helper->tau;
            electron_energy = electronEnergy(_helper->k[0], _helper->k[1], _helper->k[2], _helper->electronic_band.effective_mass);
            electron_action += electron_energy*(_helper->next->tau - tau_one);

            if(_helper->type == -2){
                phonon_index = _helper->index;
                second_helper = _external_used[i].conjugated->linked;
                tau_two = second_helper->tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_length + tau_one - tau_two);
            }
        }
        _helper = nullptr;
        second_helper = nullptr;
    }

    // compute Green function with exact estimator
    tau_length = (tau_length < 1e-9) ? 1e-9 : (tau_length >= _tau_max) ? _tau_max - 1e-9 : tau_length;
    int bin = static_cast<int>((tau_length - 0.) * 1./_points_step);

    long double prefactor = std::pow(1 + (_points[bin] - tau_length)/tau_length, current_order);
    long double exponential = std::exp(-((_points[bin] - tau_length)/tau_length)*(electron_action + phonon_action));
    long double diagrams_ratio = prefactor*exponential;

    _points_gf_exact[bin] += diagrams_ratio/(_points_step); // accumulate Green function value
};

double GreenFuncNphBands::groundStateEnergyExactEstimator(long double tau_length){
    if(tau_length <= _tau_cutoff_energy){return 0;} // reject if below cutoff
    else{

        // compute electron bare propagators action
        int current_order = _current_order_int + 2*_current_ph_ext; // total order of diagrams (number of phonon vertices)
        double electron_energy = 0;
        double electron_action = 0, phonon_action = 0;

        // phonon bare propagators action
        int phonon_index = -1;
        long double tau_one = 0;
        long double tau_two = 0;
        _helper = nullptr;
        FullVertexNode * second_helper;

        electron_energy = electronEnergy(_head->k[0], _head->k[1], _head->k[2], _head->electronic_band.effective_mass);
        electron_action = electron_energy*(_head->tau_next - _head->tau);

        // main loop
        for(int i = 0; i < std::max(_current_order_int, 2*_current_ph_ext); ++i){
        if(i < _current_order_int){
            _helper = _internal_used[i].linked;
            tau_one = _helper->tau;
            electron_energy = electronEnergy(_helper->k[0], _helper->k[1], _helper->k[2], _helper->electronic_band.effective_mass);
            electron_action += electron_energy*(_helper->tau_next - tau_one);

            if(_helper->type == 1){
                phonon_index = _helper->index;
                second_helper = _internal_used[i].conjugated->linked;
                tau_two = second_helper->tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_two - tau_one);
            }
        }
        _helper = nullptr;
        second_helper = nullptr;
        if(i < 2*_current_ph_ext){
            _helper = _external_used[i].linked;
            tau_one = _helper->tau;
            electron_energy = electronEnergy(_helper->k[0], _helper->k[1], _helper->k[2], _helper->electronic_band.effective_mass);
            electron_action += electron_energy*(_helper->tau_next - tau_one);

            if(_helper->type == -2){
                phonon_index = _helper->index;
                second_helper = _external_used[i].conjugated->linked;
                tau_two = second_helper->tau;
                phonon_action += phononEnergy(_phonon_modes, phonon_index)*(tau_length + tau_one - tau_two);
            }
        }
        _helper = nullptr;
        second_helper = nullptr;
    }

        double diagram_energy = (electron_action + phonon_action - static_cast<double>(current_order))/tau_length; // energy of current diagram

        if(_flags.blocking_analysis){groundStateEnergyBlockEstimator(diagram_energy);} // perform block analysis
        _gs_energy_count++; // update number of computed exact energy estimators
        return _current_sign*diagram_energy; // return energy of current diagram
    }
};

void GreenFuncNphBands::groundStateEnergyBlockEstimator(long double gs_energy){
    int block_number = static_cast<int>(_gs_energy_count/_block_size);
    _gs_energy_block_array[block_number] += gs_energy;

    if((_gs_energy_count+1)%_block_size == 0){
        _gs_energy_block_array[block_number] = _gs_energy_block_array[block_number]/static_cast<long double>(_block_size);
    }
};

double GreenFuncNphBands::effectiveMassExactEstimator(long double tau_length){
    if(tau_length <= _tau_cutoff_mass){return 0;}
    else{
        double mass_average_inv_x = 0; double mass_average_inv_y = 0; double mass_average_inv_z = 0;
        double electron_average_kx = 0; double electron_average_ky = 0; double electron_average_kz = 0;
        double mx, my, mz;
        double xP_inv = 1, yP_inv = 1, zP_inv = 1; 

        if(_num_bands == 1){
            mx = _m_x_el; my = _m_y_el; mz = _m_z_el;

            mass_average_inv_x = (_head->tau_next - _head->tau)/mx;
            mass_average_inv_y = (_head->tau_next - _head->tau)/my;
            mass_average_inv_z = (_head->tau_next - _head->tau)/mz;

            electron_average_kx = _head->k[0]*(_head->tau_next - _head->tau);
            electron_average_ky = _head->k[1]*(_head->tau_next - _head->tau);
            electron_average_kz = _head->k[2]*(_head->tau_next - _head->tau);

            for(int i=0; i < std::max(_current_order_int, 2*_current_ph_ext); ++i){
                if(i < _current_order_int){
                    _helper = _internal_used[i].linked;

                    mass_average_inv_x += (_helper->tau_next - _helper->tau)/mx;
                    mass_average_inv_y += (_helper->tau_next - _helper->tau)/my;
                    mass_average_inv_z += (_helper->tau_next - _helper->tau)/mz;

                    electron_average_kx += _helper->k[0]*(_helper->tau_next - _helper->tau)/mx;
                    electron_average_ky += _helper->k[1]*(_helper->tau_next - _helper->tau)/my;
                    electron_average_kz += _helper->k[2]*(_helper->tau_next - _helper->tau)/mz;
                }
                if(i < 2*_current_ph_ext){
                    _helper = _external_used[i].linked;
                    
                    mass_average_inv_x += (_helper->tau_next - _helper->tau)/mx;
                    mass_average_inv_y += (_helper->tau_next - _helper->tau)/my;
                    mass_average_inv_z += (_helper->tau_next - _helper->tau)/mz;

                    electron_average_kx += _helper->k[0]*(_helper->tau_next - _helper->tau)/mx;
                    electron_average_ky += _helper->k[1]*(_helper->tau_next - _helper->tau)/my;
                    electron_average_kz += _helper->k[2]*(_helper->tau_next - _helper->tau)/mz;
                }
            }

            mass_average_inv_x = (1./tau_length)*mass_average_inv_x;
            mass_average_inv_y = (1./tau_length)*mass_average_inv_y;
            mass_average_inv_z = (1./tau_length)*mass_average_inv_z;

            electron_average_kx = (1./tau_length)*electron_average_kx*electron_average_kx;
            electron_average_ky = (1./tau_length)*electron_average_ky*electron_average_ky;
            electron_average_kz = (1./tau_length)*electron_average_kz*electron_average_kz;

            xP_inv = static_cast<long double>(mass_average_inv_x - electron_average_kx); // x component
            yP_inv = static_cast<long double>(mass_average_inv_y - electron_average_ky); // y component
            zP_inv = static_cast<long double>(mass_average_inv_z - electron_average_kz); // z component

            _effective_masses[0] += xP_inv; // x component
            _effective_masses[1] += yP_inv; // y component
            _effective_masses[2] += zP_inv; // z component
        }
        else if (_num_bands == 3){
            // to be implemented
        }
        /*
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
        */

        // compute inverse of (average) effective mass, _D dimensionality of the system
        double inv_mass_avg = (tau_length*(std::pow(electron_average_kx,2) + std::pow(electron_average_ky,2) + std::pow(electron_average_kz,2))/_D);

        if(_flags.blocking_analysis){effectiveMassBlockEstimator(inv_mass_avg, xP_inv, yP_inv, zP_inv);}
        _effective_mass_count++;
        return inv_mass_avg;
    }
};

void GreenFuncNphBands::effectiveMassBlockEstimator(long double avg, long double xP, long double yP, long double zP){
    int block_number = static_cast<int>(_effective_mass_count/_block_size);

    _effective_mass_block_array[block_number] += avg;
    _effective_masses_block_array[3*block_number] += xP;
    _effective_masses_block_array[3*block_number+1] += yP;
    _effective_masses_block_array[3*block_number+2] += zP;
            
    if((_effective_mass_count+1)%_block_size == 0){
        _effective_mass_block_array[block_number] = _effective_mass_block_array[block_number]/static_cast<long double>(_block_size);
        _effective_masses_block_array[3*block_number] = _effective_masses_block_array[3*block_number]/static_cast<long double>(_block_size);
        _effective_masses_block_array[3*block_number+1] = _effective_masses_block_array[3*block_number+1]/static_cast<long double>(_block_size);
        _effective_masses_block_array[3*block_number+2] = _effective_masses_block_array[3*block_number+2]/static_cast<long double>(_block_size);

        _effective_mass_block_array[block_number] = (1.L)/_effective_mass_block_array[block_number];
        _effective_masses_block_array[3*block_number] = (1.L)/_effective_masses_block_array[3*block_number];
        _effective_masses_block_array[3*block_number+1] = (1.L)/_effective_masses_block_array[3*block_number+1];
        _effective_masses_block_array[3*block_number+2] = (1.L)/_effective_masses_block_array[3*block_number+2];
    }
};

void GreenFuncNphBands::ZFactorExactEstimator(long double tau_length){
    if(tau_length <= _tau_cutoff_Z){return;}
    else{
        _Z_factor_array[_current_ph_ext]++;
        _Z_factor_count++;
        return;
    }
};

void GreenFuncNphBands::computeRatioNegativeUpdates(long long int num_updates){
    long long int negative_updates_total = 0;
    for(int i = 0; i < 8; ++i){
        negative_updates_total += _num_negative_diagrams[i];
    }
    _ratio_negative_updates = static_cast<long double>(negative_updates_total)/static_cast<long double>(num_updates);
}

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

void GreenFuncNphBands::writeDiagram(std::string filename, int i, double r) const{
    std::ofstream file(filename, std::ofstream::app);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    FullVertexNode * pointer = _head;

    file << "Iteration: " << i << "\n";

    file << "index: " << 0 << " time: " << pointer->tau << " next time: " << pointer->tau_next << " wx: " << pointer->w[0] << " wy: " 
        << pointer->w[1] << " wz: " << pointer->w[2] <<" type: " << pointer->type << "\n";
        file << "\n";
        file << "propagator: " << 0 << "          kx: " << pointer->k[0] << " ky: " << pointer->k[1]
        << " kz: " << pointer->k[2] << "\n";
        file << "band number: " << pointer->electronic_band.band_number << " el eff mass: " << pointer->electronic_band.effective_mass
        << " (c1, c2, c3): (" <<  pointer->electronic_band.c1 << ", " <<  pointer->electronic_band.c2 << ", " <<  pointer->electronic_band.c3 << ")\n";
        file << "\n";

    pointer = pointer->next;
    int j = 1;
    while(pointer != _tail){
        file << "index: " << j << " time: " << pointer->tau << " next time: " << pointer->tau_next << " wx: " << pointer->w[0] << " wy: " 
        << pointer->w[1] << " wz: " << pointer->w[2] <<" type: " << pointer->type << 
        " phonon mode: " << pointer->index << " phonon energy: " << _phonon_modes[pointer->index] << "\n";
        file << "\n";
        file << "propagator: " << j << "          kx: " << pointer->k[0] << " ky: " << pointer->k[1]
        << " kz: " << pointer->k[2] << "\n";
        file << "band number: " << pointer->electronic_band.band_number << " el eff mass: " << pointer->electronic_band.effective_mass 
        << " (c1, c2, c3): (" << pointer->electronic_band.c1 << ", " << pointer->electronic_band.c2 << ", " << pointer->electronic_band.c3 << ")\n";
        file << "\n";
        pointer = pointer->next;
        ++j;
    }

    file << "index: " << j << " time: " << _tail->tau << " next time: " << _tail->tau_next << " wx: " << _tail->w[0] 
    << " wy: " << _tail->w[1] << " wz: " << _tail->w[2] 
    << " type: " << _tail->type <<  "\n";
    file << "\n";

    file << "ext phonons: " << _current_ph_ext << " int order: " << _current_order_int << " chosen update: " << r <<"\n";
    file << "\n";

    file << "\n";
    file << "\n";
    file << std::endl;
    file.close();

    checkConnections(filename);
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

void GreenFuncNphBands::checkConnections(std::string filename) const{
    std::ofstream file(filename, std::ofstream::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    int counter = 0;
    double wx = 0, wx_opp = 0;
    bool check_int = false, check_ext = false;
    for(int i = 0; i < _current_order_int; ++i){
        if(_internal_used[i].used == true){++counter;}
        if(_internal_used[i].linked->type == 1){
            wx = _internal_used[i].linked->w[0];
            wx_opp = _internal_used[i].conjugated->linked->w[0];
            if(!isEqual(wx, wx_opp)){
                check_int = true;
            }
        }
    }
    if(counter != _current_order_int){file << "Number of elements in internal array doesn't correspond to current internal order\n"; file << std::endl;}
    if(check_int){file << "Conjugated vertices do not correspond, please fix.\n"; file << std::endl;}

    counter = 0;
    for(int i = 0; i < 2*_current_ph_ext; ++i){
        if(_external_used[i].used == true){++counter;}
        if(_external_used[i].linked->type == 1){
            wx = _external_used[i].linked->w[0];
            wx_opp = _external_used[i].conjugated->linked->w[0];
            if(!isEqual(wx, wx_opp)){
                check_ext = true;
            }
        }
    }
    if(counter != 2*_current_ph_ext){file << "Number of elements in external array doesn't correspond to current external order\n"; file << std::endl;}
    if(check_ext){file << "Conjugated vertices do not correspond, please fix.\n"; file << std::endl;}

    file << "\n";
    file << "\n";

    file.close();
};

void GreenFuncNphBands::writeZFactor(const std::string& filename) const {
    std::ofstream file;

    file.open(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }
    file << "# Quasiparticle weight values (Z factor) for different number of external phonons:" << std::endl;
    file << "# Input parameters are: kx = " << _kx << ", ky = " << _ky << ", kz = " << _kz << std::endl;
    file << "# Chemical potential: " << _chem_potential << ", number of degenerate electronic bands : " << _num_bands << std::endl;
    file << "# Minimum length of diagrams for which Z factor is computed = " << _tau_cutoff_Z << "." << std::endl;
    file << "# Number of diagrams used for Z factor calculation: " << _Z_factor_count << std::endl;

    if(_num_bands == 1){
        file << "# Electronic effective masses: mx_el = " << _m_x_el << ", my_el = " 
            << _m_y_el << ", mz_el = " << _m_z_el << std::endl;
    }
    else if(_num_bands == 3){
        file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << _A_LK_el 
            << ", B_LK_el = " << _B_LK_el << ", C_LK_el = " << _C_LK_el << std::endl;
    }
    file << "# 1BZ volume: " << _V_BZ << " BvK volume: " << _V_BvK << " dielectric constant: " 
    << _dielectric_const << std::endl;
        file << std::endl;

    file << "# Number of phonon modes: " << _num_phonon_modes << std::endl;
    for(int i=0; i<_num_phonon_modes; i++){
        file << "# phonon mode (" << i << "): " << _phonon_modes[i] << ", Born effective charge (" << i << "): " 
        << _dielectric_responses[i] << std::endl;
    }
    file << std::endl;
    
    file << "# num external phonons    quasiparticle weight" << std::endl;
    for(int i=0; i<_ph_ext_max+1; i++){
        file << i << "    " << computeZFactor(i) <<std::endl;
    }
    file << std::endl;
    file.close();
};

void GreenFuncNphBands::writeMCStatistics(std::string filename) const {
    std::ofstream file(filename, std::ofstream::app);
    
    long double avg_order = static_cast<long double>(_mc_statistics.avg_order) / static_cast<long double>(_mc_statistics.num_diagrams); // average order of diagrams
    long double avg_order_sq = static_cast<long double>(_mc_statistics.avg_order_squared) / static_cast<long double>(_mc_statistics.num_diagrams); // average squared order of diagrams
    long double avg_ph_ext = static_cast<long double>(_mc_statistics.avg_ph_ext) / static_cast<long double>(_mc_statistics.num_diagrams); // average number of external phonons
    long double avg_ph_ext_sq = static_cast<long double>(_mc_statistics.avg_ph_ext_squared) / static_cast<long double>(_mc_statistics.num_diagrams); // average squared number of external phonons
    long double avg_ph_int = static_cast<long double>(_mc_statistics.avg_ph_int) / static_cast<long double>(_mc_statistics.num_diagrams); // average number of internal phonons
    long double avg_ph_int_sq = static_cast<long double>(_mc_statistics.avg_ph_int_squared)/static_cast<long double>(_mc_statistics.num_diagrams); // average squared number of internal phonons

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
    file << "Average order of diagrams: " << avg_order << "\n";
    file << "Std dev order of diagrams: " << std::sqrt(avg_order_sq - avg_order*avg_order) << "\n";
    file << "\n";
    file << "Average number of internal phonons: " << avg_ph_int << "\n";
    file << "Std dev number of internal phonons: " << std::sqrt(avg_ph_int_sq - avg_ph_int*avg_ph_int) << "\n";
    file << "\n";
    file << "Average number of external phonons: " << avg_ph_ext << "\n";
    file << "Std dev number of external phonons: " << std::sqrt(avg_ph_ext_sq - avg_ph_ext*avg_ph_ext) << "\n";
    file << "\n";
    file << "Number of zero order diagrams: " << _mc_statistics.zero_order_diagrams << "\n";
    file << "\n";

    if(_num_bands > 1){
        file << "Bold Diagrammatic Monte Carlo analysis" << "\n";
        for(int i = 0; i < 8; ++i){
            file << "Update " << i << ", negative diagram updates: " << _num_negative_diagrams[i] << "\n";
        }
        file << "\n";
        file << "Total ratio (simulation): " << _ratio_negative_updates << "." << std::endl;
        file << "\n";
    }
    file << "Probabilities: " << _p_length << ", " <<  _p_add_int << ", " << _p_rem_int << ", " << _p_add_ext << ", "
        << _p_rem_ext << ", "<< _p_swap << ", " << _p_shift << ", " << _p_stretch << ".\n";
    file << "\n";

    file.close();
    std::cout << "Values' statistics written to file " << filename << "." << std::endl;
};


/*
#####################################################################################
OLD VERSION OF STRETCH UPDATE
#####################################################################################

    initialize momentum values for first electron propagator
    double kx = _propagators[0].el_propagator_kx;
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