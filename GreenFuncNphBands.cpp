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
    _phonon_dispersions = new double[_num_phonon_modes];
    _born_effective_charges = new double[_num_phonon_modes];

    // initialize bands
    _bands = new Band[_order_int_max + 2*_ph_ext_max + 1];  // needs further development

    if(_num_bands == 1){
        for(int i = 0; i <_order_int_max + 2*_ph_ext_max + 1; i++){
            _bands[i].band_number = 1;
            _bands[i].c1 = 1;
            _bands[i].c2 = 0;
            _bands[i].c3 = 0;
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
void GreenFuncNphBands::setPhononDispersions(double * phonon_dispersions){
    _phonon_dispersions = phonon_dispersions;
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


long double GreenFuncNphBands::diagramLengthUpdate(long double tau_init){
    // initialize momentum values for last propagator
    double kx = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_kx;
    double ky = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_ky;
    double kz = _propagators[_current_order_int + 2*_current_ph_ext].el_propagator_kz;

    {// generate new time value for last vertex, reject if it goes out of bounds
    /*long double tau_fin = _last_vertex - std::log(1-drawUniformR())/(electronDispersion(kx,ky,kz, _el_eff_mass)-_chem_potential 
        + phononDispersion(_ph_dispersion)*_current_ph_ext);
    if(tau_fin <= _tau_max){
        _vertices[_current_order_int + 2*_current_ph_ext + 1].tau = tau_fin;
        return tau_fin;*/
    }
    //else{return tau_init;}
};