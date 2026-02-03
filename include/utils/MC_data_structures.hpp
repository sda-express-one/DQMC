#ifndef MC_DATA_STRUCTURES_HPP
#define MC_DATA_STRUCTURES_HPP
#include <iostream>
#include <string>


// This file contains data structures used in the Monte Carlo simulation for Green's function calculations

// parameters to set in simulation_parameters.txt
struct parameters{
    std::string type;
    unsigned long long int N_diags = 100000;
    unsigned long long int relax_steps = 100000;
    int dimensions = 3;
    long double tau_max = 100;
    double volume = 1;
    double kx = 0;
    double ky = 0;
    double kz = 0;
    double alpha = 0;
    double chem_potential = -2;
    int order_int_max = 0;
    int ph_ext_max = 0;
    int num_bands = 1;
    int num_phonon_modes = 1;
    double el_eff_mass = 1;
    double ph_dispersion = 1;
    double m_x = 1;
    double m_y = 1;
    double m_z = 1;
    double A_LK = 1;
    double B_LK = 1;
    double C_LK = 0;
    double dielectric_const = 1;
    double V_BZ = 1;
    double V_BvK = 1;
};

// parameters to set in cpu_settings.txt
struct cpu_info{
    bool parallel_mode = false;
    int num_procs = 1;
    int num_nodes = 1;
    int autocorr_steps = 0;
    bool cpu_time = false;
};

// parameters to set in simulation_settings.txt
struct settings{
    bool gf_exact = false;
    bool histo = false;
    bool gs_energy = false;
    bool effective_mass = false;
    bool Z_factor = false;
    bool write_diagrams = false;
    bool time_benchmark = false;
    bool mc_statistics = false;
    bool fix_tau_value = false;
    int num_points_exact = 100;
    int num_bins = 100;
    int selected_order = 0;
    long double tau_cutoff_energy = 10;
    long double tau_cutoff_mass = 10;
    long double tau_cutoff_statistics = 0;
};

// string to type conversion functions

inline int stringToInt(const std::string& str){
    try {
        return std::stoi(str);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid integer: " << str << std::endl;
        return 0;
    }
}

inline unsigned long long int stringToUnsignedLongLongInt(const std::string& str){
    try {
        return std::stoull(str);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid unsigned long long int: " << str << std::endl;
        return 0;
    }
}

inline bool stringToBool(const std::string& str){
    if(str == "true" || str == "True" || str == "1"){
        return true;
    } else if(str == "false" || str == "False" || str == "0"){
        return false;
    } else {
        std::cerr << "Invalid boolean: " << str << std::endl;
        return false;
    }
}

inline double stringToDouble(const std::string& str){
    try {
        return std::stod(str);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid double: " << str << std::endl;
        return 1.0;
    }
}

inline long double stringToLongDouble(const std::string& str){
    try {
        return std::stold(str);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid long double: " << str << std::endl;
        return 1.0;
    }
}

// Structure to hold flags for different types of calculations
struct Flags{
        bool gf_exact = false; // flag for exact Green function
        bool histo = false; // flag for histogram method
        bool gs_energy = false; // flag for ground state energy
        bool effective_mass = false; // flag for effective mass
        bool Z_factor = false; // flag for Z factor
        bool write_diagrams = false; // flag to write diagrams to file
        bool time_benchmark = false; // flag for time benchmarking
        bool mc_statistics = false; // flag for statistics
        bool fix_tau_value = false; // flag to fix diagram length during simulation (gs properties)
    };


// Structure to hold a vertex in the diagram
struct Vertex{
    long double tau = 0; // time coordinate of the vertex
    int type = 0; // +1 outgoing and -1 incoming (internal), +2 outgoing and -2 incoming (external), 0 unassigned (extrema)
    int linked = -1; // describes connection to other vertex of phonon propagator (-1 if not linked for extrema)
    double wx = 0.; // x component of phonon frequency
    double wy = 0.; // y component of phonon frequency
    double wz = 0.; // z component of phonon frequency
    int index = 0; // describes link to phonon mode (0 first phonon mode, 1 second phonon mode, etc.)
};

// Structure to hold a propagator in the diagram
struct Propagator{
    double el_propagator_kx = 0.; // x component of electron propagator momentum
    double el_propagator_ky = 0.; // y component of electron propagator momentum
    double el_propagator_kz = 0.; // z component of electron propagator momentum
};

struct MC_Statistics{
    unsigned long long int num_diagrams = 0; // total number of diagrams used for statistics
    long double avg_tau = 0; // average length of diagrams
    long double avg_tau_squared = 0; // average squared length of diagrams
    unsigned long long int avg_order = 0; // average order of diagrams
    unsigned long long int avg_order_squared = 0; // average squared order of diagrams
    unsigned long long int avg_ph_int = 0; // average number of internal phonons
    unsigned long long int avg_ph_int_squared = 0; // average squared number of internal phonons
    unsigned long long int avg_ph_ext = 0; // average number of external phonons
    unsigned long long int avg_ph_ext_squared = 0; // average squared number of external phonons
    unsigned long long int zero_order_diagrams = 0; // number of zero order diagrams
};

struct Band{
    int band_number = 0; // band number
    double effective_mass = 1.0; // effective mass of the band
    double c1 = 1; // coefficient for wavefunction phi_1 (band number 1)
    double c2 = 0; // coefficient for wavefunction phi_2 (band number 2)
    double c3 = 0; // coefficient for wavefunction phi_3 (band number 3)
};

#endif