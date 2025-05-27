#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <sstream>
#include "MC_data_structures.hpp"
#include "GreenFuncNph.hpp"

// read input from .txt file
double* readProbabilities(const std::string& filename);
parameters readSimParameterstxt(const std::string& filename);
settings readSimSettingstxt(const std::string& filename);

int main(){
    // read simulation parameters from file
    double* probs = readProbabilities("simulation_probabilities_MC.txt");
    parameters sim = readSimParameterstxt("simulation_parameters.txt");
    settings sets = readSimSettingstxt("simulation_settings.txt");

    // initialize GreenFuncNph object
    GreenFuncNph diagram(sim.N_diags, sim.tau_max, sim.kx, sim.ky, sim.kz, sim.chem_potential, sim.order_int_max, sim.ph_ext_max);

    // simulations settings
    diagram.setAlpha(sim.alpha);
    diagram.setVolume(sim.volume);
    diagram.setDimension(sim.dimensions);

    // Markov chain settings
    diagram.setRelaxSteps(sim.relax_steps);
    diagram.setProbabilities(probs);
    delete[] probs;

    //histogram settings
    diagram.setN_bins(sets.num_bins);

    // exact GF estimator settings
    diagram.setNumPoints(sets.num_points_exact);
    diagram.setSelectedOrder(sets.selected_order);

    // set calculations perfomed
    diagram.setCalculations(sets.gf_exact, sets.histo, sets.gs_energy, sets.effective_mass, sets.Z_factor);
    // set benchmarking
    diagram.setBenchmarking(sets.time_benchmark);
    // set MC statistics
    diagram.setMCStatistics(sets.mc_statistics);
    // print diagrams to file
    diagram.writeDiagrams(sets.write_diagrams);

    // exact estimators settings (time cutoffs)
    diagram.setTauCutoffEnergy(sets.tau_cutoff_energy);
    diagram.setTauCutoffMass(sets.tau_cutoff_mass);
    diagram.setTauCutoffStatistics(sets.tau_cutoff_statistics);

    // main simulation
    diagram.markovChainMC();

    std::cout << std::endl;
    std::cout << "Terminating the program." << std::endl;
    std::cout << std::endl;
    return 0;
}


double* readProbabilities(const std::string& filename){
    double* probs = new double[8];
    // Default probabilities
    probs[0] = 1;
    for(int i = 1; i < 8; i++){
        probs[i] = 0;
    }

    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default probabilities." << std::endl;
        return probs;
    }

    std::string line;
    while(std::getline(file, line)){
        if(line.empty() || line[0] == '#') continue; // Skip empty lines and comments

        std::istringstream iss(line);
        std::string key;
        std::string value;

        if(iss >> key){
            if(key == "length_update" || key == "prob_length"){
                iss >> value;
                probs[0] = stringToDouble(value);
            } 
            else if(key == "add_internal_update" || key == "prob_add_internal"){
                iss >> value;
                probs[1] = stringToDouble(value);
            } 
            else if(key == "remove_internal_update" || key == "prob_remove_internal"){
                iss >> value;
                probs[2] = stringToDouble(value);
            } 
            else if(key == "add_external_update" || key == "prob_add_external"){
                iss >> value;
                probs[3] = stringToDouble(value);
            } 
            else if(key == "remove_external_update" || key == "prob_remove_external"){
                iss >> value;
                probs[4] = stringToDouble(value);
            } 
            else if(key == "swap_phonon_update" || key == "prob_swap"){
                iss >> value;
                probs[5] = stringToDouble(value);
            } 
            else if(key == "shift_phonon_update" || key == "prob_shift"){
                iss >> value;
                probs[6] = stringToDouble(value);
            } 
            else if(key == "stretch_diagram_update" || key == "prob_stretch"){
                iss >> value;
                probs[7] = stringToDouble(value);
            }
        }
    }
    file.close();
    return probs;
}

parameters readSimParameterstxt(const std::string& filename){
    parameters params;
    std::ifstream file(filename);

    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default parameters." << std::endl;
        return params;
    }

    std::string line;
    while(std::getline(file, line)){
        if(line.empty() || line[0] == '#') continue; // Skip empty lines and comments

        std::istringstream iss(line);
        std::string key;

        if(iss >> key){
            if(key == "simulation_steps" || key == "steps"){
                std::string value;
                iss >> value;
                if(value == "0"){
                    std::cerr << "Warning: N_diags is set to 0. Using default value of 100000000." << std::endl;
                    params.N_diags = 100000000;
                } 
                else {
                    params.N_diags = stringToUnsignedLongLongInt(value);
                }
            } 
            else if(key == "relax_steps" || key == "thermalization_steps"){
                std::string value;
                iss >> value;
                params.relax_steps = stringToUnsignedLongLongInt(value);
            } 
            else if(key == "max_tau_value") {
                std::string value;
                iss >> value;
                if(value == "0"){
                    std::cerr << "Warning: tau_max is set to 0. Using default value of 100." << std::endl;
                    params.tau_max = 100;
                } 
                else {
                    params.tau_max = stringToLongDouble(value);
                }
            } 
            else if(key == "dimensionality"){
                std::string value;
                iss >> value;
                if(value == "2"){
                    params.dimensions = 2;
                } 
                else if(value == "3"){
                    params.dimensions = 3;
                } 
                else {
                    std::cerr << "Warning: Invalid dimensionality. Using default value of 3." << std::endl;
                    params.dimensions = 3;
                }
            }
            else if(key == "volume"){
                std::string value;
                iss >> value;
                params.volume = stringToDouble(value);
            }
            else if(key == "kx"){
                std::string value;
                iss >> value;
                params.kx = stringToDouble(value);
            }
            else if(key == "ky"){
                std::string value;
                iss >> value;
                params.ky = stringToDouble(value);
            }
            else if(key == "kz"){
                std::string value;
                iss >> value;
                params.kz = stringToDouble(value);
            } 
            else if(key == "coupling_strength"){
                std::string value;
                iss >> value;
                params.alpha = stringToDouble(value);
            }
            else if(key == "chemical_potential"){
                std::string value;
                iss >> value;
                params.chem_potential = stringToDouble(value);
            }
            else if(key == "max_internal_order"){
                std::string value;
                iss >> value;
                params.order_int_max = stringToInt(value);
            }
            else if(key == "max_num_ext_phonons"){
                std::string value;
                iss >> value;
                params.ph_ext_max = stringToInt(value);
            }
        }
    }
    file.close();
    return params;
}

settings readSimSettingstxt(const std::string& filename){
    settings sets;
    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default settings." << std::endl;
        return sets;
    }

    std::string line;
    while(std::getline(file, line)){
        if(line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        //std::cout << line << std::endl;
        std::istringstream iss(line);
        std::string key;

        if(iss >> key){
            if(key == "exact_green_function" || key == "exact_Green Function" || key == "exact_GF"){
                std::string value;
                iss >> value;
                sets.gf_exact = stringToBool(value);
            } 
            else if(key == "histogram_method" || key == "histo"){
                std::string value;
                iss >> value;
                sets.histo = stringToBool(value);
            } 
            else if(key == "gs_energy" || key == "ground_state_energy"){
                std::string value;
                iss >> value;
                sets.gs_energy = stringToBool(value);
            } 
            else if(key == "effective_mass"){
                std::string value;
                iss >> value;
                sets.effective_mass = stringToBool(value);
            } 
            else if(key == "Z_factor"){
                std::string value;
                iss >> value;
                sets.Z_factor = stringToBool(value);
            } 
            else if(key == "write_diagrams"){
                std::string value;
                iss >> value;
                sets.write_diagrams = stringToBool(value);
            }
            else if(key == "time_benchmark"){
                std::string value;
                iss >> value;
                sets.time_benchmark = stringToBool(value);
            }
            else if(key == "mc_statistics" || key == "statistics" || key == "stats"){
                std::string value;
                iss >> value;
                sets.mc_statistics = stringToBool(value);
            } 
            else if(key == "num_points_exact" || key == "num_points_exact_GF" || key == "points_exact_GF" || key == "exact_GF_points"){
                std::string value;
                iss >> value;
                sets.num_points_exact = stringToInt(value);
            }
            else if(key == "points_(exact_GF)" || key == "num_points_(exact_GF)" || key == "num_points" || key == "points"){
                std::string value;
                iss >> value;
                sets.num_points_exact = stringToInt(value);
            } 
            else if(key == "bins_(histogram)" || key == "num_bins_(histogram)" || key == "num_bins" || key == "bins"){
                std::string value;
                iss >> value;
                sets.num_bins = stringToInt(value);
            } 
            else if(key == "selected_order_(GF)"){
                std::string value;
                iss >> value;
                sets.selected_order = stringToInt(value);
            }
            else if(key == "cutoff_tau_(gs_energy)"){
                std::string value;
                iss >> value;
                sets.tau_cutoff_energy = stringToLongDouble(value);
            } 
            else if(key == "cutoff_tau_(mass)" || key == "cutoff_tau_(effective_mass)"){
                std::string value;
                iss >> value;
                sets.tau_cutoff_mass = stringToLongDouble(value);
            }
            else if(key == "cutoff_tau_stats" || key == "cutoff_tau_statistics"){
                std::string value;
                iss >> value;
                sets.tau_cutoff_statistics = stringToLongDouble(value);
            }
        }
    }

    if(sets.num_points_exact <= 0 && sets.gf_exact){
        sets.gf_exact = false;
        std::cerr << "Warning: num points for exact GF is not set. GF will not be calculated." << std::endl;
    }
    if(sets.num_bins <= 0 && sets.histo){
        sets.histo = false;
        std::cerr << "Warning: num bins for histogram is not set. Histogram will not be calculated." << std::endl;
    }
    if(sets.tau_cutoff_energy <= 0 && sets.gs_energy){
        sets.gs_energy = false;
        std::cerr << "Warning: tau cutoff for gs energy is not set. gs energy will not be calculated." << std::endl;
    }
    if(sets.tau_cutoff_mass <= 0 && sets.effective_mass){
        sets.effective_mass = false;
        std::cerr << "Warning: tau cutoff for effective mass is not set. effective mass will not be calculated." << std::endl;
    }
    file.close();
    return sets;
}