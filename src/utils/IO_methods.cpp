#include "../../include/utils/IO_methods.hpp"

void readProbabilities(const std::string& filename, double* probs, int num_updates){
    // Default probabilities
    for(int i = 0; i < num_updates; i++){
        probs[i] = (1./num_updates);
    }

    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default probabilities." << std::endl;
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
                probs[0%num_updates] = stringToDouble(value);
            } 
            else if(key == "add_internal_update" || key == "prob_add_internal"){
                iss >> value;
                probs[1%num_updates] = stringToDouble(value);
            } 
            else if(key == "remove_internal_update" || key == "prob_remove_internal"){
                iss >> value;
                probs[2%num_updates] = stringToDouble(value);
            } 
            else if(key == "add_external_update" || key == "prob_add_external"){
                iss >> value;
                probs[3%num_updates] = stringToDouble(value);
            } 
            else if(key == "remove_external_update" || key == "prob_remove_external"){
                iss >> value;
                probs[4%num_updates] = stringToDouble(value);
            } 
            else if(key == "swap_phonon_update" || key == "prob_swap"){
                iss >> value;
                probs[5%num_updates] = stringToDouble(value);
            } 
            else if(key == "shift_phonon_update" || key == "prob_shift"){
                iss >> value;
                probs[6%num_updates] = stringToDouble(value);
            } 
            else if(key == "stretch_diagram_update" || key == "prob_stretch"){
                iss >> value;
                probs[7%num_updates] = stringToDouble(value);
            }
        }
    }
    file.close();
};

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
            if(key == "type"){
                std::string value;
                iss >> value;
                if(value.empty()){
                    value = "simple";
                }
                params.type = value;
            }
            else if(key == "simulation_steps" || key == "steps"){
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
            else if(key == "dimensions"){
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
            else if(key == "electron_effective_mass"){
                std::string value;
                iss >> value;
                params.el_eff_mass = stringToDouble(value);
            }
            else if(key == "optical_phonon_dispersion"){
                std::string value;
                iss >> value;
                params.ph_dispersion = stringToDouble(value);
            }
            else if(key == "num_bands"){
                std::string value;
                iss >> value;
                params.num_bands = stringToInt(value);
            }
            else if(key == "num_phonon_modes"){
                std::string value;
                iss >> value;
                params.num_phonon_modes = stringToInt(value);
            }
            else if(key == "m_x"){
                std::string value;
                iss >> value;
                params.m_x = stringToDouble(value);
            }
            else if(key == "m_y"){
                std::string value;
                iss >> value;
                params.m_y = stringToDouble(value);
            }
            else if(key == "m_z"){
                std::string value;
                iss >> value;
                params.m_z = stringToDouble(value);
            }
            else if(key == "A_LK"){
                std::string value;
                iss >> value;
                params.A_LK = stringToDouble(value);
            }
            else if(key == "B_LK"){
                std::string value;
                iss >> value;
                params.B_LK = stringToDouble(value);
            }
            else if(key == "C_LK"){
                std::string value;
                iss >> value;
                params.C_LK = stringToDouble(value);
            }
            else if(key == "V_1BZ"){
                std::string value;
                iss >> value;
                params.V_BZ = stringToDouble(value);
            }
            else if(key == "V_BvK"){
                std::string value;
                iss >> value;
                params.V_BvK = stringToDouble(value);
            }
            else if(key == "dielectric_constant" || key == "dielectric_const" || key == "diel_const" || key == "diel_infty"){
                std::string value;
                iss >> value;
                params.dielectric_const = stringToDouble(value);
            }
        }
    }
    file.close();
    return params;
};

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
            else if(key == "blocking_method" || key == "blocking"){
                std::string value;
                iss >> value;
                sets.blocking_analysis = stringToBool(value);
            }
            else if(key == "N_blocks"){
                std::string value;
                iss >> value;
                sets.N_blocks = stringToInt(value);
            }
            else if(key == "fix_tau_value" || key == "fix_length"){
                std::string value;
                iss >> value;
                sets.fix_tau_value = stringToBool(value);
            }
        }
    }

    if(sets.num_points_exact <= 0 && sets.gf_exact){
        sets.gf_exact = false;
        std::cerr << "Warning: num points for exact GF is not correctly set. GF will not be calculated." << std::endl;
    }
    if(sets.num_bins <= 0 && sets.histo){
        sets.histo = false;
        std::cerr << "Warning: num bins for histogram is not correctly set. Histogram will not be calculated." << std::endl;
    }
    if(sets.tau_cutoff_energy <= 0 && sets.gs_energy){
        sets.gs_energy = false;
        std::cerr << "Warning: tau cutoff for gs energy is not correctly set. Ground state energy will not be calculated." << std::endl;
    }
    if(sets.tau_cutoff_mass <= 0 && sets.effective_mass){
        sets.effective_mass = false;
        std::cerr << "Warning: tau cutoff for effective mass is not correctly set. Effective mass will not be calculated." << std::endl;
    }
    file.close();
    return sets;
};

cpu_info readCPUSettingstxt(const std::string& filename){
    cpu_info cpu;
    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default CPU settings." << std::endl;
        return cpu;
    }
    std::string line;
    while(std::getline(file, line)){
        if(line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        std::istringstream iss(line);
        std::string key;
        if(iss >> key){
            if(key == "parallel_mode"){
                std::string value;
                iss >> value;
                // enable openMP parallelization
                cpu.parallel_mode = stringToBool(value);
            }
            if(key == "cpus" || key == "num_procs" || key == "number_of_processors"){
                std::string value;
                iss >> value;
                cpu.num_procs = stringToInt(value);
            } 
            else if(key == "nodes" || key == "num_nodes" || key == "number_of_nodes"){
                std::string value;
                iss >> value;
                cpu.num_nodes = stringToInt(value);
            }
            else if(key == "autocorr_steps" || key == "autocorrelation_steps"){
                std::string value;
                iss >> value;
                cpu.autocorr_steps = stringToInt(value);
            } 
            else if(key == "cpu_time"){
                std::string value;
                iss >> value;
                cpu.cpu_time = stringToBool(value);
            }
        }
    }
    file.close();
    return cpu;
};

void readPhononModes(const std::string& filename, double * phonon_modes, double * dielectric_responses, int num_phonon_modes){
    // default
    for(int i=0; i<num_phonon_modes; i++){
        phonon_modes[i] = 0.5;
        dielectric_responses[i] = 0.5;
    }

    std::ifstream file(filename);

    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        std::cerr << "Using default parameters." << std::endl;
        return;
    }

    int i = 0;
    int j = 0;
    std::string line;
    while(std::getline(file, line)){
        if(line.empty() || line[0] == '#') continue; // Skip empty lines and comments

        auto label_phonon = "phonon_mode(" + std::to_string(i) + ")";
        auto label_charge = "dielectric_response(" + std::to_string(j) + ")";

        std::istringstream iss(line);
        std::string key;

        if(iss >> key){
            if(key == label_phonon && i < num_phonon_modes){
                std::string value;
                iss >> value;
                phonon_modes[i] = stringToDouble(value);
                i++;
                
            }
            else if(key == label_charge && j < num_phonon_modes){
                std::string value;
                iss >> value;
                dielectric_responses[j] = stringToDouble(value);
                j++;
                
            }
        }
    }
    file.close();
};

void writeGS_Energy(const std::string& filename, GreenFuncNphBands * diagram, int num_threads, bool blocking, 
    long double gs_energy_mean, long double * gs_energy_threads,
    long double gs_energy_mean_var, long double * gs_energy_var_threads){
    std::cout << "Ground state energy of the system is: " << gs_energy_mean;
    if(blocking){ std::cout << " +\\- " << std::sqrt(gs_energy_mean_var);}
    std::cout << "." << std::endl;

    if(blocking){std::cout << "Number of blocks used for blocking method: " << diagram->getNumBlocks() << "." << std::endl;}

    std::cout << "Input parameters are: kx = " << diagram->getkx() << ", ky = " << diagram->getky() << ", kz = " << diagram->getkz() << "." << std::endl;

    std::cout << "Chemical potential: " << diagram->getChemPotential() << ", number of degenerate electronic bands: " << diagram->getNumBands() << "." << std::endl;
    std::cout << "Minimum length of diagrams for which gs energy is computed: " << diagram->getTauCutoffEnergy() << "." << std::endl;
    //std::cout << "Number of diagrams used for ground state energy calculation: " << _gs_energy_count << std::endl;

    if(diagram->getNumBands() == 1){
        std::cout << "Electronic effective masses: mx_el = " << diagram->get_m_x_el() << ", my_el = " 
            << diagram->get_m_y_el() << ", mz_el = " << diagram->get_m_z_el() << std::endl;
        }
    else if(diagram->getNumBands() == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << diagram->get_A_LK_el() 
            << ", B_LK_el = " << diagram->get_B_LK_el() << ", C_LK_el = " << diagram->get_C_LK_el() << std::endl;
    }

    std::cout <<"1BZ volume: " << diagram->get1BZVolume() << " BvK volume: " << diagram->getBvKVolume() << " optical dielectric constant: " 
        << diagram->getDielectricConst() << ", tau cutoff: " << diagram->getTauCutoffEnergy() << std::endl;
    std::cout << std::endl;

    std::cout << "Number of phonon modes: " << diagram->getNumPhononModes() << std::endl;

    for(int i=0; i<diagram->getNumPhononModes(); i++){
        std::cout << "phonon mode (" << i << "): " << diagram->getPhononMode(i) << ", Born effective charge (" << i << "): " 
        << diagram->getDielectricResponse(i) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Parallel process" << std::endl;
    std::cout << "Number of independent parallel processes: " << num_threads << std::endl;
    std::cout << "Final ground state energy computed as average of the independent parallel processes." << std::endl;

    std::ofstream file(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not gs_energy.txt open file " << filename << std::endl;
    }
    else{
        file << "# Ground state energy of the system is: " << gs_energy_mean;
        if(blocking){file << " +\\- " << std::sqrt(gs_energy_mean_var);}
        file << "." << std::endl;

        if(blocking){file << "# Number of blocks used for blocking method: " << diagram->getNumBlocks() << "." << std::endl;}

        file << "# Input parameters are: kx = " << diagram->getkx() << ", ky = " << diagram->getky() << ", kz = " << diagram->getkz() << "." << std::endl;
        file << "# Chemical potential: " << diagram->getChemPotential() << ", number of degenerate electronic bands : " << diagram->getNumBands() << std::endl;
        file << "# minimum length of diagrams for which gs energy is computed = " << diagram->getTauCutoffEnergy() << "." << std::endl;
        //file << "Number of diagrams used for ground state energy calculation: " << _gs_energy_count << std::endl;

        if(diagram->getNumBands() == 1){
            file << "# Electronic effective masses: mx_el = " << diagram->get_m_x_el() << ", my_el = " 
                << diagram->get_m_y_el() << ", mz_el = " << diagram->get_m_z_el() << std::endl;
        }
        else if(diagram->getNumBands() == 3){
            file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << diagram->get_A_LK_el() 
                << ", B_LK_el = " << diagram->get_B_LK_el() << ", C_LK_el = " << diagram->get_C_LK_el() << std::endl;
        }
        file << "# 1BZ volume: " << diagram->get1BZVolume() << " BvK volume: " << diagram->getBvKVolume() << " optical dielectric constant: " 
        << diagram->getDielectricConst() << ", tau cutoff: " << diagram->getTauCutoffEnergy() << std::endl;
        file << std::endl;

        file << "# Number of phonon modes: " << diagram->getNumPhononModes() << std::endl;
        for(int i=0; i<diagram->getNumPhononModes(); i++){
            file << "# phonon mode (" << i << "): " << diagram->getPhononMode(i) << ", Born effective charge (" << i << "): " 
            << diagram->getDielectricResponse(i) << std::endl;
        }
        file << std::endl;
        file << "# Parallel process" << std::endl;
        file << "# Number of independent parallel processes: " << num_threads << std::endl;
        file << "# Final ground state energy computed as average of the independent parallel processes." << std::endl;
        file << "# Computed individual values are: " << std::endl;
        file << "# mean     variance" << std::endl;
        for(int i = 0; i < num_threads; i++){
            file << gs_energy_threads[i]; 
            if(blocking){file << "  " << gs_energy_var_threads[i];}
            file << std::endl;
        }
        file << std::endl;
    }
    std::cout << std::endl;
};

void writeEffectiveMass(const std::string filename, GreenFuncNphBands * diagram, int num_threads, bool blocking, 
    long double effective_mass_avg_mean, long double effective_mass_avg_var, 
    long double * effective_masses_mean, long double ** effective_masses_threads,
    long double * effective_masses_mean_var, long double ** effective_masses_var_threads){

    if(blocking){std::cout << "Number of blocks used for blocking method: " << diagram->getNumBlocks() << "." << std::endl;}

    std::cout << "Input parameters are: chemical potential: " << diagram->getChemPotential() << ", number of degenerate electronic bands: " << diagram->getNumBands() << std::endl;
    //std::cout << "Number of diagrams used for effective mass calculation: " << _effective_mass_count << std::endl;

    if(diagram->getNumBands() == 1){
        std::cout << "Electronic effective masses: mx_el = " << diagram->get_m_x_el() << ", my_el = " 
            << diagram->get_m_y_el() << ", mz_el = " << diagram->get_m_z_el() << std::endl;
    }

    else if(diagram->getNumBands() == 3){
        std::cout << "Electronic Luttinger-Kohn parameters: A_LK_el = "  << diagram->get_A_LK_el() 
                << ", B_LK_el = " << diagram->get_B_LK_el() << ", C_LK_el = " << diagram->get_C_LK_el() << std::endl;
    }

    std::cout <<"1BZ volume: " << diagram->get1BZVolume() << " BvK volume: " << diagram->getBvKVolume() << " optical dielectric constant: " 
        << diagram->getDielectricConst() << ", tau cutoff: " << diagram->getTauCutoffMass() << std::endl;
    std::cout << std::endl;

    std::cout << "Number of phonon modes: " << diagram->getNumPhononModes() << std::endl;
    for(int i=0; i<diagram->getNumPhononModes(); i++){
        std::cout << "phonon mode (" << i << "): " << diagram->getPhononMode(i) << ", Born effective charge (" << i << "): " 
        << diagram->getDielectricResponse(i) << std::endl;
    }

    std::cout << std::endl;

    if(diagram->getNumBands() == 1){
        std::cout << "Polaronic effective masses are: mx_pol = " << effective_masses_mean[0];
        if(blocking){std::cout << " +\\- " << std::sqrt(effective_masses_mean_var[0]);}
        std::cout << ", my_pol = " << effective_masses_mean[1];
        if(blocking){std::cout << " +\\- " << std::sqrt(effective_masses_mean_var[1]);} 
        std::cout << ", mz_pol = " << effective_masses_mean[2];
        if(blocking){std::cout << " +\\- " << std::sqrt(effective_masses_mean_var[2]);}
        std::cout << "." << std::endl;
    }
    else if(diagram->getNumBands() == 3){
        // stuff
    }
    std::cout << std::endl;

    std::cout << "Average effective mass of diagrams is: " << effective_mass_avg_mean;
    if(blocking){std::cout << " +\\- " << std::sqrt(effective_mass_avg_var);}
    std::cout << "." << std::endl;

    std::cout << std::endl;

    std::cout << "Parallel process" << std::endl;
    std::cout << "Number of independent parallel processes: " << num_threads << std::endl;
    std::cout << "Final polaronic effective masses computed as average of the independent parallel processes." << std::endl;

    std::string new_filename = "effective_mass.txt";
    std::ofstream file(new_filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not effective_mass.txt open file " << filename << std::endl;
    }
    else{
        file << "# Effective mass of the system is: " << effective_mass_avg_mean;
        if(blocking){file <<" +\\- " << effective_mass_avg_var;}
        file << "." << std::endl;

        if(blocking){file << "# Number of blocks used for blocking method: " << diagram->getNumBlocks() << "." << std::endl;}

        file << "# Chemical potential: " << diagram->getChemPotential() << ", number of degenerate electronic bands: " << diagram->getNumBands() << std::endl;
        //file << "Number of diagrams used for effective mass calculation: " << _effective_mass_count << std::endl;

        if(diagram->getNumBands() == 1){
            file << "# Electronic effective masses: mx_el = " << diagram->get_m_x_el() << ", my_el = " 
                << diagram->get_m_y_el() << ", mz_el = " << diagram->get_m_z_el() << std::endl;
        }
        else if(diagram->getNumBands() == 3){
            file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << diagram->get_A_LK_el()
                    << ", B_LK_el = " << diagram->get_B_LK_el() << ", C_LK_el = " << diagram->get_C_LK_el() << std::endl;
        }
        file <<"# 1BZ volume: " << diagram->get1BZVolume() << " BvK volume: " << diagram->getBvKVolume() << " optical dielectric constant: " 
        << diagram->getDielectricConst() << ", tau cutoff: " << diagram->getTauCutoffMass() << std::endl;
        file << std::endl;

        file << "# Number of phonon modes: " << diagram->getNumPhononModes() << std::endl;
        for(int i=0; i<diagram->getNumPhononModes(); i++){
            file << "# phonon mode (" << i << "): " << diagram->getPhononMode(i) << ", Born effective charge (" << i << "): " 
                << diagram->getDielectricResponse(i) << std::endl;
        }
        file << std::endl;

        if(diagram->getNumBands() == 1){
            file << "# Polaronic effective masses are: mx_pol = " << effective_masses_mean[0];
            if(blocking){file << " +\\- " << std::sqrt(effective_masses_mean_var[0]);}
            file << ", my_pol = " << effective_masses_mean[1];
            if(blocking){file << " +\\- " << std::sqrt(effective_masses_mean_var[1]);}
            file << ", mz_pol = " << effective_masses_mean[2];
            if(blocking){file << " +\\- " << std::sqrt(effective_masses_mean_var[2]);}
            file << std::endl; 
        }
       else if(diagram->getNumBands() == 3){
            // stuff
        }

        file << std::endl;
        
        file << "# Parallel process" << std::endl;
        file << "# Number of independent parallel processes: " << num_threads << std::endl;
        file << "# Final polaronic effective masses computed as average of the independent parallel processes." << std::endl;
        
        file << std::endl;
        file << "# Computed individual values are: " << std::endl;
        file << "mx     var(mx)     my      var(my)     mz      var(mz)" << std::endl;
        for(int i = 0; i < num_threads; i++){
            file << effective_masses_threads[i][0] << " ";
            if(blocking){file << effective_masses_var_threads[i][0] << " ";} 
            file<< effective_masses_threads[i][1] << " ";
            if(blocking){file << effective_masses_var_threads[i][1] << " ";}
            file << effective_masses_threads[i][2];
            if(blocking){file << " " << effective_masses_var_threads[i][2];} 
            file << std::endl;
        }
        file << std::endl;
        file.close();
    }
    std::cout << std::endl;
};

void writeGF_Histo(const std::string filename, GreenFuncNphBands * diagram, int num_threads, 
    long double * histo_points, long double * histo_values){
    
    std::cout << "GF computed (histogram method)" << std::endl;

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }

    file << "# Green function computed using the histogram method." << std::endl;
    file << "# kx = " << diagram->getkx() << ", ky = " << diagram->getky() << ", kz = " << diagram->getkz() << ", chemical potential = " << diagram->getChemPotential() << std::endl;

    file << "# Parallel process" << std::endl;
    file << "# Number of independent parallel processes: " << num_threads << "." << std::endl;
    file << "# Final Green function computed as average of the independent parallel processes." << std::endl;

    for(int i=0; i<diagram->getNumBins(); i++){
        file << histo_points[i] << " " << histo_values[i] << "\n";
    }
    file.close();
    std::cout << "Histogram written to file " << filename << "." << std::endl;
    std::cout << std:: endl;
};

void writeGF_Exact(const std::string filename, GreenFuncNphBands * diagram, int num_threads, 
    long double * gf_points, long double * gf_values){
    
    if(diagram->getGFSelectedOrder() < 0){
        std::cout << "Exact Green's function calculated for all orders of external phonons." << std::endl;
    }
    else{
        std::cout << "Exact Green's function calculated for number of external phonons " << diagram->getGFSelectedOrder() << "." << std::endl;;
    }
    
    std::ofstream file;

    file.open(filename);

    if(!file.is_open()){
        std::cout << "Could not open file " << filename << std::endl;
        return;
    }

    if(diagram->getGFSelectedOrder() < 0){
        file << "# Exact Green's function calculated for all orders of external phonons.\n";
    }
    else{
        file << "# Exact Green's function calculated for number of external phonons " << diagram->getGFSelectedOrder() << ".\n";
    }
    file << "# kx = " << diagram->getkx() << ", ky = " << diagram->getky() << ", kz = " << diagram->getkz() << ", chemical potential = " << diagram->getChemPotential() << "\n";

    file << "# Parallel process" << std::endl;
    file << "# Number of independent parallel processes: " << num_threads << "." << std::endl;
    file << "# Final Green function computed as average of the independent parallel processes." << std::endl;

    for(int i=0; i<diagram->getNumPoints(); i++){
        file << gf_points[i] << " " << gf_values[i] << "\n";
    }
    file.close();
    std::cout << "Exact Green's function written to file " << filename << "." << std::endl;
    std::cout << std::endl;
};

void writeZFactor(const std::string filename, GreenFuncNphBands * diagram, int num_threads,
    long double * Z_Factor_array, long double * Z_Factor_array_var){
    
    std::cout << "Quasiparticle weights (Z factor) computed." << std::endl;

    std::ofstream file;
    file.open(filename, std::ofstream::app);

    if(!file.is_open()){
        std::cerr << "Could not open file " << filename << std::endl;
        return;
    }
    file << "# Quasiparticle weight values (Z factor) for different number of external phonons:" << std::endl;
    file << "# Number of independent parallel processes: " << num_threads << "." << std::endl;
    file << "# Input parameters are: kx = " << diagram->getkx() << ", ky = " << diagram->getky() << ", kz = " << diagram->getkz() << ", chemical potential = " << diagram->getChemPotential() << std::endl;
    file << "# chemical potential: " << diagram->getChemPotential() << ", number of degenerate electronic bands: " << diagram->getNumBands() << std::endl;
    file << "# minimum length of diagrams for which Z factor is computed = " << diagram->getTauCutoffZ() << "." << std::endl;
    if(diagram->getNumBands() == 1){
        file << "# Electronic effective masses: mx_el = " << diagram->get_m_x_el() << ", my_el = " 
            << diagram->get_m_y_el() << ", mz_el = " << diagram->get_m_z_el() << std::endl;
    }
    else if(diagram->getNumBands() == 3){
        file << "# Electronic Luttinger-Kohn parameters: A_LK_el = "  << diagram->get_A_LK_el()
                << ", B_LK_el = " << diagram->get_B_LK_el() << ", C_LK_el = " << diagram->get_C_LK_el() << std::endl;
    }
    file <<"# 1BZ volume: " << diagram->get1BZVolume() << " BvK volume: " << diagram->getBvKVolume() << " optical dielectric constant: " 
        << diagram->getDielectricConst() << std::endl;

    file << "# Number of phonon modes: " << diagram->getNumPhononModes() << std::endl;
    for(int i=0; i<diagram->getNumPhononModes(); i++){
        file << "# phonon mode (" << i << "): " << diagram->getPhononMode(i) << ", Born effective charge (" << i << "): " 
        << diagram->getDielectricResponse(i) << std::endl;
    }
    file << std::endl;

    file << "# number of external phonons     Z factor     variance(Z factor)" << std::endl;
    for(int i=0; i<diagram->getPhExtMax(); i++){
        file << i << "  " << Z_Factor_array[i] << " " << Z_Factor_array_var[i] << std::endl;
    }
    file << std::endl;
    file.close();

    std::cout << "Z factor values written to file " << filename << "." << std::endl;
    std::cout << std::endl;
};