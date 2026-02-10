#include <iostream>
#include <string>
#include <omp.h>
#include "../include/utils/IO_methods.hpp"
#include "../include/utils/computational_methods.hpp"
#include "../include/utils/MC_data_structures.hpp"
#include "../include/GreenFuncNph.hpp"
#include "../include/GreenFuncNphBands.hpp"

uint64_t getClockTime();

int main(){
    // read simulation parameters from file
    double probs[8];
    readProbabilities("simulation_probabilities_MC.txt", probs, 8);
    parameters sim = readSimParameterstxt("simulation_parameters.txt");
    settings sets = readSimSettingstxt("simulation_settings.txt");
    cpu_info cpu = readCPUSettingstxt("simulation_cpu_settings.txt");

    // initialize GreenFuncNph object
    uint64_t seed = getClockTime();
    Diagram::setSeed(seed);

    // type == simple standard DMC computation
    if(sim.type == "simple"){
        GreenFuncNph diagram(sim.N_diags, sim.tau_max, sim.kx, sim.ky, sim.kz, sim.chem_potential, 
            sim.order_int_max, sim.ph_ext_max, sim.el_eff_mass, sim.ph_dispersion);
        // simulations settings
        diagram.setAlpha(sim.alpha);
        diagram.setVolume(sim.volume);
        diagram.setDimension(sim.dimensions);

        // Markov chain settings
        diagram.setRelaxSteps(sim.relax_steps);
        diagram.setProbabilities(probs);

        //histogram settings
        diagram.setN_bins(sets.num_bins);

        // exact GF estimator settings
        diagram.setNumPoints(sets.num_points_exact);
        diagram.setSelectedOrder(sets.selected_order);

        // set calculations perfomed
        diagram.setCalculations(sets.gf_exact, sets.histo, sets.gs_energy, sets.effective_mass, sets.Z_factor, sets.fix_tau_value);
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
    }

    // type == bands anisotropic 1 band model or 3 band LK model
    else if(sim.type == "bands"){

        double phonon_modes[sim.num_phonon_modes];
        double dielectric_responses[sim.num_phonon_modes];
        readPhononModes("simulation_parameters.txt", phonon_modes, dielectric_responses, sim.num_phonon_modes);

        if(cpu.parallel_mode == false){
            GreenFuncNphBands diagram(sim.N_diags, sim.tau_max, sim.kx, sim.ky, sim.kz, sim.chem_potential, sim.order_int_max,
                sim.ph_ext_max, sim.num_bands, sim.num_phonon_modes);
        
            diagram.setPhononModes(phonon_modes);
            diagram.setDielectricResponses(dielectric_responses);
            diagram.set1BZVolume(sim.V_BZ);
            diagram.setBvKVolume(sim.V_BvK);
            diagram.setDielectricConstant(sim.dielectric_const);

            if(sim.num_bands == 1){
                diagram.setEffectiveMasses(sim.m_x, sim.m_y, sim.m_z);
            }
            else if(sim.num_bands == 3){
                diagram.setLuttingerKohnParameters(sim.A_LK, sim.B_LK, sim.C_LK);
            }

            // Markov chain settings
            diagram.setRelaxSteps(sim.relax_steps);
            diagram.setProbabilities(probs);

            //histogram settings
            diagram.setN_bins(sets.num_bins);

            // exact GF estimator settings
            diagram.setNumPoints(sets.num_points_exact);
            diagram.setSelectedOrder(sets.selected_order);

            // set calculations perfomed
            diagram.setCalculations(sets.gf_exact, sets.histo, sets.gs_energy, sets.effective_mass, sets.Z_factor, sets.blocking_analysis, sets.fix_tau_value);
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
        }
        else{
            GreenFuncNphBands diagram_relax(sim.N_diags, sim.tau_max, sim.kx, sim.ky, sim.kz, sim.chem_potential, sim.order_int_max,
                sim.ph_ext_max, sim.num_bands, sim.num_phonon_modes);
            
            if(omp_get_max_threads() < cpu.num_procs){
                std::cerr << "Warning! Number of cpus per nodes set exceeds current architecture capabilities." << std::endl;
                std::cerr << "Setting number of parallel processes per node to " << omp_get_max_threads() << "." << std::endl;
                std::cerr << std::endl;
                cpu.num_procs = omp_get_max_threads();
            }
            if(cpu.num_procs <= 0){
                cpu.num_procs = omp_get_num_procs();
            }
            
            diagram_relax.setMaster(cpu.parallel_mode);
            diagram_relax.setNumNodes(cpu.num_nodes);
            diagram_relax.setNumProcs(cpu.num_procs);
            
            diagram_relax.setPhononModes(phonon_modes);
            diagram_relax.setDielectricResponses(dielectric_responses);
            diagram_relax.set1BZVolume(sim.V_BZ);
            diagram_relax.setBvKVolume(sim.V_BvK);
            diagram_relax.setDielectricConstant(sim.dielectric_const);

            if(sim.num_bands == 1){
                diagram_relax.setEffectiveMasses(sim.m_x, sim.m_y, sim.m_z);
            }
            else if(sim.num_bands == 3){
                diagram_relax.setLuttingerKohnParameters(sim.A_LK, sim.B_LK, sim.C_LK);
            }

            // Markov chain settings
            diagram_relax.setRelaxSteps(sim.relax_steps);
            diagram_relax.setAutcorrSteps(cpu.autocorr_steps);
            diagram_relax.setProbabilities(probs);

            //histogram settings
            diagram_relax.setN_bins(sets.num_bins);

            // exact GF estimator settings
            diagram_relax.setNumPoints(sets.num_points_exact);
            diagram_relax.setSelectedOrder(sets.selected_order);

            // blocking method
            diagram_relax.setNumBlocks(sets.N_blocks);

            // set calculations perfomed
            diagram_relax.setCalculations(sets.gf_exact, sets.histo, sets.gs_energy, sets.effective_mass, sets.Z_factor, sets.blocking_analysis, sets.fix_tau_value);
            // set benchmarking
            diagram_relax.setBenchmarking(sets.time_benchmark);
            // set MC statistics
            diagram_relax.setMCStatistics(sets.mc_statistics);
            // print diagrams to file
            diagram_relax.writeDiagrams(sets.write_diagrams);
            
            // exact estimators settings (time cutoffs)
            diagram_relax.setTauCutoffEnergy(sets.tau_cutoff_energy);
            diagram_relax.setTauCutoffMass(sets.tau_cutoff_mass);
            diagram_relax.setTauCutoffStatistics(sets.tau_cutoff_statistics);

            // main simulation
            diagram_relax.markovChainMCOnlyRelax();

            long double * gs_energy = nullptr;
            long double * gs_energy_var = nullptr;
            if(sets.gs_energy){
                gs_energy = new long double[cpu.num_procs];
                gs_energy_var = new long double[cpu.num_procs];
            }

            long double * effective_mass = nullptr;
            long double * effective_mass_var = nullptr;
            long double ** effective_masses = nullptr;
            long double ** effective_masses_var = nullptr;
            if(sets.effective_mass){
                effective_mass = new long double[cpu.num_procs];
                effective_mass_var = new long double[cpu.num_procs];
                effective_masses = new long double * [cpu.num_procs];
                effective_masses_var = new long double * [cpu.num_procs];
                for (int i = 0; i < cpu.num_procs; i++){
                    effective_masses[i] = new long double[3];
                    effective_masses_var[i] = new long double[3];
                }
            }

            long double ** points_histogram = nullptr;
            long double ** gf_histo = nullptr;
            int num_bins_histo = diagram_relax.getNumBins();
            int num_points_exact_gf = diagram_relax.getNumPoints();

            if(sets.histo){
                points_histogram = new long double * [cpu.num_procs];
                gf_histo = new long double * [cpu.num_procs];
                for(int i = 0; i < cpu.num_procs; i++){    
                    points_histogram[i] = new long double[num_bins_histo];
                    gf_histo[i] = new long double[num_bins_histo];
                }
            }

            long double ** points_gf_exact = nullptr;
            long double ** gf_values_exact = nullptr;
            if(sets.gf_exact){
                points_gf_exact = new long double * [cpu.num_procs];
                gf_values_exact = new long double * [cpu.num_procs];
                for(int i = 0; i < cpu.num_procs; i++){
                    points_gf_exact[i] = new long double[num_points_exact_gf];
                    gf_values_exact[i] = new long double[num_points_exact_gf];
                }
            }

            
            /*Propagator * propagators_thermalized = new Propagator[sim.order_int_max + 2*sim.ph_ext_max + 1];
            Band * bands_thermalized = new Band[sim.order_int_max + 2*sim.ph_ext_max + 1];
            Vertex * vertices_thermalized = new Vertex[sim.order_int_max + 2*sim.ph_ext_max + 2];

            for(int i = 0; i < sim.order_int_max + 2*sim.ph_ext_max + 1; i++){
                propagators_thermalized[i] = diagram_relax.getPropagator(i);
                bands_thermalized[i] = diagram_relax.getBand(i);
                vertices_thermalized[i] = diagram_relax.getVertex(i);
            }
            vertices_thermalized[sim.order_int_max + 2*sim.ph_ext_max + 1] = diagram_relax.getVertex(sim.order_int_max + 2*sim.ph_ext_max + 1);*/

            int current_order_int = diagram_relax.getCurrentOrderInt();
            int current_ph_ext = diagram_relax.getCurrentPhExt();

            int num_threads = 1;
            omp_set_num_threads(cpu.num_procs);

            seed = getClockTime();

            #pragma omp parallel
            {
                int ID = omp_get_thread_num();

                Propagator * propagators_thermalized = nullptr;
                Band * bands_thermalized = nullptr;
                Vertex * vertices_thermalized = nullptr;

                #pragma omp critical

                {
                    propagators_thermalized = new Propagator[sim.order_int_max + 2*sim.ph_ext_max + 1];
                    bands_thermalized = new Band[sim.order_int_max + 2*sim.ph_ext_max + 1];
                    vertices_thermalized = new Vertex[sim.order_int_max + 2*sim.ph_ext_max + 2];

                    for(int i = 0; i < sim.order_int_max + 2*sim.ph_ext_max + 1; i++){
                        propagators_thermalized[i] = diagram_relax.getPropagator(i);
                        bands_thermalized[i] = diagram_relax.getBand(i);
                        vertices_thermalized[i] = diagram_relax.getVertex(i);
                    }

                    
                }

                GreenFuncNphBands diagram_simulate(propagators_thermalized, vertices_thermalized, bands_thermalized,
                    sim.N_diags, sim.tau_max, sim.kx, sim.ky, sim.kz, sim.chem_potential, sim.order_int_max,
                    sim.ph_ext_max, sim.num_bands, sim.num_phonon_modes);
                #pragma omp critical
                {
                    Diagram::setSeed(seed, ID*2654435761U);
                    delete[] propagators_thermalized;
                    delete[] bands_thermalized;
                    delete[] vertices_thermalized;
                }

                #pragma omp barrier
                #pragma omp master
                {
                    diagram_simulate.setMaster(true);
                    num_threads = omp_get_num_threads();
                    /*delete[] propagators_thermalized;
                    delete[] bands_thermalized;
                    delete[] vertices_thermalized;*/
                    
                }
                #pragma omp barrier

                diagram_simulate.setPhononModes(phonon_modes);
                diagram_simulate.setDielectricResponses(dielectric_responses);
                diagram_simulate.set1BZVolume(sim.V_BZ);
                diagram_simulate.setBvKVolume(sim.V_BvK);
                diagram_simulate.setDielectricConstant(sim.dielectric_const);
                diagram_simulate.setCurrentOrderInt(current_order_int);
                diagram_simulate.setCurrentPhExt(current_ph_ext);

                if(sim.num_bands == 1){
                    diagram_simulate.setEffectiveMasses(sim.m_x, sim.m_y, sim.m_z);
                }
                else if(sim.num_bands == 3){
                    diagram_simulate.setLuttingerKohnParameters(sim.A_LK, sim.B_LK, sim.C_LK);
                }

                // Markov chain settings
                diagram_simulate.setAutcorrSteps(cpu.autocorr_steps);
                diagram_simulate.setProbabilities(probs);

                //histogram settings
                diagram_simulate.setN_bins(sets.num_bins);

                // exact GF estimator settings
                diagram_simulate.setNumPoints(sets.num_points_exact);
                diagram_simulate.setSelectedOrder(sets.selected_order);

                // blocking method
                diagram_simulate.setNumBlocks(sets.N_blocks);

                // set calculations perfomed
                diagram_simulate.setCalculations(sets.gf_exact, sets.histo, sets.gs_energy, sets.effective_mass, sets.Z_factor, sets.blocking_analysis, sets.fix_tau_value);
                // set benchmarking
                diagram_simulate.setBenchmarking(sets.time_benchmark);
                // set MC statistics
                diagram_simulate.setMCStatistics(sets.mc_statistics);
                // print diagrams to file
                diagram_simulate.writeDiagrams(sets.write_diagrams);

                // exact estimators settings (time cutoffs)
                diagram_simulate.setTauCutoffEnergy(sets.tau_cutoff_energy);
                diagram_simulate.setTauCutoffMass(sets.tau_cutoff_mass);
                diagram_simulate.setTauCutoffStatistics(sets.tau_cutoff_statistics);

                // main simulation
                if(cpu.cpu_time == true){

                    std::string filename = "cpu_times.txt";
                    std::ofstream file;

                    #pragma omp master 
                    {
                    file.open(filename);
                    file.close();
                    }

                    auto start = std::chrono::steady_clock::now();
                    diagram_simulate.markovChainMCOnlySample();
                    auto end = std::chrono::steady_clock::now();

                    #pragma omp critical
                    {
                    file.open(filename, std::ofstream::app);
                    if(!file.is_open()){std::cout << "Could not open file " << filename << std::endl;}
                    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                    file << ID << " finished the computation." << std::endl;
                    file << "Elapsed time: " << duration_s.count() << " seconds" << std::endl;
                    file.close();
                    }
                }
                else{
                    diagram_simulate.markovChainMCOnlySample();
                }

                # pragma omp critical 
                {   
                    if(sets.gs_energy){
                        gs_energy[ID] = diagram_simulate.getGSEnergy();
                        gs_energy_var[ID] = diagram_simulate.getGSEnergyVar();
                    }
                    if(sets.effective_mass){
                        effective_mass[ID] = diagram_simulate.getEffectiveMass();
                        effective_mass_var[ID] = diagram_simulate.getEffectiveMassVar();
                        diagram_simulate.getEffectiveMasses(effective_masses[ID]);
                        diagram_simulate.getEffectiveMassesVar(effective_masses_var[ID]);
                    }
                    if(sets.histo){diagram_simulate.getHistogram(points_histogram[ID], gf_histo[ID]);}
                    if(sets.gf_exact){diagram_simulate.getGFExactPoints(points_gf_exact[ID], gf_values_exact[ID]);}
                    /*delete[] propagators_thermalized;
                    delete[] bands_thermalized;
                    delete[] vertices_thermalized;*/
                }

                
            }

            if(sets.gs_energy){
                long double gs_energy_mean = computeMean(gs_energy, num_threads);
                long double gs_energy_mean_var = computeMean(gs_energy_var, num_threads)*(num_threads)/(num_threads-1);
    
                writeGS_Energy("gs_energy.txt", &diagram_relax, num_threads, sets.blocking_analysis,
                    gs_energy_mean, gs_energy,
                    gs_energy_mean_var, gs_energy_var
                );
                
                delete[] gs_energy;
                delete[] gs_energy_var;
            }

            if(sets.effective_mass){
                long double effective_mass_mean = computeMean(effective_mass, num_threads);
                long double effective_mass_mean_var = computeMean(effective_mass_var, num_threads)*num_threads/(num_threads-1);
                long double * effective_masses_values = new long double[num_threads];
                long double * effective_masses_values_var = new long double[num_threads];
                long double effective_masses_mean[3];
                long double effective_masses_mean_var[3];

                for(int i=0; i<3; i++){
                    for(int j=0; j < num_threads; j++){
                        effective_masses_values[j] = effective_masses[j][i];
                        effective_masses_values_var[j] = effective_masses_var[j][i];
                    }
                    effective_masses_mean[i] = computeMean(effective_masses_values, num_threads);
                    effective_masses_mean_var[i] = computeMean(effective_masses_values_var, num_threads)*num_threads/(num_threads-1);
                }

                writeEffectiveMass("effective_mass.txt", &diagram_relax, num_threads, sets.blocking_analysis,
                    effective_mass_mean, effective_mass_mean_var,
                    effective_masses_mean, effective_masses,
                    effective_masses_mean_var, effective_masses_var
                );

                delete[] effective_masses_values;
                delete[] effective_masses_values_var;
                delete[] effective_mass;
                delete[] effective_mass_var;
                for(int i=0; i < cpu.num_procs; i++){
                    delete[] effective_masses[i];
                    delete[] effective_masses_var[i];
                }
                delete[] effective_masses;
                delete[] effective_masses_var;
            }

            if(sets.histo){
                long double * gf_histo_threads = new long double[num_threads];
                long double * gf_histo_mean = new long double[num_bins_histo];
                for(int i = 0; i < num_bins_histo; i++){
                    for(int j = 0; j < num_threads; j++){
                        gf_histo_threads[j] = gf_histo[j][i];
                    }
                    gf_histo_mean[i] = computeMean(gf_histo_threads, num_threads);
                }
                
                writeGF_Histo("histo.txt", &diagram_relax, num_threads, points_histogram[0], gf_histo_mean);

                delete[] gf_histo_threads;
                delete[] gf_histo_mean;
                for(int i=0;i < cpu.num_procs; i++){
                    delete[] points_histogram[i];
                    delete[] gf_histo[i];
                }
                delete[] points_histogram;
                delete[] gf_histo;
            }

            if(sets.gf_exact){
                long double * gf_exact_threads = new long double[num_threads];
                long double * gf_exact_mean = new long double[num_points_exact_gf];
                for(int i=0; i < num_points_exact_gf; i++){
                    for(int j=0; j < num_threads; j++){
                        gf_exact_threads[j] = gf_values_exact[j][i];
                    }
                    gf_exact_mean[i] = computeMean(gf_exact_threads, num_threads);
                }
                std::string a = "GF_";
                auto b = std::to_string(diagram_relax.getGFSelectedOrder());
                if(diagram_relax.getGFSelectedOrder() < 0){
                    b = "total";
                }
                std::string c = "_exact.txt";
                
                writeGF_Exact(a+b+c, &diagram_relax, num_threads, points_gf_exact[0], gf_exact_mean);

                delete[] gf_exact_threads;
                delete[] gf_exact_mean;
                for(int i=0; i < cpu.num_procs; i++){
                    delete[] points_gf_exact[i];
                    delete[] gf_values_exact[i];
                }
                delete[] points_gf_exact;
                delete[] gf_values_exact;
            }
        }
    }
    
    //std::cout << std::endl;
    std::cout << "Terminating the program." << std::endl;
    std::cout << std::endl;
    return 0;
};

uint64_t getClockTime(){
    uint64_t seed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() ^
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    return seed;
};