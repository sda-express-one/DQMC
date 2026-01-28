#ifndef IO_METHODS_HPP
#define IO_METHODS_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include "MC_data_structures.hpp"
#include "../GreenFuncNphBands.hpp"

// read from .txt file
void readProbabilities(const std::string& filename, double * probs, int num_updates);

void readPhononModes(const std::string& filename, double * phonon_modes, double * dielectric_responses, int num_phonon_modes);

parameters readSimParameterstxt(const std::string& filename);

settings readSimSettingstxt(const std::string& filename);

cpu_info readCPUSettingstxt(const std::string& filename);

// write to .txt file (multithread process)
void writeGS_Energy(const std::string& filename, GreenFuncNphBands * diagram, int num_threads, 
    long double gs_energy_mean, long double * gs_energy_threads);

void writeEffectiveMass(const std::string filename, GreenFuncNphBands * diagram, int num_threads, 
    long double effective_mass_avg_mean, long double * effective_mass_mean, long double ** effective_mass_threads);

void writeGF_Histo(const std::string filename, GreenFuncNphBands * diagram, int num_threads,
    long double * histo_points, long double * histo_values);

void writeGF_Exact(const std::string filename, GreenFuncNphBands * diagram, int num_threads, 
    long double * gf_points, long double * gf_values);
#endif