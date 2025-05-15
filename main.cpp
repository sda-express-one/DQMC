#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
//#include <nlohmann/json.hpp> // Requires JSON for Modern C++
#include "GreenFuncNph.hpp"

int main(){
    unsigned long long int N = 500000000;
    GreenFuncNph diagram(N, 50, 0, 0, 0, -2.2, 500, 50);
    
    // simulations settings
    diagram.setAlpha(2.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);

    // Markov chain settings
    diagram.setRelaxSteps(100000000);
    //diagram.setProbabilities(0.2, 0.18, 0.18, 0.12, 0.12, 0.1, 0.08, 0.02);
    diagram.setProbabilities(0.3, 0.3, 0.3, 0, 0, 0, 0.1, 0);

    //histogram settings
    diagram.setN_bins(500);

    // exact GF estimator settings
    diagram.setNumPoints(500);
    diagram.setSelectedOrder(0);

    // print diagrams to file
    diagram.writeDiagrams(false);

    // exact estimators settings (time cutoffs)
    diagram.setTauCutoffEnergy(15.0);
    diagram.setTauCutoffMass(15.0);

    // main simulation
    diagram.markovChainMC(0, true, true, false, false, false);
    std::cout << "Terminating the program." << std::endl;
    return 0;
}



struct parameters
{   
    unsigned long long int N_diags = 0;
    double kx = 0;
    double ky = 0;
    double kz = 0;
    double alpha = 0;
    double chem_potential = -10;
    int order_int_max = 0;
    int ph_ext_max = 0;
};


struct settings
{
    bool gf_exact = false;
    bool histo = false;
    bool gs_energy = false;
    bool effective_mass = false;
    bool Z_factor = false;
};
