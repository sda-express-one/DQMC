#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFuncNph.hpp"

int main(){
    long long int unsigned N = 5000000000;
    GreenFuncNph diagram(N, 50, 0, 0, 0, -8.00, 500, 50);
    
    // simulations settings
    diagram.setAlpha(6.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);

    // Markov chain settings
    diagram.setRelaxSteps(100000000);
    diagram.setProbabilities(0.2, 0.18, 0.18, 0.12, 0.12, 0.1, 0.08, 0.02);
    //diagram.setProbabilities(0.3, 0.3, 0.3, 0, 0, 0, 0.1, 0);

    //histogram settings
    diagram.setN_bins(500);

    // exact GF estimator settings
    diagram.setNumPoints(500);
    diagram.setWidthEval(5.0);
    diagram.setSelectedOrder(0);

    // print diagrams to file
    diagram.writeDiagrams(false);

    // exact estimators settings (time cutoffs)
    diagram.setTauCutoffEnergy(10.0);
    diagram.setTauCutoffMass(10.0);

    // main simulation
    diagram.markovChainMC(0, false, true, true, true, false);
    std::cout << "Terminating the program." << std::endl;
    return 0;
}