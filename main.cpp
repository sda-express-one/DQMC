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
    unsigned long long int N = 5000000000;
    GreenFuncNph diagram(N, 50, 0, 0, 0, -1.25, 500, 51);
    
    // simulations settings
    diagram.setAlpha(1.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);

    // Markov chain settings
    diagram.setRelaxSteps(100000000);
    diagram.setProbabilities(0.2, 0.18, 0.18, 0.12, 0.12, 0.1, 0.08, 0.02);
    //diagram.setProbabilities(0.3, 0.3, 0.3, 0, 0, 0, 0.1, 0);

    //histogram settings
    diagram.setN_bins(500);

    // exact GF estimator settings
    diagram.setNumPoints(200);
    diagram.setWidthEval(5);
    diagram.setSelectedOrder(0);

    // print diagrams to file
    diagram.writeDiagrams(false);

    // exact estimators settings (time cutoffs)
    diagram.setTauCutoffEnergy(15.0);
    diagram.setTauCutoffMass(15.0);

    // main simulation
    diagram.markovChainMC(0, false, true, true, true, true);
    std::cout << "Terminating the program." << std::endl;
    return 0;
}