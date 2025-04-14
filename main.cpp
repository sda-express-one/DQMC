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
    long long int unsigned N = 1000000000;
    GreenFuncNph diagram(N, 50, 0, 0, 0, -2.20, 500, 100);
    
    // settings
    diagram.setAlpha(2.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);
    diagram.setN_bins(500);
    diagram.setRelaxSteps(100000000);
    diagram.setProbabilities(0.4, 0.3, 0.3, 0, 0, 0, 0, 0);
    diagram.writeDiagrams(true);
    diagram.setTauCutoffEnergy(15.0);
    diagram.setTauCutoffMass(15.0);

    // main simulation
    diagram.markovChainMC(0, false, true, false, false, false);
    std::cout << "Terminating the program." << std::endl;
    return 0;
}