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
    long long int unsigned N = 10000;
    GreenFuncNph diagram(N, 50, 0, 0, 0, -5.75, 500, 100);
    
    // settings
    diagram.setAlpha(5.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);
    diagram.setN_bins(500);
    diagram.setRelaxSteps(100000000);
    diagram.setProbabilities(0.2, 0.2, 0.2, 0.2, 0.2);
    diagram.writeDiagrams(true);

    // main simulation
    diagram.markovChainMC(0, false, false, false, false, false);
    std::cout << "Terminating the program." << std::endl;
    return 0;
}