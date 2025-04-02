#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "GreenFunc.hpp"

int main(){
    long long int N = 1000000000;
    GreenFunc diagram(N, 50.0, 0., 0., 0., -2.2, 500);
    
    // settings
    diagram.setAlpha(2.0);
    diagram.setVolume(1.0);
    diagram.setDimension(3);
    diagram.setN_bins(200);
    diagram.setRelaxSteps(10000000);
    //diagram.setProbabilities(1/3, 1/3);

    // main simulation
    diagram.markovChain(-1);
    std::cout << "Computation finished." << std::endl;
    diagram.writeHistogram("histo.txt");
    std::cout << "Terminating the program." << std::endl;
    return 0;
}