#ifndef MC_BENCHMARKING_HPP
#define MC_BENCHMARKING_HPP

#include <iostream>
#include <fstream>
#include <chrono>


class MC_Benchmarking {
    public:

    // constructor
    MC_Benchmarking(unsigned long long int total_iterations , int num_updates = 1);

    // copy constructor
    MC_Benchmarking(const MC_Benchmarking& other);

    // copy assignment operator
    MC_Benchmarking& operator=(const MC_Benchmarking& other);

    // destructor 
    ~MC_Benchmarking(){
        if(_num_updates > 0){
            delete[] _updates_time;
            delete[] _updates_iterations;
        }
    };

    // start and stop timers
    // for total simulation and each update
    void startTimer();
    void startUpdateTimer();
    void stopTimer();
    void stopUpdateTimer(int update_index);
    
    // print results to console and write to file
    void printResults();
    void writeResultsToFile(const std::string& filename);

    private:
    // store time taken for total simulation
    std::chrono::high_resolution_clock::time_point _start_time;
    std::chrono::high_resolution_clock::time_point _end_time;

    // store time taken for each update
    std::chrono::high_resolution_clock::time_point _start_time_update;
    std::chrono::high_resolution_clock::time_point _end_time_update;

    // full simulation time
    const unsigned long long int _total_iterations = 0;
    long double _total_time = 0;
    long double _total_time_avg = 0;

    // update time
    const int _num_updates = 0;
    unsigned long long int* _updates_iterations = nullptr;
    long double* _updates_time = nullptr; 
};

#endif 