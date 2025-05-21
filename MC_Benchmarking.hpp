#ifndef MC_BENCHMARKING_HPP
#define MC_BENCHMARKING_HPP

#include <iostream>
#include <fstream>
#include <chrono>


class MC_Benchmarking {
    public:
    MC_Benchmarking(int num_updates = 9);
    ~MC_Benchmarking();
    void startTimer();
    void startUpdateTimer();
    void stopTimer();
    void stopUpdateTimer(int update_index);
    void calculateAverageTime();
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
    long double _total_time = 0;
    long long int _total_iterations = 0;
    long double _total_time_avg = 0;
    // update time
    const int _num_updates = 0;
    long double* _updates_time; 
    long long int* _updates_iterations;
};

#endif 