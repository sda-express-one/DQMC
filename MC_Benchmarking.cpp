#include <iostream>
#include <fstream>
#include <chrono>
#include "MC_Benchmarking.hpp"

MC_Benchmarking::MC_Benchmarking(unsigned long long int total_iterations, int num_updates) : _total_iterations(total_iterations), 
    _num_updates(num_updates) {
    if(_num_updates > 0){   
        _updates_time = new long double[num_updates];
        _updates_iterations = new unsigned long long int[num_updates];
        for(int i = 0; i < num_updates; ++i) {
            _updates_time[i] = 0;
            _updates_iterations[i] = 0;
        }
    }
};

MC_Benchmarking::MC_Benchmarking(const MC_Benchmarking& other) : _total_iterations(other._total_iterations), 
    _num_updates(other._num_updates) {
    if(_num_updates > 0){
        _updates_time = new long double[_num_updates];
        _updates_iterations = new unsigned long long int[_num_updates];
        for(int i = 0; i < _num_updates; ++i) {
            _updates_time[i] = other._updates_time[i];
            _updates_iterations[i] = other._updates_iterations[i];
        }
    }
};

MC_Benchmarking& MC_Benchmarking::operator=(const MC_Benchmarking&  other){
    if(this != &other) {
        if(_num_updates > 0) {
            delete[] _updates_time;
            delete[] _updates_iterations;
            _updates_time = new long double[_num_updates];
            _updates_iterations = new unsigned long long int[_num_updates];
            for(int i = 0; i < _num_updates; ++i) {
                _updates_time[i] = other._updates_time[i];
                _updates_iterations[i] = other._updates_iterations[i];
            }
        }
    }
    return *this;
}

void MC_Benchmarking::startTimer() {
    _start_time = std::chrono::high_resolution_clock::now();
};

void MC_Benchmarking::startUpdateTimer() {
    _start_time_update = std::chrono::high_resolution_clock::now();
};

void MC_Benchmarking::stopTimer() {
    _end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(_end_time - _start_time);
    _total_time += static_cast<long double>(duration.count()) / 1e6; // Convert to seconds
};

void MC_Benchmarking::stopUpdateTimer(int update_index) {
    _end_time_update = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(_end_time_update - _start_time_update);
    _updates_time[update_index] += static_cast<long double>(duration.count()) / 1e6; // Convert to seconds
    _updates_iterations[update_index]++;
};

void MC_Benchmarking::printResults(){
    std::cout << "Total time taken: " << _total_time << " seconds." << std::endl;
    std::cout << "Average time per iteration: " << _total_time / _total_iterations << " seconds." << std::endl;
    std::cout << "Total iterations: " << _total_iterations << "." << std::endl;
    for(int i = 0; i < _num_updates; ++i) {
        std::cout << "Average time for update " << i << ": " << _updates_time[i] / _updates_iterations[i] << " seconds" << std::endl;
        std::cout << "Number of iterations for update " << i << ": " << _updates_iterations[i] << "." << std::endl;
    }
    std::cout << std::endl;
};

void MC_Benchmarking::writeResultsToFile(const std::string& filename) {
    std::ofstream file;
    file.open(filename, std::ios_base::app);
    if(!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    file << "Total time taken: " << _total_time << " seconds." << std::endl;
        file << "Average time per iteration: " << _total_time / _total_iterations << " seconds." << std::endl;
    file << "Total iterations: " << _total_iterations << "." << std::endl;
    for(int i = 0; i < _num_updates; ++i) {
        file << "Average time for update " << i << ": " << _updates_time[i] / _updates_iterations[i] << " seconds" << std::endl;
        file << "Number of iterations for update " << i << ": " << _updates_iterations[i] << std::endl;
    }
    file << std::endl;
    file.close();
};