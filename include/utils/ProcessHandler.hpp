#ifndef PROCESS_HANDLER_HPP
#define PROCESS_HANDLER_HPP
#include <iostream>
#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include "ElectronBand.hpp"

// ====================================================================
// Configuration handler for templates
// built once at startup
// ====================================================================

typedef std::function<double(double)> ProcessHandler;

struct PipelineConfig {
    ProcessHandler electronBandProcess;
    // ProcessHandler phononModeProcess;
    // ProcessHandler couplingProcess;

    // Runs a full buffer through all the inputs in sequence
    void processBuffer(std::vector<double>& buffer) const noexcept {
        for (double& s : buffer)
            s = electronBandProcess(s);
    }
};

// ====================================================================
// Registries, one map per subsystem, mapping string identifiers to processing functions
// ====================================================================

typedef std::function<ProcessHandler()> ElectronBandFactory;
// typedef std::function<ProcessHandler()> PhononModeFactory;
// typedef std::function<ProcessHandler()> CouplingFactory;

ProcessHandler createElectronBandProcess(const std::string& identifier) {
    static const std::unordered_map<std::string, ElectronBandFactory> registry = {
        {"parab", []() { return ProcessHandler(ElectronBandEnergy<ParabolicBandSimple>()); }},
        //{"IBZ", []() { return 0; }},
        // Add more electron band processing functions here
    };
    auto it = registry.find(identifier);
    if(it == registry.end()) throw std::invalid_argument("Unknown electron band identifier: " + identifier);
    return it->second();
};

#endif