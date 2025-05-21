#include "progressbar.hpp"
#include <iostream>
#include <iomanip>

ProgressBar::ProgressBar(int total, int width) : _total(total), _width(width) {
    _is_terminal = checkTerminal();
}



bool ProgressBar::checkTerminal() const {
    // Check if the output is a terminal
        #ifdef _WIN32
            return _isatty(_fileno(stdout));
        #else
            return isatty(fileno(stdout));
        #endif
};

void ProgressBar::update(int current){
    if(_is_terminal){
        if(static_cast<int>(current%(_total/100) == 0)){
            // Calculate the percentage completed
            double percentage = static_cast<double>(current) / _total;
            int pos = static_cast<int>(_width * percentage);

            // Print the progress bar
            std::cout << "\rProgress: [";
            for (int i = 0; i < _width; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << percentage * 100. << "%";
            std::cout.flush();
        }
    }
    else{
        if(static_cast<int>(current%(_total/100) == 0)){
            std::cout << "Progress: " << static_cast<double>(current)/static_cast<double>(_total)*100. << "%" <<std::endl;
        }
    }
}

void ProgressBar::finish(){
    if(_is_terminal){
        int pos = static_cast<int>(_width * 1.);
        std::cout << "\rProgress: [";
        for(int i = 0; i < _width; ++i){
            if(i < pos) std::cout << "=";
            else if(i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] 100%";
        std::cout << std::endl;
    }
    else{
        std::cout << "Progress: 100%" << std::endl;
    }
}