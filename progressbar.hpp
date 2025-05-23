#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
#include <iomanip>

// Platform-specific includes
#ifdef _WIN32
    #include <io.h> // windows-specific headers
#else
    #include <unistd.h>   // UNIX/Linux headers
#endif

class ProgressBar{
    public:

    ProgressBar(unsigned long long int total, int width = 50);
    ~ProgressBar() = default;
    void update(unsigned long long int current);
    void finish();
    inline void setTotal(unsigned long long int total) {_total = total;}

    private:

    unsigned long long int _total; // total number of iterations
    int _width; // width of the progress bar
    bool _is_terminal; // flag to check if output is a terminal, for backend compatibility

    bool checkTerminal() const; // check if standard output is a terminal
};

#endif