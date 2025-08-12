#include "SimulationParams.h"
#include <iostream>

void SimulationParameters::print() const {
    std::cout << "nmode: " << nmode << "\n";
    std::cout << "nsurfz: " << nsurfz << "\n";
    std::cout << "nstep: " << nstep << "\n";
    std::cout << "nwrite: " << nwrite << "\n";
    std::cout << "dt: " << dt << "\n";
    std::cout << "zeta: " << zeta << "\n";
    std::cout << "kc1: " << kc1 << ", kc2: " << kc2 << "\n";
    std::cout << "ncount: " << ncount << "\n";
    std::cout << "freqFile: " << freqFile << "\n";
    std::cout << "modeFile: " << modeFile << "\n";
    std::cout << "surfFile: " << surfFile << "\n";
    std::cout << "inputDir: " << inputDir << "\n";
    std::cout << "resultDir: " << resultDir << "\n";
    std::cout << "iforce: " << iforce << "\n";
    if (iforce == 1) {
        std::cout << "forcef: " << forcef << "\n";
        std::cout << "famp: " << famp << "\n";
    } else {
        std::cout << "ps: " << ps << "\n";
        std::cout << "rho: " << rho << "\n";
        std::cout << "mu: " << mu << "\n";
    }
    std::cout << "mass: " << mass << "\n";
}
