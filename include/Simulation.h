#include "ModeShape.h"
#include "SimulationParams.h"
#include "ForceCalculator.h"
#include "TimeIntegrator.h"
#include <iostream>

class Simulation {
public:
    SimulationParameters params;
    ModeShape modeShape;
    ForceCalculator forceCalc;
    TimeIntegrator integrator;

    void initialize();
    void run();
    void outputResults();
};
