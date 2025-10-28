
#include "SimulationParams.h"
#include "ForceCalculator.h"
#include "TimeIntegrator.h"
#include <iostream>

class Simulation {
public:
    SimulationParams params;
    State state; 
    Geometry geom;
    ModeData mdata;
    ForceCalculator fCalc;
    TimeIntegrator integrator;

    Simulation()
        : fCalc(geom, mdata, state, params )            // 必須引数を渡して初期化// TimeIntegrator も同様
    {}


    void initialize();
    void run();
    void writeVTK(int step, const Geometry& geom, const State& state, const std::string& rdir, int nwrite);
};
