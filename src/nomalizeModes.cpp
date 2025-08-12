#include "ModeShape.h"
#include "SimulationParams.h"
#include <iostream>
#include <cmath>


void ModeShape::normalizeModes(double mass){
    
    SimulationParameters params;

    double nMass = mass / nPoints;
    double ci;

    for (int imode=0; imode<nModes; imode++){
        ci =0.0;
        for (int j=0; j<nPoints; j++){
            ci += std::pow(modes[imode][j][0],2) 
                + std::pow(modes[imode][j][1],2) 
                + std::pow(modes[imode][j][2],2);
        }
        ci = 1.0 / std::sqrt( nMass * ci );
    
        for (int j=0; j<nPoints; j++){
            modes[imode][j][0] *= ci;
            modes[imode][j][1] *= ci;
            modes[imode][j][2] *= ci;
        }
    }
    
}