#include "ModeShape.h"
#include <iostream>
#include <fstream>
#include <vector>

class ForceCalculator {
public:
    
    ModeShape& modeShape;

    std::vector<std::vector<double>> fx, fy, fz; // 力の一般座標
    std::vector<double> fi; // モード力

    void calculateExternalForce(double timeStep);
    void convertForceToMode();
};
