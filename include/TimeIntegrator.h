#include <iostream>
#include <fstream>
#include <vector>

class TimeIntegrator {
public:
    double dt;

    void rungeKuttaStep(double f, double q, double qdot, double omg, double zeta,
                        double& qf, double& qfdot);
};
