#include <iostream>
#include <fstream>
#include <vector>

class TimeIntegrator {
public:

    void rungeStep(double f, double q, double qdot, double dt,
                   double omg, double zeta,
                   double &qf, double &qfdot);
    
    void newmarkStep( double f, double q, double qdot, double qddot,
    double dt, double omg, double zeta,
    double beta, double gamma,
    double &qf, double &qfdot, double &qfddot);
};
