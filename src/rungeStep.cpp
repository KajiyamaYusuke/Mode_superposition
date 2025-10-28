#include "TimeIntegrator.h"
#include <vector>

// 1モード分をRunge-Kuttaで時間発展させる
void TimeIntegrator::rungeStep(double f, double q, double qdot, double dt,
               double omg, double zeta,
               double &qf, double &qfdot)
{   
    double L1 = dt * qdot;
    double K1 = (f - 2.0 * zeta * omg * qdot - omg * omg * q) * dt;

    double L2 = dt * (qdot + K1 / 2.0);
    double K2 = (f - 2.0 * zeta * omg * (qdot + K1 / 2.0) - omg * omg * (q + L1 / 2.0)) * dt;

    double L3 = dt * (qdot + K2 / 2.0);
    double K3 = (f - 2.0 * zeta * omg * (qdot + K2 / 2.0) - omg * omg * (q + L2 / 2.0)) * dt;

    double L4 = dt * (qdot + K3);
    double K4 = (f - 2.0 * zeta * omg * (qdot + K3) - omg * omg * (q + L3)) * dt;

    qf    = q    + (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;
    qfdot = qdot + (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
}
