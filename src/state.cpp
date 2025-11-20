#include "State.h"
#include <iostream>
#include <cmath>

void State::initialize(int nPoints_, int nModes_, int nSteps_, const Geometry& geom) {
    nPoints = nPoints_;
    nModes  = nModes_;
    nSteps  = nSteps_;

    // モード座標の初期化
    q.assign(nModes, 0.0);
    qdot.assign(nModes, 0.0);
    qddot.assign(nModes, 0.0);

    // 節点変位
    disp.assign(nPoints, Displacement());
    predictedDisp.assign(nPoints, Displacement());
    vel.assign(nPoints, Displacement());

    for(int i = 0; i < nPoints; i++){
        disp[i].ux = 0.0 + geom.points[i].x;
        disp[i].uy = 0.0 + geom.points[i].y;
        disp[i].uz = 0.0 + geom.points[i].z;
    }

    for(int i = 0; i < nPoints; i++){
        vel[i].ux = 0.0 ;
        vel[i].uy = 0.0 ;
        vel[i].uz = 0.0 ;
    }

    // 面積・角度は後で geom.nxsup / nsurfz に基づき初期化
    harea.clear();
    degree.clear();
}

void State::mode2uf(const Geometry& geom, const ModeData& modeData, int step) {
    if (step < 0 || step >= nSteps) return;

    for (int i = 0; i < nPoints; ++i) {
        predictedDisp[i].ufx = 0.0;
        predictedDisp[i].ufy = 0.0;
        predictedDisp[i].ufz = 0.0;

        vel[i].ux = 0.0;
        vel[i].uy = 0.0;
        vel[i].uz = 0.0;

        //std::cout<<"after_mode2uf|disp[1]= "<<disp[1].ux<<std::endl;

        for (int m = 0; m < nModes; ++m) {
            double qi = qf[m];
            predictedDisp[i].ufx += modeData.modes[m][i].ux * qi* 1.0e3;
            predictedDisp[i].ufy += modeData.modes[m][i].uy * qi* 1.0e3;
            predictedDisp[i].ufz += modeData.modes[m][i].uz * qi* 1.0e3;

            vel[i].ux += modeData.modes[m][i].ux * qfdot[m] * 1.0e3;
            vel[i].uy += modeData.modes[m][i].uy * qfdot[m]* 1.0e3;
            vel[i].uz += modeData.modes[m][i].uz * qfdot[m]* 1.0e3;
        }
        predictedDisp[i].ufx += geom.points[i].x ;
        predictedDisp[i].ufy += geom.points[i].y ;
        predictedDisp[i].ufz += geom.points[i].z ;

    }

}


void State::calcArea(const Geometry& geom) {
    if (geom.nxsup < 2 || geom.nsurfz < 2) return;


    // harea
    harea.assign(geom.nxsup, 0.0);
    for (int i = 0; i < geom.nxsup; ++i) {
        double hi = 0.0;
        for (int j = 0; j < geom.nsurfz-1; ++j) {
            int pid1 = geom.surfp[i][j];
            int pid2 = geom.surfp[i][j+1];

            if (pid1 < 0 || pid2 < 0) continue;

            double ztmp = std::abs(disp[pid2].uz - disp[pid1].uz);
            double ytmp = geom.ymid[j] - 0.5 * (disp[pid2].uy + disp[pid1].uy);

            if (ytmp < 0.0) ytmp = 0.0;
            if (!std::isfinite(ytmp)) ytmp = 0.0;
            if (!std::isfinite(ztmp)) ztmp = 0.0;

            hi += 2.0 * ytmp * ztmp;
        }

        harea[i] = hi;
        //std::cout<<"harea["<< 25 << "] = "<< harea[25]<<std::endl;
    }

    // degree
    degree.assign(2, std::vector<std::vector<double>>(geom.nxsup, std::vector<double>(geom.nsurfz, 0.0)));

    for (int i = 1; i < geom.nxsup-1; ++i) {
        for (int j = 1; j < geom.nsurfz-1; ++j) {
            int pid_left  = geom.surfp[i-1][j];
            int pid_right = geom.surfp[i+1][j];
            int pid_down  = geom.surfp[i][j-1];
            int pid_up    = geom.surfp[i][j+1];

            if (pid_left<0 || pid_right<0 || pid_down<0 || pid_up<0) continue;

            double dx  = 0.5*(disp[pid_right].ux - disp[pid_left].ux);
            double dy1 = 0.5*(disp[pid_right].uy - disp[pid_left].uy);
            double dy2 = 0.5*(disp[pid_up].uy - disp[pid_down].uy);
            double dz  = 0.5*(disp[pid_up].uz - disp[pid_down].uz);
            
            

            if (dx != 0.0) degree[0][i][j] = std::atan(dy1/dx);
            if (dz != 0.0) degree[1][i][j] = std::atan(dy2/dz);
        }
    }

}


void State::uf2u() {
    // モード座標を現在値に上書き

    //std::cout<<"before"<<disp[1].ufx<<std::endl;

    for(int i = 0; i < nModes; ++i){
        q[i] = qf[i];
        qdot[i] = qfdot[i];
        qddot[i] = qfddot[i];
    }

    // 節点座標を上書き
    for(int i = 0; i < nPoints; ++i) {
        disp[i].ux = predictedDisp[i].ufx;
        disp[i].uy = predictedDisp[i].ufy;
        disp[i].uz = predictedDisp[i].ufz;
    }
    //std::cout<<"after"<<disp[1].ux<<std::endl;
}

