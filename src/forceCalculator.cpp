#include "ForceCalculator.h"
#include "Displacement.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

auto checkNaN = [](double val, const std::string& name) {
    if (std::isnan(val) || std::isinf(val)) {
        std::cerr << "[NaN DETECTED] " << name << " = " << val << std::endl;
    }
};


ForceCalculator::ForceCalculator(const Geometry& geom_, const ModeData& md_, State& st_, const SimulationParams& sp_)
    : geom(geom_), modeData(md_), state(st_), sp(sp_) {}




void ForceCalculator::initialize() {
    int nPoints = geom.nPoints;
    int nsurfl  = geom.nsurfl;
    int nsurfz  = geom.nsurfz;
    int nModes  = modeData.nModes;

    // 節点ごとの外力
    fx.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fy.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fz.assign(nsurfl, std::vector<double>(nsurfz, 0.0));

    // モード力
    fi.assign(nModes, 0.0);

    // 内部バッファ
    nxsup = geom.nxsup;  // 断面数 (x方向分割)
    fdis.assign(nxsup, std::vector<double>(geom.nsurfz - 1, 0.0));

    psurf.assign(nxsup, 0.0);
    Ug.assign(state.nSteps, 0.0);
    minHarea.assign(state.nSteps, 0.0);

    std::cout << "[ForceCalculator] initialized: "
              << "nPoints=" << nPoints
              << ", nModes=" << nModes
              << ", nxsup=" << nxsup << std::endl;
}

void ForceCalculator::calcForce(double t, int n) {
    int nsurfl = geom.nsurfl;
    int nsurfz = geom.nsurfz;
    int nxsup  = geom.nxsup;

    // まず全てゼロクリア
    fx.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fy.assign(nsurfl, std::vector<double>(nsurfz, 0.0));
    fz.assign(nsurfl, std::vector<double>(nsurfz, 0.0));


    for (int i = 0; i < nxsup; i++) {
        for (int j = 0; j < nsurfz - 1; j++) {
            fdis[i][j] = 0.0;
        }
    }


    if (sp.iforce == 1) {
        // ==== sin波加振 ====

        minHarea[n] = *std::min_element(state.harea.begin(), state.harea.end());
        
        
        for (int i = 1; i < 6; i++) {
            for (int j = 1; j < nsurfz-1; j++) {
                int pid = geom.surfp[i][j];
                int pid_ip1 = geom.surfp[i+1][j];
                int pid_im1 = geom.surfp[i-1][j];
                int pid_jp1 = geom.surfp[i][j+1];
                int pid_jm1 = geom.surfp[i][j-1];

                double dx = 0.5 * (state.disp[pid_ip1].ux - state.disp[pid_im1].ux);
                double dy = 0.5 * (state.disp[pid_ip1].uy - state.disp[pid_im1].uy);
                double ds = std::sqrt(dx*dx + dy*dy);
                double dz = 0.5 * (state.disp[pid_jp1].uz - state.disp[pid_jm1].uz);
                if (pid < 0) continue;
                fx[i][j] = 2900 * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::sin(state.degree[0][i][j]);
                fy[i][j] = -2900 * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::cos(state.degree[0][i][j]);
                fz[i][j] = 2900 * ds * dz * 1.0e-6 * std::sin(state.degree[1][i][j]);
            }
        }
    } else if (sp.iforce == 0) {

        // ==== 1D flow model ====
        minHarea[n] = *std::min_element(state.harea.begin(), state.harea.end());


        // separation point
        int nsep = findNsep(minHarea[n]) ;


        // 流量 Ug
        if (minHarea[n] > 0.0) {
            Ug[n] = std::sqrt(2.0 * sp.ps / sp.rho) * minHarea[n];
        } else {
            Ug[n] = 0.0;
        }


        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = sp.ps;

        if (minHarea[n] > 0.0) {
            for (int i = 1; i < nsep; i++) {
                double dx = std::abs(geom.points[geom.surfp[i][ nsurfz-2]].x - geom.points[geom.surfp[i-1][ nsurfz-2]].x);
                double h  = (state.harea[i] + state.harea[i-1]) / (2.0 * geom.zmax);

                psurf[i] = psurf[i-1] + 0.5 * sp.rho * Ug[n] * Ug[n] * (1.0 / (state.harea[i-1]*state.harea[i-1]) - 1.0 / (state.harea[i]*state.harea[i]))
                                      - 12.0 * sp.mu * dx / (geom.zmax * h * h * h) * Ug[n] * 1e3;
            }
        } else {
            for (int i = 1; i < nsep-1; i++) {
                psurf[i] = sp.ps;
            }
        }


        // 力 fx, fy, fz
        for (int i = 1; i < nsep-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {
                int pid_ip1 = geom.surfp[i+1][j];
                int pid_im1 = geom.surfp[i-1][j];
                int pid_jp1 = geom.surfp[i][j+1];
                int pid_jm1 = geom.surfp[i][j-1];

                double dx = 0.5 * (state.disp[pid_ip1].ux - state.disp[pid_im1].ux);
                double dy = 0.5 * (state.disp[pid_ip1].uy - state.disp[pid_im1].uy);
                double ds = std::sqrt(dx*dx + dy*dy);
                double dz = 0.5 * (state.disp[pid_jp1].uz - state.disp[pid_jm1].uz);

                fx[i][j] = psurf[i] * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::sin(state.degree[0][i][j]);
                fy[i][j] = -psurf[i] * ds * dz * 1.0e-6 * std::cos(state.degree[1][i][j]) * std::cos(state.degree[0][i][j]);
                fz[i][j] = psurf[i] * ds * dz * 1.0e-6 * std::sin(state.degree[1][i][j]);

            }
        }


        
    }
    

}

void ForceCalculator::f2mode() {



    for (int imode = 0; imode < modeData.nModes; imode++) {
        fi[imode] = 0.0;
        for (int i = 0; i < geom.nsurfl; i++) {
            for (int j = 0; j < geom.nsurfz; j++) {
                int pid = geom.surfp[i][j];
                if (pid < 0) continue; // 節点が存在しない場合スキップ

                fi[imode] += fx[i][j] * modeData.modes[imode][pid].ux
                           + fy[i][j] * modeData.modes[imode][pid].uy
                           + fz[i][j] * modeData.modes[imode][pid].uz;

            }
        }
    }

     //std::cout<<"|disp[1]= "<<state.disp[11].ux<<std::endl;
}

void ForceCalculator::contactForce() {

    contactFlag = false;
    double omg1 = 2.0 * M_PI * modeData.frequencies[0]; // 1次固有振動数
    double omg2 = omg1 * omg1;

    for (int i = 0; i < geom.nxsup; ++i) {           // nxsup は計算範囲
        for (int j = 1; j < geom.nsurfz - 1; ++j) {  // 2..nsurfz-1 (0-index)
            double yval = state.disp[geom.surfp[i][j]].uy;  // 節点変位 v
            if (yval > geom.ymid[j]) {
                double ytmp = (geom.ymid[j] - yval) * 1e-3;  // mm->m
                double ydot = state.vel[geom.surfp[i][j]].uy;
                fy[i][j] += (sp.kc1 * omg2 * ytmp * (1.0 + sp.kc2 * omg2 * ytmp * ytmp) - sp.kc3 * ydot) * geom.sarea[i][j] * 1e-6 ;
                contactFlag = true;
            }
        }
    }
}

void ForceCalculator::calcDis() {

    contactFlag = false;

    // 単位質量の仮定
    double mass = sp.mass;

    for (int i = 1; i < nxsup; ++i) {            // 2..nxsup (0-indexなので1スタート)
        for (int j = 1; j < geom.nsurfz - 1; ++j) { 
            fdis[i][j] = 0.0;

            int pid = geom.surfp[i][j];
            if (pid < 0) continue;

            double v_now = state.disp[pid].uy;
            double v_next = state.predictedDisp[pid].ufy;           // Runge-Kuttaで予測した変位
            

            if (v_now <= geom.ymid[j] && v_next > geom.ymid[j]) {
                double vc = (v_next - v_now) / sp.dt; // mm/s -> m/s は必要なら換算
                vc *= 1e-3;

                fdis[i][j] = -mass * vc / (sp.dt * geom.nPoints);
                fy[i][j] += fdis[i][j];

                contactFlag = true;
            }
        }
    }
}

double ForceCalculator::findMinHarea() {
    return *std::min_element(state.harea.begin(), state.harea.end());
}

int ForceCalculator::findNsep(double minH) {
    for (int i = 1; i < geom.nxsup; i++) {
        if (std::fabs(state.harea[i] - minH) < 1e-16 || state.harea[i] <= 0.0) {
            return i+1;
        }
    }
    return geom.nxsup;
}


void ForceCalculator::outputForceVectors(int step) const {
    std::ostringstream stepStr;
    stepStr << std::setw(4) << std::setfill('0') << step;

    std::ofstream fout("../output2/force_" + stepStr.str() + ".csv");
    fout << "x,y,z,Fx,Fy,Fz\n";  // CSVヘッダー

    for (int i = 1; i < geom.nxsup - 1; ++i) {
        for (int j = 1; j < geom.nsurfz - 1; ++j) {
            int pid = geom.surfp[i][j];
            if (pid < 0) continue;

            const auto &p = geom.points[pid];
            fout << p.x << "," << p.y << "," << p.z << ","
                 << fx[i][j] << "," << fy[i][j] << "," << fz[i][j] << "\n";
        }
    }

    //std::cout << "[Output] force vectors written for step " << step << std::endl;
}
