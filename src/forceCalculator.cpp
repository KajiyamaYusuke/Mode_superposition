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

    // 実験条件の設定例
    double flow_L_min = 160.0; 
    Qin = (flow_L_min / 60.0) * 1.0e-3; // L/min -> m^3/s 
    
    // チューブの体積 
    double r = 0.0125; 
    double L = 0.25;
    Volume = M_PI * r * r * L; 
    
    K_air = 1.13 * 340.0 * 340.0; // rho * c^2
    
    current_psub = 0.0; // 初期圧力

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
        
        
        for (int i = 1; i < nxsup-1; i++) {
            for (int j = 1; j < nsurfz-1; j++) {

                fx[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
                //fy[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
                //fz[i][j] = sp.famp * std::sin(2.0* M_PI *sp.forcef* t);
            }
        }
    } else if (sp.iforce == 0) {

        // ==== 1D flow model ====
        minHarea[n] = *std::min_element(state.harea.begin(), state.harea.end());


        // separation point
        int nsep = findNsep(minHarea[n]) ;

        double Qout = 0.0;

        // 流量 Ug
        if (minHarea[n] > 0.0 && current_psub > 0.0) {
            Qout = std::sqrt(2.0 * current_psub / sp.rho) * minHarea[n];
        } else {
            Qout = 0.0;
        }

        Ug[n] = Qout;

        double dP = (K_air / Volume) * (Qin - Qout) * sp.dt;
        current_psub += dP;

        // psurf 計算
        std::fill(psurf.begin(), psurf.end(), 0.0);
        psurf[0] = current_psub;

        if (minHarea[n] - 0.0 > 1e-3) {
            for (int i = 1; i < nsep; i++) {
                double dx = std::abs(geom.points[geom.surfp[i][ int(nsurfz/2)]].x - geom.points[geom.surfp[i-1][int(nsurfz/2)]].x);
                double h  = (state.harea[i] + state.harea[i-1]) / (2.0 * geom.zmax);

                double ha1 = std::max(state.harea[i-1], 1e-6);
                double ha2 = std::max(state.harea[i], 1e-6);

                psurf[i] = psurf[i-1] + 0.5 * sp.rho * Ug[n] * Ug[n] * (1.0 / (state.harea[i-1]*state.harea[i-1]) - 1.0 / (state.harea[i]*state.harea[i]))
                                      - 12.0 * sp.mu * dx / (geom.zmax * h * h * h) * Ug[n] * 1e3;

            }
        } else {
            for (int i = 1; i < nsep-1; i++) {
                psurf[i] = current_psub;
            }
            //std::cout<<"cloze, n = "<< n <<std::endl;
        }

        if ( n%100 == 0){
            //std::cout<<" n = "<<n<<" minarea = "<<minHarea[n]<<" psurf[20]="  <<psurf[20]<<std::endl;
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
        if ( n%200 == 0){
            //std::cout<<" n = "<<std::setw(4)<<n<<" fx[10][15] = "<<fx[10][15]<<" fy[10][15]="  <<fy[10][15]<<" nsep = "<<nsep<<std::endl;
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
            int node = geom.surfp[i][j];

            double y     = state.disp[node].uy;        // 現在変位
            double ydot  = state.vel[node].uy;         // 現在速度
            double ymid  = geom.ymid[j];
            double yhat  = y + sp.dt * ydot;           // 予測位置


            // 接触状態を判定
            bool contact_now    = (y > ymid);          // 現時点で接触
            bool contact_future = (yhat > ymid);       // 次ステップで接触

            if (!contact_now && !contact_future) {

                continue;
            }

            double pen = (ymid - y) * 1e-3; 

            double f_contact = sp.kc1 * omg2 * pen * (1.0 + sp.kc2 * omg2 * pen * pen);

            double kc3eff = contact_now ? sp.kc3 : 0.0;

            double f_damp = - kc3eff * ydot;

            double f_total = (f_contact + f_damp) * geom.sarea[i][j] * 1e-6;

            fy[i][j] += f_total;
            contactFlag = true;
        }
    } 

        // パラメータ設定（調整可）
    /* const double v_th   = 0.05;     // [m/s] 以下なら「密着」とみなす速度閾値
    const double eps_adh = 0.2;     // [–] 粘着の強さ（0〜1）
    const double delta_adh = 5e-2;  // [m] 接触面からどれくらい離れても引き合うか

    for (int i = 0; i < geom.nxsup; ++i) {
        for (int j = 1; j < geom.nsurfz - 1; ++j) {
            int idx = geom.surfp[i][j];
            double yval = state.disp[idx].uy;   // 節点変位
            double ydot = state.vel[idx].uy;    // y速度
            double gap = geom.ymid[j] - yval;   // +: 離れ, -: 食い込み
            double fy_tmp = 0.0;

            bool inContact = (gap < 0.0);  // ペナルティ判定

            // ---- 通常接触力（ペナルティ法）----
            if (inContact) {
                double ytmp = gap * 1e-3; // mm→m
                fy_tmp += (sp.kc1 * omg2 * ytmp * (1.0 + sp.kc2 * omg2 * ytmp * ytmp)
                          - sp.kc3 * ydot);
                contactFlag = true;
            }

            // ---- 粘着的な保持条件（離れかけでも保持）----
            else if (contactFlag && fabs(ydot) < v_th && gap < delta_adh) {
                // 距離に応じて弱まる“粘着力”
                double Fadh = -eps_adh * sp.kc1 * omg2 * (delta_adh - gap);
                fy_tmp += Fadh;
                contactFlag = true;
            } else {
                contactFlag = false;
            }

            // ---- 力を加算 ----
            if (contactFlag) {
                fy[i][j] += fy_tmp * geom.sarea[i][j] * 1e-6;
                contactFlag = true;
            }
        }
    } */
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
    int min_idx = geom.nxsup;

for (int i = 1; i < geom.nxsup; i++) {
        if (std::fabs(state.harea[i] - minH) < 1e-16 || state.harea[i] <= 0.0) {
            min_idx = i; 
            break; // 見つけたらループを抜ける
        }
    }

    int extended_nsep = std::min(min_idx + 1, geom.nxsup); 
    return extended_nsep;
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
