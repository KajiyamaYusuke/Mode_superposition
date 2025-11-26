#pragma once
#include "Geometry.h"
#include "ModeData.h"
#include "State.h"
#include "SimulationParams.h"
#include <vector>

class ForceCalculator {
public:
    ForceCalculator(const Geometry& geom_, const ModeData& md_, State& st_, const SimulationParams& sp_);

    // 外力ベクトル
    std::vector<std::vector<double>> fx, fy, fz; // [nPoints][]
    std::vector<double> fi; // モード力 [nModes]

    double current_psub; // 現在の声門下圧
    double Qin;          // 供給流量 (m^3/s)
    double Volume;       // 声門下容積 (m^3)
    double K_air;        // 体積弾性率 (Pa)

    // 出力用
    std::vector<double> Ug, minHarea; 
    std::vector<double> psurf; 
    bool contactFlag = false;

    // 初期化（メモリ確保など）
    void initialize();

    // 時刻 t, ステップ n での力計算
    void calcForce(double t, int step);
    void f2mode();
    void contactForce();
    void calcDis();
    void outputForceVectors(int step) const;

private:
    const Geometry& geom;
    const ModeData& modeData;
    State& state;
    const SimulationParams& sp;

    // 内部バッファ
    std::vector<std::vector<double>> fdis;
    int nxsup = 0;


    // ヘルパー
    double findMinHarea();
    int findNsep(double minH);
    double firstOmega() const;
};

