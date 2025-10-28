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

