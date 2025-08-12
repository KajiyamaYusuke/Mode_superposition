#include <iostream>
#include <fstream>
#include <vector>
#include <array>


class ModeShape {
public:
    int nPoints = 0;
    int nCells = 0;
    int nModes = 20;

    std::vector<double> x, y, z;                  // 座標
    std::vector<std::vector<int>> connect;        // 8点接続のセル（各セルは8点）
    std::vector<int> offsets;
    std::vector<int> types;

    // mode[modeIndex][pointIndex][xyz(0,1,2)]
    std::vector<std::vector<std::array<double,3>>> modes; 

    std::vector<double> frequencies; // 固有振動数
    std::vector<double> dampingRatios;

    // モードの変位・速度・加速度を保存
    std::vector<std::vector<double>> q, qdot, qddot;

    void normalizeModes(double mass);
    void loadFromVTK(const std::string& filename);
    void printVTK() const;
};
