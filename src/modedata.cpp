#include "ModeData.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <regex>


void ModeData::initialize(int nModes_, const Geometry& geom) {
    nModes  = nModes_;
    int nPoints = geom.nPoints;

    // モード形状の初期化 [nModes][nPoints]
    modes.assign(nModes, std::vector<Displacement>(nPoints, Displacement()));

    // 固有振動数・減衰比の初期化
    frequencies.assign(nModes, 0.0);
    dampingRatios.assign(nModes, 0.0);
}

void ModeData::loadFromVTU(const std::string& filename, const Geometry& geom) {
    int nPoints = geom.nPoints;
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open VTU file: " + filename);

    std::string line;
    // 周波数ごとに {X, Y, Z} データをまとめる辞書
    std::map<std::string, std::vector<double>> xData, yData, zData;

    // X/Y/Z 成分を判定する正規表現
    std::regex rxX("Name=\"Displacement_field,_X-?component_@_([0-9.]+_Hz)\"");
    std::regex rxY("Name=\"Displacement_field,_Y-?component_@_([0-9.]+_Hz)\"");
    std::regex rxZ("Name=\"Displacement_field,_Z-?component_@_([0-9.]+_Hz)\"");

    std::string currentKey;
    std::vector<double>* currentBuffer = nullptr;

    while (std::getline(file, line)) {
        std::smatch match;

        if (std::regex_search(line, match, rxX)) {
            currentKey = match[1]; // 例: "87.739_Hz"
            xData[currentKey] = {};
            currentBuffer = &xData[currentKey];
            continue;
        } else if (std::regex_search(line, match, rxY)) {
            currentKey = match[1];
            yData[currentKey] = {};
            currentBuffer = &yData[currentKey];
            continue;
        } else if (std::regex_search(line, match, rxZ)) {
            currentKey = match[1];
            zData[currentKey] = {};
            currentBuffer = &zData[currentKey];
            continue;
        } else if (line.find("</DataArray>") != std::string::npos) {
            currentBuffer = nullptr;
            continue;
        }

        // DataArray内部の数値読み込み
        if (currentBuffer) {
            std::istringstream ss(line);
            double v;
            while (ss >> v) currentBuffer->push_back(v);
        }
    }

    // --- モード名（周波数）の一覧を抽出 ---
    std::vector<std::string> modeKeys;
    for (auto& [key, _] : xData) modeKeys.push_back(key);

    std::sort(modeKeys.begin(), modeKeys.end(), [](const std::string& a, const std::string& b) {
    auto parseHz = [](const std::string& s) {
        return std::stod(s.substr(0, s.find("_Hz"))); // "_Hz"前までをdoubleに変換
    };
        return parseHz(a) < parseHz(b);
    });

    nModes = modeKeys.size();
    modes.resize(nModes);

    for (int m = 0; m < nModes; ++m) {
        const std::string& key = modeKeys[m];

        const auto& x = xData[key];
        const auto& y = yData[key];
        const auto& z = zData[key];

        if (x.size() != nPoints || y.size() != nPoints || z.size() != nPoints) {
            throw std::runtime_error("Mode " + key + ": size mismatch with geometry points");
        }

        std::vector<Displacement> disp(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            disp[i] = { z[i], x[i], y[i] }; // 座標入れ替え
        }
        modes[m] = disp;

    }

    /* std::cout << "\n=== ModeData Debug Output ===\n";
    for (int m = 0; m < std::min(nModes, 3); ++m) {  // 最初の3モードだけ
        std::cout << "Mode " << m 
                << " (" << modeKeys[m] << "), " 
                << modes[m].size() << " points\n";

        for (int i = 0; i < std::min(nPoints, 10); ++i) {  // 各モードの最初の5点だけ
            const auto& d = modes[m][i];
            std::cout << "  Point " << i 
                    << ": (" << d.ux << ", " << d.uy << ", " << d.uz << ")\n";
        }
        std::cout << std::endl;
    } */


    std::cout << "ModeData: Loaded " << nModes 
              << " modes (" << nPoints << " points each)." << std::endl;
}


void ModeData::loadFreqDamping(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open frequency file: " + filename);

    frequencies.clear();
    dampingRatios.clear();

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream ss(line);
        double freqHz, omega, damping, Q;
        ss >> freqHz >> omega >> damping >> Q;

        frequencies.push_back(freqHz);
        dampingRatios.push_back(damping);
    }

    if (frequencies.size() != nModes) {
        throw std::runtime_error("Mode count mismatch between VTU and frequency file");
    }

    nModes = frequencies.size();

    std::cout << "ModeData: Loaded frequencies for " << nModes << " modes." << std::endl;
}

void ModeData::normalizeModes(double mass, const Geometry& geom){
    
    int nPoints = geom.nPoints;
    double nMass = mass / nPoints ;
    double ci;
    //std::cout<<"nMass = "<<nMass<<"\n";

    for (int imode=0; imode<nModes; imode++){
        ci = 0.0;
        for (int j=0; j<nPoints; j++){
            ci += modes[imode][j].ux*modes[imode][j].ux
                + modes[imode][j].uy*modes[imode][j].uy
                + modes[imode][j].uz*modes[imode][j].uz;
        }
        ci = 1.0 / std::sqrt( nMass * ci );
    
        for (int j=0; j<nPoints; j++){
            modes[imode][j].ux *= ci;
            modes[imode][j].uy *= ci;
            modes[imode][j].uz *= ci;
        }
    }
    
}


void ModeData::loadFromVTU_old(const std::string& filename, const Geometry& geom) {
    int nPoints = geom.nPoints;
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open VTU file: " + filename);

    std::string line;

    // 各モード番号ごとに変位ベクトルを格納
    std::map<std::string, std::vector<Displacement>> modeData;

    // 正規表現：<DataArray ... Name="mode01" ...>
    std::regex rxMode("Name=\"mode(\\d+)\"");

    std::string currentKey;
    std::vector<Displacement>* currentBuffer = nullptr;

    while (std::getline(file, line)) {
        std::smatch match;

        // 新しいモードの開始
        if (std::regex_search(line, match, rxMode)) {
            currentKey = "mode" + match[1].str();  // 例: "mode01"
            modeData[currentKey] = {};
            currentBuffer = &modeData[currentKey];
            continue;
        }

        // DataArray終了タグで終了
        if (line.find("</DataArray>") != std::string::npos) {
            currentBuffer = nullptr;
            continue;
        }

        // DataArray内部の数値読み込み
        if (currentBuffer) {
            std::istringstream ss(line);
            double x, y, z;
            while (ss >> x >> y >> z) {
                currentBuffer->push_back({z, x, y});
            }
        }
    }

    // --- モード名の一覧を抽出してソート ---
    std::vector<std::string> modeKeys;
    for (auto& [key, _] : modeData) modeKeys.push_back(key);

    std::sort(modeKeys.begin(), modeKeys.end(), [](const std::string& a, const std::string& b) {
        auto getNum = [](const std::string& s) {
            return std::stoi(s.substr(4)); // "mode"以降を整数に変換
        };
        return getNum(a) < getNum(b);
    });

    nModes = modeKeys.size();
    modes.resize(nModes);

    for (int m = 0; m < nModes; ++m) {
        const std::string& key = modeKeys[m];
        const auto& disp = modeData[key];

        if (disp.size() != nPoints) {
            throw std::runtime_error("Mode " + key + ": size mismatch with geometry points");
        }

        modes[m] = disp;
    }
}
