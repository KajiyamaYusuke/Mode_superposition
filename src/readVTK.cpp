#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <cctype>

#include "ModeShape.h"

void ModeShape::loadFromVTK(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Failed to open VTU file: " + filename);
    }

    std::string line;
    int pointCount = 0, cellCount = 0;

    // 一時バッファ
    std::vector<double> pointsBuffer;
    std::vector<int> connectivityBuffer;
    std::vector<int> offsetsBuffer;
    std::vector<int> typesBuffer;

    // モードデータを一時保存
    std::vector<std::vector<std::array<double,3>>> modesTemp;

    while (std::getline(file, line)) {
        // Points
        if (line.find("<Points>") != std::string::npos) {
            while (std::getline(file, line)) {
                if (line.find("<DataArray") != std::string::npos) {
                    // 座標値読み込み
                    pointsBuffer.clear();
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        double px, py, pz;
                        while (ss >> px >> py >> pz) {
                            pointsBuffer.push_back(px);
                            pointsBuffer.push_back(py);
                            pointsBuffer.push_back(pz);
                        }
                    }
                }
                if (line.find("</Points>") != std::string::npos) break;
            }
        }

        // Cells
        else if (line.find("<Cells>") != std::string::npos) {
            while (std::getline(file, line)) {
                if (line.find("connectivity") != std::string::npos) {
                    connectivityBuffer.clear();
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) connectivityBuffer.push_back(id);
                    }
                }
                else if (line.find("offsets") != std::string::npos) {
                    offsetsBuffer.clear();
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) offsetsBuffer.push_back(id);
                    }
                }
                else if (line.find("types") != std::string::npos) {
                    typesBuffer.clear();
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        int id;
                        while (ss >> id) typesBuffer.push_back(id);
                    }
                }
                if (line.find("</Cells>") != std::string::npos) break;
            }
        }

        // PointData
        else if (line.find("<PointData>") != std::string::npos) {
            while (std::getline(file, line)) {
                if (line.find("DataArray") != std::string::npos && line.find("mode") != std::string::npos) {
                    std::vector<std::array<double,3>> modePoints;
                    while (std::getline(file, line) && line.find("</DataArray>") == std::string::npos) {
                        std::istringstream ss(line);
                        double vx, vy, vz;
                        while (ss >> vx >> vy >> vz) {
                            modePoints.push_back({vx, vy, vz});
                        }
                    }
                    modesTemp.push_back(modePoints);
                }
                if (line.find("</PointData>") != std::string::npos) break;
            }
        }
    }

    // 格納
    nPoints = pointsBuffer.size() / 3;
    x.resize(nPoints);
    y.resize(nPoints);
    z.resize(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        x[i] = pointsBuffer[i*3];
        y[i] = pointsBuffer[i*3+1];
        z[i] = pointsBuffer[i*3+2];
    }

    nCells = offsetsBuffer.size();
    connect.clear();
    connect.resize(nCells);
    int start = 0;
    for (int i = 0; i < nCells; ++i) {
        int end = offsetsBuffer[i];
        connect[i].assign(connectivityBuffer.begin() + start,
                          connectivityBuffer.begin() + end);
        start = end;
    }
    offsets = offsetsBuffer;
    types = typesBuffer;

    nModes = modesTemp.size();
    modes = modesTemp;
}
