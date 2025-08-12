#include <iostream>
#include <iomanip>  
#include <vector>
#include <array>


#include "ModeShape.h"

void ModeShape::printVTK() const {
    std::cout << "Number of Points: " << nPoints << std::endl;
    std::cout << "Coordinates (first 10 points):" << std::endl;
    for (int i = 0; i < std::min(nPoints, 10); ++i) {
        std::cout << "  Point " << i << ": ("
                  << std::setprecision(6) << x[i] << ", "
                  << y[i] << ", " << z[i] << ")" << std::endl;
    }

    std::cout << "\nNumber of Cells: " << nCells << std::endl;
    std::cout << "Cell connectivity (first 20 cells):" << std::endl;
    for (int i = 0; i < std::min(nCells, 20); ++i) {
        std::cout << "  Cell " << i << ": ";
        for (int idx : connect[i]) {
            std::cout << idx << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nOffsets (first 10): ";
    for (int i = 0; i < std::min(nCells, 10); ++i) {
        std::cout << offsets[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "\nTypes (first 10): ";
    for (int i = 0; i < std::min(nCells, 10); ++i) {
        std::cout << types[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "\nNumber of Modes: " << nModes << std::endl;
    for (int m = 0; m < std::min(nModes, 5); ++m) {
        std::cout << "Mode " << m+1 << " (first 15 points):" << std::endl;
        for (int i = 0; i < std::min(nPoints, 15); ++i) {
            const auto& vec = modes[m][i];
            std::cout << "  Point " << i << ": (" 
                      << vec[0] << ", " << vec[1] << ", " << vec[2] << ")" << std::endl;
        }
        std::cout << std::endl;
    }
}
