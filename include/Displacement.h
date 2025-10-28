#pragma once

struct Displacement {
    double ux = 0.0;   // 現在のx変位
    double uy = 0.0;   // 現在のy変位
    double uz = 0.0;   // 現在のz変位

    double ufx = 0.0;  // 予測x変位
    double ufy = 0.0;  // 予測y変位
    double ufz = 0.0;  // 予測z変位
};
