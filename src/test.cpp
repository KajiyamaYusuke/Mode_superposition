#include "Simulation.h"
#include <iostream>
#include <iomanip>
#include <iostream>

int main() {
    // 1. Simulation パラメータ設定
    Simulation sim;


    // 2. 初期化

    sim.initialize();
    //sim.params.print(std::cout);
    sim.params.iforce = 0;
    //sim.geom.print();

     

    // 4. シミュレーション実行
    sim.run();


    // 6. 簡単に q の最後の値を表示
    for (int i = 0; i < sim.mdata.nModes; i++) {
        //std::cout << "Mode " << i << " final q: " << sim.state.q[i] << std::endl;
    } 

    return 0;
}



