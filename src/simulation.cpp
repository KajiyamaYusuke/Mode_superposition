
#define _USE_MATH_DEFINES  
#include "Simulation.h"
#include <iostream>
#include <algorithm>
#include <cmath>

void Simulation::initialize() {
    std::cout << "[Simulation] Initializing..." << std::endl;
    std::string err ="error";

    params.loadFromFile("/home/kajiyama/code/simulation/input/param.txt", err );

    geom.loadFromVTK("/home/kajiyama/code/simulation/input/mode_renewal3mm.vtu");
    geom.surfExtractFromNAS("/home/kajiyama/code/simulation/input/surface_data_renewal.nas",12,21);
    //geom.surfExtract("/home/kajiyama/code/simulation/input/surface.txt", 20);
    geom.surfArea();
    //geom.print();
 
    geom.jtypes[5] = 3;   // 三角形
    geom.jtypes[9] = 4;   // 四角形
    geom.jtypes[10] = 4;
    geom.jtypes[13] = 6;  // 六面体

    mdata.initialize(params.nmode, geom);
    mdata.loadFromVTU("/home/kajiyama/code/simulation/input/mode_renewal3mm.vtu", geom);
    mdata.loadFreqDamping("/home/kajiyama/code/simulation/input/frequency_renewal3mm.txt");


    mdata.normalizeModes( params.mass, geom);

    state.initialize(geom.nPoints, params.nmode, params.nstep, geom);

    // ForceCalculator 初期化
    fCalc.initialize(); 


    std::cout << "[Simulation] Initialization complete." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Running..." << std::endl;

    // nSteps+1 に対応
    state.qf.resize(mdata.nModes, 0.0);
    state.qfdot.resize(mdata.nModes, 0.0);


    double P = 1;
    int num = 1;

    std::ofstream fa("../output/area.dat");
    std::ofstream fu("../output/velocity.dat");
    std::ofstream fp("../output/pressure.dat");

    fa << "# x[m]  area[m^2]\n";
    fu << "# x[m]  velocity[m/s]\n";
    fp << "# x[m]  pressure[Pa]\n"; 

    std::vector<double> zeta(mdata.nModes, 0);
    double omega1 = 200 * 2 * M_PI;
    double omega2 = 400 * 2 * M_PI;
    double alpha = 2*omega1*omega2*((0.0015*omega2 - 0.0025*omega1)/(omega2*omega2 - omega1*omega1));
    double beta = 2*(omega2*0.0025 - omega1* 0.0015)/(omega2*omega2 - omega1*omega1);

    for ( int i = 0; i < mdata.nModes; ++i){
        zeta[i] = 1/2*(alpha/(2.0 * M_PI * mdata.frequencies[i]) + beta * 2.0 * M_PI * mdata.frequencies[i]);
    }

    for (int n = 0
        ; n < params.nstep; n++) {
        double t = n * params.dt;


        // 2. 断面積や角度を更新
        state.calcArea(geom);


        fCalc.calcForce(t, n);

        if (n % 100 == 0) {
            fCalc.outputForceVectors(n);
        }

        fCalc.contactForce();



        if ( n%20 == 0){
            fa <<std::setw(3)<< n;
            fp <<std::setw(3)<< n;
            for (int i = 0; i < geom.nxsup; ++i) {
                
                fa << " " <<std::setw(8)<< state.harea[i] << " ";
                fp << " " <<std::setw(8)<< fCalc.psurf[i] << " ";
            }
            fa << "\n";
            fp << "\n";
            fu << n << " " << fCalc.Ug[n] << "\n";
        }


        int icont = 0;

        for (int icont = 1; icont <= params.ncont; ++icont) {

            // 4. モード力への変換
            fCalc.f2mode();



            // 5. 時間積分（RK4）
            for (int i = 0; i < mdata.nModes; i++) {
                double f = fCalc.fi[i];
                double q = state.q[i];
                double qdot = state.qdot[i];
                double omg = 2.0 * M_PI * mdata.frequencies[i];
                double qf, qfdot;

                integrator.rungeStep(f, q, qdot, params.dt, omg, zeta[i], qf, qfdot);

                state.qf[i]    = qf;
                state.qfdot[i] = qfdot;

            }


            // 6. モード変位 → 節点変位
            state.mode2uf(geom, mdata, n+1);

            // calculate dissipation force for contact
            fCalc.contactFlag = false;
            fCalc.calcDis();

            if (!fCalc.contactFlag) break;  // contactFlg == false の場合はループを抜ける
        }
    
        state.uf2u();

        if( n % 20 == 0){
            writeVTK(num, geom, state, "../result", 200);
            num++;
        }

    }

    std::cout << "[Simulation] Run complete." << std::endl;
} 

void Simulation::writeVTK(int step, const Geometry& geom, const State& state, const std::string& rdir, int nwrite) {
    // ファイル名
    std::ostringstream num;
    num << std::setw(4) << std::setfill('0') << step;
    std::string filename = rdir + "/deform" + num.str() + ".vtu";

    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

/*     std::cout << "step: " << step * nwrite << std::endl;
    std::cout << "output: " << filename << std::endl;  */

    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    fout << "  <UnstructuredGrid>\n";
    fout << "    <Piece NumberOfPoints=\"" << geom.nPoints 
         << "\" NumberOfCells=\"" << geom.nCells << "\">\n";

    // Points
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nPoints; i++) {
        fout << std::scientific << std::setprecision(6)
             << state.disp[i].ux << " " << state.disp[i].uy << " " << state.disp[i].uz << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Points>\n";

    // Cells
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        int nverts = geom.jtypes[geom.types[i]];
        for (int j = 0; j < nverts; j++) {
            fout << geom.connect[i][j] << " ";
        }
        fout << "\n";
    }
    fout << "        </DataArray>\n";

    fout << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        fout << geom.offsets[i] << "\n";
    }
    fout << "        </DataArray>\n";

    fout << "        <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < geom.nCells; i++) {
        fout << geom.types[i] << "\n";
    }
    fout << "        </DataArray>\n";
    fout << "      </Cells>\n";

    fout << "    </Piece>\n";
    fout << "  </UnstructuredGrid>\n";
    fout << "</VTKFile>\n";

    fout.close();
}
