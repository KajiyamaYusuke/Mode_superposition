
#define _USE_MATH_DEFINES  
#include "Simulation.h"
#include "wavwrite.h"
#include <iostream>
#include <algorithm>
#include <cmath>

void Simulation::initialize() {
    std::cout << "[Simulation] Initializing..." << std::endl;
    std::string err ="error";

    params.loadFromFile("../input/param.txt", err );

    geom.loadFromVTK("../input/M5/M5_mode_T2_d2_b15c4.vtu");

    geom.surfExtractFromNAS("../input/M5/M5_surface_T2_d2.nas",68,70);

    geom.surfArea();
    geom.print();
 
    geom.jtypes[5] = 3;   // 三角形
    geom.jtypes[9] = 4;   // 四角形
    geom.jtypes[10] = 4;
    geom.jtypes[13] = 6;  // 六面体

    mdata.initialize(params.nmode, geom);

    mdata.loadFromVTU("../input/M5/M5_mode_T2_d2_b15c4.vtu", geom);

    mdata.loadFreqDamping("../input/M5/M5_freq_T2_d2_b15c4.txt");



    mdata.normalizeModes( params.mass, geom);
    

    state.initialize(geom.nPoints, params.nmode, params.nstep, geom);

    // ForceCalculator 初期化
    fCalc.initialize(); 

    std::ofstream initDbg("../output/debug_step20_30.txt", std::ios::trunc);
    if (initDbg) {
        initDbg << "=== Detailed Parameter Log (Steps 20-30) ===\n";
        initDbg.close();
    }


    std::cout << "[Simulation] Initialization complete." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Running..." << std::endl;

    // nSteps+1 に対応
    state.qf.resize(mdata.nModes, 0.0);
    state.qfdot.resize(mdata.nModes, 0.0);
    state.qfddot.resize(mdata.nModes, 0.0);

    double P = 1;
    int num = 0;

    std::ofstream fa("../output/area.dat");
    std::ofstream fu("../output/displace.dat");
    std::ofstream fu2("../output/displace2.dat");
    std::ofstream fp("../output/pressure.dat");
    std::ofstream fpv("../output/pressure_vt.dat");
    std::ofstream fuv("../output/airflow_vt.dat");

    fa << "# x[m]  area[m^2]\n";
    fu << "# x[m]  displace\n";
    fu2 << "# x[m] displace\n";
    fp << "# x[m]  pressure[Pa]\n"; 
    fpv << "# x[m]  pressure[Pa]\n";
    fuv << "# x[m]  airflow[l/s]\n";

    std::vector<double> zeta(mdata.nModes, 0);
    double omega1 = 50 * 2 * M_PI;
    double omega2 = 250 * 2 * M_PI;
    double alpha = 2*omega1*omega2*((0.0015*omega2 - 0.01*omega1)/(omega2*omega2 - omega1*omega1));
    double beta = 2*(omega2*0.01 - omega1* 0.0015)/(omega2*omega2 - omega1*omega1);

    double minDist2 = 1e2;
    int nearestIdx = -1;
    for (int i = 0; i < geom.nsurfl; ++i) {
        for (int j = 0; j < geom.nsurfz; ++j) {
            int idx = geom.surfp[i][j];
            double dx = geom.points[idx].x - 10;
            double dz = geom.points[idx].z - 8.6;
            double dist2 = dx*dx + dz*dz ;


            if (dist2 < minDist2) {
                minDist2 = dist2;
                nearestIdx = idx;
            }
        }
    }
    minDist2 = 1e2;
    int nearestIdx2 = -1;
    for (int i = 0; i < geom.nsurfl; ++i) {
        for (int j = 0; j < geom.nsurfz; ++j) {
            int idx = geom.surfp[i][j];
            double dx = geom.points[idx].x - 7.2;
            double dz = geom.points[idx].z - 8.6;
            double dist2 = dx*dx + dz*dz ;

            if (dist2 < minDist2) {
                minDist2 = dist2;
                nearestIdx2 = idx;
            }
        }
    }
    std::cout<<"idx="<<geom.points[nearestIdx].x<<", "<<geom.points[nearestIdx].y<<", "<<geom.points[nearestIdx].z<<"\n";
    std::cout<<"idx2="<<geom.points[nearestIdx2].x<<", "<<geom.points[nearestIdx2].y<<", "<<geom.points[nearestIdx2].z<<"\n";

    for ( int i = 0; i < mdata.nModes; ++i){
        zeta[i] = 1/2*(alpha/(2.0 * M_PI * mdata.frequencies[i]) + beta * 2.0 * M_PI * mdata.frequencies[i]);
    }

    state.mode2uf(geom, mdata, 0); 
    state.uf2u(); // dispを更新
    writeVTK(num, geom, state, "../result", 1); // step 0 を出力
    num++;
    std::cout << "[Simulation] Output step 0 (Initial State). check this if bumpy." << std::endl;

    std::vector<double> soundSignal;
        soundSignal.reserve(params.nstep);

    for (int n = 0
        ; n < params.nstep; n++) {
        double t = n * params.dt;


        // 2. 断面積や角度を更新
        state.calcArea(geom);


        fCalc.calcForce(t, n);



        //fCalc.contactForce();

        static int debug_step = 0;
        if (debug_step < 10) {
            // 1ループ目は新規作成(上書き)、2ループ目以降は追記(app)モードで開く
            std::ios_base::openmode mode = (debug_step == 0) ? std::ios::out : std::ios::app;
            std::ofstream fsDeg("../output/Degree_debug.dat", mode);
            
            if (fsDeg) {
                fsDeg << "=== Step: " << debug_step << " ===" << std::endl;
                fsDeg << "j\tdegL\t\tdegR" << std::endl; // 見出し
                
                int test_j = 30; // 観察したいi座標
                // j全体のループ
                for (int i = 1; i < geom.nxsup - 1; ++i) {
                    fsDeg << i << "\t" 
                        << state.degree[0][i][test_j] << "\n" ;
                }
                fsDeg << std::endl;
                fsDeg.close();
            }
            debug_step++;
        }



        if ( n%5 == 0){
            fa <<std::setw(4)<< n;
            fp <<std::setw(4)<< n;
            for (int i = 0; i < geom.nxsup; ++i) {
                
                fa << " " <<std::setw(8)<< state.harea[i] << " ";
                fp << " " <<std::setw(8)<< fCalc.psurf[i] << " ";
            }
            fa << "\n";
            fp << "\n";
            
        }


        for (int icont = 1; icont <= params.ncont; ++icont) {

            // 4. モード力への変換
            fCalc.f2mode();

        // double rampDuration = 0.025; // 0.1秒かけて立ち上げ（状況により0.5など長くする）
        // double rampFactor = 1.0;
        
        // if (t < rampDuration) {
        //     rampFactor = t / rampDuration; 
        //     // 例: t=0なら0倍, t=0.05なら0.5倍, t=0.1以上なら1.0倍
        // }

        // // 計算されたモード力すべてに係数をかける
        // for(int i=0; i<mdata.nModes; ++i) {
        //     fCalc.fi[i] *= rampFactor;
        // }

            // 5. 時間積分（RK4）
/*             for (int i = 0; i < mdata.nModes; i++) {
                double f = fCalc.fi[i];
                double q = state.q[i];
                double qdot = state.qdot[i];
                double omg = 2.0 * M_PI * mdata.frequencies[i];
                double qf, qfdot;

                integrator.rungeStep(f, q, qdot, params.dt, omg, params.zeta, qf, qfdot);

                state.qf[i]    = qf;
                state.qfdot[i] = qfdot;

            } */
            
            // Newmark parameters (average acceleration)
            const double beta  = 0.275625;
            const double gamma = 0.55;


            for (int i = 0; i < mdata.nModes; ++i) {
                double f    = fCalc.fi[i];                      // モード力 (tilde f_i)
                double q    = state.q[i];                       // q_n
                double qdot = state.qdot[i];                    // qdot_n
                double qdd  = state.qddot[i];                   // qddot_n (保持しておく)
                double omega = 2.0 * M_PI * mdata.frequencies[i];

                double qf, qfdot, qfddot;
                integrator.newmarkStep(f, q, qdot, qdd, params.dt, omega, params.zeta,
                                    beta, gamma, qf, qfdot, qfddot);

                state.qf[i]    = qf;
                state.qfdot[i] = qfdot;
                state.qfddot[i] = qfddot;
            }




            // 6. モード変位 → 節点変位
            state.mode2uf(geom, mdata, n+1);



            // calculate dissipation force for contact
            fCalc.contactFlag = false;
            fCalc.calcDis();
            if (fCalc.contactFlag && fCalc.max_force_diff < 1.0e-6) { 
                break; 
            }

            if (!fCalc.contactFlag) break;  // contactFlg == false の場合はループを抜ける
        }

        if (n%20 ==0){
            //fu << n *params.dt << " "<<state.predictedDisp[nearestIdx].ufy - geom.points[nearestIdx].y<<" "<<state.predictedDisp[nearestIdx].ufx - geom.points[nearestIdx].x<< "\n";
            //fu2 << n *params.dt << " "<<state.predictedDisp[nearestIdx2].ufy - geom.points[nearestIdx2].y<<" "<<state.predictedDisp[nearestIdx2].ufx - geom.points[nearestIdx2].x<< "\n";
            fu << n *params.dt << " "<<state.predictedDisp[nearestIdx].ufy <<" "<<state.predictedDisp[nearestIdx].ufx << "\n";
            fu2 << n *params.dt << " "<<state.predictedDisp[nearestIdx2].ufy <<" "<<state.predictedDisp[nearestIdx2].ufx << "\n";
        }
    
        state.uf2u();

        if (n >= 200 && n <= 400) {
            std::ofstream dbgFile("../output/debug_step20_30.txt", std::ios::app);
            if (dbgFile) {
                // 小数点以下12桁まで高精度で出力し、微小なズレを逃さない
                dbgFile << std::scientific << std::setprecision(12);
                dbgFile << "================ Step " << n << " ================\n";
                
                // 1. 流体力学の主要パラメータ
                dbgFile << "[Fluid]\n";
                dbgFile << "  minHarea  : " << fCalc.minHarea[n] << "\n";
                dbgFile << "  currentUg : " << fCalc.currentUg << "\n";
                dbgFile << "  currentPg : " << fCalc.currentPg << "\n";
                // 圧力分布の代表点（メッシュ中央付近）
                int mid_i = geom.nxsup / 2;
                dbgFile << "  psurf[mid]: " << fCalc.psurf[mid_i] << "\n";

                // 2. モード力と構造力学（1次モード代表）
                dbgFile << "[Structure - Mode 0]\n";
                dbgFile << "  fiL[0]    : " << fCalc.fi[0] << "\n";
                dbgFile << "  qL[0]     : " << state.q[0] << "\n";
                dbgFile << "  qdotL[0]  : " << state.qdot[0] << "\n";
                dbgFile << "  qddotL[0] : " << state.qddot[0] << "\n";

                // 3. 実際の節点変位（モニター用の nearestIdxL を使用）
                dbgFile << "[Nodal Displacement (nearestIdxL)]\n";
                dbgFile << "  dispL.uy  : " << state.disp[nearestIdx].uy << "\n";
                // 速度や予測変位も確認
                dbgFile << "  velL.uy   : " << state.vel[nearestIdx].uy << "\n";
                
                // 4. 接触とループ制御
                dbgFile << "[System]\n";
                dbgFile << "  contact   : " << (fCalc.contactFlag ? "TRUE" : "FALSE") << "\n";
                
                dbgFile << "---------------------------------------\n";
                dbgFile.close();
            }
        }

        if( n % 20 == 0){
            writeVTK(num, geom, state, "../result", 20);
            fCalc.outputForceVectors(n);
            num++;
        }
        if ( n%5== 0){
            fpv <<std::setw(4)<< n;
            fpv << " " <<std::setw(8)<< fCalc.Pd[9] << " ";
            fpv << "\n";
            fuv <<std::setw(4)<< n;
            fuv << " " <<std::setw(8)<< fCalc.currentUg << " ";
            fuv << "\n";
        }
        soundSignal.push_back(fCalc.Pd[9]);

        if (t <= 0.5 && n%20 == 0) {
            static int step = 0;
            std::ostringstream filename;
            filename << "../output_force/contact_force_" << std::setw(6) << std::setfill('0') << step << ".txt";
            step++;
            // そのステップ用のファイルを開く
            std::ofstream fcOut(filename.str());
            if (fcOut) {
                fcOut << "Time: " << t << "\n";
                
                fcOut << "[Left]\n";
                for (int i = 0; i < geom.nxsup; ++i) {
                    for (int j = 0; j < geom.nsurfz; ++j) {
                        fcOut << std::scientific << fCalc.contactForce_ij[i][j] << " ";
                    }
                    fcOut << "\n";
                }
                fcOut.close(); // 書き終わったら閉じる
            }
        }
    }

    WavWriter::save(soundSignal, params.dt, "../output/output_sound.wav");

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

    std::cout << "step: " << step * nwrite << std::endl;
    std::cout << "output: " << filename << std::endl;   

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
