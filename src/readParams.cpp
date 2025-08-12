#include "SimulationParams.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <limits>


void SimulationParameters::loadFromFile(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string tmp;

    // nmode
    std::getline(fin, tmp); // コメント行
    std::getline(fin, tmp); // コメント行
    fin >> nmode;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // nsurfz
    std::getline(fin, tmp);
    fin >> nsurfz;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // nstep
    std::getline(fin, tmp);
    fin >> nstep;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // nwrite
    std::getline(fin, tmp);
    fin >> nwrite;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // dt
    std::getline(fin, tmp);
    fin >> dt;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // zeta
    std::getline(fin, tmp);
    fin >> zeta;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // kc1, kc2
    std::getline(fin, tmp);
    fin >> kc1 >> kc2;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // ncount
    std::getline(fin, tmp);
    fin >> ncount;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // freqFile
    std::getline(fin, tmp);
    std::getline(fin, freqFile);

    // modeFile
    std::getline(fin, tmp);
    std::getline(fin, modeFile);

    // surfFile
    std::getline(fin, tmp);
    std::getline(fin, surfFile);

    // inputDir
    std::getline(fin, tmp);
    std::getline(fin, inputDir);

    // resultDir
    std::getline(fin, tmp);
    std::getline(fin, resultDir);

    // iforce
    std::getline(fin, tmp);
    fin >> iforce;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (iforce == 1) {
        // forcef
        std::getline(fin, tmp);
        fin >> forcef;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // famp
        std::getline(fin, tmp);
        fin >> famp;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    } else if (iforce == 0) {
        // ps
        std::getline(fin, tmp);
        fin >> ps;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // rho
        std::getline(fin, tmp);
        fin >> rho;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // mu
        std::getline(fin, tmp);
        fin >> mu;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // mass
    std::getline(fin, tmp);
    fin >> mass;
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}
