#include <iostream>
#include <vector>



class SimulationParameters {
public:
    int nmode;
    int nsurfz;
    int nstep;
    int nwrite;
    double dt;
    double zeta;
    double kc1, kc2;

    int ncount;

    std::string inputDir, resultDir;
    std::string freqFile, modeFile, surfFile;
    

    int iforce;
    double forcef;
    double famp;
    double ps;
    double rho;
    double mu;

    double mass;

    void loadFromFile(const std::string& filename);
    void print() const;
};

