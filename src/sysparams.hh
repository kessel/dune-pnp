#ifndef SYSPARAMS_H
#define SYSPARAMS_H


#include <string>
#include <vector>

class Surface;

class Sysparams {
  public:
    Sysparams();
    ~Sysparams();
    void readConfigFile(std::string filename);
    double K;
    std::vector<Surface> surfaces;
    std::string meshfile;
    int n_surfaces;
    double dt;
    double tend;
    int verbosity;
};

class Surface {
  public:
    Surface(); 
    
    int coulombBtype;
    double coulombFlux;
    double coulombPotential;
    double coulombSigma;
    double coulombEpsilon;
    double coulombChargeability;

    int plusDiffusionBtype;
    double plusDiffusionFlux;
    double plusDiffusionConcentration;

    int minusDiffusionBtype;
    double minusDiffusionFlux;
    double minusDiffusionConcentration;
};

#endif
