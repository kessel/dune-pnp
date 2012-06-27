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
    bool cylindrical;
    double l_b;
    int linearSolverIterations;
    double newtonReassembleThreshold;
    double newtonReduction;
    double newtonMinLinearReduction;
    double newtonMaxIterations;
    double newtonLineSearchMaxIteration;
    double c0;
    double tau;
    int outputFreq;
    int nSteps;
    int potentialUpdateFreq;

    int printStiffnessMatrix;
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
