
#include "sysparams.hh"


#include <iostream>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>


Sysparams::Sysparams() {
}

Sysparams::~Sysparams() {
}

void Sysparams::readConfigFile(std::string config_file) {

  Dune::ParameterTree configuration;
  Dune::ParameterTreeParser parser;
  
  try{
      parser.readINITree( config_file, configuration );
  }
  catch(...){
      std::cerr << "Could not read config file \"" 
                << config_file << "\"!" << std::endl;
      exit(1);
  }

  // Mesh file must be specified
  meshfile = configuration.get<std::string>("mesh.filename");

  // Create particles
  n_surfaces = configuration.get<int>("system.n_surfaces");
  
  Surface temp;
  for (int i = 0; i < n_surfaces; i++)
    surfaces.push_back(temp);

  for(int i = 0; i < n_surfaces; i++)
  {
    std::string p_name = "surface_";
    std::ostringstream s;
    s << p_name << i; // we don't need to set 0 and 1
                        // remember that vector starts with 0, so access the boundaries with -2
                        // TODO find a clever way!
    p_name = s.str();
    surfaces[i].coulombBtype=configuration.get<int>(p_name+".coulombBtype");
    switch ( surfaces[i].coulombBtype ) {
      case 0:
        surfaces[i].coulombPotential =  configuration.get<double>(p_name+".coulombPotential");
        break;
      case 1:
        surfaces[i].coulombFlux =  configuration.get<double>(p_name+".coulombFlux");
        break;
    }
    surfaces[i].plusDiffusionBtype =configuration.get<int>(p_name+".plusDiffusionBtype");
    switch ( surfaces[i].plusDiffusionBtype ) {
      case 0:
        surfaces[i].plusDiffusionConcentration =  configuration.get<double>(p_name+".plusDiffusionConcentration");
        break;
      case 1:
        surfaces[i].plusDiffusionFlux =  configuration.get<double>(p_name+".plusDiffusionFlux");
        break;
    }
    surfaces[i].minusDiffusionBtype =configuration.get<int>(p_name+".minusDiffusionBtype");
    switch ( surfaces[i].minusDiffusionBtype ) { 
      case 0:
        surfaces[i].minusDiffusionConcentration =  configuration.get<double>(p_name+".minusDiffusionConcentration");
        break;
      case 1:
        surfaces[i].minusDiffusionFlux =  configuration.get<double>(p_name+".minusDiffusionFlux");
        break;
    }

  }
}


Surface::Surface() {
  coulombBtype = 1;
  coulombPotential = 0;
  coulombFlux = 0;
  coulombSigma = 0;
  coulombEpsilon = 1;
  coulombChargeability = 0;
  
  plusDiffusionBtype = 1;
  plusDiffusionConcentration = 0;
  plusDiffusionFlux = 0;
  
  minusDiffusionBtype = 1;
  minusDiffusionConcentration = 0;
  minusDiffusionFlux = 0;
}

