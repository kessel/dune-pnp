#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"dune/common/exceptions.hh" // We use exceptions
#include"pnp_solver_main.hh"

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-pnp." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
        <<" processes!"<<std::endl;

    std::string f("hallo");

    PnpSolverMain p;
    p.run(f);


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
