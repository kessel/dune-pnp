
#include <string>
#include"dune/common/mpihelper.hh"

class PnpSolverMain {
  public:
    PnpSolverMain(Dune::MPIHelper &helper_);
    ~PnpSolverMain();
    void run(std::string configfile);
  private:
    Dune::MPIHelper &helper;
};
