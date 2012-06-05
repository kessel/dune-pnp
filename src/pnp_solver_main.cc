

#include "../config.h"
#include "pnp_solver_main.hh"


#include <iostream>
#include <dune/grid/sgrid.hh>
#include<dune/grid/common/gridinfo.hh>
#include "GridLearning.hh"

#include "sysparams.hh"
#include "potential.hh"

#include<dune/grid/io/file/gmshreader.hh>


#include<dune/grid/uggrid.hh>
#include<dune/grid/uggrid/uggridfactory.hh>


#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/pk2dfem.hh>



#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// instationary
#include<dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>

#include<dune/pdelab/instationary/onestep.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
// /instationary


#include <dune/grid/io/file/gnuplot.hh>

#include "btype.hh"
#include "dirichlet_bc.hh"
#include "diff_operator.hh"
#include "example03_toperator.hh"
#include "example03_operator.hh"


#include "stationary_diffusion.hh"
#include "stationary_pnp.hh"


PnpSolverMain::PnpSolverMain(Dune::MPIHelper &helper_) : helper(helper_) {
}

PnpSolverMain::~PnpSolverMain() {
}

void PnpSolverMain::run(std::string configfile) {

  Sysparams s;
  s.readConfigFile(configfile);

  std::string meshfile = s.meshfile;

  double dt = s.dt;
  double tend = s.tend;

  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  const int dimgrid=2;
 
  typedef Dune::UGGrid<dimgrid> GridType;
  Dune::GridFactory<GridType> factory;

   
  if(helper.rank() == 0)
  {
    // read a gmsh file
    Dune::GmshReader<GridType> gmshreader;
    gmshreader.read(factory, s.meshfile, boundaryIndexToEntity, elementIndexToEntity, true, true);
  }


    // read a gmsh file

  Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> colCom(helper.getCommunicator());

  // Here go the main settings what we want to use:
  const int dim=2;

    // Communicate boundary vector
  int size = boundaryIndexToEntity.size();
  colCom.broadcast (&size, 1, 0);
  if (helper.rank() > 0)
    boundaryIndexToEntity.resize(size);
  colCom.broadcast(&boundaryIndexToEntity[0],size,0);

  std::cout << helper.rank() << " " << boundaryIndexToEntity.size() << std::endl;



  GridType* grid = factory.createGrid();

  std::cout << "Grid has been modified by load balancing: " << grid->loadBalance() << std::endl;

//  Dune::gridinfo(grid);

  typedef GridType::LeafGridView GV;

  GV gv = grid->leafView();

  stationary_pnp<GV>(s, gv, boundaryIndexToEntity, elementIndexToEntity, helper);

  /*

  typedef GV::Grid::ctype Coord;
  typedef double Real;
  
  Real time = 0;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  typedef BCType<GV, std::vector<int> > B;                                         // boundary condition type
  B b(gv, boundaryIndexToEntity);

  typedef GFS::ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); 

  typedef BCExtension<GV, Real,  std::vector<int> > BCE;
  BCE bce(gv, boundaryIndexToEntity);

  typedef GFS::VectorContainer<Real>::Type U;
  U uold(gfs ,0.0);

//  Dune::PDELab::interpolate(bce,gfs ,uold);
  

  // This is example now!
    typedef Example03LocalOperator<B> LOP; 
  LOP lop(b,4);                                                 // local operator r
  typedef Example03TimeLocalOperator TLOP; 
  TLOP tlop(4);                                                 // local operator m
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,U,GFS,GFS,LOP,TLOP,CC,CC,MBE> IGOS;
  IGOS igos(gfs,cc,gfs,cc,lop,tlop);                            // new grid operator space

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  // <<<6>>> Solver for linear problem per stage
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGOS,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igos,ls,1e-10);

  // <<<7>>> time-stepper
  Dune::PDELab::Alexander2Parameter<Real> method;               // defines coefficients
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,U,U> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(2);                                     // time stepping scheme

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example03_Q2");              // append number to file name
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF udgf(gfs,uold);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
    fn.increment();                                             // increase file number
  }

  // <<<9>>> time loop
  U unew(gfs,0.0);                                              // solution to be computed
  while (time<tend-1e-8) {
      // do time step
     // b.setTime(time+dt);                                       // compute constraints
      cc.clear();                                               // for this time step
      Dune::PDELab::constraints(b,gfs,cc);
      osm.apply(time,dt,uold,bce,unew);                           // do one time step

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
      DGF udgf(gfs,unew);
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
      fn.increment();

      uold = unew;                                              // advance time step
      time += dt;
    }
*/


}

