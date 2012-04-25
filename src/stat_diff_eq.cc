

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

#include <dune/grid/io/file/gnuplot.hh>

#include "btype.hh"
#include "dirichlet_bc.hh"
#include "diff_operator.hh"



PnpSolverMain::PnpSolverMain() {
}

PnpSolverMain::~PnpSolverMain() {
}

void PnpSolverMain::run(std::string configfile) {
  std::cout << "I want to solve!!" << std::endl;


  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  const int dimgrid=2;
 
  typedef Dune::UGGrid<dimgrid> GridType;
  Dune::GridFactory<GridType> factory;

    // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(factory, "mesh.msh", boundaryIndexToEntity, elementIndexToEntity, true, true);

  for (int i =0; i<10; i++) {
    std::cout << i << " " << boundaryIndexToEntity[i] << " " << elementIndexToEntity[i] << std::endl;
  }


  // Here go the main settings what we want to use:
  const int dim=2;

  GridType* grid = factory.createGrid();

//  Dune::gridinfo(grid);

  typedef GridType::LeafGridView GV;

  GV gv = grid->leafView();

  typedef GV::Grid::ctype Coord;
  typedef double Real;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;     // constraints class
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
  U u(gfs ,0.0);

  Dune::PDELab::interpolate(bce,gfs ,u);
  
  Sysparams s;
  
  typedef Force<Sysparams, Dune::FieldVector<GridType::ctype,dim> > F;
  F f(s);

  typedef DiffOperator<B, F> LOP;
  LOP lop(b, f, 5);



  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> assemble and solve linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  SLP slp(gos,u,ls,1e-10);
  slp.apply();

  // <<<7>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("yeah",Dune::VTKOptions::binaryappended);

  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write("yeah.dat"); 
  

  
//  GridLearning<GV> testclass = GridLearning<GV>(gv);


//  Potential<Sysparams, Dune::FieldVector<GridType::ctype,dim> > p(s);

//  for (double x = 0; x<2; x+=0.1) {
//    Dune::FieldVector<GridType::ctype,dim> X(x); 
//    std::cout << X << " " << p(X) << std::endl;
//  }


//  testclass.coutElementPositions();
//  testclass.coutCornerPositions();
//  testclass.coutCornerByElement();

//  typedef GV::Codim<0>::Iterator asdf;

//  LGV::Codim<0>::Iterator it;

//  std::cout << "Now lets print the centers of all elements" << std::endl;
  
//  int counter = 0;



}

