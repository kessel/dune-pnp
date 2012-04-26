
#include "datawriter.hh"
#include <dune/common/mpihelper.hh>
#include <string>

template<class GV>
void stationary_diffusion(Sysparams s, GV gv, std::vector<int> boundaryIndexToEntity, std::vector<int> elementIndexToEntity, Dune::MPIHelper& helper) {
  typedef typename GV::Grid::ctype Coord;
  typedef typename GV::Grid GridType;
  typedef double Real;

  const int dim = 2;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> FEM;
  FEM fem(gv);
  //typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real, dim> FEM;
  //FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  typedef BCType<GV, std::vector<int>, Sysparams, 0 > B;                                         // boundary condition type
  B b(gv, boundaryIndexToEntity, s);

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); 

  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 0 > BCE;
  BCE bce(gv, boundaryIndexToEntity, s);


  typedef typename GFS::template VectorContainer<Real>::Type U;
  U u(gfs ,0.0);

  Dune::PDELab::interpolate(bce,gfs ,u);
      
//  typedef Force<Sysparams, template Dune::FieldVector<GridType::ctype,dim> > F;
//  F f(s);
//
  typedef std::vector<double> FluxContainer;
  FluxContainer fluxContainer;
  for (int i = 0; i<gv.size(1); i++) {
      fluxContainer.push_back(-1000);
  }

  typedef typename GV::template Codim<0>::Iterator GVIT;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  
  GVIT it = gv.template begin<0>();

  for (; it != gv.template end<0>(); ++it) {
    for (IntersectionIterator iit = gv.ibegin(*it); iit!=  gv.iend(*it); ++iit) {
      if (iit->boundary()) {
        int physgroup_index=boundaryIndexToEntity[iit->boundarySegmentIndex()];
        fluxContainer[iit->boundarySegmentIndex()]=s.surfaces[physgroup_index].coulombFlux;
      }
    }
  }


  int f = 5;
  typedef DiffOperator<B, int, Sysparams, FluxContainer> LOP;
  LOP lop(b, f, s, fluxContainer);



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
  
  int argc = 0;
//  char**& argv = 0;
  
  DataWriter<GV, double, 0> dw(gv, helper);
  std::string s1("solution");
  std::string s2("solution.dat");
  dw.writeIpbsCellData(gfs, u, s1, s2);

}

