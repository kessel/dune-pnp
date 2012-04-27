
#include "datawriter.hh"
#include <dune/common/mpihelper.hh>
#include <string>
#include"pnp_operator.hh"
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>


template<class GV>
void stationary_pnp(Sysparams s, GV gv, std::vector<int> boundaryIndexToEntity, std::vector<int> elementIndexToEntity, Dune::MPIHelper& helper) {
  
  typedef typename GV::Grid::ctype Coord;
  typedef typename GV::Grid GridType;
  typedef double Real;

  const int dim = 2;

  
  typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;

  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> PhiFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CpFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CmFEM;
  PhiFEM phiFem(gv);
  CpFEM cpFEM(gv);
  CmFEM cmFEM(gv);

  
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  typedef Dune::PDELab::SimpleGridFunctionStaticSize GFSSize; // what is that?

  typedef Dune::PDELab::GridFunctionSpace
   <GV, PhiFEM, CON, VectorBackend, GFSSize> PhiGFS;
  PhiGFS phiGFS(gv, phiFem);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CpFEM, CON, VectorBackend, GFSSize> CpGFS;
  CpGFS cpGFS(gv, cpFEM);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CmFEM, CON, VectorBackend, GFSSize> CmGFS;
  CmGFS cmGFS(gv, cmFEM);

  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;

  typedef Dune::PDELab::CompositeGridFunctionSpace<GFMapper,PhiGFS, CpGFS, CmGFS> GFS;
  GFS gfs(phiGFS, cpGFS, cmGFS);

  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 0 > PhiBC;
  PhiBC phiB(gv, boundaryIndexToEntity, s);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 1 > CpBC;
  CpBC cpB(gv, boundaryIndexToEntity, s);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 2 > CmBC;
  CmBC cmB(gv, boundaryIndexToEntity, s);
  typedef Dune::PDELab::CompositeGridFunction<PhiBC, CpBC, CmBC> BCE;
  BCE bce(phiB, cpB, cmB);

  typedef BCType<GV, std::vector<int>, Sysparams, 0 > PhiBC_T;
  PhiBC_T phiB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 1 > CpBC_T;
  CpBC_T cpB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 2 > CmBC_T;
  CmBC_T cmB_t(gv, boundaryIndexToEntity, s);
  typedef Dune::PDELab::CompositeGridFunction<PhiBC_T, CpBC_T, CmBC_T> BT;
  BT bt(phiB_t, cpB_t, cmB_t);

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bt,gfs,cc); 




/*  typedef BCType<GV, std::vector<int>, Sysparams, 0 > B;              // boundary condition type
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
//  F f(s);*/

  
  typedef std::vector<std::vector<double > > FluxContainer;
  FluxContainer fluxContainer;
  
  std::vector<double> temp;
  temp.push_back(-1000);
  temp.push_back(-1000);
  temp.push_back(-1000);
  for (int i = 0; i<gv.size(1); i++) {
      fluxContainer.push_back(temp);
  } 

  typedef typename GV::template Codim<0>::Iterator GVIT;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  
  GVIT it = gv.template begin<0>();

  for (; it != gv.template end<0>(); ++it) {
    for (IntersectionIterator iit = gv.ibegin(*it); iit!=  gv.iend(*it); ++iit) {
      if (iit->boundary()) {
        int physgroup_index=boundaryIndexToEntity[iit->boundarySegmentIndex()];
        fluxContainer[iit->boundarySegmentIndex()][0]=s.surfaces[physgroup_index].coulombFlux;
        fluxContainer[iit->boundarySegmentIndex()][1]=s.surfaces[physgroup_index].plusDiffusionFlux;
        fluxContainer[iit->boundarySegmentIndex()][2]=s.surfaces[physgroup_index].minusDiffusionFlux;
      }
    }
  }


  int f = 5; 
  typedef PnpOperator<PhiBC, CpBC, CmBC, int, Sysparams, FluxContainer> LOP;
  LOP lop(phiB, cpB, cmB, f, s, fluxContainer);

  typedef typename GFS::template VectorContainer<Real>::Type U;
  U u(gfs ,0.0);

  Dune::PDELab::interpolate(bce,gfs ,u);

  
  
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
//  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> PhiGFSS;
  typedef PhiGFS PhiGFSS;
  PhiGFSS phiGFSS(gfs);
  typedef Dune::PDELab::VectorDiscreteGridFunction<PhiGFSS, U> PhiDGF;
  PhiDGF phiDGF(phiGFSS, u);

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phiDGF,"Phi"));
  vtkwriter.write("solution",Dune::VTKOptions::binaryappended);

  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>,"solution");
  gnuplotwriter.write("yeah.dat"); 
  
/*  int argc = 0;
//  char**& argv = 0;
  
  DataWriter<GV, double, 0> dw(gv, helper);
  std::string s1("solution");
  std::string s2("solution.dat");
  dw.writeIpbsCellData(gfs, u, s1, s2); */

}

