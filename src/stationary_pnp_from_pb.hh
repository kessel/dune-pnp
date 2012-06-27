
#include "datawriter.hh"
#include <dune/common/mpihelper.hh>
#include <string>
#include"pnp_operator.hh"
#include"pb_operator.hh"
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<ionFlux.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/backend/novlpistlsolverbackend.hh>
//#include<dune/istl/superlu.hh>
//#include<dune/pdelab/backend/seqistlsolverbackend.hh>

template<typename GV, typename RF>
class Phi0Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, Phi0Initial<GV,RF> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  Phi0Initial (const GV& gv_) 
    : gv(gv_)  
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    y=0;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

template<typename GV, typename RF>
class CpInitial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, CpInitial<GV,RF> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  CpInitial (const GV& gv_) 
    : gv(gv_)  
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    y=0.06;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

template<typename GV, typename RF>
class CmInitial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, CmInitial<GV,RF> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  CmInitial (const GV& gv_) 
    : gv(gv_)  
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    y=0.06;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

template<class GV>
void stationary_pnp(Sysparams s, GV gv, std::vector<int> boundaryIndexToEntity, std::vector<int> elementIndexToEntity, Dune::MPIHelper& helper) {
  
  typedef typename GV::Grid::ctype Coord;
  typedef typename GV::Grid GridType;
  typedef double Real;

  typedef std::vector<int> PG;
  
  typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;



  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> PbFEM;
  PbFEM pbFem(gv);
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints PbCON;     // constraints class
  PbCON pbCon;
  typedef Dune::PDELab::SimpleGridFunctionStaticSize GFSSize; // what is that?
  typedef Dune::PDELab::GridFunctionSpace
   <GV, PbFEM, PbCON, VectorBackend, GFSSize> PbGFS;
  PbGFS pbGFS(gv, pbFem, pbCon);

  pbCon.compute_ghosts(pbGFS);
  
  typedef BCType<GV, std::vector<int>, Sysparams, 0 > PbBC_T;
  PbBC_T pbB_t(gv, boundaryIndexToEntity, s);
 
//  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, -1 > PbBC;
//  PbBC pbB(gv, boundaryIndexToEntity, s);
  
  typedef typename PbGFS::template ConstraintsContainer<Real>::Type PbCC;
  PbCC pbcc;
  Dune::PDELab::constraints(pbB_t,pbGFS,pbcc,false); 

  typedef typename Dune::PDELab::BackendVectorSelector<PbGFS,double>::Type PbU;
  PbU pbu(pbGFS,0.0);



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


  typedef PBOperator<PbBC_T, int, Sysparams, FluxContainer> PbLOP;
  PbLOP pblop(pbB_t, 1, s, fluxContainer);


  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  
  typedef Dune::PDELab::GridOperator     <PbGFS,PbGFS,PbLOP,MBE,double, double, double,PbCC,PbCC,true> PbGO;
  PbGO pbgo(pbGFS,pbcc,pbGFS,pbcc,pblop);
  
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<PbGO> PbLS; 
  PbLS pbls( pbGFS, s.linearSolverIterations, 1, s.verbosity );

  typedef Dune::PDELab::Newton<PbGO,PbLS,PbU> PbNEWTON;
  PbNEWTON pbnewton(pbgo,pbu,pbls);
  pbnewton.setLineSearchStrategy(pbnewton.hackbuschReuskenAcceptBest);
//pb    newton.setLineSearchStrategy(newton.noLineSearch);
  pbnewton.setReassembleThreshold(s.newtonReassembleThreshold);
  pbnewton.setVerbosityLevel(s.verbosity);
  pbnewton.setReduction(s.newtonReduction);
  pbnewton.setMinLinearReduction(s.newtonMinLinearReduction);
  pbnewton.setMaxIterations(s.newtonMaxIterations);
  pbnewton.setLineSearchMaxIterations(s.newtonLineSearchMaxIteration);
  try {
    pbnewton.apply();
  } catch (...) {
      std::cout << "Something has happened" << std::endl;
  }

  typedef Dune::PDELab::DiscreteGridFunction<PbGFS, PbU> PbDGF;
  PbDGF pbDGF(pbGFS, pbu);

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PbDGF>(pbDGF,"pb"));
  vtkwriter.write("pb",Dune::VTKOptions::binaryappended);

//  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
//  gnuplotwriter.addVertexData(pbu,"pb");
//  gnuplotwriter.write("pb.dat"); 

//  DataWriter<GV, double, 0> dw(gv, helper);
//  std::string s1("pb");
//  std::string s2("pb");
//  dw.writeIpbsCellData(pbDGF, pbu, s1, s2);


  // Here comes the PNP Part:
  
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> PhiFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CpFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CmFEM;
  PhiFEM phiFem(gv);
  CpFEM cpFEM(gv);
  CmFEM cmFEM(gv);

  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;     // constraints class
  CON con;

  typedef Dune::PDELab::GridFunctionSpace
   <GV, PhiFEM, CON, VectorBackend, GFSSize> PhiGFS;
  PhiGFS phiGFS(gv, phiFem, con);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CpFEM, CON, VectorBackend, GFSSize> CpGFS;
  CpGFS cpGFS(gv, cpFEM, con);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CmFEM, CON, VectorBackend, GFSSize> CmGFS;
  CmGFS cmGFS(gv, cmFEM, con);

  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;

  typedef Dune::PDELab::CompositeGridFunctionSpace<GFMapper,PhiGFS, CpGFS, CmGFS> GFS;
  GFS gfs(phiGFS, cpGFS, cmGFS);

  con.compute_ghosts(gfs);

  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 0, PbDGF > PhiBC;
  PhiBC phiB(gv, boundaryIndexToEntity, s, pbDGF);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 1, PbDGF > CpBC;
  CpBC cpB(gv, boundaryIndexToEntity, s, pbDGF);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 2, PbDGF > CmBC;
  CmBC cmB(gv, boundaryIndexToEntity, s, pbDGF);
  typedef Dune::PDELab::CompositeGridFunction<PhiBC, CpBC, CmBC> BCE;
  BCE bce(phiB, cpB, cmB);

  typedef BCType<GV, std::vector<int>, Sysparams, 0 > PhiBC_T;
  PhiBC_T phiB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 1 > CpBC_T;
  CpBC_T cpB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 2 > CmBC_T;
  CmBC_T cmB_t(gv, boundaryIndexToEntity, s);

  typedef Dune::PDELab::CompositeConstraintsParameters<PhiBC_T, CpBC_T, CmBC_T> BT;
//  typedef Dune::PDELab::CompositeGridFunction<PhiBC_T, CpBC_T, CmBC_T> BT;
  BT bt(phiB_t, cpB_t, cmB_t);

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bt,gfs,cc,false); 

  std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cc.size()
    << " of " << gfs.globalSize() << std::endl;



  int f = 5; 

//  typedef typename GFS::template VectorContainer<Real>::Type U;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,double>::Type U;
  U u(gfs ,0.0);
//  Dune::PDELab::set_shifted_dofs(cc,0.0,u);

  typedef Phi0Initial<GV, double> Phi;
  Phi phi(gv);
  typedef CpInitial<GV, double> Cp;
  Cp cp(gv);
  typedef CmInitial<GV, double> Cm;
  Cm cm(gv);

  typedef Dune::PDELab::CompositeGridFunction<Phi, Cp, Cm> InitialValues;
  InitialValues initial(phi, cp, cm);

//  Dune::PDELab::interpolate(initial,gfs ,u);
  Dune::PDELab::interpolate(bce,gfs ,u);
//  Dune::PDELab::set_nonconstrained_dofs(cc,0.1,u);

  
  DataWriter<GV, double, 0> dw(gv, helper);
  std::string s1;
  std::string s2;
  s1 = s2 = "after_constraint_interpolation";

  dw.writePNPCellData(gfs, u, s1, s2); 
  
  typedef PnpOperator<PhiBC_T, CpBC_T, CmBC_T, int, Sysparams, FluxContainer> LOP;
  LOP lop(phiB_t, cpB_t, cmB_t, f, s, fluxContainer);

  
  
//  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
//  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
//  GOS gos(gfs,cc,gfs,cc,lop);
//
//  // <<<5>>> Select a linear solver backend
//  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
//  LS ls(s.linearSolverIterations,true);
//
//  // <<<6>>> assemble and solve linear problem
//  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
//  SLP slp(gos,u,ls,1e-10);
//  slp.apply();


      typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
   // typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
   // GOS gos(gfs,cc,gfs,cc,lop);
    typedef Dune::PDELab::GridOperator     <GFS,GFS,LOP,MBE,double, double, double,CC,CC,true> GO;
    GO go(gfs,cc,gfs,cc,lop);
    
    typedef typename GO::template MatrixContainer<double>::Type M;
    M m(go);
    m = 0.0;
    go.jacobian(u,m);
//    Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);


    // <<<5a>>> Select a linear solver backend
//    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS; 
//    LS ls( gfs, s.linearSolverIterations, 1, s.verbosity );
    
    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
//    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
    LS ls( gfs, s.linearSolverIterations, s.verbosity );

//    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_Jacobi<GFS> LS;
//    LS ls( gfs, s.linearSolverIterations, s.verbosity );
//
//    
//    ISTLBackend_AMG_NOVLP
//    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<GOS,double> LS;
//    LS ls( gfs, 2, s.linearSolverIterations, s.verbosity );
//    
//    
//    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
//    LS ls( gfs, s.linearSolverIterations, s.verbosity );


// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
///typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
///   LS(gfs);
//typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
//typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SOR<GFS> LS;
//LS ls(gfs);


    // <<<5b>>> Solve nonlinear problem
    typedef Dune::PDELab::Newton<GO,LS,U> NEWTON;
    NEWTON newton(go,u,ls);
    newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
//    newton.setLineSearchStrategy(newton.noLineSearch);
    newton.setReassembleThreshold(s.newtonReassembleThreshold);
    newton.setVerbosityLevel(s.verbosity);
    newton.setReduction(s.newtonReduction);
    newton.setMinLinearReduction(s.newtonMinLinearReduction);
    newton.setMaxIterations(s.newtonMaxIterations);
    newton.setLineSearchMaxIterations(s.newtonLineSearchMaxIteration);
    try {
      newton.apply();
    } catch (...) {
        std::cout << "Something has happened" << std::endl;
    }


  // <<<7>>> graphical output
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> PhiGFSS;
  PhiGFSS phiGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<PhiGFSS, U> PhiDGF;
  PhiDGF phiDGF(phiGFSS, u);
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,1> CpGFSS;
  CpGFSS cpGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CpGFSS, U> CpDGF;
  CpDGF cpDGF(cpGFSS, u);
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,2> CmGFSS;
  CmGFSS cmGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CmGFSS, U> CmDGF;
  CmDGF cmDGF(cmGFSS, u);

//  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phiDGF,"Phi"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<CpDGF>(cpDGF,"Cp"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<CmDGF>(cmDGF,"Cm"));
  vtkwriter.write("solution",Dune::VTKOptions::binaryappended);

//  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
//  gnuplotwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>,"solution");
//  gnuplotwriter.write("yeah.dat"); 
  s1 = s2 = "sol";

  dw.writePNPCellData(gfs, u, s1, s2); 

  

 

  typedef typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFSS, U>::Traits::RangeType I;
  std::vector<I> ip;
  std::vector<I> im;
  I temp2(0);
  for (int i = 0; i<s.n_surfaces; i++) {
    ip.push_back(temp2);
    im.push_back(temp2);
  }


//    typedef typename GV::template Codim<0>::template Partition
//       <Dune::Interior_Partition>::Iterator LeafIterator;
//
//
//
//
//    for (LeafIterator it = 
//            gv.template begin<0,Dune::Interior_Partition>();
//            it!=gv.template end<0,Dune::Interior_Partition>(); ++it) {
//
//    std::cout << "yeah" << std::endl;
//    for (typename GV::IntersectionIterator ii = gv.ibegin(*it); ii!=gv.ibegin(*it); ++ii) {
//      std::cout << "intersectino" << std::endl;
//      if (ii->boundary()) {
//        std::cout << "boundary" << std::endl;
//      }
//    }
//    }
 

  //std::cout << "------- ION CURRENTS --------" << std::endl;
  //calcIonFlux<GV, GFS, U, PG, std::vector<I>, Sysparams >(gv, gfs, u, boundaryIndexToEntity, ip, im, s);
  //for (unsigned int i = 0; i<s.n_surfaces; i++) {
  //  std::cout << i << " " << ip[i] << " " << im[i] << std::endl;
  //}
  

}

