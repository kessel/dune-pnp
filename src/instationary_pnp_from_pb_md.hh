
#include "datawriter.hh"
#include <dune/common/mpihelper.hh>
#include <string>
#include"pnp_operator.hh"
#include"pnp_toperator.hh"
#include"pb_operator.hh"
#include "poisson_operator.hh"
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<ionFlux.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/backend/novlpistlsolverbackend.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include"datawriter.hh"
#include "diffusion_operator.hh"
#include "diffusion_toperator.hh"
//#include<dune/istl/superlu.hh>
//#include<dune/pdelab/backend/seqistlsolverbackend.hh>
//
#define BCGS_SSORk    1
#define BCGS_NOPREC   2
#define CG_NOPREC     3
#define CG_Jacobi     4
#define CG_AMG_SSOR   5

#ifndef PDEGREE 
#define PDEGREE 1
#endif

#ifndef LINEARSOLVER
#define LINEARSOLVER BCGS_SSORk
#endif

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
void instationary_pnp_md(Sysparams s, GV gv, std::vector<int> boundaryIndexToEntity, std::vector<int> elementIndexToEntity, Dune::MPIHelper& helper) {
  
  typedef typename GV::Grid::ctype Coord;
  typedef typename GV::Grid GridType;
  typedef double Real;

  typedef std::vector<int> PG;
  
  typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;



  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, PDEGREE> PbFEM;
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
  
#if LINEARSOLVER == BCGS_SSORk
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<PbGO> PbLS; 
  PbLS pbls( pbGFS, s.linearSolverIterations, 1, s.verbosity );
#endif

#if LINEARSOLVER == BCGS_NOPREC
    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<PbGFS> PbLS;
    PbLS pbls( pbGFS, s.linearSolverIterations, s.verbosity );
#endif

#if LINEARSOLVER == CG_NOPREC
    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<PbGFS> PbLS;
    PbLS pbls( pbGFS, s.linearSolverIterations, s.verbosity );
#endif

#if LINEARSOLVER == CG_Jacobi
    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_Jacobi<PbGFS> PbLS;
    PbLS pbls( pbGFS, s.linearSolverIterations, s.verbosity );
#endif

#if LINEARSOLVER == CG_AMG_SSOR
    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<PbGOS,double> PbLS;
    PbLS pbls( pbGFS, 2, s.linearSolverIterations, s.verbosity );
#endif


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

//  DataWriter<GV, double, 0> dw(gv, helper);
//  std::string s1("pb");
//  std::string s2("pb");
//  dw.writeIpbsCellData(pbDGF, pbu, s1, s2);


  // Here comes the PNP Part:
  
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, PDEGREE> PhiFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, PDEGREE> CpFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, PDEGREE> CmFEM;
  PhiFEM phiFem(gv);
  CpFEM cpFEM(gv);
  CmFEM cmFEM(gv);

  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;     // constraints class
  CON phicon;
  CON cpcon;
  CON cmcon;

  typedef Dune::PDELab::GridFunctionSpace
   <GV, PhiFEM, CON, VectorBackend, GFSSize> PhiGFS;
  PhiGFS phiGFS(gv, phiFem, phicon);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CpFEM, CON, VectorBackend, GFSSize> CpGFS;
  CpGFS cpGFS(gv, cpFEM, cpcon);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CmFEM, CON, VectorBackend, GFSSize> CmGFS;
  CmGFS cmGFS(gv, cmFEM, cmcon);

  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;

//  typedef Dune::PDELab::CompositeGridFunctionSpace<GFMapper,PhiGFS, CpGFS, CmGFS> GFS;
//  GFS gfs(phiGFS, cpGFS, cmGFS);

  phicon.compute_ghosts(phiGFS);
  cpcon.compute_ghosts(cpGFS);
  cmcon.compute_ghosts(cmGFS);

  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 0, PbDGF > PhiBC;
  PhiBC phiB(gv, boundaryIndexToEntity, s, pbDGF);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 1, PbDGF > CpBC;
  CpBC cpB(gv, boundaryIndexToEntity, s, pbDGF);
  typedef BCExtension<GV, Real,  std::vector<int>, Sysparams, 2, PbDGF > CmBC;
  CmBC cmB(gv, boundaryIndexToEntity, s, pbDGF);
//  typedef Dune::PDELab::CompositeGridFunction<PhiBC, CpBC, CmBC> BCE;
//  BCE bce(phiB, cpB, cmB);

  typedef BCType<GV, std::vector<int>, Sysparams, 0 > PhiBC_T;
  PhiBC_T phiB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 1 > CpBC_T;
  CpBC_T cpB_t(gv, boundaryIndexToEntity, s);
  typedef BCType<GV, std::vector<int>, Sysparams, 2 > CmBC_T;
  CmBC_T cmB_t(gv, boundaryIndexToEntity, s);

//  typedef Dune::PDELab::CompositeConstraintsParameters<PhiBC_T, CpBC_T, CmBC_T> BT;
//  typedef Dune::PDELab::CompositeGridFunction<PhiBC_T, CpBC_T, CmBC_T> BT;
//  BT bt(phiB_t, cpB_t, cmB_t);

  typedef typename PhiGFS::template ConstraintsContainer<Real>::Type PhiCC;
  typedef typename CpGFS::template ConstraintsContainer<Real>::Type  CpCC;
  typedef typename CpGFS::template ConstraintsContainer<Real>::Type  CmCC;
  PhiCC phicc;
  CpCC cpcc;
  CmCC cmcc;
  Dune::PDELab::constraints(phiB_t,phiGFS,phicc,false); 
  Dune::PDELab::constraints(cpB_t,cmGFS,cpcc,false); 
  Dune::PDELab::constraints(cmB_t,cpGFS,cmcc,false); 

//  std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cc.size()
//    << " of " << gfs.globalSize() << std::endl;

  int f = 5; 

//  typedef typename GFS::template VectorContainer<Real>::Type U;
  typedef typename Dune::PDELab::BackendVectorSelector<PhiGFS,double>::Type UPhi;
  typedef typename Dune::PDELab::BackendVectorSelector<PhiGFS,double>::Type UCp ;
  typedef typename Dune::PDELab::BackendVectorSelector<PhiGFS,double>::Type UCm ;
  UPhi uphi(phiGFS ,0.0);
  UCp ucp(cpGFS ,0.0);
  UCm ucm(cmGFS ,0.0);
//  Dune::PDELab::set_shifted_dofs(cc,0.0,u);
  typedef Dune::PDELab::DiscreteGridFunction<PhiGFS,UPhi> PhiDGF;
  PhiDGF phiDGF(phiGFS,uphi);
  typedef Dune::PDELab::DiscreteGridFunction<CpGFS,UCp> CpDGF;
  CpDGF cpDGF(cpGFS,ucp);
  typedef Dune::PDELab::DiscreteGridFunction<CmGFS,UCm> CmDGF;
  CmDGF cmDGF(cmGFS,ucm);


  Dune::PDELab::interpolate(phiB, phiGFS ,uphi);
  Dune::PDELab::interpolate(cpB, cpGFS ,ucp);
  Dune::PDELab::interpolate(cmB, cmGFS ,ucm);
  DataWriter<GV> dw(gv);
  dw.writeData(phiGFS, uphi, "phi.dat");
  dw.writeData(cpGFS, ucp, "cp.dat");
  dw.writeData(cmGFS, ucm, "cm.dat");

  typedef PoissonOperator<PhiBC_T, int, Sysparams, FluxContainer, CpDGF, CmDGF> PoissonLOP;
  PoissonLOP poissonlop(pbB_t, 1, s, fluxContainer, cpDGF, cmDGF);

  typedef Dune::PDELab::GridOperator<PhiGFS,PhiGFS,PoissonLOP,MBE, Real , Real , Real ,PhiCC,PhiCC> PhiGO;
  PhiGO phigo(phiGFS, phicc, phiGFS, phicc, poissonlop);

  typedef Dune::PDELab::StationaryLinearProblemSolver<PhiGO,PbLS,UPhi> SLP;
  SLP slp(phigo,uphi,pbls,1e-10);



  dw.writeData(phiGFS, uphi, "phi1.dat");


  typedef DiffusionOperator<CpBC_T, int, Sysparams, PhiGFS, UPhi> CpOP;
  CpOP cpop(cpB_t, 5, s, phiGFS, uphi, 1);
  typedef DiffusionOperator<CpBC_T, int, Sysparams, PhiGFS, UPhi> CmOP;
  CmOP cmop(cpB_t, 5, s, phiGFS, uphi, -1);

  typedef DiffusionTOperator CpTop;
  CpTop cptop(5);
  
  typedef DiffusionTOperator CmTop;
  CmTop cmtop(5);

  typedef Dune::PDELab::GridOperator<CpGFS,CpGFS,CpOP,MBE, Real , Real , Real ,CpCC,CpCC> CpGO0;
  CpGO0 cpgo0(cpGFS,cpcc,cpGFS,cpcc,cpop);
  typedef Dune::PDELab::GridOperator<CpGFS,CpGFS,CpTop,MBE, Real , Real , Real ,CpCC,CpCC> CpGO1;
  CpGO1 cpgo1(cpGFS,cpcc,cpGFS,cpcc,cptop);
  typedef Dune::PDELab::OneStepGridOperator<CpGO0,CpGO1> CpIGO; 
  CpIGO cpigo(cpgo0,cpgo1);
  
  typedef Dune::PDELab::GridOperator<CmGFS,CmGFS,CmOP,MBE, Real , Real , Real ,CmCC,CmCC> CmGO0;
  CmGO0 cmgo0(cmGFS,cmcc,cmGFS,cmcc,cmop);
  typedef Dune::PDELab::GridOperator<CmGFS,CmGFS,CmTop,MBE, Real , Real , Real ,CmCC,CmCC> CmGO1;
  CmGO1 cmgo1(cmGFS,cmcc,cmGFS,cmcc,cmtop);
  typedef Dune::PDELab::OneStepGridOperator<CmGO0,CmGO1> CmIGO; 
  CmIGO cmigo(cmgo0,cmgo1);


  typedef Dune::PDELab:: StationaryLinearProblemSolver<CpIGO,PbLS,UCp> CpPDESOLVER;
  CpPDESOLVER cppdesolver(cpigo, pbls, 1e-5);
  typedef Dune::PDELab:: StationaryLinearProblemSolver<CmIGO,PbLS,UCp> CmPDESOLVER;
  CmPDESOLVER cmpdesolver(cmigo, pbls, 1e-5);


  Dune::PDELab:: Alexander2Parameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,CpIGO,CpPDESOLVER,UCp,UCp> cposm(method,cpigo,cppdesolver);
  Dune::PDELab::OneStepMethod<Real,CmIGO,CmPDESOLVER,UCm,UCm> cmosm(method,cmigo,cmpdesolver);

  UCp ucpNew(ucp);
  UCm ucmNew(ucm);



      typedef typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFS, UPhi>::Traits::RangeType I;

      std::vector<I> ip;
      std::vector<I> im;
      I temp2(0);
      for (int i = 0; i<s.n_surfaces; i++) {
        ip.push_back(temp2);
        im.push_back(temp2);
      }
      std::ofstream currentfile;
      if (gv.comm().rank() == 0) {
        currentfile.open("current.dat");
      }



  double time=0;

  double dt=s.tau;
  int output_counter = 0;
  int output_freq=10;
  int potentialUpdateFrequency=10;
  char filename[20];
  for (int i = 0; i<s.nSteps; i++) { 
    cposm.apply(time,dt,ucp,cpB,ucpNew);
    ucp=ucpNew;
    cmosm.apply(time,dt,ucm,cmB,ucmNew);
    ucm=ucmNew;
    time+=dt;
    if (i % s.potentialUpdateFreq == 0) {
        slp.apply();
    }
    if (i % s.outputFreq == 0) {
      output_counter++;
      std::cout << "Writing output" << std::endl;
      sprintf(filename, "phi%03d.dat", output_counter);
      dw.writeData(phiGFS, uphi, std::string(filename));
      sprintf(filename, "cp%03d.dat", output_counter);
      dw.writeData(cpGFS, ucp, std::string(filename));
      sprintf(filename, "cm%03d.dat", output_counter);
      dw.writeData(cmGFS, ucm, std::string(filename));
  

      calcIonFlux<GV, PhiGFS, CpGFS, CmGFS, UPhi, UCp, UCm, PG, std::vector<I>, Sysparams >(gv, phiGFS, cpGFS, cmGFS, uphi, ucp, ucm, boundaryIndexToEntity, ip, im, s);
      if (gv.comm().rank()==0) {
        currentfile << time;
        for (unsigned int i = 0; i<s.n_surfaces; i++) {
           currentfile << " " << ip[i] << " " << im[i];
        }
        currentfile << std::endl;
        currentfile.flush();
      }
    }
  }
  slp.apply();





//  Dune::PDELab::set_nonconstrained_dofs(cc,0.1,u);

//  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
//  gnuplotwriter.addVertexData<template UPhi>(uphi,"phi");
//  gnuplotwriter.write("phi.dat"); 


 /* 
  std::string s1;
  std::string s2;
  s1 = s2 = "after_constraint_interpolation";

  dw.writePNPCellData(gfs, u, s1, s2); 
  
  typedef PnpOperator<PhiBC_T, CpBC_T, CmBC_T, int, Sysparams, FluxContainer> LOP;
  LOP lop(phiB_t, cpB_t, cmB_t, f, s, fluxContainer);

  typedef PnpTOperator<Sysparams> TLOP;
  TLOP tlop(s, s.tau);

  
  
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
    typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
    GO1 go1(gfs, gfs, tlop);
    typedef Dune::PDELab::OneStepGridOperator<GO,GO1, false> IGO;
    IGO igo(go,go1);



    


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


//  Dune::PDELab::Alexander2Parameter<Real> method;
  Dune::PDELab::ExplicitEulerParameter<Real> method;
//  Dune::PDELab::OneStepMethod<Real,IGO,NEWTON,U,U> osm(method,igo,newton);
//
  typedef Dune::PDELab::CFLTimeController<Real,IGO> TC;
  TC tc(0.001,igo);
//
  Dune::PDELab::ExplicitOneStepMethod<Real,IGO,LS,U,U,TC> osm(method,igo,ls,tc);

  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
//  Dune::PDELab::FilenameHelper fn("example05_Q1Q1");
//  {
//    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
//    U0DGF u0dgf(u0sub,u);
//    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
//    U1DGF u1dgf(u1sub,u);
//    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
//    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
//    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
//    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
//    fn.increment();
//  }
    typedef typename GO::template MatrixContainer<double>::Type M;
    M m(go);
    m = 0.0;
//    igo.jacobian(u,m);
//    if (s.printStiffnessMatrix)
//      Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);


  // <<<9>>> time loop
  U unew(gfs,0.0);
  unew = u;
  double dtstart=s.tau;
  double dt = dtstart;
  for (int i = 0; i<100; i++) //(time<tend-1e-8)
    {
      // do time step
      double time=i*dt;
      osm.apply(time,dt,u,unew);

//      // graphics
//      typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
//      U0DGF u0dgf(u0sub,unew);
//      typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
//      U1DGF u1dgf(u1sub,unew);
//      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
//      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
//      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
//      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
//      fn.increment();
//
      u = unew;
//      time += dt;
//      if (dt<dtmax-1e-8) dt = std::min(dt*1.1,dtmax);           // time step adaption
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
  
*/
}

