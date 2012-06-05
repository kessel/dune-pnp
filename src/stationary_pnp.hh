
#include "datawriter.hh"
#include <dune/common/mpihelper.hh>
#include <string>
#include"pnp_operator.hh"
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<ionFlux.hh>
#include<dune/pdelab/newton/newton.hh>

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

  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> PhiFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CpFEM;
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, Coord,Real, 1> CmFEM;
  PhiFEM phiFem(gv);
  CpFEM cpFEM(gv);
  CmFEM cmFEM(gv);

  
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;     // constraints class
  CON phiCon;
  CON cplusCon;
  CON cminusCon;
  typedef Dune::PDELab::SimpleGridFunctionStaticSize GFSSize; // what is that?

  typedef Dune::PDELab::GridFunctionSpace
   <GV, PhiFEM, CON, VectorBackend, GFSSize> PhiGFS;
  PhiGFS phiGFS(gv, phiFem, phiCon);
  phiCon.compute_ghosts(phiGFS);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CpFEM, CON, VectorBackend, GFSSize> CpGFS;
  CpGFS cpGFS(gv, cpFEM, cplusCon);
  cplusCon.compute_ghosts(cpGFS);

  typedef Dune::PDELab::GridFunctionSpace
   <GV, CmFEM, CON, VectorBackend, GFSSize> CmGFS;
  CmGFS cmGFS(gv, cmFEM, cminusCon);
  cminusCon.compute_ghosts(cmGFS);

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

  typedef Dune::PDELab::CompositeConstraintsParameters<PhiBC_T, CpBC_T, CmBC_T> BT;
//  typedef Dune::PDELab::CompositeGridFunction<PhiBC_T, CpBC_T, CmBC_T> BT;
  BT bt(phiB_t, cpB_t, cmB_t);

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bt,gfs,cc); 

  std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cc.size()
    << " of " << gfs.globalSize() << std::endl;


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

//  typedef typename GFS::template VectorContainer<Real>::Type U;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,double>::Type U;
  U u(gfs ,0.0);
  Dune::PDELab::set_shifted_dofs(cc,0.0,u);

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
  
  typedef PnpOperator<PhiBC_T, CpBC_T, CmBC_T, int, Sysparams, FluxContainer> LOP;
  LOP lop(phiB_t, cpB_t, cmB_t, f, s, fluxContainer);

  
  
//  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
//  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
//  GOS gos(gfs,cc,gfs,cc,lop);
//
//  // <<<5>>> Select a linear solver backend
//  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
//  LS ls(5000,true);
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

    // <<<5a>>> Select a linear solver backend
//    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS; 
//    LS ls( gfs, 5000, 5, s.verbosity );
//    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
//    LS ls( gfs, 5000, s.verbosity );

    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_Jacobi<GFS> LS;
    LS ls( gfs, 5000, s.verbosity );
//
//    
//    ISTLBackend_AMG_NOVLP
//    typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<GOS,double> LS;
//    LS ls( gfs, 2, 5000, s.verbosity );
//    
//    
//    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
//    LS ls( gfs, 5000, s.verbosity );


// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
//typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
//typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
//typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SOR<GFS> LS;


    // <<<5b>>> Solve nonlinear problem
    typedef Dune::PDELab::Newton<GO,LS,U> NEWTON;
    NEWTON newton(go,u,ls);
    newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
//    newton.setLineSearchStrategy(newton.noLineSearch);
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(s.verbosity);
    newton.setReduction(1e-5);
    newton.setMinLinearReduction(1e-2);
    newton.setMaxIterations(5000);
    newton.setLineSearchMaxIterations(50);
    newton.apply();


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

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>(phiDGF,"Phi"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<CpDGF>(cpDGF,"Cp"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<CmDGF>(cmDGF,"Cm"));
  vtkwriter.write("solution",Dune::VTKOptions::binaryappended);

//  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
//  gnuplotwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PhiDGF>,"solution");
//  gnuplotwriter.write("yeah.dat"); 

  
  DataWriter<GV, double, 0> dw(gv, helper);
  std::string s1("solution");
  std::string s2("solution");
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

