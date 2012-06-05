


#define PI 3.1415
typedef double Real;

template<class GV, class GFS, class U, class PG, class C, class S>
void calcIonFlux(const GV& gv, const GFS& gfs, const U& u, const PG& pg, C& ip, C& im, S& s) {

  const int dim = GV::dimensionworld;
  
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> PhiGFSS;
  PhiGFSS phiGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<PhiGFSS, U> PhiDGF;
  PhiDGF phiDGF(phiGFSS, u);
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFSS, U> gradPhiF(phiGFSS,u);

  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,1> CpGFSS;
  CpGFSS cpGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CpGFSS, U> CpDGF;
  CpDGF cpDGF(cpGFSS, u);
  typename Dune::PDELab::DiscreteGridFunctionGradient < CpGFSS, U> gradCpF(cpGFSS,u);
  
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,2> CmGFSS;
  CmGFSS cmGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CmGFSS, U> CmDGF;
  CmDGF cmDGF(cmGFSS, u);
  typename Dune::PDELab::DiscreteGridFunctionGradient < CmGFSS, U> gradCmF(cmGFSS,u);

  
  typedef typename PhiDGF::Traits::RangeType phiRT;
  phiRT phi;
  typedef typename PhiDGF::Traits::RangeType cpRT;
  phiRT cp;
  typedef typename PhiDGF::Traits::RangeType cmRT;
  phiRT cm;


//  Dune::FieldVector<typename PhiDGF::Traits::RangeType,dim> gradphi;
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFSS, U>::Traits::RangeType gradphi;
  typename Dune::PDELab::DiscreteGridFunctionGradient < CpGFSS, U>::Traits::RangeType gradCp;
  typename Dune::PDELab::DiscreteGridFunctionGradient < CmGFSS, U>::Traits::RangeType gradCm;
  
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFSS, U>::Traits::RangeType i;

  const int codim=0;
  typedef typename GV::template Codim<0>::template Partition
       <Dune::Interior_Partition>::Iterator LeafIterator;

//    // loop over quadrature points
  for (LeafIterator it = 
            gv.template begin<0,Dune::Interior_Partition>();
            it!=gv.template end<0,Dune::Interior_Partition>(); ++it) {

    for (typename GV::IntersectionIterator ii = gv.ibegin(*it); ii!=gv.iend(*it); ++ii) {
            Dune::FieldVector<Real, dim> evalPos = 
              ii->geometry().center();
            Dune::FieldVector<Real, dim> local = 
              it->geometry().local(evalPos);

            phiDGF.evaluate(*it, local, phi);
            cpDGF.evaluate(*it, local, cp);
            cmDGF.evaluate(*it, local, cm);

            gradPhiF.evaluate(*it, local, gradphi);
            gradCpF.evaluate(*it, local, gradCp);
            gradCmF.evaluate(*it, local, gradCm);


            Real factor = ii->geometry().volume();
            if (s.cylindrical)
                factor *= 2*PI*it->geometry().global(local)[1]; 
            gradCp *= -factor;
            gradCm *= -factor;
            gradphi *= factor;

            gradphi*=cp;

      if (ii->boundary()) {
            int physgroup_index = pg[ii->boundarySegmentIndex()]; 
            ip[physgroup_index][0] += ( gradCp + gradphi )*ii->unitOuterNormal(ii->geometry().local(evalPos));
            std::cout  << it->geometry().global(local)  << " "
              << phi << " " << cp << " " << cm << " " << gradphi << " " << gradCp << " " << gradCm << " flux" << std::endl;

            gradphi*=cm/cp;
            im[physgroup_index][0] += ( gradCm - gradphi )*ii->unitOuterNormal(ii->geometry().local(evalPos));

 //           std::cout << "flux " << it->geometry().global(local) << " " << ( gradCp ) * ii->outerNormal(ii->geometry().local(evalPos)) << " " << ( gradCm ) * ii->outerNormal(ii->geometry().local(evalPos))  << " " <<  ( gradphi ) * ii->outerNormal(ii->geometry().local(evalPos)) << " " << cp << " " << cm << " " << physgroup_index << std::endl;


//            im[physgroup_index] += ( -1* gradCm + cm*gradphi )*factor;
      }


    }
  }
}

 
