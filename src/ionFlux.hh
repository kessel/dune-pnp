


#define PI 3.1415
typedef double Real;

template<class GV, class PhiGFS, class CpGFS, class CmGFS, class UPhi, class UCp, class UCm, class PG, class C, class S>
void calcIonFlux(const GV& gv, const PhiGFS& phiGFS, const CmGFS & cmGFS, const CpGFS & cpGFS, const UPhi& uphi, const UCp& ucp, const UCm & ucm,  const PG& pg, C& ip, C& im, S& s) {

  const int dim = GV::dimensionworld;
  
  typedef Dune::PDELab::DiscreteGridFunction<PhiGFS, UPhi> PhiDGF;
  PhiDGF phiDGF(phiGFS, uphi);
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFS, UPhi> gradPhiF(phiGFS,uphi);

  typedef Dune::PDELab::DiscreteGridFunction<CpGFS, UCp> CpDGF;
  CpDGF cpDGF(cpGFS, ucp);
  typename Dune::PDELab::DiscreteGridFunctionGradient < CpGFS, UCp> gradCpF(cpGFS,ucp);
  
  typedef Dune::PDELab::DiscreteGridFunction<CmGFS, UCm> CmDGF;
  CmDGF cmDGF(cmGFS, ucm);
  typename Dune::PDELab::DiscreteGridFunctionGradient < CmGFS, UCm> gradCmF(cmGFS,ucm);

  
  typedef typename PhiDGF::Traits::RangeType phiRT;
  phiRT phi;
  typedef typename PhiDGF::Traits::RangeType cpRT;
  phiRT cp;
  typedef typename PhiDGF::Traits::RangeType cmRT;
  phiRT cm;


//  Dune::FieldVector<typename PhiDGF::Traits::RangeType,dim> gradphi;
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFS, UPhi>::Traits::RangeType gradphi;
  typename Dune::PDELab::DiscreteGridFunctionGradient < CpGFS, UCp>::Traits::RangeType gradCp;
  typename Dune::PDELab::DiscreteGridFunctionGradient < CmGFS, UCm>::Traits::RangeType gradCm;
  
  typename Dune::PDELab::DiscreteGridFunctionGradient < PhiGFS, UPhi>::Traits::RangeType i;

  const int codim=0;
  typedef typename GV::template Codim<0>::template Partition
       <Dune::Interior_Partition>::Iterator LeafIterator;

  for (int i = 0; i<im.size(); i++) {
      im[i][0]=0;
      ip[i][0]=0;
  }

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

 
