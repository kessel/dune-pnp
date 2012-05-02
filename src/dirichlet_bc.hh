


template<typename GV, typename RF, typename PGMap, class S, const int component>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension<GV,RF,PGMap, S, component> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  //! construct from grid view
  BCExtension(const GV& gv_, const PGMap& pg_, const S& s_) : gv(gv_), pg(pg_), s(s_) {
    t = 0;
  }

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //! Set value for potential at outer domain boundaries
    //
    typedef typename GV::ctype ctype;
    const int dim=GV::dimensionworld;

    Dune::FieldVector<ctype,dim> integrationPointGlobal =  e.geometry().global(xlocal);
    int counter = 0;
    ctype dist;
    ctype this_dist;
    int physgroup_index = -1;
    Dune::FieldVector<ctype,dim> xglobal;
    Dune::FieldVector<ctype,dim> dist_vec;
    for (IntersectionIterator ii = gv.ibegin(e); ii != gv.iend(e) ; ++ii) {
      if (ii->boundary()) {
      xglobal = ii->geometry().center();
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii->boundarySegmentIndex()];
      xglobal = ii->geometry().corner(0);
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii->boundarySegmentIndex()];
      xglobal = ii->geometry().corner(1);
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii->boundarySegmentIndex()];
      } else {


      typename GV::Traits::Grid::template Codim<0>::EntityPointer o( ii->outside() );
      counter++;
      for  (IntersectionIterator ii2 = gv.ibegin(*o); ii2 != gv.iend(*o) ; ++ii2) {
        if (ii2->boundary()) {
      xglobal = ii2->geometry().center();
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii2->boundarySegmentIndex()];
      xglobal = ii2->geometry().corner(0);
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii2->boundarySegmentIndex()];
      xglobal = ii2->geometry().corner(1);
      dist_vec = xglobal - integrationPointGlobal;
      dist=dist_vec * dist_vec;
      if (dist < 1e-3) 
        physgroup_index = pg[ii2->boundarySegmentIndex()];
        }
      }
      }

      }
    if (physgroup_index >= 0) {
        switch (component){
          case 0:
            if (s.surfaces[physgroup_index].coulombBtype == 0) {
              y = s.surfaces[physgroup_index].coulombPotential; break;
            }
          case 1:
            if (s.surfaces[physgroup_index].plusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].plusDiffusionConcentration; break;
            }
          case 2:
            if (s.surfaces[physgroup_index].minusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].minusDiffusionConcentration; break;
            }
      }
    }
    if (s.verbosity > 1) {
       std::cout << "boundary " << integrationPointGlobal << " type " << physgroup_index << std::endl;
    }
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

  void setTime(double t_) {
    t=t_;
  }

private :

  const GV& gv;
  const PGMap& pg;
  const S& s;
  double t;
};





