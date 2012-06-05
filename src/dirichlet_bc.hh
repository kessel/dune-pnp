


template<typename GV, typename RF, typename PGMap, class S, const int component>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension<GV,RF,PGMap, S, component> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::ctype ctype;

  //! construct from grid view
  BCExtension(const GV& gv_, const PGMap& pg_, const S& s_) : gv(gv_), pg(pg_), s(s_) {
    t = 0;
  }

  const inline bool global_on_intersection(Dune::FieldVector<ctype, GV::dimensionworld> integrationPointGlobal, IntersectionIterator& ii ) const 
  {
    ctype dist;
    ctype this_dist;


    if (GV::dimensionworld == 2) {
      Dune::FieldVector<ctype,GV::dimensionworld> p_vec = ii->geometry().corner(1) - ii->geometry().corner(0);
      p_vec/=p_vec.two_norm();
      Dune::FieldVector<ctype,GV::dimensionworld> p2 = p_vec;
      p2 *= ((integrationPointGlobal - ii->geometry().corner(0))*p_vec);
      Dune::FieldVector<ctype,GV::dimensionworld> dist_vec = p2 - (integrationPointGlobal - ii->geometry().corner(0));

      if (dist_vec.two_norm() < 1e-9)
        return true;
      return false;
    } else {
//      std::cout << "implement me in " << __FILE__ << __LINE__ << std:endl;
//      DUNE_THROW("Not implemented!");
    }

  }
  const inline int bctype(int physgroup_index) const {
    int bctype;
    switch (component){
    case 0:
      bctype=s.surfaces[physgroup_index].coulombBtype;
    case 1:
      bctype=s.surfaces[physgroup_index].plusDiffusionBtype;
    case 2:
      bctype=s.surfaces[physgroup_index].minusDiffusionBtype;;
    }
    return bctype;
  }

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //! Set value for potential at outer domain boundaries
    //

    Dune::FieldVector<ctype,GV::dimensionworld> integrationPointGlobal =  e.geometry().global(xlocal);
    int counter = 0;
    int physgroup_index = -1;
    for (IntersectionIterator ii = gv.ibegin(e); ii != gv.iend(e) ; ++ii) {
      if (ii->boundary()) {
        if (global_on_intersection(integrationPointGlobal, ii)) {
          if (physgroup_index==-1 || bctype(physgroup_index) != 0) {
            physgroup_index = pg[ii->boundarySegmentIndex()];
          }
        }
      } else if (ii->neighbor()) {
        typename GV::Traits::Grid::template Codim<0>::EntityPointer o( ii->outside() );
        for  (IntersectionIterator ii2 = gv.ibegin(*o); ii2 != gv.iend(*o) ; ++ii2) {
          if (ii2->boundary()) {
            if (global_on_intersection(integrationPointGlobal, ii2)) {
              if (physgroup_index==-1 || bctype(physgroup_index) != 0) {
                physgroup_index = pg[ii2->boundarySegmentIndex()];
              } else {
                int bctype;
                switch (component){
                  case 0:
                    bctype=s.surfaces[physgroup_index].coulombBtype;
                  case 1:
                    bctype=s.surfaces[physgroup_index].plusDiffusionBtype;
                  case 2:
                    bctype=s.surfaces[physgroup_index].minusDiffusionBtype;;
                }
                if (bctype != 0)
                  physgroup_index = pg[ii2->boundarySegmentIndex()];
              }
            }
          }
        }
      }
    }
        switch (component){
          case 0:
            if (physgroup_index > -1 && s.surfaces[physgroup_index].coulombBtype == 0) {
              y = s.surfaces[physgroup_index].coulombPotential;
            } else {
              y = 0;
            }
 break;
          case 1:
            if (physgroup_index > -1 && s.surfaces[physgroup_index].plusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].plusDiffusionConcentration; 
            } else {
              y = 0.06;
            }
 break;
          case 2:
            if (physgroup_index > -1 && s.surfaces[physgroup_index].minusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].minusDiffusionConcentration; 
            } else {
              y = 0.06;
            }
 break;
      }
    if (s.verbosity > 1) {
       //std::cout << "boundary " << integrationPointGlobal << " type " << physgroup_index << " " << y << std::endl;
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





