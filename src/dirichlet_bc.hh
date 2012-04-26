


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

    int physgroup_index =-1;
    int counter = 0;
    for (IntersectionIterator ii = gv.ibegin(e); ii != gv.iend(e) ; ++ii) {
      counter++;
      std::cout << ii->geometry().center() << " " << (int) ii->boundary() << std::endl;
      if (ii->boundary() ) {
        physgroup_index = pg[ii->boundarySegmentIndex()];
        switch (component){
          case 0:
            if (s.surfaces[physgroup_index].coulombBtype == 0) {
              std::cout << ii->geometry().center() << " " << e.geometry().center() << " " << physgroup_index << " asdf" << std::endl;
              y = s.surfaces[physgroup_index].coulombPotential; return;
            }
          case 1:
            if (s.surfaces[physgroup_index].plusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].plusDiffusionConcentration; return;
            }
          case 2:
            if (s.surfaces[physgroup_index].minusDiffusionBtype == 0) {
              y = s.surfaces[physgroup_index].minusDiffusionConcentration; return;
            }
        }
      }
    }
    std::cout << e.geometry().global(xlocal) << " " << e.geometry().center() << " " << counter << " inner" << std::endl;
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

