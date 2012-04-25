



template<typename GV, typename PGMap, class S, const int component>
class BCType : public Dune::PDELab::BoundaryGridFunctionBase<
        Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
        Dune::FieldVector<int,1> >,BCType<GV,PGMap, S, component> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view
  BCType (const GV& gv_, const PGMap& pg_, const S& s_) : gv(gv_), pg(pg_), s(s_) {
    t=0;
  }

  //! return bc type at point on intersection
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    /** use physical index to determine B.C.
    *   values are specified in .geo file
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for iPBS iterated boundaries  */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    switch (component) {
      case 0:
        if (s.surfaces[physgroup_index].coulombBtype==0) {
          y = 1; return;
        } else {
          y = 0; return;
        }
      case 1:
        if (s.surfaces[physgroup_index].plusDiffusionBtype==0) {
          y = 1; return;
        } else {
          y = 0; return;
        }
      case 2:
        if (s.surfaces[physgroup_index].plusDiffusionBtype==0) {
          y = 1; return;
        } else {
          y = 0; return;
        }
    return;
    }
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

    void setTime(double t_) {
    t=t_;
  }



private:

  const GV&    gv;
  const PGMap& pg;
  const S& s;
  
  double t;
};

